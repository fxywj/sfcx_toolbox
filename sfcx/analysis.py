# -*- coding: utf-8 -*-
"""Analyses that can be run at any point of the coordinate system."""

from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.cluster import hierarchy as sch
from scipy.spatial import cKDTree

from .data import vec_to_mat


# ------------------------------------------------------------- null models --

def _random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


class NullEngine:
    """Surrogate maps that respect spatial autocorrelation.

    kind='spin'        spherical rotation, needs parcel coordinates
    kind='msr'         Moran spectral randomisation (eigen sign flips)
    kind='permutation' plain shuffling, anti-conservative, for smoke tests
    """

    def __init__(self, cohort, kind="msr", n_perm=1000, seed=0, k=15):
        self.cohort = cohort
        self.n_perm = int(n_perm)
        self.seed = int(seed)
        self.kind = kind
        if kind == "spin" and cohort.coords is None:
            raise ValueError("spin null requires cohort.coords (n_nodes, 3)")
        self._prep(k)

    def _prep(self, k):
        n = self.cohort.n_nodes
        self._base = None
        if self.kind == "spin":
            xyz = self.cohort.coords
            xyz = xyz - xyz.mean(axis=0)
            r = np.linalg.norm(xyz, axis=1, keepdims=True)
            r[r < 1e-12] = 1.0
            self._xyz = xyz / r
        elif self.kind == "msr":
            a = self.cohort.adjacency()
            d = np.sqrt(((self.cohort.hierarchy[:, None] - self.cohort.hierarchy[None, :]) ** 2))
            w = a * np.exp(-d)
            np.fill_diagonal(w, 0.0)
            m = np.eye(n) - np.ones((n, n)) / n
            wc = m @ w @ m
            wc = 0.5 * (wc + wc.T)
            lam, v = np.linalg.eigh(wc)
            self._lam, self._v = lam, v

    def surrogates(self, v):
        v = np.asarray(v, dtype=float).ravel()
        rng = np.random.default_rng(self.seed)
        out = np.empty((self.n_perm, v.size))
        if self.kind == "permutation":
            for i in range(self.n_perm):
                out[i] = rng.permutation(v)
        elif self.kind == "spin":
            tree = cKDTree(self._xyz)
            for i in range(self.n_perm):
                rot = _random_rotation(rng)
                _, idx = tree.query(self._xyz @ rot.T)
                out[i] = v[idx]
        else:  # msr: flip eigenvector signs inside near-degenerate blocks
            proj = self._v.T @ v
            lam = self._lam
            tol = 1e-6 * (np.abs(lam).max() + 1e-12)
            blocks, cur = [], [0]
            for i in range(1, lam.size):
                if abs(abs(lam[i]) - abs(lam[i - 1])) <= tol:
                    cur.append(i)
                else:
                    blocks.append(cur)
                    cur = [i]
            blocks.append(cur)
            for i in range(self.n_perm):
                signs = np.ones_like(proj)
                for b in blocks:
                    if len(b) > 1:
                        s = rng.choice([-1.0, 1.0], size=len(b))
                        signs[b] = s
                out[i] = self._v @ (proj * signs)
        return out


def _two_sided_p(obs, null):
    null = np.asarray(null, dtype=float)
    ge = np.count_nonzero(null >= obs)
    le = np.count_nonzero(null <= obs)
    return (1.0 + min(ge, le)) / (1.0 + null.size)


# ------------------------------------------------------------ estimators ---

def _nanstd(x, axis=0):
    """NaN-aware SD that returns NaN (silently) for all-NaN slices."""
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    n = ok.sum(axis=axis)
    s = np.where(ok, x, 0.0)
    mean = np.divide(s.sum(axis=axis), n, out=np.zeros_like(n, dtype=float), where=n > 0)
    dev = np.where(ok, x - np.expand_dims(mean, axis), 0.0)
    var = np.divide((dev ** 2).sum(axis=axis), n, out=np.full(n.shape, np.nan), where=n > 1)
    return np.sqrt(var)


def variability(coupling, cohort, estimator="regional"):
    """Dispersion map for one estimator.

    regional : SD across subjects of the regional read-out (handled by caller)
    node     : per region, the mean of its edges' SD
    edge     : SD across subjects of each edge
    """
    if estimator == "edge":
        return _nanstd(coupling, axis=0)
    edge_sd = _nanstd(coupling, axis=0)
    if estimator == "node":
        out = np.full(cohort.n_nodes, np.nan)
        for i, idx in enumerate(cohort._region_profile):
            if idx.size:
                out[i] = np.nanmean(edge_sd[idx])
        return out
    raise ValueError("use regional_variability() for the regional estimator")


def regional_variability(regional, covariates=None):
    """SD across subjects of regional coupling, optionally adjusted."""
    x = np.asarray(regional, dtype=float)
    if covariates is not None and covariates.size:
        x = _residualise(x, covariates)
    return _nanstd(x, axis=0)


def _residualise(y, cov):
    cov = np.asarray(cov, dtype=float)
    if cov.ndim == 1:
        cov = cov[:, None]
    x = np.column_stack([np.ones(cov.shape[0]), cov])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


# ------------------------------------------------------------------ axis ---

def axis(v, hierarchy, null_engine=None, n_perm=1000):
    """Spearman r between a dispersion map and the cortical hierarchy.

    Regional maps are compared with the regional hierarchy.  Edge maps are
    compared with the endpoint-mean hierarchy, and the null is built by
    transforming the *node* hierarchy and re-expanding it, so that the
    spatial structure of the null is inherited from the parcel layout.
    """
    v = np.asarray(v, dtype=float).ravel()
    h = np.asarray(hierarchy, dtype=float).ravel()
    if h.size != v.size:
        raise ValueError("hierarchy length %d does not match map length %d" % (h.size, v.size))
    ok = np.isfinite(v) & np.isfinite(h)
    if ok.sum() < 2:
        # Too few valid regions to correlate (e.g. a dispersion map that is
        # undefined for every region).  Return NaN instead of letting scipy
        # raise "`x` and `y` must have length at least 2".
        return dict(r=float("nan"), p=float("nan"), null="none", n_perm=0,
                    n=int(ok.sum()), note="fewer than 2 valid entries")
    r = stats.spearmanr(v[ok], h[ok]).statistic
    res = dict(r=float(r), p=float("nan"), null="none", n_perm=0, n=int(ok.sum()))

    if null_engine is None:
        return res

    coh = getattr(null_engine, "cohort", None)
    n_nodes = coh.n_nodes if coh is not None else h.size
    is_edge = (v.size != n_nodes)

    v_filled = np.where(ok, v, np.nanmean(v[ok]))
    if not is_edge:
        surr = null_engine.surrogates(v_filled)
        null_r = np.array([stats.spearmanr(s[ok], h[ok]).statistic for s in surr], dtype=float)
    else:
        base = coh.hierarchy
        surr_h = null_engine.surrogates(base)
        null_r = []
        for sh in surr_h:
            eh = 0.5 * (sh[coh.iu] + sh[coh.iv])
            null_r.append(stats.spearmanr(v[ok], eh[ok]).statistic)
        null_r = np.array(null_r, dtype=float)

    null_r = null_r[np.isfinite(null_r)]
    if null_r.size:
        res["p"] = float(_two_sided_p(r, null_r))
        res["null_mean"] = float(np.mean(null_r))
        res["null_sd"] = float(np.std(null_r))
    res["null"] = null_engine.kind
    res["n_perm"] = int(null_engine.n_perm)
    res["n_null_valid"] = int(null_r.size)
    return res


def axis_controlled(v, hierarchy, control, null_engine=None, n_perm=1000):
    """Partial Spearman of v and hierarchy given a control map (rank residual)."""
    v = np.asarray(v, dtype=float).ravel()
    h = np.asarray(hierarchy, dtype=float).ravel()
    c = np.asarray(control, dtype=float).ravel()
    ok = np.isfinite(v) & np.isfinite(h) & np.isfinite(c)
    if ok.sum() < 2:
        return dict(r=float("nan"), p=float("nan"), null="none", n_perm=0,
                    n=int(ok.sum()), note="fewer than 2 valid entries")
    rv, rh = stats.rankdata(v[ok]), stats.rankdata(h[ok])
    rc = stats.rankdata(c[ok])[:, None]
    e_v = _residualise(rv, rc).ravel()
    e_h = _residualise(rh, rc).ravel()
    r = stats.pearsonr(e_v, e_h).statistic
    res = dict(r=float(r), p=float("nan"), null="none", n_perm=0, n=int(ok.sum()))
    if null_engine is not None:
        surr = null_engine.surrogates(np.where(ok, v, np.nanmean(v[ok])))
        null_r = []
        for s in surr:
            ev = _residualise(stats.rankdata(s[ok]), rc).ravel()
            null_r.append(stats.pearsonr(ev, e_h).statistic)
        res["p"] = float(_two_sided_p(r, np.array(null_r)))
        res["null"] = null_engine.kind
        res["n_perm"] = int(null_engine.n_perm)
    return res


# ---------------------------------------------------------------- layers ---

def layers(edge_v, cohort, n_perm=1000, seed=0):
    """Within- versus between-network dispersion of an edge-level map."""
    v = np.asarray(edge_v, dtype=float)
    same = cohort.networks[cohort.iu] == cohort.networks[cohort.iv]
    w = np.nanmean(v[same])
    b = np.nanmean(v[~same])
    delta = w - b

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    base = cohort.networks.copy()
    for i in range(n_perm):
        lab = rng.permutation(base)
        s = lab[cohort.iu] == lab[cohort.iv]
        null[i] = np.nanmean(v[s]) - np.nanmean(v[~s])
    return dict(within=float(w), between=float(b), delta=float(delta),
                p=float(_two_sided_p(delta, null)), n_edges=int(v.size),
                same_net_edges=int(same.sum()), null="network-label permutation",
                n_perm=int(n_perm))


# ------------------------------------------------------------- epicentres --

def epicentres(v, cohort, k=2, seed=0, n_perm=1000):
    v = np.asarray(v, dtype=float)
    a = cohort.adjacency() > 0
    n = cohort.n_nodes
    med = np.nanmedian(v)
    focal = []
    for i in range(n):
        nb = np.flatnonzero(a[i])
        if nb.size == 0:
            continue
        if v[i] >= np.nanmax(v[nb]) and np.nanmean(v[nb]) > med:
            focal.append(i)
    focal = np.array(sorted(focal), dtype=int)

    sep = 0
    for x in range(focal.size):
        for y in range(x + 1, focal.size):
            if a[focal[x], focal[y]]:
                sep += 1
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.choice(n, size=max(focal.size, 1), replace=False)
        c = 0
        for x in range(idx.size):
            for y in range(x + 1, idx.size):
                if a[idx[x], idx[y]]:
                    c += 1
        null[i] = c
    z = (sep - null.mean()) / (null.std() + 1e-12)

    prof = np.nan_to_num(v, nan=np.nanmedian(v))
    prof = (prof - prof.mean()) / (prof.std() + 1e-12)
    link = sch.linkage(prof[:, None], method="ward")
    lab = sch.fcluster(link, t=k, criterion="maxclust")

    counts = {}
    for i in focal:
        counts[int(cohort.networks[i])] = counts.get(int(cohort.networks[i]), 0) + 1
    return dict(regions=[int(i) for i in focal], n=int(focal.size),
                internal_edges=int(sep), z=float(z),
                p=float((1 + np.count_nonzero(null <= sep)) / (1 + n_perm)),
                null_mean=float(null.mean()),
                cluster_labels=lab.astype(int).tolist(),
                network_counts=counts)


# ------------------------------------------------------------- individual --

def deviation_spectrum(coupling, cohort, scales=("regional", "network", "edge"),
                       max_edges=20000, seed=0):
    """Per-subject deviation-spectrum slope at several observation scales."""
    rng = np.random.default_rng(seed)
    out = {}
    if "regional" in scales:
        reg = cohort.regional_coupling()
        out["regional"] = _slope(reg)
    if "network" in scales:
        reg = cohort.regional_coupling()
        net = np.full((reg.shape[0], cohort.n_networks), np.nan)
        for g in range(cohort.n_networks):
            m = cohort.networks == g
            if m.any():
                net[:, g] = np.nanmean(reg[:, m], axis=1)
        out["network"] = _slope(net)
    if "edge" in scales:
        e = coupling
        if e.shape[1] > max_edges:
            idx = rng.choice(e.shape[1], size=max_edges, replace=False)
            e = e[:, idx]
        out["edge"] = _slope(e)
    if "gap" in scales:
        same = cohort.networks[cohort.iu] == cohort.networks[cohort.iv]
        w = np.nanmean(coupling[:, same], axis=1)
        b = np.nanmean(coupling[:, ~same], axis=1)
        out["gap"] = w - b
    return out


def _slope(x):
    """Slope of the deviation profile against rank (sorted by signed value)."""
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    s = np.where(ok, x, 0.0).sum(axis=0)
    c = ok.sum(axis=0)
    mean = np.divide(s, c, out=np.zeros_like(s), where=c > 0)
    dev = np.where(ok, x - mean, 0.0)
    order = np.argsort(-dev, axis=1)
    ranked = np.take_along_axis(dev, order, axis=1)
    m = ranked.shape[1]
    r = np.arange(1, m + 1, dtype=float)
    r = (r - r.mean()) / r.std()
    return (ranked * r[None, :]).mean(axis=1) / (r ** 2).mean()


def partial_spearman(x, y, cov=None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if cov is not None:
        cov = np.asarray(cov, dtype=float)
        if cov.ndim == 1:
            cov = cov[:, None]
        ok &= np.all(np.isfinite(cov), axis=1)
    if ok.sum() < 2:
        return dict(r=float("nan"), p=float("nan"), n=int(ok.sum()),
                    note="fewer than 2 valid observations")
    rx, ry = stats.rankdata(x[ok]), stats.rankdata(y[ok])
    if cov is not None:
        c = stats.rankdata(cov[ok], axis=0)
        rx = _residualise(rx, c).ravel()
        ry = _residualise(ry, c).ravel()
    res = stats.pearsonr(rx, ry)
    return dict(r=float(res.statistic), p=float(res.pvalue), n=int(ok.sum()))


def reliability(coupling, cohort, seed=0):
    a, b = cohort.split_half(seed)
    va = _nanstd(coupling[a], axis=0)
    vb = _nanstd(coupling[b], axis=0)
    ok = np.isfinite(va) & np.isfinite(vb)
    if ok.sum() < 2 or a.size < 2 or b.size < 2:
        # Split-half reliability needs at least 2 subjects per half; with
        # smaller cohorts the dispersion in each half is undefined.
        return dict(split_half_r=float("nan"), n_edges=0, n_half=int(a.size),
                    note="cohort too small for split-half reliability")
    r = stats.spearmanr(va[ok], vb[ok]).statistic
    return dict(split_half_r=float(r), n_edges=int(ok.sum()), n_half=int(a.size))
