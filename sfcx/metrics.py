# -*- coding: utf-8 -*-
"""Coupling definitions, layer 1 of the coordinate system.

Each definition returns an (n_subjects, n_edges) array for edge-native measures,
or (n_subjects, n_nodes) for regional measures.  Matrix definitions first build
a fixed structural measure from the group-mean connectome and then couple it to
each subject's functional edges.

The full set of 23 definitions from the manuscript is implemented:
- 9 static communication measures (group-mean matrix -> edge coupling)
- 4 reliance measures (per-subject BOLD -> regional coupling)
- 5 hybrid models (per-subject FC + group-mean SC -> regional coupling)
- 4 frontier measures (per-subject BOLD + SC -> regional coupling)
- 2 innovative measures (per-subject BOLD + SC -> regional coupling)
"""

from __future__ import annotations

import numpy as np
from scipy import linalg, sparse
from scipy.sparse import csgraph
from scipy.stats import pearsonr, rankdata

from .data import mat_to_vec, vec_to_mat, _zscore_by_group


def _standardise_edges(x, groups, n_groups):
    return _zscore_by_group(np.asarray(x, dtype=np.float64), groups, n_groups)


def _couple(measure_vec, fc, iu, n_nodes):
    """z(measure) * z(FC) per subject."""
    mg = _zscore_by_group(measure_vec[None, :], iu, n_nodes)[0]
    fg = _zscore_by_group(fc.astype(np.float64), iu, n_nodes)
    return fg * mg[None, :]


def _incidence(cohort):
    """Sparse region-edge incidence, cached on the cohort."""
    cached = getattr(cohort, "_incidence_cache", None)
    if cached is not None:
        return cached
    from scipy import sparse

    e = cohort.n_edges
    cols = np.arange(e)
    m = sparse.csr_matrix(
        (np.ones(2 * e), (np.concatenate([cohort.iu, cohort.iv]),
                          np.concatenate([cols, cols])),
         ),
        shape=(cohort.n_nodes, e))
    mi = sparse.csr_matrix((np.ones(e), (cohort.iu, cols)), shape=(cohort.n_nodes, e))
    mj = sparse.csr_matrix((np.ones(e), (cohort.iv, cols)), shape=(cohort.n_nodes, e))
    counts = np.asarray(m.sum(axis=1)).ravel()
    counts[counts < 1] = 1.0
    cohort._incidence_cache = (m, mi, mj, counts)
    return cohort._incidence_cache


def _profile_z(x, cohort):
    """Standardise each edge inside the profile of each of its endpoints."""
    m, _, _, counts = _incidence(cohort)
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x[None, :]
    s1 = x @ m.T
    s2 = (x * x) @ m.T
    mean = s1 / counts[None, :]
    var = np.maximum(s2 / counts[None, :] - mean * mean, 0.0)
    sd = np.sqrt(var)
    sd[sd < 1e-12] = 1.0
    zi = (x - mean[:, cohort.iu]) / sd[:, cohort.iu]
    zj = (x - mean[:, cohort.iv]) / sd[:, cohort.iv]
    return zi, zj


def _couple_profile(measure_vec, fc, cohort):
    """Symmetric coupling field from per-endpoint standardisation."""
    mi, mj = _profile_z(measure_vec, cohort)
    fi, fj = _profile_z(fc, cohort)
    return 0.5 * (mi * fi + mj * fj)


def _pair_scale(m):
    """Scale a non-negative pair matrix to a comparable range."""
    m = np.asarray(m, dtype=np.float64)
    m = np.nan_to_num(m, nan=0.0, posinf=0.0, neginf=0.0)
    m = np.log1p(np.clip(m, 0, None))
    sd = m.std()
    return m / sd if sd > 1e-12 else m


def _col_normalize(W):
    """Column-normalise a weight matrix (used by reliance / hybrid / frontier)."""
    W = np.asarray(W, dtype=np.float64)
    W = W / (W.sum(axis=0, keepdims=True) + 1e-10)
    np.fill_diagonal(W, 0.0)
    W = W / (W.sum(axis=0, keepdims=True) + 1e-10)
    return W


def _robust_pearsonr(x, y):
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 3:
        return 0.0, 1.0
    return pearsonr(x[mask], y[mask])


def _fisher_z(r):
    return np.arctanh(np.clip(r, -0.9999, 0.9999))


# ------------------------------------------------------------- matrix ops --

def _structural(cohort, binarise=True):
    a = cohort.group_matrix("sc")
    if binarise:
        a = (a > 1e-9).astype(np.float64)
    return a


def _cosine(a):
    nrm = np.sqrt((a * a).sum(axis=1))
    nrm[nrm < 1e-12] = 1.0
    an = a / nrm[:, None]
    return an @ an.T


def _matching(a):
    a2 = a @ a
    k = a.sum(axis=1)
    den = k[:, None] + k[None, :] - 2.0 * a
    with np.errstate(divide="ignore", invalid="ignore"):
        m = np.where(den > 0, a2 / den, 0.0)
    return np.nan_to_num(m)


def _communicability(a):
    w, v = linalg.eigh(a)
    w = np.clip(w, 0, None)
    return (v * np.exp(w)) @ v.T


def _degree_product(a):
    k = a.sum(axis=1)
    return k[:, None] * k[None, :]


def _diffusion_efficiency(a):
    w = np.where(a > 0, 1.0 / np.maximum(a, 1e-12), 0.0)
    sp = sparse.csr_matrix(w)
    d = csgraph.dijkstra(sp, directed=False)
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = 1.0 / d
    eff[~np.isfinite(eff)] = 0.0
    return eff


def _flow(a):
    """Two-step random-walk reachability, (D^-1 A)^2, symmetrised."""
    k = a.sum(axis=1)
    k[k < 1e-12] = 1.0
    p = a / k[:, None]
    p2 = p @ p
    return 0.5 * (p2 + p2.T)


def _mfpt(a):
    n = a.shape[0]
    k = a.sum(axis=1)
    k[k < 1e-12] = 1.0
    p = a / k[:, None]
    pi = k / k.sum()
    z = linalg.inv(np.eye(n) - p + np.ones((n, 1)) @ pi[None, :])
    d = np.diag(z)
    m = (d[None, :] - z) / np.clip(pi[None, :], 1e-12, None)
    return np.clip(m, 0, None)


def _path_ease(a):
    w = np.where(a > 0, 1.0 / np.maximum(a, 1e-12), 0.0)
    sp = sparse.csr_matrix(w)
    d = csgraph.dijkstra(sp, directed=False)
    d[~np.isfinite(d)] = np.nan
    with np.errstate(invalid="ignore"):
        ease = 1.0 / (1.0 + d)
    return np.nan_to_num(ease)


def _search_info(a):
    w = np.where(a > 0, 1.0 / np.maximum(a, 1e-12), 0.0)
    sp = sparse.csr_matrix(w)
    d, pred = csgraph.dijkstra(sp, directed=False, return_predecessors=True)
    n = a.shape[0]
    k = a.sum(axis=1)
    k[k < 1e-12] = 1.0
    logp = -np.log2(1.0 / k)
    out = np.zeros((n, n))
    for i in range(n):
        p = pred[i]
        for j in range(n):
            if i == j or not np.isfinite(d[i, j]):
                continue
            s = 0.0
            node = j
            guard = 0
            while node != i and node >= 0 and guard < n:
                s += logp[node]
                node = p[node]
                guard += 1
            out[i, j] = s
    return out


# --------------------------------------------------- per-subject regional --

def _reliance(cohort):
    """Per-subject, per-region slope of FC on SC, re-expressed at edge level."""
    n, e = cohort.n, cohort.n_edges
    out = np.zeros((n, e), dtype=np.float64)
    for i, idx in enumerate(cohort._region_edges):
        if idx.size < 3:
            continue
        sc = cohort.sc[:, idx].astype(np.float64)
        fc = cohort.fc[:, idx].astype(np.float64)
        scm = sc - sc.mean(axis=1, keepdims=True)
        fcm = fc - fc.mean(axis=1, keepdims=True)
        var = (scm * scm).sum(axis=1, keepdims=True)
        var[var < 1e-12] = np.inf
        b = (scm * fcm).sum(axis=1, keepdims=True) / var
        pred = b * scm
        out[:, idx] = scm * pred
    return _zscore_by_group(out, cohort.iu, cohort.n_nodes)


def _reliance_raw(cohort):
    """Classical reliance: corr(ts_{t-1} W, ts_t) per region."""
    if cohort.bold is None:
        raise ValueError("reliance_raw requires BOLD time series in the cohort")
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        y_pred = ts[:-1] @ W
        y_actual = ts[1:]
        for i in range(p):
            out[s, i], _ = _robust_pearsonr(y_pred[:, i], y_actual[:, i])
    return out


def _reliance_deriv(cohort):
    """Derivative reliance: corr(-ts + ts @ W, dts/dt) per region."""
    if cohort.bold is None:
        raise ValueError("reliance_deriv requires BOLD time series in the cohort")
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        dy = ts[1:] - ts[:-1]
        y_curr = ts[:-1]
        dy_pred = -y_curr + y_curr @ W
        for i in range(p):
            out[s, i], _ = _robust_pearsonr(dy_pred[:, i], dy[:, i])
    return out


def _reliance_raw_nonlinear(cohort, alpha=1.0):
    """Nonlinear reliance (leak + tanh): corr((1-a)ts + a*tanh(ts@W), ts_t)."""
    if cohort.bold is None:
        raise ValueError("reliance_raw_nonlinear requires BOLD time series")
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        drive = ts[:-1] @ W
        y_pred = (1 - alpha) * ts[:-1] + alpha * np.tanh(drive)
        y_actual = ts[1:]
        for i in range(p):
            out[s, i], _ = _robust_pearsonr(y_pred[:, i], y_actual[:, i])
    return out


def _reliance_deriv_nonlinear(cohort, alpha=1.0):
    """Nonlinear derivative reliance: corr(a*(tanh(ts@W)-ts), dts/dt)."""
    if cohort.bold is None:
        raise ValueError("reliance_deriv_nonlinear requires BOLD time series")
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        dy = ts[1:] - ts[:-1]
        y_curr = ts[:-1]
        drive = y_curr @ W
        dy_pred = alpha * (np.tanh(drive) - y_curr)
        for i in range(p):
            out[s, i], _ = _robust_pearsonr(dy_pred[:, i], dy[:, i])
    return out


# --------------------------------------------------------------- hybrid ----

def _zscore_cols(M):
    M = np.asarray(M, dtype=np.float64)
    out = np.zeros_like(M)
    for i in range(M.shape[1]):
        col = M[:, i]
        s = col.std()
        if s > 1e-10:
            out[:, i] = (col - col.mean()) / s
    return out


def _hybrid_matrices(cohort):
    """Z-scored group-mean SC + communication matrices for hybrid models."""
    a = cohort.group_matrix("sc")
    a_log = np.log1p(np.clip(a, 0, None))
    mats = [a_log]
    for fn in (_cosine, _matching, _communicability, _mfpt, _path_ease,
               _search_info, _flow):
        m = fn((a_log > 1e-9).astype(np.float64))
        mats.append(m)
    return [_zscore_cols(m) for m in mats]


def _loo_col_corr(C, F):
    """Per-column Pearson r between matrices C and F, leaving out each column's
    own row.  Fully vectorised: O(p^2), no Python loop over regions.
    """
    C = np.asarray(C, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    n = C.shape[0]
    idx = np.arange(n)
    sC, qC = C.sum(0), (C * C).sum(0)
    sF, qF = F.sum(0), (F * F).sum(0)
    sCF = (C * F).sum(0)
    dC, dF = C[idx, idx], F[idx, idx]
    m = max(n - 1, 1)
    muC = (sC - dC) / m
    muF = (sF - dF) / m
    cov = (sCF - dC * dF) / m - muC * muF
    varC = (qC - dC * dC) / m - muC * muC
    varF = (qF - dF * dF) / m - muF * muF
    den = np.sqrt(np.maximum(varC, 0.0) * np.maximum(varF, 0.0))
    return np.where(den > 1e-12, cov / den, 0.0)


def _safe_weights(w, k):
    """Non-negative weights summing to one; never divides by zero."""
    w = np.clip(np.asarray(w, dtype=float), 0.0, 1.0)
    s = float(np.sum(w))
    if not np.isfinite(s) or s < 1e-9:
        return np.full(k, 1.0 / max(k, 1))
    return w / s


def _combine(M_all_z, w):
    out = np.zeros_like(M_all_z[0])
    for wk, mk in zip(w, M_all_z):
        out = out + wk * mk
    return out


def _slsqp_weights(obj, k):
    """Solve min obj(w) s.t. w >= 0, sum(w) = 1, crash-safe."""
    from scipy.optimize import minimize
    res = minimize(obj, np.full(k, 1.0 / max(k, 1)), method="SLSQP",
                   bounds=[(0.0, 1.0)] * k,
                   constraints=[{"type": "eq", "fun": lambda w: float(np.sum(w) - 1.0)}],
                   options={"maxiter": 200, "ftol": 1e-8})
    return _safe_weights(res.x, k)


def _hybrid_global(cohort):
    """Global hybrid: one set of weights across all regions."""
    M_all_z = _hybrid_matrices(cohort)
    k = len(M_all_z)
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, p)
        w = _slsqp_weights(
            lambda w: -float(np.mean(_loo_col_corr(_combine(M_all_z, w), fc))), k)
        out[s] = _loo_col_corr(_combine(M_all_z, w), fc)
    return out


def _hybrid_per_region(cohort):
    """Per-region hybrid: one weight vector per region."""
    M_all_z = _hybrid_matrices(cohort)
    k = len(M_all_z)
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)
    stack = np.stack(M_all_z, axis=0)               # (k, p, p)
    for s in range(n):
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, p)
        for i in range(p):
            cols = stack[:, :, i]                    # (k, p)
            mask = np.ones(p, dtype=bool)
            mask[i] = False
            y = fc[mask, i]

            def obj(w, cols=cols, mask=mask, y=y):
                x = np.tensordot(_safe_weights(w, k), cols, axes=1)[mask]
                r, _ = _robust_pearsonr(x, y)
                return -r

            w = _slsqp_weights(obj, k)
            x = np.tensordot(w, cols, axes=1)[mask]
            out[s, i], _ = _robust_pearsonr(x, y)
    return out


def _hybrid_cv(cohort, n_splits=5):
    """Cross-validated global hybrid.

    Weights are fitted on the training regions and evaluated on the held-out
    regions; the per-column correlation always leaves out the column's own row
    (vectorised leave-one-out, see ``_loo_col_corr``).
    """
    from sklearn.model_selection import KFold
    M_all_z = _hybrid_matrices(cohort)
    k = len(M_all_z)
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    regions = np.arange(p)
    for s in range(n):
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, p)
        for train_idx, test_idx in kf.split(regions):
            w = _slsqp_weights(
                lambda w: -float(np.mean(
                    _loo_col_corr(_combine(M_all_z, w), fc)[train_idx])), k)
            r = _loo_col_corr(_combine(M_all_z, w), fc)
            out[s, test_idx] = r[test_idx]
    return out


def _hybrid_ridge(cohort):
    """Ridge regression hybrid."""
    from sklearn.linear_model import RidgeCV
    M_all_z = _hybrid_matrices(cohort)
    alphas = np.logspace(-3, 3, 20)
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)

    for s in range(n):
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, cohort.n_nodes)
        for i in range(p):
            mask = np.ones(p, dtype=bool)
            mask[i] = False
            X = np.column_stack([M[mask, i] for M in M_all_z])
            y = fc[mask, i]
            model = RidgeCV(alphas=alphas, cv=5).fit(X, y)
            out[s, i], _ = _robust_pearsonr(model.predict(X), y)
    return out


def _hybrid_network(cohort, n_clusters=7, n_splits=3):
    """Network-constrained hybrid."""
    from sklearn.cluster import KMeans
    from sklearn.model_selection import KFold
    from scipy.optimize import minimize
    M_all_z = _hybrid_matrices(cohort)
    n_methods = len(M_all_z)
    bounds = [(0, 1)] * n_methods
    constraints = [{'type': 'eq', 'fun': lambda w: np.sum(w) - 1}]
    x0 = np.ones(n_methods) / n_methods
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)
    region_indices = np.arange(p)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    for s in range(n):
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, cohort.n_nodes)
        for train_idx, test_idx in kf.split(region_indices):
            fc_train = fc[train_idx][:, train_idx]
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            train_labels = kmeans.fit_predict(fc_train)
            unique_labels = sorted(set(train_labels))
            label_map = {orig: new for new, orig in enumerate(unique_labels)}
            train_labels = np.array([label_map[l] for l in train_labels])
            test_labels = kmeans.predict(fc[test_idx][:, train_idx])
            test_labels = np.array([label_map.get(l, -1) for l in test_labels])

            for net in range(len(unique_labels)):
                net_train = train_idx[train_labels == net]
                net_test = test_idx[test_labels == net]
                if len(net_train) < 2 or len(net_test) == 0:
                    continue
                def net_neg(w):
                    w = np.maximum(np.asarray(w, dtype=float), 0)
                    w = w / w.sum()
                    combined = sum(w[k] * M_all_z[k] for k in range(n_methods))
                    r_sum = 0.0
                    for i in net_train:
                        mask = np.ones(p, dtype=bool)
                        mask[i] = False
                        mask[~np.isin(np.arange(p), net_train)] = False
                        r, _ = _robust_pearsonr(combined[:, i][mask], fc[:, i][mask])
                        r_sum += r
                    return -(r_sum / max(len(net_train), 1))
                res = minimize(net_neg, x0, method='SLSQP', bounds=bounds,
                               constraints=constraints,
                               options={'maxiter': 500, 'ftol': 1e-10, 'disp': False})
                w = res.x / res.x.sum()
                combined = sum(w[k] * M_all_z[k] for k in range(n_methods))
                for i in net_test:
                    mask = np.ones(p, dtype=bool)
                    mask[i] = False
                    out[s, i], _ = _robust_pearsonr(combined[:, i][mask], fc[:, i][mask])
    return out


# ------------------------------------------------------------- frontier ----

def _information_geometry(cohort):
    """Log-Euclidean geodesic coupling on the SPD manifold."""
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        sc = vec_to_mat(cohort.sc[s], cohort.iu, cohort.iv, cohort.n_nodes)
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, cohort.n_nodes)
        min_eig_sc = np.min(np.linalg.eigvalsh(sc))
        if min_eig_sc <= 0:
            sc = sc + (abs(min_eig_sc) + 1e-6) * np.eye(p)
        min_eig_fc = np.min(np.linalg.eigvalsh(fc))
        if min_eig_fc <= 0:
            fc = fc + (abs(min_eig_fc) + 1e-6) * np.eye(p)
        log_sc = linalg.logm(sc)
        log_fc = linalg.logm(fc)
        for i in range(p):
            mask = np.ones(p, dtype=bool)
            mask[i] = False
            sc_col = np.real(log_sc[:, i][mask])
            fc_col = np.real(log_fc[:, i][mask])
            out[s, i], _ = _robust_pearsonr(sc_col, fc_col)
    return out


def _multi_scale_transfer_entropy(cohort, scales=(1, 2, 4, 8), n_cond_bins=5):
    """Multi-scale transfer entropy from SC to BOLD."""
    if cohort.bold is None:
        raise ValueError("multi_scale_transfer_entropy requires BOLD time series")
    from sklearn.metrics import mutual_info_score as _mi
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    out = np.zeros((n, p), dtype=np.float64)

    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        te_scores = []
        for scale in scales:
            ts_down = ts[::scale]
            T = ts_down.shape[0]
            if T < 20:
                te_scores.append(np.zeros(p))
                continue
            sc_pred = ts_down[:-1] @ W
            te = np.zeros(p)
            for i in range(p):
                sc_vals = sc_pred[:, i]
                fc_future = ts_down[1:, i]
                fc_past = ts_down[:-1, i]
                n_groups = min(n_cond_bins, len(fc_past) // 20)
                if n_groups < 2:
                    te[i] = 0
                    continue
                quantiles = np.percentile(fc_past, np.linspace(0, 100, n_groups + 1))
                weighted_mi = 0.0
                total_weight = 0
                for g in range(n_groups):
                    if g == n_groups - 1:
                        mask = (fc_past >= quantiles[g]) & (fc_past <= quantiles[g + 1])
                    else:
                        mask = (fc_past >= quantiles[g]) & (fc_past < quantiles[g + 1])
                    if mask.sum() < 10:
                        continue
                    sc_g = sc_vals[mask]
                    fc_f_g = fc_future[mask]
                    sc_bin = np.digitize(sc_g, bins=np.linspace(sc_g.min(), sc_g.max(), 15))
                    fc_f_bin = np.digitize(fc_f_g, bins=np.linspace(fc_f_g.min(), fc_f_g.max(), 15))
                    mi_g = _mi(sc_bin, fc_f_bin)
                    weighted_mi += mi_g * mask.sum()
                    total_weight += mask.sum()
                te[i] = weighted_mi / max(total_weight, 1)
            te_scores.append(te)
        out[s] = np.mean(te_scores, axis=0)
    return out


def _criticality(cohort):
    """Criticality index: distance of SC eigenvalue spectrum from lambda=1."""
    n, p = cohort.n, cohort.n_nodes
    out = np.zeros((n, p), dtype=np.float64)
    for s in range(n):
        sc = vec_to_mat(cohort.sc[s], cohort.iu, cohort.iv, cohort.n_nodes)
        fc = vec_to_mat(cohort.fc[s], cohort.iu, cohort.iv, cohort.n_nodes)
        sc_coupling = sc / (sc.sum(axis=1, keepdims=True) + 1e-10)
        eigvals = np.linalg.eigvalsh(sc_coupling)
        lambda_max = np.max(np.abs(eigvals))
        ci = 1.0 / (1.0 + abs(lambda_max - 1.0))
        for i in range(p):
            mask = np.ones(p, dtype=bool)
            mask[i] = False
            sc_row = sc_coupling[i, :][mask]
            fc_col = fc[:, i][mask]
            r, _ = _robust_pearsonr(sc_row, fc_col)
            out[s, i] = r * ci
    return out


def _granger(cohort, max_lag=3):
    """Granger causality F-statistic: SC-weighted past predicts BOLD future."""
    if cohort.bold is None:
        raise ValueError("granger requires BOLD time series")
    from sklearn.linear_model import LinearRegression
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    out = np.zeros((n, p), dtype=np.float64)

    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        T = ts.shape[0]
        if T <= max_lag + 2:
            continue
        sc_weighted_full = ts[:-1] @ W
        for i in range(p):
            y = ts[max_lag:, i]
            X_restricted = np.column_stack([
                ts[max_lag - lag:T - lag, i] for lag in range(1, max_lag + 1)
            ])
            sc_cols = []
            for lag in range(1, max_lag + 1):
                start = max_lag - lag
                end = T - lag - 1
                if start < end:
                    sc_cols.append(sc_weighted_full[start:end, i:i+1])
            if not sc_cols:
                continue
            min_len = min(c.shape[0] for c in sc_cols)
            sc_cols = [c[-min_len:] for c in sc_cols]
            X_restricted_adj = X_restricted[-min_len:]
            y_adj = y[-min_len:]
            X_full = np.column_stack([X_restricted_adj] + sc_cols)
            try:
                reg_r = LinearRegression().fit(X_restricted_adj, y_adj)
                reg_f = LinearRegression().fit(X_full, y_adj)
                rss_r = np.sum((y_adj - reg_r.predict(X_restricted_adj)) ** 2)
                rss_f = np.sum((y_adj - reg_f.predict(X_full)) ** 2)
                n_obs = len(y_adj)
                q = X_full.shape[1]
                if rss_f > 1e-10:
                    F = ((rss_r - rss_f) / max_lag) / (rss_f / (n_obs - q - 1))
                    out[s, i] = max(F, 0)
            except Exception:
                out[s, i] = 0
    return out


# ----------------------------------------------------------- innovative ----

def _xorshift32(seed):
    """Deterministic xorshift32 for cross-language reproducibility."""
    x = np.uint32(seed)
    while True:
        x ^= x << 13
        x ^= x >> 17
        x ^= x << 5
        yield int(x)


def _weight_shuffle_null(a_log, seed):
    """Weight-multiset permutation null model, deterministic across languages."""
    n = a_log.shape[0]
    iu, iv = np.triu_indices(n, 1)
    w = a_log[iu]
    wv = w[w > 0]
    if wv.size == 0:
        return np.zeros_like(a_log)
    gen = _xorshift32(seed)
    perm = np.argsort([next(gen) for _ in range(wv.size)])
    wv_shuffled = wv[perm]
    tmp = np.zeros_like(w)
    tmp[w > 0] = wv_shuffled
    a_null = np.zeros_like(a_log)
    a_null[iu] = tmp
    return a_null + a_null.T


def _rnc(cohort):
    """Reliance-Null Coupling: reliance minus weight-shuffle null."""
    if cohort.bold is None:
        raise ValueError("rnc requires BOLD time series")
    n, p, t = cohort.bold.shape
    a = cohort.group_matrix("sc")
    a_log = np.log1p(np.clip(a, 0, None))
    W = _col_normalize(a_log)
    out = np.zeros((n, p), dtype=np.float64)

    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        y_pred = ts[:-1] @ W
        y_actual = ts[1:]
        r_obs = np.zeros(p)
        for i in range(p):
            r_obs[i], _ = _robust_pearsonr(y_pred[:, i], y_actual[:, i])
        seed = hash((cohort.state, s)) % (2 ** 31)
        a_null = _weight_shuffle_null(a_log, seed)
        W_null = _col_normalize(a_null)
        y_pred_null = ts[:-1] @ W_null
        r_null = np.zeros(p)
        for i in range(p):
            r_null[i], _ = _robust_pearsonr(y_pred_null[:, i], y_actual[:, i])
        out[s] = r_obs - r_null
    return out


def _trr(cohort, lags=(1, 2, 3, 4, 5)):
    """Temporal-Reach Reliance: weighted sum of multi-lag reliance."""
    if cohort.bold is None:
        raise ValueError("trr requires BOLD time series")
    n, p, t = cohort.bold.shape
    W = _col_normalize(cohort.group_matrix("sc"))
    w = 1.0 / np.array(lags, dtype=float)
    w = w / w.sum()
    out = np.zeros((n, p), dtype=np.float64)

    for s in range(n):
        ts = cohort.bold[s].T.astype(np.float64)  # (T, n_nodes)
        P = ts @ W
        R = np.full((len(lags), p), np.nan)
        for row, k in enumerate(lags):
            if ts.shape[0] - k < 10:
                continue
            for i in range(p):
                R[row, i], _ = _robust_pearsonr(P[:-k, i], ts[k:, i])
        out[s] = np.nansum(w[:, None] * R, axis=0)
    return out


# --------------------------------------------------------------- public ----

_EDGE_DEFS = {"direct", "rank_direct", "reliance"}
_MATRIX_DEFS = {
    "cosine": _cosine,
    "matching": _matching,
    "communicability": _communicability,
    "degree_product": _degree_product,
    "mfpt": _mfpt,
    "path_ease": _path_ease,
    "search_info": _search_info,
    "diffusion_efficiency": _diffusion_efficiency,
    "flow": _flow,
}
_REGIONAL_DEFS = {
    "reliance_raw": _reliance_raw,
    "reliance_deriv": _reliance_deriv,
    "reliance_raw_nonlinear": _reliance_raw_nonlinear,
    "reliance_deriv_nonlinear": _reliance_deriv_nonlinear,
    "hybrid_global": _hybrid_global,
    "hybrid_per_region": _hybrid_per_region,
    "hybrid_cv": _hybrid_cv,
    "hybrid_ridge": _hybrid_ridge,
    "hybrid_network": _hybrid_network,
    "information_geometry": _information_geometry,
    "multi_scale_transfer_entropy": _multi_scale_transfer_entropy,
    "criticality": _criticality,
    "granger": _granger,
    "rnc": _rnc,
    "trr": _trr,
}


def coupling(cohort, metric="communicability"):
    """Return ``(kind, array)`` for one definition.

    kind is ``"edge"`` for the canonical (n_subjects, n_edges) field, or
    ``"regional"`` for (n_subjects, n_nodes) read-outs.
    """
    if cohort.precomputed is not None:
        raise ValueError("this cohort carries a precomputed variability map only")
    from .custom import CUSTOM_METRICS, evaluate

    if metric in CUSTOM_METRICS:
        return evaluate(CUSTOM_METRICS[metric]["fn"], cohort)
    if metric == "direct":
        sc = _standardise_edges(cohort.sc, cohort.iu, cohort.n_nodes)
        fc = _standardise_edges(cohort.fc, cohort.iu, cohort.n_nodes)
        return "edge", sc * fc
    if metric == "rank_direct":
        sc = np.apply_along_axis(rankdata, 1, cohort.sc.astype(np.float64))
        fc = np.apply_along_axis(rankdata, 1, cohort.fc.astype(np.float64))
        return "edge", (_standardise_edges(sc, cohort.iu, cohort.n_nodes)
                        * _standardise_edges(fc, cohort.iu, cohort.n_nodes))
    if metric == "reliance":
        return "edge", _reliance(cohort)
    if metric in _MATRIX_DEFS:
        a = _structural(cohort)
        m = _pair_scale(_MATRIX_DEFS[metric](a))
        v = mat_to_vec(m, cohort.iu, cohort.iv)
        return "edge", _couple(v, cohort.fc, cohort.iu, cohort.n_nodes)
    if metric in _REGIONAL_DEFS:
        return "regional", _REGIONAL_DEFS[metric](cohort)
    raise ValueError("unknown metric %r" % metric)


def edge_coupling(cohort, metric="communicability"):
    """Return the canonical (n_subjects, n_edges) coupling for one definition."""
    kind, arr = coupling(cohort, metric)
    if kind != "edge":
        raise ValueError(
            "metric %r returns a regional read-out; use coupling() instead" % metric)
    return arr


def regional_readout(cohort, metric, coupling_field):
    """Regional read-out of a coupling field.

    For the direct definition this is the classical c_i = corr(SC_i., FC_i.).
    For every other definition the edge field is averaged within each region's
    connection profile, which is the natural regional summary of that field.
    """
    if metric == "direct":
        return cohort.regional_coupling()
    arr = np.asarray(coupling_field, dtype=np.float64)
    if arr.shape == (cohort.n, cohort.n_nodes):
        return arr
    out = np.full((arr.shape[0], cohort.n_nodes), np.nan)
    for i, idx in enumerate(cohort._region_profile):
        if idx.size == 0:
            continue
        out[:, i] = arr[:, idx].mean(axis=1)
    return out


def edge_projection(regional, cohort):
    """Endpoint mean of a regional map, giving an edge-compatible field."""
    reg = np.asarray(regional, dtype=np.float64)
    return 0.5 * (reg[:, cohort.iu] + reg[:, cohort.iv])


def metric_table():
    from .coords import METRICS
    rows = []
    for k, v in METRICS.items():
        rows.append((k, v["family"], v["needs"], v["description"]))
    return rows
