# -*- coding: utf-8 -*-
"""Cohort container and IO.

Everything downstream is edge-first: structural and functional connectivity are
held as (n_subjects, n_edges) vectors plus an (iu, iv) index pair.  This keeps
memory flat and makes the three dispersion estimators cheap to derive from one
another.
"""

from __future__ import annotations

import numpy as np


def upper_tri_index(n_nodes: int):
    """Row/column indices of the strict upper triangle, C order."""
    iu, iv = np.triu_indices(n_nodes, k=1)
    return iu.astype(np.int32), iv.astype(np.int32)


def vec_to_mat(v, iu, iv, n_nodes, sym=True, diag=0.0):
    """Scatter a length-E vector into an (n_nodes, n_nodes) matrix."""
    m = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    m[iu, iv] = v
    if sym:
        m[iv, iu] = v
    np.fill_diagonal(m, diag)
    return m


def mat_to_vec(m, iu, iv):
    return np.asarray(m)[iu, iv]


def _zscore_by_group(x, groups, n_groups):
    """Standardise each row of x within each group of columns.

    x      : (n, E)
    groups : (E,) integer group label, values in [0, n_groups)
    """
    out = np.empty_like(x, dtype=np.float64)
    order = np.argsort(groups, kind="stable")
    bounds = np.searchsorted(groups[order], np.arange(n_groups + 1))
    for g in range(n_groups):
        idx = order[bounds[g]:bounds[g + 1]]
        if idx.size < 2:
            out[:, idx] = 0.0
            continue
        block = x[:, idx].astype(np.float64)
        mu = block.mean(axis=1, keepdims=True)
        sd = block.std(axis=1, keepdims=True)
        sd[sd < 1e-12] = 1.0
        out[:, idx] = (block - mu) / sd
    return out


class Cohort:
    """A group of subjects described on a fixed parcellation.

    Parameters
    ----------
    sc, fc : (n_subjects, n_edges) array_like
        Structural and functional connectivity, upper-triangle order.
    n_nodes : int
    networks : (n_nodes,) int
        Community label per region (e.g. Yeo-7).
    hierarchy : (n_nodes,) float
        Position on the cortical hierarchy, low = unimodal.
    coords : (n_nodes, 3) float, optional
        Sphere coordinates for spin tests.
    behaviour : (n_subjects,) float, optional
        Target variable for individual-difference analyses.
    covariates : (n_subjects, k) float, optional
    state : str
    subject_ids : (n_subjects,) array_like, optional
    """

    def __init__(self, sc, fc, n_nodes, networks, hierarchy,
                 coords=None, behaviour=None, covariates=None,
                 state="rest", subject_ids=None, bold=None):
        self.sc = np.ascontiguousarray(sc, dtype=np.float32)
        self.fc = np.ascontiguousarray(fc, dtype=np.float32)
        if self.sc.shape != self.fc.shape:
            raise ValueError("sc and fc must have the same shape, got %s and %s"
                             % (self.sc.shape, self.fc.shape))
        if self.sc.ndim != 2:
            raise ValueError("sc and fc must be 2-D (n_subjects, n_edges)")
        self.n, self.n_edges = self.sc.shape
        self.n_nodes = int(n_nodes)
        if self.n_edges != self.n_nodes * (self.n_nodes - 1) // 2:
            raise ValueError(
                "n_edges=%d is not the upper triangle of n_nodes=%d"
                % (self.n_edges, self.n_nodes))
        self.iu, self.iv = upper_tri_index(self.n_nodes)
        self.networks = np.asarray(networks, dtype=int)
        self.hierarchy = np.asarray(hierarchy, dtype=float)
        if self.networks.shape != (self.n_nodes,):
            raise ValueError("networks must have length n_nodes")
        if self.hierarchy.shape != (self.n_nodes,):
            raise ValueError("hierarchy must have length n_nodes")
        self.coords = None if coords is None else np.asarray(coords, dtype=float)
        self.behaviour = None if behaviour is None else np.asarray(behaviour, dtype=float)
        self.covariates = None if covariates is None else np.asarray(covariates, dtype=float)
        self.state = state
        self.subject_ids = (np.arange(self.n) if subject_ids is None
                            else np.asarray(subject_ids))
        if bold is not None:
            bold = np.asarray(bold, dtype=np.float32)
            if bold.ndim == 2:
                bold = bold[None, :, :]
            if bold.ndim != 3:
                raise ValueError("bold must be 3-D (n_subjects, n_nodes, n_time)")
            if bold.shape[0] != self.n or bold.shape[1] != self.n_nodes:
                raise ValueError(
                    "bold must be (n_subjects=%d, n_nodes=%d, n_time), got %s"
                    % (self.n, self.n_nodes, bold.shape))
            self.bold = np.ascontiguousarray(bold)
        else:
            self.bold = None

        self._region_edges = [np.flatnonzero(self.iu == i) for i in range(self.n_nodes)]
        self._region_profile = [np.flatnonzero((self.iu == i) | (self.iv == i))
                                for i in range(self.n_nodes)]
        self._adjacency = None

    # ------------------------------------------------------------ basics --

    def __repr__(self):
        return ("Cohort(n=%d, n_nodes=%d, n_edges=%d, state=%r, networks=%d, behaviour=%s)"
                % (self.n, self.n_nodes, self.n_edges, self.state,
                   len(np.unique(self.networks)), self.behaviour is not None))

    @property
    def n_networks(self):
        return int(self.networks.max()) + 1

    def adjacency(self, threshold=0.0):
        """Group-mean structural matrix, binarised above `threshold`."""
        if self._adjacency is None:
            mean_sc = self.sc.mean(axis=0).astype(np.float64)
            m = vec_to_mat(mean_sc, self.iu, self.iv, self.n_nodes)
            self._adjacency = (m > max(threshold, 1e-9)).astype(np.float64)
        return self._adjacency

    def group_matrix(self, which="sc"):
        v = (self.sc if which == "sc" else self.fc).mean(axis=0).astype(np.float64)
        return vec_to_mat(v, self.iu, self.iv, self.n_nodes)

    # ------------------------------------------------------- coupling ----

    def regional_coupling(self, standardise=True):
        """c_i = corr(SC_i., FC_i.) for each subject and region -> (n, n_nodes).

        Regions with no valid edges (or zero variance) come back as NaN.
        """
        out = np.full((self.n, self.n_nodes), np.nan, dtype=np.float64)
        for i, idx in enumerate(self._region_profile):
            if idx.size < 3:
                continue
            a = self.sc[:, idx].astype(np.float64)
            b = self.fc[:, idx].astype(np.float64)
            a = a - a.mean(axis=1, keepdims=True)
            b = b - b.mean(axis=1, keepdims=True)
            sa = np.sqrt((a * a).sum(axis=1))
            sb = np.sqrt((b * b).sum(axis=1))
            ok = (sa > 1e-12) & (sb > 1e-12)
            c = np.full(a.shape[0], np.nan)
            c[ok] = (a[ok] * b[ok]).sum(axis=1) / (sa[ok] * sb[ok])
            out[:, i] = c
        return out

    def split_half(self, seed=0):
        """Random split of subjects into two halves (indices)."""
        rng = np.random.default_rng(seed)
        perm = rng.permutation(self.n)
        h = self.n // 2
        return perm[:h], perm[h:2 * h]

    # ------------------------------------------------------------- io ----

    @classmethod
    def from_arrays(cls, sc, fc, networks, hierarchy, **kw):
        sc = np.asarray(sc)
        n_nodes = int(round((1 + np.sqrt(1 + 8 * sc.shape[1])) / 2))
        return cls(sc, fc, n_nodes, networks, hierarchy, **kw)

    @classmethod
    def from_npz(cls, path, sc_key="sc", fc_key="fc", bold_key="bold",
                 iu_key="iu_i", iv_key="iv_j",
                 networks_key="networks", hierarchy_key="hierarchy",
                 n_nodes_key="n_nodes", **kw):
        """Load an .npz with SC/FC edge vectors and parcel annotation."""
        z = np.load(path, allow_pickle=True)
        keys = set(z.files)
        sc = z[sc_key]
        fc = z[fc_key]
        if sc.ndim == 1:
            sc = sc[None, :]
            fc = fc[None, :]
        if iu_key in keys and iv_key in keys:
            iu, iv = z[iu_key].astype(np.int64), z[iv_key].astype(np.int64)
            n_nodes = int(max(iu.max(), iv.max()) + 1)
        else:
            n_nodes = int(z[n_nodes_key]) if n_nodes_key in keys else \
                int(round((1 + np.sqrt(1 + 8 * sc.shape[1])) / 2))
        networks = z[networks_key] if networks_key in keys else np.zeros(n_nodes, dtype=int)
        hierarchy = z[hierarchy_key] if hierarchy_key in keys else np.linspace(1, 7, n_nodes)
        bold = z[bold_key] if bold_key in keys else None
        return cls(sc, fc, n_nodes, networks, hierarchy, bold=bold, **kw)

    @classmethod
    def from_precomputed_variability(cls, v_regional, networks, hierarchy, **kw):
        """Build a light cohort whose SC/FC placeholders carry a ready map.

        Used by the reproduce route, where per-subject coupling matrices are not
        available but the variability map is.
        """
        v = np.asarray(v_regional, dtype=float)
        n_nodes = v.size
        obj = cls.__new__(cls)
        obj.sc = np.zeros((1, n_nodes * (n_nodes - 1) // 2), dtype=np.float32)
        obj.fc = obj.sc.copy()
        obj.n, obj.n_edges = 1, obj.sc.shape[1]
        obj.n_nodes = n_nodes
        obj.iu, obj.iv = upper_tri_index(n_nodes)
        obj.networks = np.asarray(networks, dtype=int)
        obj.hierarchy = np.asarray(hierarchy, dtype=float)
        obj.coords = None
        obj.behaviour = None
        obj.covariates = None
        obj.state = kw.get("state", "precomputed")
        obj.subject_ids = np.array([0])
        obj._region_edges = [np.flatnonzero(obj.iu == i) for i in range(n_nodes)]
        obj._region_profile = [np.flatnonzero((obj.iu == i) | (obj.iv == i))
                               for i in range(n_nodes)]
        obj._adjacency = None
        obj._precomputed = v
        return obj

    @property
    def precomputed(self):
        return getattr(self, "_precomputed", None)
