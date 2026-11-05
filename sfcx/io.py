# -*- coding: utf-8 -*-
"""Data intake: matrices, BOLD series, annotation, and format conversion.

Three shapes are accepted for connectivity:

    (n_nodes, n_nodes)              one subject
    (n_subjects, n_nodes, n_nodes)  full matrices
    (n_subjects, n_edges)           already vectorised, upper triangle

BOLD is (n_subjects, n_nodes, n_time), (n_nodes, n_time) or (n_time, n_nodes).
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np

from .data import Cohort, upper_tri_index
from .functional import bold_to_fc, bold_to_matrices


# --------------------------------------------------------------- helpers ----

def n_nodes_from_edges(n_edges):
    k = int(round((1 + np.sqrt(1 + 8 * int(n_edges))) / 2))
    if k * (k - 1) // 2 != int(n_edges):
        raise ValueError("%d edges is not an upper triangle" % int(n_edges))
    return k


def _is_square(x):
    return x.ndim == 2 and x.shape[0] == x.shape[1]


def to_edge_vectors(x, name="connectivity"):
    """Normalise a connectivity array to (n_subjects, n_edges) plus n_nodes."""
    arr = np.asarray(x, dtype=np.float64)
    if arr.ndim == 2:
        if _is_square(arr):
            p = arr.shape[0]
            iu, iv = upper_tri_index(p)
            return arr[iu, iv][None, :].astype(np.float32), p
        p = n_nodes_from_edges(arr.shape[1])
        return arr.astype(np.float32), p
    if arr.ndim == 3:
        n, p, q = arr.shape
        if p != q:
            raise ValueError("%s: last two dimensions must be square, got %s" % (name, arr.shape))
        iu, iv = upper_tri_index(p)
        return arr[:, iu, iv].astype(np.float32), p
    raise ValueError("%s must be 2-D or 3-D, got shape %s" % (name, arr.shape))


def _symmetrise_and_zero_diag(m):
    m = np.asarray(m, dtype=np.float64)
    m = np.nan_to_num(m, nan=0.0, posinf=0.0, neginf=0.0)
    m = 0.5 * (m + m.T)
    np.fill_diagonal(m, 0.0)
    return m


def default_hierarchy(cohort):
    """Node strength rank in [1, 7] when no cortical hierarchy is supplied."""
    s = np.abs(cohort.group_matrix("sc")).sum(axis=1)
    if np.allclose(s, s[0]):
        return np.linspace(1.0, 7.0, cohort.n_nodes)
    order = np.argsort(np.argsort(s)).astype(float)
    return 1.0 + 6.0 * order / max(cohort.n_nodes - 1, 1)


def default_networks(n_nodes, n_networks=1):
    return np.zeros(int(n_nodes), dtype=int) if n_networks <= 1 else \
        np.arange(int(n_nodes)) % int(n_networks)


# ---------------------------------------------------------------- build ----

def build_cohort(sc=None, fc=None, bold=None, networks=None, hierarchy=None,
                 coords=None, behaviour=None, covariates=None, state="rest",
                 subject_ids=None, fc_kind="pearson", time_axis="auto",
                 n_networks=None):
    """Assemble a :class:`Cohort` from any mix of connectivity / BOLD inputs.

    Returns ``(cohort, report)``; ``report`` records what was inferred and what
    was defaulted, so the caller can surface honest warnings.
    """
    report = {"inferred": [], "warnings": [], "sources": {}}
    if sc is None and bold is None:
        raise ValueError("need at least a structural matrix (sc) or BOLD series (bold)")

    n_nodes = None
    if sc is not None:
        sc_v, p = to_edge_vectors(sc, "sc")
        n_nodes = p
        report["sources"]["sc"] = "matrix" if np.asarray(sc).ndim else "vector"
        report["inferred"].append("SC: %d subjects x %d edges (P=%d)" % (sc_v.shape[0], sc_v.shape[1], p))
    if bold is not None:
        mats = bold_to_matrices(bold, kind=fc_kind, time_axis=time_axis)
        p_bold = mats.shape[1]
        if n_nodes is not None and p_bold != n_nodes:
            raise ValueError("SC has %d nodes but BOLD has %d" % (n_nodes, p_bold))
        n_nodes = p_bold
        report["inferred"].append(
            "BOLD -> FC (%s): %d subjects x %d time points" % (fc_kind, mats.shape[0], np.asarray(bold).shape[-1]))
    if fc is not None:
        fc_v, p = to_edge_vectors(fc, "fc")
        if n_nodes is not None and p != n_nodes:
            raise ValueError("SC has %d nodes but FC has %d" % (n_nodes, p))
        n_nodes = p
        report["sources"]["fc"] = "matrix"
    elif bold is not None:
        fc_v = bold_to_fc(bold, kind=fc_kind, time_axis=time_axis)
        report["sources"]["fc"] = "derived from BOLD (%s)" % fc_kind
    else:
        raise ValueError("need an FC matrix or a BOLD series")

    if sc_v.shape[0] != fc_v.shape[0]:
        raise ValueError("SC has %d subjects but FC has %d"
                         % (sc_v.shape[0], fc_v.shape[0]))

    if networks is None:
        n_networks = 1 if n_networks is None else int(n_networks)
        networks = default_networks(n_nodes, n_networks)
        report["warnings"].append({"zh":
            "未提供网络标签：已生成 %d 个社区；网络内/网络间对比只有使用真实标签才有意义" % n_networks,
            "en":
            "no network labels supplied: %d community(ies) generated; "
            "the within/between layer contrast is only meaningful with real labels" % n_networks})
    else:
        networks = np.asarray(networks).ravel().astype(int)
        if networks.size != n_nodes:
            raise ValueError("networks must have %d entries, got %d" % (n_nodes, networks.size))

    hierarchy_provided = hierarchy is not None
    if hierarchy is None:
        s = np.asarray(sc_v, dtype=np.float64)
        strength = np.zeros(n_nodes)
        iu, iv = upper_tri_index(n_nodes)
        np.add.at(strength, iu, np.abs(s).mean(axis=0))
        np.add.at(strength, iv, np.abs(s).mean(axis=0))
        if np.allclose(strength, strength[0]):
            hierarchy = np.linspace(1.0, 7.0, n_nodes)
        else:
            order = np.argsort(np.argsort(strength)).astype(float)
            hierarchy = 1.0 + 6.0 * order / max(n_nodes - 1, 1)
        report["warnings"].append({"zh":
            "未提供皮层层级：以节点结构强度作为代理；此时变异轴是强度轴，而非层状轴",
            "en":
            "no cortical hierarchy supplied: node structural strength used as a proxy; "
            "the variability axis is then a strength axis, not a laminar one"})
    else:
        hierarchy = np.asarray(hierarchy, dtype=float).ravel()
        if hierarchy.size != n_nodes:
            raise ValueError("hierarchy must have %d entries, got %d" % (n_nodes, hierarchy.size))

    if coords is not None:
        coords = np.asarray(coords, dtype=float)
        if coords.shape != (n_nodes, 3):
            raise ValueError("coords must be (n_nodes, 3), got %s" % (coords.shape,))

    if behaviour is not None:
        behaviour = np.asarray(behaviour, dtype=float).ravel()
        if behaviour.size != sc_v.shape[0]:
            raise ValueError("behaviour must have %d entries" % sc_v.shape[0])
    if covariates is not None:
        covariates = np.asarray(covariates, dtype=float)
        if covariates.ndim == 1:
            covariates = covariates[:, None]
        if covariates.shape[0] != sc_v.shape[0]:
            raise ValueError("covariates must have %d rows" % sc_v.shape[0])

    coh = Cohort(sc=sc_v, fc=fc_v, n_nodes=n_nodes, networks=networks,
                 hierarchy=hierarchy, coords=coords, behaviour=behaviour,
                 covariates=covariates, state=state, subject_ids=subject_ids,
                 bold=bold)
    coh.meta = dict(report=report, hierarchy_provided=hierarchy_provided,
                    n_networks_provided=networks is not None)
    return coh, report


# ----------------------------------------------------------------- load ----

def load_any(path, key=None):
    """Load an array or a dict of arrays from .npz/.npy/.mat/.csv/.tsv/.json."""
    ext = os.path.splitext(str(path))[1].lower()
    if ext == ".npz":
        z = np.load(path, allow_pickle=True)
        if key:
            return z[key]
        return {k: z[k] for k in z.files}
    if ext == ".npy":
        return np.load(path, allow_pickle=True)
    if ext in (".mat", ".mat73"):
        return _load_mat(path)
    if ext in (".csv", ".tsv", ".txt", ".dat"):
        return _load_table(path, delimiter="\t" if ext == ".tsv" else ",")
    if ext == ".json":
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    raise ValueError("unsupported file type %r" % ext)


def _load_table(path, delimiter=","):
    rows = []
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if not row:
                continue
            try:
                rows.append([float(x) for x in row])
            except ValueError:
                continue  # header line
    if not rows:
        raise ValueError("no numeric rows found in %s" % path)
    width = max(len(r) for r in rows)
    rows = [r for r in rows if len(r) == width]
    return np.asarray(rows, dtype=np.float64)


def _load_mat(path):
    try:
        from scipy.io import loadmat
    except ImportError as exc:  # pragma: no cover
        raise ImportError("scipy is required to read .mat files") from exc
    try:
        raw = loadmat(path, squeeze_me=True, struct_as_record=False)
        return {k: v for k, v in raw.items() if not k.startswith("__")}
    except NotImplementedError:
        return _load_mat73(path)


def _load_mat73(path):
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "this .mat is v7.3 (HDF5); install h5py or re-save as v7") from exc
    out = {}
    with h5py.File(path, "r") as f:
        def walk(node, prefix=""):
            for k, v in node.items():
                if isinstance(v, h5py.Dataset):
                    out[prefix + k] = np.asarray(v).T
                else:
                    walk(v, prefix + k + ".")
        walk(f)
    return out


def _pick(container, names, default=None):
    if not isinstance(container, dict):
        return default
    for n in names:
        if n in container:
            return container[n]
    return default


def load_cohort(path, sc_key="sc", fc_key="fc", bold_key="bold",
                networks_key="networks", hierarchy_key="hierarchy",
                coords_key="coords", behaviour_key="behaviour",
                covariates_key="covariates", state="rest", fc_kind="pearson",
                time_axis="auto"):
    """Load a cohort from a single file holding any of the supported keys."""
    obj = load_any(path)
    if isinstance(obj, np.ndarray):
        if obj.ndim == 2 and _is_square(obj):
            return build_cohort(sc=obj, fc=obj.copy(), state=state)
        return build_cohort(sc=obj, fc=obj.copy(), state=state)
    sc = _pick(obj, [sc_key, "SC", "sc_matrix", "struct"])
    fc = _pick(obj, [fc_key, "FC", "fc_matrix", "func"])
    bold = _pick(obj, [bold_key, "BOLD", "bold", "timeseries", "ts"])
    networks = _pick(obj, [networks_key, "network", "yeo", "communities"])
    hierarchy = _pick(obj, [hierarchy_key, "hier", "level"])
    coords = _pick(obj, [coords_key, "xyz", "sphere"])
    behaviour = _pick(obj, [behaviour_key, "behavior", "behav", "y"])
    covariates = _pick(obj, [covariates_key, "cov", "covariates"])
    if sc is None and bold is None:
        arrays = {k: v for k, v in obj.items() if isinstance(v, np.ndarray)}
        if len(arrays) == 2:
            keys = sorted(arrays)
            sc, fc = arrays[keys[0]], arrays[keys[1]]
        else:
            raise ValueError("could not find sc/fc/bold keys in %s; found: %s"
                             % (path, ", ".join(sorted(obj))))
    return build_cohort(sc=sc, fc=fc, bold=bold, networks=networks,
                        hierarchy=hierarchy, coords=coords, behaviour=behaviour,
                        covariates=covariates, state=state, fc_kind=fc_kind,
                        time_axis=time_axis)


# --------------------------------------------------------------- export ----

def export_npz(path, cohort, extra=None):
    payload = dict(sc=cohort.sc, fc=cohort.fc, iu_i=cohort.iu, iv_j=cohort.iv,
                   networks=cohort.networks, hierarchy=cohort.hierarchy,
                   n_nodes=cohort.n_nodes, state=cohort.state)
    if cohort.bold is not None:
        payload["bold"] = cohort.bold
    if cohort.coords is not None:
        payload["coords"] = cohort.coords
    if cohort.behaviour is not None:
        payload["behaviour"] = cohort.behaviour
    if cohort.covariates is not None:
        payload["covariates"] = cohort.covariates
    if extra:
        payload.update(extra)
    np.savez_compressed(path, **payload)
    return path


def export_mat(path, cohort, extra=None):
    """Write a MATLAB-friendly v7 .mat bundle (readable by the MATLAB toolbox)."""
    from scipy.io import savemat
    payload = dict(sc=cohort.sc, fc=cohort.fc, iu=cohort.iu + 1, iv=cohort.iv + 1,
                   n_nodes=cohort.n_nodes, networks=cohort.networks,
                   hierarchy=cohort.hierarchy, state=cohort.state)
    if cohort.bold is not None:
        payload["bold"] = cohort.bold
    if cohort.coords is not None:
        payload["coords"] = cohort.coords
    if cohort.behaviour is not None:
        payload["behaviour"] = cohort.behaviour
    if cohort.covariates is not None:
        payload["covariates"] = cohort.covariates
    if extra:
        payload.update(extra)
    savemat(path, payload)
    return path
