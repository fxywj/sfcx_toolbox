# -*- coding: utf-8 -*-
"""Functional representation: BOLD time series -> edge-wise FC.

The two pipelines follow the paper's layer 2:

    pearson    full correlation between regional mean series
    precision  partial correlation, Ledoit-Wolf shrinkage on the covariance

Everything returns the canonical (n_subjects, n_edges) upper-triangle vector.
"""

from __future__ import annotations

import numpy as np

from .data import upper_tri_index


def as_bold3d(bold, time_axis="auto"):
    """Normalise BOLD input to (n_subjects, n_nodes, n_time).

    Accepted shapes: (nodes, time), (subjects, nodes, time), (time, nodes).
    With ``time_axis='auto'`` a 2-D array is read as (nodes, time) when the
    second dimension is the longer one, and transposed otherwise.
    """
    x = np.asarray(bold, dtype=np.float64)
    if x.ndim == 2:
        if time_axis == "time":
            x = x.T
        elif time_axis == "nodes":
            pass
        else:
            x = x if x.shape[1] >= x.shape[0] else x.T
        x = x[None, :, :]
    elif x.ndim == 3:
        if time_axis == "time":
            x = np.transpose(x, (0, 2, 1))
        elif time_axis == "auto" and x.shape[1] < x.shape[2] and x.shape[0] == x.shape[2]:
            # ambiguous (n_time, n_nodes, n_subjects) is not supported; keep as is
            pass
    else:
        raise ValueError("BOLD must be 2-D or 3-D, got shape %s" % (x.shape,))
    if x.shape[2] < 3:
        raise ValueError("BOLD needs at least 3 time points, got %d" % x.shape[2])
    return x


def _interpolate_missing(x):
    """Linear interpolation along time for each (subject, node) series."""
    n, p, t = x.shape
    if np.all(np.isfinite(x)):
        return x
    out = x.copy()
    tt = np.arange(t, dtype=float)
    for s in range(n):
        block = out[s]
        bad = ~np.isfinite(block)
        if not bad.any():
            continue
        rows = np.flatnonzero(bad.any(axis=1))
        for i in rows:
            ok = ~bad[i]
            if ok.sum() < 2:
                block[i] = 0.0
            else:
                block[i, ~ok] = np.interp(tt[~ok], tt[ok], block[i, ok])
    return out


def _standardise(x):
    """Demean and scale each series to unit SD; flat series become zero."""
    x = x - x.mean(axis=2, keepdims=True)
    sd = np.sqrt((x * x).sum(axis=2, keepdims=True))
    flat = sd[..., 0] < 1e-12
    sd[sd < 1e-12] = 1.0
    z = x / sd
    z[flat] = 0.0
    return z


def _corr_from_standardised(z):
    n, p, _ = z.shape
    out = np.empty((n, p, p), dtype=np.float64)
    for s in range(n):
        out[s] = (z[s] @ z[s].T) / z.shape[2]
    return out


def _ledoit_wolf_cov(x):
    """Ledoit-Wolf shrinkage of the sample covariance towards mu * I.

    x : (n_time, n_nodes).  Returns the shrunk (n_nodes, n_nodes) covariance.
    """
    t, p = x.shape
    xc = x - x.mean(axis=0, keepdims=True)
    s = (xc.T @ xc) / t
    mu = np.trace(s) / p
    target = mu * np.eye(p)
    d2 = ((s - target) ** 2).sum()
    if d2 < 1e-18:
        return target
    b2 = 0.0
    for k in range(t):
        outer = np.outer(xc[k], xc[k])
        b2 += ((outer - s) ** 2).sum()
    b2 /= t * t
    lam = min(b2 / d2, 1.0)
    return lam * target + (1.0 - lam) * s


def _precision_from_cov(cov):
    """Partial correlation matrix from a covariance matrix."""
    d = np.sqrt(np.clip(np.diag(cov), 1e-12, None))
    cor = cov / d[:, None] / d[None, :]
    np.fill_diagonal(cor, 1.0)
    try:
        inv = np.linalg.inv(cor)
    except np.linalg.LinAlgError:
        inv = np.linalg.pinv(cor)
    dd = np.sqrt(np.clip(np.abs(np.diag(inv)), 1e-12, None))
    pc = -inv / dd[:, None] / dd[None, :]
    np.fill_diagonal(pc, 0.0)
    return pc


def bold_to_matrices(bold, kind="pearson", time_axis="auto", interpolate=True):
    """Return per-subject FC matrices, shape (n_subjects, n_nodes, n_nodes)."""
    x = as_bold3d(bold, time_axis=time_axis)
    if interpolate:
        x = _interpolate_missing(x)
    n, p, _ = x.shape
    if kind == "pearson":
        z = _standardise(x)
        return _corr_from_standardised(z)
    if kind in ("precision", "partial"):
        out = np.empty((n, p, p), dtype=np.float64)
        for s in range(n):
            out[s] = _precision_from_cov(_ledoit_wolf_cov(x[s].T))
        return out
    raise ValueError("unknown functional pipeline %r; use 'pearson' or 'precision'" % kind)


def bold_to_fc(bold, kind="pearson", time_axis="auto", interpolate=True):
    """Return the upper-triangle FC edge vector, (n_subjects, n_edges)."""
    mats = bold_to_matrices(bold, kind=kind, time_axis=time_axis, interpolate=interpolate)
    n, p, _ = mats.shape
    iu, iv = upper_tri_index(p)
    return mats[:, iu, iv].astype(np.float32)
