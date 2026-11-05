# -*- coding: utf-8 -*-
"""Generate the 5-subject teaching dataset shipped with the toolbox.

    python examples/make_example_data.py

Everything is deterministic (SEED below).  BOLD is generated first and FC is
computed from it with the toolbox's own ``pearson`` pipeline, so a user who
imports the BOLD series and a user who imports the FC matrices get the same
numbers -- that is the point of the dataset.

How the data are made
---------------------
1. Six network caps are placed on a sphere; a cortical-hierarchy gradient runs
   from network 0 (low) to network 5 (high).
2. SC is a distance kernel with a within-network boost, thresholded to ~25 %
   density, then perturbed per subject.
3. A diffusion kernel ``exp(-tau * L)`` on the group SC gives a base
   correlation matrix ``C0`` -- the functional pattern that structure predicts.
4. Each subject's FC mixes ``C0`` with a subject-specific random correlation
   matrix, and the mixing weight rises with the hierarchy.  The mixture
   ``diag(a) C0 diag(a) + diag(b) R diag(b)`` is still a valid correlation
   matrix, so the data are well formed while higher-order regions carry more
   inter-individual scatter -- which is what produces a positive variability
   axis.
5. BOLD is drawn with exactly that sample correlation (AR(1) spectrum), so
   ``pearson FC`` recovers the stored FC to numerical precision.

Output (examples/data/):
    demo5_group.npz / demo5_group.mat   five subjects in one bundle
    demo5_networks.csv                  per-region network label (0..5)
    demo5_hierarchy.csv                 per-region cortical hierarchy (1..7)
    demo5_coords.csv                    per-region unit-sphere coordinates
    demo5_behaviour.csv                 one illustrative score per subject
    demo5_covariates.csv                two illustrative nuisance variables
    sub-01_sc.csv ... sub-05_bold.csv   per-subject matrices
    demo5_expected.json                 reference numbers to check an import
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy import linalg

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
sys.path.insert(0, os.path.dirname(HERE))

from sfcx.functional import bold_to_matrices  # noqa: E402

SEED = 20261105
N_SUBJECTS = 5
N_NODES = 64
N_TIME = 120
N_NETWORKS = 6
TR = 2.0

K_NN = 14              # neighbours kept per row -> ~25 % density
TAU = 2.5              # diffusion time of the structure-predicted FC
G_LOW, G_HIGH = 0.30, 0.50   # subject-specific share at hierarchy 0 / 1
AR1 = 0.45


# ------------------------------------------------------------- geometry ----

def _geometry(rng):
    """Network caps on a sphere, a hierarchy gradient, and a thresholded SC."""
    p, k = N_NODES, N_NETWORKS
    centres = rng.normal(size=(k, 3))
    centres /= np.linalg.norm(centres, axis=1, keepdims=True)

    networks = np.repeat(np.arange(k), p // k)
    if networks.size < p:
        networks = np.concatenate([networks, np.full(p - networks.size, k - 1)])

    coords = np.zeros((p, 3))
    for g in range(k):
        idx = np.flatnonzero(networks == g)
        d = centres[g] + 0.55 * rng.normal(size=(idx.size, 3))
        coords[idx] = d / np.linalg.norm(d, axis=1, keepdims=True)

    hierarchy = networks / (k - 1) + 0.08 * rng.normal(size=p)
    hierarchy = np.clip(hierarchy, 0, 1)

    dist = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    same = networks[:, None] == networks[None, :]
    sc = np.exp(-dist / 0.55) * (1.0 + 0.9 * same)
    thr = np.partition(sc, -K_NN, axis=1)[:, -K_NN][:, None]
    sc = np.where(sc >= thr, sc, 0.0)
    sc = np.maximum(sc, sc.T)
    np.fill_diagonal(sc, 0.0)
    return sc, networks.astype(int), 1.0 + 6.0 * hierarchy, coords


def _subject_sc(rng, sc_group):
    """Keep the wiring, let the weights wobble a little per subject."""
    mask = sc_group > 0
    out = np.empty((N_SUBJECTS,) + sc_group.shape, dtype=np.float64)
    for s in range(N_SUBJECTS):
        w = np.where(mask, sc_group * np.exp(0.08 * rng.normal(size=sc_group.shape)), 0.0)
        w = np.maximum(w, w.T)
        np.fill_diagonal(w, 0.0)
        out[s] = w
    return out


# ------------------------------------------------------- functional model ---

def _diffusion_corr(sc, tau=TAU):
    """exp(-tau * L) normalised to a correlation matrix."""
    deg = sc.sum(axis=1)
    deg[deg < 1e-12] = 1.0
    lap = np.diag(deg) - sc
    k = linalg.expm(-tau * lap)
    dg = np.sqrt(np.diag(k))
    c = k / dg[:, None] / dg[None, :]
    np.fill_diagonal(c, 1.0)
    return c


def _ar1(x, phi=AR1):
    """Slow drift along the last axis, variance restored afterwards."""
    eps = np.sqrt(1.0 - phi ** 2) * x
    out = np.empty_like(eps)
    out[..., 0] = eps[..., 0]
    for t in range(1, eps.shape[-1]):
        out[..., t] = phi * out[..., t - 1] + eps[..., t]
    return out


def _bold_with_corr(rng, corr):
    """BOLD (P, T) whose sample correlation matrix equals ``corr``."""
    p, t = corr.shape[0], N_TIME
    z = _ar1(rng.normal(size=(t, p)))
    z = z - z.mean(axis=0, keepdims=True)
    s = (z.T @ z) / t
    w, v = np.linalg.eigh(s)
    w = np.clip(w, 1e-10, None)
    z = z @ ((v / np.sqrt(w)) @ v.T)
    chol = np.linalg.cholesky(corr + 1e-8 * np.eye(p))
    return (z @ chol.T).T


def _subject_fc(rng, sc_group, hierarchy):
    """Per-subject correlation matrix from the structure prior plus scatter."""
    h = (hierarchy - 1.0) / 6.0
    base = _diffusion_corr(sc_group)
    g = G_LOW + (G_HIGH - G_LOW) * h
    a, b = np.sqrt(1.0 - g), np.sqrt(g)

    fc = np.empty((N_SUBJECTS, N_NODES, N_NODES), dtype=np.float64)
    bold = np.empty((N_SUBJECTS, N_NODES, N_TIME), dtype=np.float64)
    for s in range(N_SUBJECTS):
        random_corr = np.corrcoef(rng.normal(size=(N_NODES, N_TIME)))
        target = (a[:, None] * a[None, :]) * base + (b[:, None] * b[None, :]) * random_corr
        bold[s] = _bold_with_corr(rng, target)
        fc[s] = np.corrcoef(bold[s])
    return bold, fc


# ---------------------------------------------------------------- writing ---

def _write_csv(path, arr, fmt="%.6g"):
    np.savetxt(path, np.atleast_2d(arr), delimiter=",", fmt=fmt)


def _write_labelled(path, header, values):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(",".join(header) + "\n")
        for row in values:
            f.write(",".join(str(v) for v in row) + "\n")


def main():
    rng = np.random.default_rng(SEED)
    os.makedirs(OUT, exist_ok=True)

    sc_group, networks, hierarchy, coords = _geometry(rng)
    sc = _subject_sc(rng, sc_group)
    bold, fc = _subject_fc(rng, sc_group, hierarchy)
    bold = bold.astype(np.float32)
    fc = bold_to_matrices(bold, kind="pearson").astype(np.float32)

    subject_ids = np.array(["sub-%02d" % (i + 1) for i in range(N_SUBJECTS)])

    # Illustrative behaviour: a subject-level coupling summary plus scatter.
    # With five subjects this is a teaching prop, not a result.
    regional = np.array([np.nanmean(np.diag(fc[s] @ sc[s]) / (np.abs(sc[s]).sum(1) + 1e-12))
                         for s in range(N_SUBJECTS)])
    z = (regional - regional.mean()) / (regional.std() + 1e-12)
    behaviour = 0.9 * z + 0.5 * rng.normal(size=N_SUBJECTS)
    covariates = np.column_stack([rng.normal(size=N_SUBJECTS),
                                  rng.normal(size=N_SUBJECTS)])

    bundle = dict(
        sc=sc.astype(np.float32), fc=fc, bold=bold,
        networks=networks.astype(np.int32),
        hierarchy=hierarchy.astype(np.float32), coords=coords.astype(np.float32),
        behaviour=behaviour.astype(np.float32),
        covariates=covariates.astype(np.float32),
        subject_ids=subject_ids, n_nodes=np.int32(N_NODES), tr=np.float32(TR),
    )
    np.savez_compressed(os.path.join(OUT, "demo5_group.npz"), **bundle)

    from scipy.io import savemat
    savemat(os.path.join(OUT, "demo5_group.mat"), bundle, do_compression=True)

    for s, sid in enumerate(subject_ids):
        _write_csv(os.path.join(OUT, "%s_sc.csv" % sid), sc[s])
        _write_csv(os.path.join(OUT, "%s_fc.csv" % sid), fc[s])
        _write_csv(os.path.join(OUT, "%s_bold.csv" % sid), bold[s])

    _write_labelled(os.path.join(OUT, "demo5_networks.csv"),
                    ["region", "network"], enumerate(networks))
    _write_labelled(os.path.join(OUT, "demo5_hierarchy.csv"),
                    ["region", "hierarchy"],
                    ((i, round(float(hierarchy[i]), 4)) for i in range(N_NODES)))
    _write_labelled(os.path.join(OUT, "demo5_coords.csv"),
                    ["region", "x", "y", "z"],
                    ((i,) + tuple(np.round(coords[i], 5)) for i in range(N_NODES)))
    _write_labelled(os.path.join(OUT, "demo5_behaviour.csv"),
                    ["subject", "score"],
                    ((subject_ids[i], round(float(behaviour[i]), 4))
                     for i in range(N_SUBJECTS)))
    _write_labelled(os.path.join(OUT, "demo5_covariates.csv"),
                    ["subject", "age_like", "motion_like"],
                    ((subject_ids[i], round(float(covariates[i, 0]), 4),
                      round(float(covariates[i, 1]), 4)) for i in range(N_SUBJECTS)))

    print("written to %s" % OUT)
    print("  sc %s  fc %s  bold %s" % (sc.shape, fc.shape, bold.shape))
    print("  %d subjects, %d nodes, %d time points, TR=%.1f s"
          % (N_SUBJECTS, N_NODES, N_TIME, TR))
    print("  sc density %.1f%%" % (100.0 * (sc[0] > 0).sum() / (N_NODES * (N_NODES - 1)))
          )


if __name__ == "__main__":
    main()
