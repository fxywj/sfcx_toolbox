# -*- coding: utf-8 -*-
"""A small generative cohort.

Nothing about the coordinate system is hard-coded here: structure is modular,
function is a noisy, individually weighted copy of structure, and the noise load
rises towards association cortex.  Whether the axis flips across estimators is
an outcome of those assumptions, not a scripted result.
"""

from __future__ import annotations

import numpy as np

from .data import Cohort, upper_tri_index


def make_cohort(n_subjects=320, n_nodes=180, n_networks=6, seed=0,
                behaviour=True, noise_slope=0.9):
    rng = np.random.default_rng(seed)

    # ---- geometry: networks are caps on a sphere --------------------------
    centres = rng.normal(size=(n_networks, 3))
    centres /= np.linalg.norm(centres, axis=1, keepdims=True)
    net = np.repeat(np.arange(n_networks), n_nodes // n_networks)
    if net.size < n_nodes:
        net = np.concatenate([net, np.full(n_nodes - net.size, n_networks - 1)])

    xyz = np.zeros((n_nodes, 3))
    for g in range(n_networks):
        idx = np.flatnonzero(net == g)
        dirs = centres[g] + 0.55 * rng.normal(size=(idx.size, 3))
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        xyz[idx] = dirs

    # hierarchy: unimodal (network 0) low, association (last) high
    hierarchy = net / (n_networks - 1) + 0.08 * rng.normal(size=n_nodes)
    hierarchy = np.clip(hierarchy, 0, 1)
    hierarchy_order = 1 + hierarchy * 6.0

    # ---- structure: distance decay plus a within-network boost ------------
    d = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=2))
    same = net[:, None] == net[None, :]
    sc = np.exp(-d / 0.55) * (1.0 + 0.9 * same)
    # keep the k strongest links per region, as a thresholded connectome would
    k_nn = 14
    thr = np.partition(sc, -k_nn, axis=1)[:, -k_nn][:, None]
    sc = np.where(sc >= thr, sc, 0.0)
    sc = np.maximum(sc, sc.T)
    np.fill_diagonal(sc, 0.0)

    iu, iv = upper_tri_index(n_nodes)
    sc_e = sc[iu, iv].astype(np.float32)
    keep = sc_e > 0
    sc_e[~keep] = 0.0

    # ---- function: weighted copy of structure plus subject noise ----------
    base = sc_e / (np.abs(sc_e).max() + 1e-9)

    # smooth spatial jitter plus a few hotspots: real dispersion maps carry
    # local peaks rather than a pure gradient.
    jitter = rng.normal(size=n_nodes)
    adj = (sc > 0).astype(float)
    deg = np.maximum(adj.sum(axis=1), 1.0)
    for _ in range(2):
        jitter = 0.35 * jitter + 0.65 * (adj @ jitter) / deg
    jitter = (jitter - jitter.mean()) / (jitter.std() + 1e-12)
    for c in rng.choice(n_nodes, size=max(6, n_nodes // 12), replace=False):
        jitter += 1.4 * np.exp(-((np.arange(n_nodes) - c) ** 2) / (2 * 3.0 ** 2))

    noise_load = (0.26 + 1.35 * hierarchy) * np.exp(0.7 * jitter)
    # regions with looser structural constraint vary more across people: this
    # is what turns a smooth gradient into local dispersion peaks
    sd_map = 0.05 + 0.30 * (noise_load / (noise_load.max() + 1e-12))
    subj_scale = np.exp(rng.normal(size=(n_subjects, n_nodes)) * sd_map[None, :])

    same_e = same[iu, iv].astype(np.float32)
    amp_e = (0.5 * (noise_load[iu] + noise_load[iv])).astype(np.float32)
    # functional noise rides on the strength of the structural link, which is
    # what makes within-network edges the more dispersed ones
    amp_e = amp_e * (0.15 + 0.85 * np.abs(base)).astype(np.float32)

    fc = np.empty((n_subjects, sc_e.size), dtype=np.float32)
    for s in range(n_subjects):
        gain = 1.0 + rng.normal(0, 0.20)
        eps = rng.normal(size=sc_e.size).astype(np.float32)
        fc[s] = gain * base + noise_slope * 0.22 * amp_e * eps
    fc = np.clip(fc, -1.4, 1.4)

    cov = np.column_stack([rng.normal(size=n_subjects), rng.normal(size=n_subjects)])

    # ---- synthetic BOLD time series (n_subjects, n_nodes, n_time) ---------
    n_time = 200
    bold = np.empty((n_subjects, n_nodes, n_time), dtype=np.float32)
    W = sc / (sc.sum(axis=0, keepdims=True) + 1e-10)
    for s in range(n_subjects):
        ts = rng.normal(size=(n_nodes, n_time)).astype(np.float32)
        for t in range(1, n_time):
            ts[:, t] = 0.7 * ts[:, t - 1] + 0.3 * ts[:, t] + 0.2 * (W @ ts[:, t - 1])
        bold[s] = ts

    coh = Cohort(sc=sc_e[None, :].repeat(n_subjects, axis=0), fc=fc,
                 n_nodes=n_nodes, networks=net, hierarchy=hierarchy_order,
                 coords=xyz, behaviour=None, covariates=cov, state="synthetic",
                 bold=bold)

    if behaviour:
        # a weak individual signal planted at the edge scale, as the paper
        # argues: the regional summary averages it away.
        from .metrics import edge_coupling
        coup = edge_coupling(coh, "communicability")
        dev = coup - coup.mean(axis=0, keepdims=True)
        order = np.argsort(-dev, axis=1)
        ranked = np.take_along_axis(dev, order, axis=1)
        r = np.arange(1, ranked.shape[1] + 1, dtype=float)
        r = (r - r.mean()) / r.std()
        slope = (ranked * r[None, :]).mean(axis=1)
        z = (slope - slope.mean()) / (slope.std() + 1e-12)
        coh.behaviour = (0.35 * z + rng.normal(0, 1.0, size=n_subjects)).astype(float)

    return coh


def write_npz(path, cohort):
    np.savez_compressed(path, sc=cohort.sc, fc=cohort.fc,
                        iu_i=cohort.iu, iv_j=cohort.iv,
                        networks=cohort.networks, hierarchy=cohort.hierarchy,
                        coords=cohort.coords,
                        behaviour=cohort.behaviour, covariates=cohort.covariates,
                        n_nodes=cohort.n_nodes)
    return path
