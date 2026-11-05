# -*- coding: utf-8 -*-
"""Route execution: one coordinate, one result bundle."""

from __future__ import annotations

import json
import os
import time
import numpy as np

from . import analysis, metrics, viz
from .coords import Coordinate, METRICS


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _null_engine(cohort, kind, n_perm, seed):
    if kind == "auto":
        kind = "spin" if cohort.coords is not None else "msr"
    return analysis.NullEngine(cohort, kind=kind, n_perm=n_perm, seed=seed)


def _hierarchy_for(cohort, v):
    """Regional hierarchy for regional maps, endpoint-mean hierarchy for edges."""
    v = np.asarray(v)
    if v.size == cohort.n_nodes:
        return cohort.hierarchy
    return 0.5 * (cohort.hierarchy[cohort.iu] + cohort.hierarchy[cohort.iv])


def run(cohort, coord=None, outdir="sfcx_out", n_perm=1000, null="auto", seed=0,
        make_figures=True, behaviour_label="behaviour"):
    """Execute one route and write results + figures into `outdir`."""
    coord = coord or Coordinate()
    t0 = time.time()
    os.makedirs(outdir, exist_ok=True)
    figdir = os.path.join(outdir, "figures")
    tabdir = os.path.join(outdir, "tables")
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(tabdir, exist_ok=True)

    engine = _null_engine(cohort, null, n_perm, seed)
    kind, field = metrics.coupling(cohort, coord.metric)
    if kind == "regional":
        regional = field
        coupling = metrics.edge_projection(regional, cohort)
    else:
        coupling = field
        regional = metrics.regional_readout(cohort, coord.metric, coupling)
    cov = cohort.covariates

    if coord.estimator == "regional":
        v = analysis.regional_variability(regional, cov)
        level = "regional"
    elif coord.estimator == "node":
        v = analysis.variability(coupling, cohort, "node")
        level = "regional"
    else:
        v = analysis.variability(coupling, cohort, "edge")
        level = "edge"

    res = {
        "coordinate": coord.as_dict(),
        "data": {"n_subjects": cohort.n, "n_nodes": cohort.n_nodes,
                 "n_edges": cohort.n_edges, "state": cohort.state,
                 "n_networks": cohort.n_networks},
        "variability": {"level": level, "n": int(v.size),
                        "min": float(np.nanmin(v)), "max": float(np.nanmax(v)),
                        "median": float(np.nanmedian(v))},
        "null_model": engine.kind,
    }

    fig_paths = {}
    hier = _hierarchy_for(cohort, v)
    res["axis"] = analysis.axis(v, hier, engine, n_perm)
    res["axis_level"] = level
    if level == "regional":
        if make_figures:
            vv = np.sort(v[np.isfinite(v)])[::-1]
            fig_paths["spectrum"] = viz.spectrum(
                v, os.path.join(figdir, "spectrum.png"),
                slope=float(np.polyfit(np.arange(vv.size), vv, 1)[0]) if vv.size > 1 else None)
        epi = analysis.epicentres(v, cohort, seed=seed, n_perm=min(n_perm, 1000))
        res["epicentres"] = epi
        if make_figures:
            fig_paths["epicentres"] = viz.epicentre_strip(
                v, cohort.networks, epi["regions"], os.path.join(figdir, "epicentres.png"))
        if cov is not None and cov.size:
            ok = np.isfinite(regional)
            mag = np.divide(np.where(ok, regional, 0.0).sum(0), ok.sum(0),
                            out=np.zeros(regional.shape[1]), where=ok.sum(0) > 0)
            res["axis_controlled_for_magnitude"] = analysis.axis_controlled(
                v, cohort.hierarchy, mag, engine, n_perm)

    edge_v = v if level == "edge" else analysis.variability(coupling, cohort, "edge")
    res["layers"] = analysis.layers(edge_v, cohort, n_perm=min(n_perm, 1000), seed=seed)
    if make_figures:
        fig_paths["layers"] = viz.layers_bar(
            res["layers"]["within"], res["layers"]["between"],
            os.path.join(figdir, "layers.png"), p=res["layers"]["p"])

    res["reliability"] = analysis.reliability(coupling, cohort, seed)
    if make_figures:
        a, b = cohort.split_half(seed)
        va = np.nanstd(coupling[a], axis=0)
        vb = np.nanstd(coupling[b], axis=0)
        fig_paths["reliability"] = viz.reliability_scatter(
            va, vb, os.path.join(figdir, "reliability.png"),
            r=res["reliability"]["split_half_r"])

    if cohort.behaviour is not None:
        dev = analysis.deviation_spectrum(coupling, cohort,
                                          scales=("regional", "network", "edge", "gap"),
                                          seed=seed)
        indiv = {}
        for k, x in dev.items():
            st = analysis.partial_spearman(x, cohort.behaviour, cov)
            st["values"] = x
            indiv[k] = st
        res["individual"] = {k: {kk: vv for kk, vv in st.items() if kk != "values"}
                             for k, st in indiv.items()}
        best = min(indiv, key=lambda k: indiv[k]["p"])
        res["individual"]["best_scale"] = best
        if make_figures:
            st = indiv[best]
            fig_paths["individual"] = viz.individual_scatter(
                st["values"], cohort.behaviour, os.path.join(figdir, "individual.png"),
                r=st["r"], p=st["p"], xlabel="%s-scale slope" % best,
                ylabel=behaviour_label)

    if make_figures:
        card = ("%s\nmetric     : %s (%s)\npipeline   : %s\nestimator  : %s\nscale      : %s"
                "\nnull       : %s (%d draws)"
                % (coord.label, coord.metric, METRICS[coord.metric]["family"],
                   coord.pipeline, coord.estimator, coord.scale, engine.kind, n_perm))
        fig_paths["route_card"] = viz.route_card(card, os.path.join(figdir, "route_card.png"))

    res["figures"] = fig_paths
    res["runtime_s"] = round(time.time() - t0, 2)

    # tables
    if level == "regional":
        np.savetxt(os.path.join(tabdir, "variability_regional.csv"), v, delimiter=",",
                   header="dispersion", comments="")
    else:
        np.savetxt(os.path.join(tabdir, "variability_edge.csv"), v, delimiter=",",
                   header="dispersion", comments="")
    with open(os.path.join(outdir, "results.json"), "w", encoding="utf-8") as f:
        json.dump(_jsonable(res), f, indent=2)
    return res


def sweep(cohort, outdir="sfcx_sweep", metric_list=None, estimator_list=None,
          n_perm=500, null="auto", seed=0, make_figures=True):
    """The signature analysis: axis sign as a function of the coordinate.

    Returns a dict with the r / p matrices over metric x estimator.
    """
    from .coords import METRICS, fast_metrics
    metric_list = metric_list or fast_metrics()
    estimator_list = estimator_list or ["regional", "node", "edge"]
    os.makedirs(outdir, exist_ok=True)
    figdir = os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    engine = _null_engine(cohort, null, n_perm, seed)

    r = np.full((len(metric_list), len(estimator_list)), np.nan)
    p = np.full_like(r, np.nan)
    cache = {}
    t0 = time.time()
    for i, m in enumerate(metric_list):
        kind, field = metrics.coupling(cohort, m)
        if kind == "regional":
            regional = field
            coupling = metrics.edge_projection(regional, cohort)
        else:
            coupling = field
            regional = metrics.regional_readout(cohort, m, coupling)
        cache[m] = coupling
        for j, est in enumerate(estimator_list):
            if est == "regional":
                v = analysis.regional_variability(regional, cohort.covariates)
            else:
                v = analysis.variability(coupling, cohort, est)
            st = analysis.axis(v, _hierarchy_for(cohort, v), engine, n_perm)
            r[i, j] = st["r"]
            p[i, j] = st["p"]

    out = {
        "data": {"n_subjects": cohort.n, "n_nodes": cohort.n_nodes,
                 "n_edges": cohort.n_edges, "state": cohort.state},
        "rows": metric_list,
        "cols": estimator_list,
        "r": r, "p": p,
        "null_model": engine.kind,
        "runtime_s": round(time.time() - t0, 2),
    }
    if make_figures:
        out["figure"] = viz.coordinate_matrix(
            r, metric_list, estimator_list, os.path.join(figdir, "coordinate_matrix.png"),
            pmat=p,
            title="Axis strength across the coordinate system (%s null)" % engine.kind)

    np.savetxt(os.path.join(outdir, "axis_r.csv"), r, delimiter=",",
               header=",".join(estimator_list), comments="")
    np.savetxt(os.path.join(outdir, "axis_p.csv"), p, delimiter=",",
               header=",".join(estimator_list), comments="")
    with open(os.path.join(outdir, "sweep.json"), "w", encoding="utf-8") as f:
        json.dump({"meta": {k: v for k, v in out.items() if k not in ("r", "p")},
                   "r": r.tolist(), "p": p.tolist()}, f, indent=2)
    return out
