# -*- coding: utf-8 -*-
"""Unified analysis engine.

One entry point used by the command line, the web service and the MATLAB
bridge: give it a cohort and a small settings dict, get back a serialisable
bundle with coupling, variability, the variability axis, epicentres, layers,
reliability and individual deviation-spectrum slopes.
"""

from __future__ import annotations

import os
import time

import numpy as np

from . import analysis, metrics, viz
from .coords import METRICS, Coordinate


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    return o


def _null_engine(cohort, kind, n_perm, seed):
    if kind == "auto":
        kind = "spin" if cohort.coords is not None else "msr"
    return analysis.NullEngine(cohort, kind=kind, n_perm=n_perm, seed=seed)


def _hierarchy_for(cohort, v):
    v = np.asarray(v)
    if v.size == cohort.n_nodes:
        return cohort.hierarchy
    return 0.5 * (cohort.hierarchy[cohort.iu] + cohort.hierarchy[cohort.iv])


def _axis_block(cohort, v, engine, n_perm, control=None):
    hier = _hierarchy_for(cohort, v)
    res = analysis.axis(v, hier, engine, n_perm)
    if v.size == cohort.n_nodes and control is not None:
        res["controlled_for_magnitude"] = analysis.axis_controlled(
            v, cohort.hierarchy, control, engine, n_perm)
    res["magnitude"] = float(np.nanmean(v))
    return res


def _nanmean_cols(x):
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    total = np.where(ok, x, 0.0).sum(axis=0)
    count = ok.sum(axis=0)
    return np.divide(total, count, out=np.full(total.shape, np.nan), where=count > 0)


def _magnitude_map(cohort, estimator, regional, edge_field):
    if estimator == "regional":
        return _nanmean_cols(regional)
    if estimator == "node":
        out = np.full(cohort.n_nodes, np.nan)
        mean_edge = _nanmean_cols(edge_field)
        for i, idx in enumerate(cohort._region_profile):
            if idx.size:
                out[i] = np.nanmean(mean_edge[idx])
        return out
    return _nanmean_cols(edge_field)


def _variability_for(cohort, estimator, regional, edge_field):
    if estimator == "regional":
        return analysis.regional_variability(regional, cohort.covariates), "regional"
    if estimator == "node":
        return analysis.variability(edge_field, cohort, "node"), "regional"
    return analysis.variability(edge_field, cohort, "edge"), "edge"


def magnitude(coh, metrics_list, progress=None):
    """Coupling-value mode: per-subject coupling magnitudes, no dispersion.

    Works for a single subject (n = 1) and any cohort size.  Returns a
    bundle with, per metric: kind ('edge'/'regional'), the full field
    (n_subjects x n_edges or n_subjects x P), the regional read-out, the
    group-mean edge vector / matrix and subject ids.
    """
    t0 = time.time()
    metrics_list = list(metrics_list or ["communicability"])
    if isinstance(metrics_list, str):
        metrics_list = [metrics_list]

    def tick(frac, msg):
        if progress is not None:
            progress(float(np.clip(frac, 0, 1)), msg)

    iu = coh.iu
    iv = coh.iv
    p = coh.n_nodes
    subject_ids = getattr(coh, "subject_ids", None)
    if subject_ids is None:
        subject_ids = ["sub%03d" % (s + 1) for s in range(coh.n)]
    subject_ids = [str(s) for s in subject_ids]
    out = {
        "mode": "magnitude",
        "meta": {"n_subjects": coh.n, "n_nodes": p, "n_edges": coh.n_edges,
                 "state": coh.state},
        "metrics": {},
    }
    total = max(len(metrics_list), 1)
    for i, metric in enumerate(metrics_list):
        tick(i / total, "coupling: %s" % metric)
        kind, field = metrics.coupling(coh, metric)
        field = np.asarray(field, dtype=np.float64)
        entry = {"kind": kind, "subject_ids": list(subject_ids)}
        if kind == "edge":
            regional = np.full((coh.n, p), np.nan)
            for j, idx in enumerate(coh._region_profile):
                if idx.size:
                    regional[:, j] = np.nanmean(field[:, idx], axis=1)
            mean_vec = np.nanmean(field, axis=0)
            M = np.zeros((p, p))
            M[iu, iv] = mean_vec
            M = M + M.T
            entry["field"] = field
            entry["regional"] = regional
            entry["mean_matrix"] = M
        else:
            entry["field"] = None
            entry["regional"] = np.asarray(field, dtype=np.float64)
            entry["mean_matrix"] = None
        out["metrics"][metric] = entry
    out["runtime_s"] = round(time.time() - t0, 3)
    tick(1.0, "Done · 完成")
    return out


def export_coupling(bundle, out_dir, fmt=("csv", "npz")):
    """Write coupling values from :func:`magnitude` to disk.

    CSV: ``<metric>_coupling_edges.csv`` (subject_id + every edge) and
    ``<metric>_coupling_regional.csv`` (subject_id + every parcel).
    NPZ: fields, regional read-outs and mean matrices.
    Returns the list of files written.
    """
    os.makedirs(out_dir, exist_ok=True)
    written = []
    npz_payload = {}
    for name, m in bundle["metrics"].items():
        ids = m["subject_ids"]
        if "csv" in fmt:
            reg = np.asarray(m["regional"], dtype=np.float64)
            path = os.path.join(out_dir, "%s_coupling_regional.csv" % name)
            with open(path, "w", encoding="utf-8", newline="") as fh:
                fh.write("subject_id," + ",".join(
                    "P%03d" % k for k in range(reg.shape[1])) + "\n")
                for s, sid in enumerate(ids):
                    fh.write(sid + "," + ",".join("%.6g" % v for v in reg[s]) + "\n")
            written.append(path)
            if m.get("field") is not None:
                fld = np.asarray(m["field"], dtype=np.float64)
                path = os.path.join(out_dir, "%s_coupling_edges.csv" % name)
                with open(path, "w", encoding="utf-8", newline="") as fh:
                    fh.write("subject_id," + ",".join(
                        "e%d" % k for k in range(fld.shape[1])) + "\n")
                    for s, sid in enumerate(ids):
                        fh.write(sid + "," + ",".join("%.6g" % v for v in fld[s]) + "\n")
                written.append(path)
        npz_payload["%s_regional" % name] = np.asarray(m["regional"], dtype=np.float64)
        if m.get("field") is not None:
            npz_payload["%s_field" % name] = np.asarray(m["field"], dtype=np.float64)
        if m.get("mean_matrix") is not None:
            npz_payload["%s_mean_matrix" % name] = np.asarray(
                m["mean_matrix"], dtype=np.float64)
    if "npz" in fmt:
        path = os.path.join(out_dir, "coupling_values.npz")
        np.savez_compressed(path, **npz_payload)
        written.append(path)
    return written


def analyze(cohort, metrics_list=None, estimators=("regional", "node", "edge"),
            n_perm=500, null="auto", seed=0, outdir=None, make_figures=True,
            behaviour_label="behaviour", progress=None, do_sweep=True,
            state=None):
    """Run the full analysis battery and return a serialisable bundle."""
    t0 = time.time()
    if cohort.n < 2:
        raise ValueError(
            "变异分析需要 ≥2 名被试（当前 %d）：单个被试没有个体间差异，SD 在数学上未定义。"
            "请改用「耦合值计算」模式。 / Variability analysis needs ≥2 subjects "
            "(got %d): a single subject has no inter-individual differences, so the SD "
            "is undefined. Use the coupling-values mode instead." % (cohort.n, cohort.n))
    if cohort.n < 8:
        warnings_min = {"zh": "被试数 %d < 8：变异度、split-half 可靠性与个体信号估计在这一样本量下"
                              "噪声很大，结果仅作演示用途" % cohort.n,
                        "en": "n = %d < 8: dispersion, split-half reliability and individual-signal "
                              "estimates are very noisy at this sample size; treat results as a demo." % cohort.n}
    else:
        warnings_min = None
    metrics_list = list(metrics_list or ["communicability"])
    estimators = tuple(estimators)
    outdir = outdir or os.path.join(os.getcwd(), "sfcx_out")
    figdir = os.path.join(outdir, "figures")
    tabdir = os.path.join(outdir, "tables")
    if make_figures:
        os.makedirs(figdir, exist_ok=True)
    os.makedirs(tabdir, exist_ok=True)

    engine = _null_engine(cohort, null, n_perm, seed)
    warnings = list(getattr(cohort, "meta", {}).get("report", {}).get("warnings", []))

    def tick(frac, msg):
        if progress is not None:
            progress(float(np.clip(frac, 0, 1)), msg)

    total = max(len(metrics_list), 1)
    bundle = {
        "meta": {
            "n_subjects": cohort.n, "n_nodes": cohort.n_nodes,
            "n_edges": cohort.n_edges, "n_networks": cohort.n_networks,
            "state": state or cohort.state,
            "has_behaviour": cohort.behaviour is not None,
            "has_covariates": cohort.covariates is not None,
            "has_coords": cohort.coords is not None,
            "null_model": engine.kind, "n_perm": int(n_perm),
            "hierarchy_supplied": bool(getattr(cohort, "meta", {}).get(
                "hierarchy_provided", True)),
        },
        "settings": {
            "metrics": metrics_list, "estimators": list(estimators),
            "n_perm": int(n_perm), "null": engine.kind, "seed": int(seed),
            "behaviour_label": behaviour_label,
        },
        "metrics": {},
        "warnings": (warnings + [warnings_min]) if warnings_min else warnings,
    }

    sweep_r = np.full((len(metrics_list), len(estimators)), np.nan)
    sweep_p = np.full_like(sweep_r, np.nan)
    cache = {}

    for i, metric in enumerate(metrics_list):
        tick(i / total, "coupling: %s" % metric)
        kind, field = metrics.coupling(cohort, metric)
        if kind == "regional":
            regional = field
            edge_field = metrics.edge_projection(regional, cohort)
        else:
            regional = metrics.regional_readout(cohort, metric, field)
            edge_field = field
        cache[metric] = dict(kind=kind, edge=edge_field, regional=regional,
                             coupling=field)

        block = {
            "family": METRICS.get(metric, {}).get("family", "custom"),
            "label": METRICS.get(metric, {}).get("label", metric),
            "description": METRICS.get(metric, {}).get("description", ""),
            "field_kind": kind,
            "coupling_summary": {
                "mean": float(np.nanmean(field)), "sd": float(np.nanstd(field)),
                "min": float(np.nanmin(field)), "max": float(np.nanmax(field)),
                "nan_fraction": float(np.mean(~np.isfinite(field))),
            },
            "estimators": {},
        }

        primary_est = "regional" if "regional" in estimators else estimators[0]
        for est in estimators:
            v, level = _variability_for(cohort, est, regional, edge_field)
            control = _magnitude_map(cohort, est, regional, edge_field)
            ax = _axis_block(cohort, v, engine, n_perm, control=control)
            sweep_r[i, estimators.index(est)] = ax["r"]
            sweep_p[i, estimators.index(est)] = ax["p"]
            item = {
                "level": level, "n": int(v.size),
                "min": float(np.nanmin(v)), "max": float(np.nanmax(v)),
                "median": float(np.nanmedian(v)), "mean": float(np.nanmean(v)),
                "axis": ax, "values": v.tolist(),
            }
            if est == primary_est:
                block["primary_estimator"] = est
            block["estimators"][est] = item

        prim = block["estimators"][primary_est]
        prim_v = np.asarray(prim["values"], dtype=float)

        if prim["level"] == "regional":
            prim["epicentres"] = analysis.epicentres(
                prim_v, cohort, seed=seed, n_perm=min(n_perm, 1000))
            if make_figures:
                prim["figures"] = {
                    "spectrum": viz.spectrum(
                        prim_v, os.path.join(figdir, "%s_spectrum.png" % metric),
                        slope=_descending_slope(prim_v)),
                    "epicentres": viz.epicentre_strip(
                        prim_v, cohort.networks,
                        prim["epicentres"]["regions"],
                        os.path.join(figdir, "%s_epicentres.png" % metric)),
                }

        lay = analysis.layers(analysis.variability(edge_field, cohort, "edge"),
                              cohort, n_perm=min(n_perm, 1000), seed=seed)
        block["layers"] = lay
        block["reliability"] = analysis.reliability(edge_field, cohort, seed)

        dev = analysis.deviation_spectrum(
            edge_field, cohort, scales=("regional", "network", "edge", "gap"), seed=seed)
        indiv = {}
        for scale, values in dev.items():
            entry = {"summary": {"mean": float(np.nanmean(values)),
                                 "sd": float(np.nanstd(values))}}
            if cohort.behaviour is not None:
                entry.update(analysis.partial_spearman(
                    values, cohort.behaviour, cohort.covariates))
            indiv[scale] = entry
        if cohort.behaviour is not None and indiv:
            best = min(indiv, key=lambda k: (indiv[k].get("p", np.inf)
                                             if np.isfinite(indiv[k].get("p", np.inf)) else np.inf))
            block["individual"] = indiv
            block["individual"]["best_scale"] = best
        else:
            block["individual"] = indiv

        if make_figures:
            figs = dict(prim.get("figures", {}))
            figs["layers"] = viz.layers_bar(
                lay["within"], lay["between"],
                os.path.join(figdir, "%s_layers.png" % metric), p=lay["p"])
            a, b = cohort.split_half(seed)
            va = np.nanstd(edge_field[a], axis=0)
            vb = np.nanstd(edge_field[b], axis=0)
            figs["reliability"] = viz.reliability_scatter(
                va, vb, os.path.join(figdir, "%s_reliability.png" % metric),
                r=block["reliability"]["split_half_r"])
            if cohort.behaviour is not None and "best_scale" in block["individual"]:
                bs = block["individual"]["best_scale"]
                st = indiv[bs]
                if "r" in st:
                    figs["individual"] = viz.individual_scatter(
                        dev[bs], cohort.behaviour,
                        os.path.join(figdir, "%s_individual.png" % metric),
                        r=st.get("r"), p=st.get("p"),
                        xlabel="%s-scale slope" % bs, ylabel=behaviour_label)
            card = ("%s\nfamily     : %s\nestimator  : %s\nnull       : %s (%d draws)"
                    % (metric, block["family"], primary_est, engine.kind, n_perm))
            figs["route_card"] = viz.route_card(
                card, os.path.join(figdir, "%s_card.png" % metric))
            block["figures"] = figs

        np.savetxt(os.path.join(tabdir, "%s_variability.csv" % metric), prim_v,
                   delimiter=",", header="dispersion", comments="")
        bundle["metrics"][metric] = block
        tick((i + 1) / total, "done: %s" % metric)

    if do_sweep:
        tick(0.98, "coordinate sweep")
        bundle["sweep"] = {
            "rows": metrics_list, "cols": list(estimators),
            "r": sweep_r.tolist(), "p": sweep_p.tolist(),
        }
        if make_figures:
            bundle["sweep"]["figure"] = viz.coordinate_matrix(
                sweep_r, metrics_list, list(estimators),
                os.path.join(figdir, "coordinate_matrix.png"), pmat=sweep_p,
                title="Variability axis across the coordinate system (%s null)" % engine.kind)
        np.savetxt(os.path.join(tabdir, "axis_r.csv"), sweep_r, delimiter=",",
                   header=",".join(estimators), comments="")
        np.savetxt(os.path.join(tabdir, "axis_p.csv"), sweep_p, delimiter=",",
                   header=",".join(estimators), comments="")

    bundle["tables"] = _tables(bundle, tabdir)
    bundle["runtime_s"] = round(time.time() - t0, 2)
    return _jsonable(bundle)


def _descending_slope(v):
    v = np.asarray(v, dtype=float)
    ok = np.isfinite(v)
    if ok.sum() < 2:
        return None
    vv = np.sort(v[ok])[::-1]
    return float(np.polyfit(np.arange(vv.size), vv, 1)[0])


def _tables(bundle, tabdir):
    out = {}
    for name in sorted(os.listdir(tabdir)):
        if not name.endswith(".csv"):
            continue
        with open(os.path.join(tabdir, name), "r", encoding="utf-8") as f:
            out[name] = f.read()
    rows = ["metric,estimator,axis_r,axis_p"]
    sw = bundle.get("sweep")
    if sw:
        for i, m in enumerate(sw["rows"]):
            for j, e in enumerate(sw["cols"]):
                r = sw["r"][i][j]
                p = sw["p"][i][j]
                rows.append("%s,%s,%s,%s" % (m, e,
                                             "" if r is None else "%.6f" % r,
                                             "" if p is None else "%.6g" % p))
    out["axis_summary.csv"] = "\n".join(rows) + "\n"
    return out


def run_route(cohort, coord=None, outdir="sfcx_out", n_perm=500, null="auto",
              seed=0, make_figures=True, behaviour_label="behaviour"):
    """Single-coordinate convenience wrapper returning one metric block."""
    coord = coord or Coordinate()
    bundle = analyze(cohort, [coord.metric], estimators=(coord.estimator,),
                     n_perm=n_perm, null=null, seed=seed, outdir=outdir,
                     make_figures=make_figures, behaviour_label=behaviour_label,
                     do_sweep=False)
    block = bundle["metrics"][coord.metric]["estimators"][coord.estimator]
    return {
        "coordinate": coord.as_dict(),
        "data": bundle["meta"],
        "null_model": bundle["meta"]["null_model"],
        "variability": {k: block[k] for k in ("level", "n", "min", "max", "median")},
        "axis": block["axis"],
        "layers": bundle["metrics"][coord.metric]["layers"],
        "reliability": bundle["metrics"][coord.metric]["reliability"],
        "epicentres": block.get("epicentres"),
        "individual": bundle["metrics"][coord.metric].get("individual"),
        "figures": bundle["metrics"][coord.metric].get("figures", {}),
        "runtime_s": bundle["runtime_s"],
        "warnings": bundle["warnings"],
    }
