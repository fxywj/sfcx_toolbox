# -*- coding: utf-8 -*-
"""Exercise the v0.2 surface: BOLD intake, custom coupling, engine bundle."""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import sfcx
from sfcx import demo, io, engine, custom, functional


def main():
    out = tempfile.mkdtemp()
    coh = demo.make_cohort(n_subjects=40, n_nodes=48, seed=5)

    # ---- BOLD -> FC round trip -----------------------------------------
    rng = np.random.default_rng(0)
    bold = rng.normal(size=(40, 48, 60)).astype(np.float32)
    fc_from_bold = functional.bold_to_fc(bold, kind="pearson")
    assert fc_from_bold.shape == (40, 48 * 47 // 2), fc_from_bold.shape
    assert np.all(np.isfinite(fc_from_bold))
    fc_prec = functional.bold_to_fc(bold, kind="precision")
    assert fc_prec.shape == fc_from_bold.shape
    print("PASS  BOLD -> FC (pearson, precision)")

    # ---- build_cohort with BOLD ----------------------------------------
    sc = coh.group_matrix("sc")[None, :, :].repeat(40, axis=0)
    built, report = io.build_cohort(sc=sc, bold=bold, networks=coh.networks,
                                    hierarchy=coh.hierarchy)
    assert built.n == 40 and built.n_nodes == 48
    print("PASS  build_cohort from SC matrices + BOLD")

    # ---- custom coupling ------------------------------------------------
    src = """
def couple(sc, fc, iu, iv, n_nodes):
    zs = zscore_edges(sc, iu, n_nodes)
    zf = zscore_edges(fc, iu, n_nodes)
    return zs * zf * (1.0 + 0.0 * iu)
"""
    name = custom.register_custom("demo_product", src, "standardised product")
    assert name in sfcx.METRICS
    kind, arr = custom.evaluate(custom.CUSTOM_METRICS[name]["fn"], coh)
    assert kind == "edge" and arr.shape == (40, 48 * 47 // 2)
    print("PASS  custom coupling registered and evaluated")

    bad = "def couple(sc, fc, iu, iv, n_nodes):\n    import os\n    return sc"
    try:
        custom.validate_source(bad)
    except custom.CustomCodeError as exc:
        print("PASS  sandbox rejects import: %s" % exc)
    else:
        raise AssertionError("sandbox should reject import os")

    bad2 = "def couple(sc, fc, iu, iv, n_nodes):\n    return open('/etc/passwd')"
    try:
        custom.validate_source(bad2)
    except custom.CustomCodeError as exc:
        print("PASS  sandbox rejects open(): %s" % exc)
    else:
        raise AssertionError("sandbox should reject open")

    dr = custom.dry_run(src)
    assert dr["kind"] == "edge" and dr["finite"]
    print("PASS  custom dry_run -> %s" % dr["shape"])

    # ---- engine bundle -------------------------------------------------
    bundle = engine.analyze(coh, ["communicability", "direct", "demo_product"],
                            n_perm=60, outdir=out, make_figures=True, seed=1)
    assert "metrics" in bundle and len(bundle["metrics"]) == 3
    ax = bundle["metrics"]["communicability"]["estimators"]["regional"]["axis"]
    assert np.isfinite(ax["r"])
    assert "controlled_for_magnitude" in ax
    assert bundle["sweep"]["r"][0][0] is not None
    assert "axis_summary.csv" in bundle["tables"]
    assert os.path.exists(os.path.join(out, "figures", "coordinate_matrix.png"))
    print("PASS  engine bundle: axis r=%.3f p=%.4f, figures+tables"
          % (ax["r"], ax["p"] or float("nan")))

    # ---- exports -------------------------------------------------------
    npz = io.export_npz(os.path.join(out, "cohort.npz"), coh)
    mat = io.export_mat(os.path.join(out, "cohort.mat"), coh)
    assert os.path.exists(npz) and os.path.exists(mat)
    reloaded, _ = io.load_cohort(npz)
    assert reloaded.n == coh.n and reloaded.n_nodes == coh.n_nodes
    print("PASS  npz/mat export and reload")

    print("\noutputs:", out)


if __name__ == "__main__":
    main()
