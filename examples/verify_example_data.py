# -*- coding: utf-8 -*-
"""Check the 5-subject teaching dataset and write its reference numbers.

    python examples/verify_example_data.py

Runs the same pipeline a learner would run, asserts that the data are well
formed, and writes ``examples/data/demo5_expected.json`` so that anyone can
confirm their import produced the same numbers.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
sys.path.insert(0, os.path.dirname(HERE))

import sfcx  # noqa: E402
from sfcx import io  # noqa: E402

METRICS = ("communicability", "direct", "reliance")
ESTIMATORS = ("regional", "node", "edge")
N_PERM = 500


def main():
    z = np.load(os.path.join(DATA, "demo5_group.npz"), allow_pickle=True)
    print("bundle keys :", ", ".join(sorted(z.files)))
    print("shapes      : sc %s  fc %s  bold %s"
          % (z["sc"].shape, z["fc"].shape, z["bold"].shape))

    coh, _ = io.build_cohort(sc=z["sc"], fc=z["fc"], networks=z["networks"],
                             hierarchy=z["hierarchy"], coords=z["coords"],
                             behaviour=z["behaviour"], covariates=z["covariates"],
                             subject_ids=z["subject_ids"])
    print("cohort      :", coh)

    # 1. BOLD must reproduce the stored FC exactly
    coh_b, _ = io.build_cohort(sc=z["sc"], bold=z["bold"], networks=z["networks"],
                              hierarchy=z["hierarchy"], coords=z["coords"])
    delta = float(np.abs(coh_b.fc - coh.fc).max())
    print("BOLD -> FC  : max |difference| = %.3g" % delta)
    assert delta < 1e-5, "BOLD and stored FC disagree"

    # 2. matrices must be well formed
    sc, fc = np.asarray(z["sc"]), np.asarray(z["fc"])
    assert np.allclose(sc, sc.transpose(0, 2, 1), atol=1e-6)
    assert np.allclose(fc, fc.transpose(0, 2, 1), atol=1e-6)
    assert np.allclose(np.diagonal(sc, axis1=1, axis2=2), 0.0, atol=1e-6)
    density = 100.0 * (sc[0] > 0).sum() / (sc.shape[1] * (sc.shape[1] - 1))
    print("matrices    : symmetric, zero diagonal, SC density %.1f %%" % density)

    # 3. the analysis a learner would run
    bundle = sfcx.analyze(coh, list(METRICS), estimators=ESTIMATORS,
                          n_perm=N_PERM, seed=0, make_figures=False,
                          outdir=os.path.join(HERE, "_verify_out"))

    print("\n%-16s %-9s %8s %8s %10s" % ("metric", "estimator", "axis r", "spin p", "variability"))
    print("-" * 58)
    for metric in METRICS:
        block = bundle["metrics"][metric]
        for est in ESTIMATORS:
            item = block["estimators"][est]
            ax = item["axis"]
            print("%-16s %-9s %+8.3f %8.3f %10.3f"
                  % (metric, est, ax["r"], ax["p"], item["mean"]))

    print("\nreliability (split-half r): %s"
          % {m: round(float(bundle["metrics"][m]["reliability"]["split_half_r"]), 3)
             for m in METRICS})

    expected = {
        "dataset": {
            "n_subjects": int(coh.n),
            "n_nodes": int(coh.n_nodes),
            "n_edges": int(coh.n_edges),
            "n_time": int(np.asarray(z["bold"]).shape[-1]),
            "tr": float(z["tr"]),
            "sc_density_percent": round(density, 2),
            "bold_vs_fc_max_abs_diff": delta,
        },
        "settings": {"metrics": list(METRICS), "estimators": list(ESTIMATORS),
                     "n_perm": N_PERM, "null": bundle["meta"]["null_model"], "seed": 0},
        "axis": {
            m: {e: {"r": round(float(bundle["metrics"][m]["estimators"][e]["axis"]["r"]), 4),
                    "p": round(float(bundle["metrics"][m]["estimators"][e]["axis"]["p"]), 4)}
                for e in ESTIMATORS}
            for m in METRICS
        },
        "note": ("Synthetic teaching data. Numbers are reproducible with seed 0 "
                 "and %d permutations; five subjects make every p-value a "
                 "demonstration of the machinery, not a result." % N_PERM),
    }
    path = os.path.join(DATA, "demo5_expected.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(expected, f, indent=2, ensure_ascii=False)
    print("\nreference numbers -> %s" % path)


if __name__ == "__main__":
    main()
