# -*- coding: utf-8 -*-
"""Quickstart: one dataset, four coordinates, four different answers.

Run:  python examples/quickstart.py
"""

import os
import numpy as np

import sfcx
from sfcx import Cohort, Coordinate, demo

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "quickstart_out")

# ------------------------------------------------------------------ data ----
coh = demo.make_cohort(n_subjects=200, n_nodes=120, seed=1)
print(coh)

# ------------------------------------------------- one coordinate at a time --
routes = {
    "default (regional)": Coordinate.from_preset("default"),
    "node contrast": Coordinate.from_preset("node-contrast"),
    "edge signal": Coordinate.from_preset("edge-signal"),
    "direct definition": Coordinate(metric="direct", estimator="regional"),
}

print("\n%-22s %8s %8s   %s" % ("coordinate", "axis r", "spin p", "null"))
print("-" * 62)
for name, coord in routes.items():
    res = sfcx.run(coh, coord, outdir=os.path.join(OUT, coord.slug),
                   n_perm=200, make_figures=False)
    ax = res.get("axis", {})
    print("%-22s %8.3f %8.3f   %s" % (name, ax.get("r", np.nan),
                                      ax.get("p", np.nan), res.get("null_model")))

# ---------------------------------------------------------------- the sweep --
sw = sfcx.sweep(coh, outdir=os.path.join(OUT, "sweep"), n_perm=200)
print("\nsweep: axis r over %d definitions x %d estimators -> %.2f .. %.2f"
      % (len(sw["rows"]), len(sw["cols"]), float(sw["r"].min()), float(sw["r"].max())))

# ------------------------------------------------------------- full report --
res = sfcx.run(coh, Coordinate.from_preset("default"),
               outdir=os.path.join(OUT, "report"), n_perm=300)
html = sfcx.build_report(res, os.path.join(OUT, "report", "report.html"),
                         title="SFCX quickstart", sweep=sw)
print("\nreport -> %s" % os.path.abspath(html))
