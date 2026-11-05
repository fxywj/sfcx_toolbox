# -*- coding: utf-8 -*-
"""Write reference.mat so the MATLAB toolbox can be checked against Python.

    python matlab/tests/make_reference.py

The cohort is the Python demo cohort; the expected values are the axis r of
the definition x estimator matrix.  r is deterministic for every metric, so
cells can be compared directly -- with two documented exceptions:

  * search_info depends on which shortest path is chosen among ties, and
    scipy and MATLAB break ties differently; it is excluded from parity.
  * reliance differs at the 1e-2 level through floating-point accumulation
    in the per-region profile regression; it gets a looser tolerance.

p values are never compared: the two runtimes use different random streams.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), ".."))

import sfcx  # noqa: E402
from sfcx import demo  # noqa: E402
from scipy.io import savemat  # noqa: E402

METRICS = ["direct", "rank_direct", "cosine", "matching", "communicability",
           "degree_product", "mfpt", "path_ease", "diffusion_efficiency",
           "search_info", "flow", "reliance"]

PARITY_EXCLUDE = {"search_info"}
PARITY_TOL = {m: 1e-5 for m in METRICS}
PARITY_TOL["reliance"] = 5e-2


def main():
    coh = demo.make_cohort(n_subjects=60, n_nodes=40, seed=0, behaviour=False)
    bundle = sfcx.engine.analyze(
        coh, METRICS,
        estimators=("regional", "node", "edge"),
        n_perm=50, null="auto", seed=0, make_figures=False, do_sweep=True,
        outdir=os.path.join(HERE, "_reference_out"))
    expected = np.asarray(bundle["sweep"]["r"], dtype=np.float64)
    strict = np.array([[m not in PARITY_EXCLUDE] * 3 for m in METRICS], dtype=bool)
    tol = np.array([[PARITY_TOL[m]] * 3 for m in METRICS], dtype=np.float64)
    names = np.empty((len(METRICS), 1), dtype=object)   # -> MATLAB cell array
    for i, m in enumerate(METRICS):
        names[i, 0] = m
    savemat(os.path.join(HERE, "reference.mat"), dict(
        sc=coh.sc, fc=coh.fc, networks=coh.networks, hierarchy=coh.hierarchy,
        coords=coh.coords, covariates=coh.covariates,
        expected_r=expected, metric_names=names,
        strict=strict, tol=tol))
    print("expected axis r (rows: %s; cols: regional, node, edge)" % ", ".join(METRICS))
    print(np.round(expected, 4))
    print("strict parity cells: %d/%d (search_info excluded, reliance tol 5e-2)"
          % (strict.sum(), strict.size))


if __name__ == "__main__":
    main()
