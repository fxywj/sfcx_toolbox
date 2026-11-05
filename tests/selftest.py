# -*- coding: utf-8 -*-
"""Dependency-free self test.  Run: python tests/selftest.py"""

import os
import sys
import tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import sfcx
from sfcx import Cohort, Coordinate, demo, analysis, metrics
from sfcx.coords import CoordinateError

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


@check("cohort shapes")
def _c1():
    c = demo.make_cohort(n_subjects=60, n_nodes=60, seed=3)
    assert c.n == 60 and c.n_nodes == 60 and c.n_edges == 60 * 59 // 2
    assert c.regional_coupling().shape == (60, 60)


@check("coordinate validation")
def _c2():
    Coordinate(metric="direct", estimator="edge")
    for bad in ({"estimator": "region"}, {"metric": "nope"}, {"scale": "voxel"}):
        try:
            Coordinate(**bad)
        except CoordinateError:
            pass
        else:
            raise AssertionError("should have raised for %r" % bad)


@check("presets are immutable")
def _c3():
    c = Coordinate.from_preset("default")
    c2 = c.with_(estimator="edge")
    assert c2.estimator == "edge" and c.estimator == "regional"


@check("route runs and returns axis")
def _c4():
    c = demo.make_cohort(n_subjects=60, n_nodes=60, seed=3)
    res = sfcx.run(c, Coordinate(metric="direct"), outdir=tempfile.mkdtemp(),
                   n_perm=40, make_figures=False)
    assert np.isfinite(res["axis"]["r"])
    assert res["variability"]["n"] == c.n_nodes
    assert 0 <= res["axis"]["p"] <= 1


@check("layers: within > between")
def _c5():
    c = demo.make_cohort(n_subjects=80, n_nodes=60, seed=3)
    coup = metrics.edge_coupling(c, "direct")
    lay = analysis.layers(np.nanstd(coup, axis=0), c, n_perm=80)
    assert lay["within"] > lay["between"], lay


@check("edge axis is computable")
def _c6():
    c = demo.make_cohort(n_subjects=60, n_nodes=60, seed=3)
    coup = metrics.edge_coupling(c, "direct")
    v = analysis.variability(coup, c, "edge")
    hier = 0.5 * (c.hierarchy[c.iu] + c.hierarchy[c.iv])
    assert np.isfinite(analysis.axis(v, hier, None)["r"])


@check("sweep spans a range of r")
def _c7():
    c = demo.make_cohort(n_subjects=80, n_nodes=60, seed=3)
    sw = sfcx.sweep(c, outdir=tempfile.mkdtemp(), n_perm=30, make_figures=False)
    r = sw["r"]
    assert r.shape == (len(sw["rows"]), len(sw["cols"]))
    assert np.nanmax(r) - np.nanmin(r) > 0.1


@check("report is written with figures")
def _c8():
    c = demo.make_cohort(n_subjects=60, n_nodes=60, seed=3)
    d = tempfile.mkdtemp()
    res = sfcx.run(c, Coordinate(metric="direct"), outdir=d, n_perm=30, make_figures=True)
    html = sfcx.build_report(res, os.path.join(d, "r.html"))
    txt = open(html, encoding="utf-8").read()
    assert "<html" in txt and "base64" in txt


def main():
    ok = 0
    for name, fn in CHECKS:
        try:
            fn()
            print("PASS  %s" % name)
            ok += 1
        except Exception as e:
            print("FAIL  %s  ->  %r" % (name, e))
    print("\n%d/%d passed" % (ok, len(CHECKS)))
    return 0 if ok == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
