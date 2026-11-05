# -*- coding: utf-8 -*-
"""Smoke tests.  Run: pytest -q"""

import numpy as np
import pytest

import sfcx
from sfcx import Cohort, Coordinate, demo
from sfcx.coords import CoordinateError


@pytest.fixture(scope="module")
def cohort():
    return demo.make_cohort(n_subjects=60, n_nodes=60, seed=3)


def test_cohort_shapes(cohort):
    assert cohort.n == 60
    assert cohort.n_nodes == 60
    assert cohort.n_edges == 60 * 59 // 2
    reg = cohort.regional_coupling()
    assert reg.shape == (60, 60)


def test_coordinate_validation():
    Coordinate(metric="direct", estimator="edge")
    with pytest.raises(CoordinateError):
        Coordinate(estimator="region")
    with pytest.raises(CoordinateError):
        Coordinate(metric="nope")


def test_presets_and_with():
    c = Coordinate.from_preset("default")
    c2 = c.with_(estimator="edge")
    assert c2.estimator == "edge"
    assert c.estimator == "regional"


def test_route_runs(cohort):
    res = sfcx.run(cohort, Coordinate(metric="direct"), outdir="test_out",
                   n_perm=50, make_figures=False)
    assert "axis" in res
    assert np.isfinite(res["axis"]["r"])
    assert res["variability"]["n"] == cohort.n_nodes


def test_sweep_sign_spread(cohort):
    sw = sfcx.sweep(cohort, outdir="test_sweep", n_perm=40, make_figures=False)
    r = sw["r"]
    assert r.shape == (len(sw["rows"]), len(sw["cols"]))
    assert np.nanmin(r) < np.nanmax(r)


def test_layers_positive(cohort):
    from sfcx import analysis, metrics
    coup = metrics.edge_coupling(cohort, "direct")
    lay = analysis.layers(np.nanstd(coup, axis=0), cohort, n_perm=100)
    assert lay["within"] > lay["between"]
    assert 0 <= lay["p"] <= 1


def test_axis_edge_level(cohort):
    from sfcx import analysis, metrics
    coup = metrics.edge_coupling(cohort, "direct")
    v = analysis.variability(coup, cohort, "edge")
    hier = 0.5 * (cohort.hierarchy[cohort.iu] + cohort.hierarchy[cohort.iv])
    st = analysis.axis(v, hier, None)
    assert np.isfinite(st["r"])


def test_report_writes(cohort, tmp_path):
    res = sfcx.run(cohort, Coordinate(metric="direct"), outdir=str(tmp_path),
                   n_perm=30, make_figures=True)
    html = sfcx.build_report(res, str(tmp_path / "r.html"))
    text = open(html, encoding="utf-8").read()
    assert "<html" in text
    assert "base64" in text
