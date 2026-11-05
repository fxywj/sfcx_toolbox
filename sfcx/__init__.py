# -*- coding: utf-8 -*-
"""SFCX - structure-function coupling explorer."""

from .coords import (Coordinate, METRICS, PIPELINES, ESTIMATORS, SCALES,
                     PRESETS, describe, catalog, register_metric,
                     fast_metrics, SLOW_METRICS)  # noqa: F401
from .data import Cohort, upper_tri_index, vec_to_mat, mat_to_vec  # noqa: F401
from .pipeline import run, sweep  # noqa: F401
from .demo import make_cohort, write_npz  # noqa: F401
from .report import build as build_report  # noqa: F401
from .functional import bold_to_fc, bold_to_matrices  # noqa: F401
from .custom import (register_custom, unregister_custom, list_custom,
                     validate_source, dry_run, CustomCodeError)  # noqa: F401
from .io import (build_cohort, load_cohort, load_any, export_npz, export_mat,
                 to_edge_vectors)  # noqa: F401
from .engine import analyze, magnitude, export_coupling, run_route  # noqa: F401

__version__ = "0.2.0"
