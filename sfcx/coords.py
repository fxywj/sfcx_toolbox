# -*- coding: utf-8 -*-
"""The readout coordinate system.

An SFC result is not a property of the data alone: it is a property of
(data, coordinate).  The coordinate has four layers, and any claim is only
meaningful when all four are stated.

    definition : which structural-functional correspondence is measured
    pipeline   : which functional representation is used
    estimator  : of what the dispersion is taken
    scale      : at which granularity the map is summarised
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

# ---------------------------------------------------------------- layers ----

METRICS = {
    # -- 9 static communication (edge) --
    "direct": dict(
        family="static-communication",
        label="Direct correspondence",
        label_zh="直接对应",
        description="Edge-wise product of profile-standardised SC and FC values.",
        needs="edges",
    ),
    "rank_direct": dict(
        family="static-communication",
        label="Rank-product correspondence",
        label_zh="秩乘积对应",
        description="Direct correspondence on per-subject edge ranks (Spearman-type).",
        needs="edges",
    ),
    "cosine": dict(
        family="static-communication",
        label="Cosine similarity",
        label_zh="余弦相似度",
        description="Cosine between the binarised structural profiles of two regions.",
        needs="matrix",
    ),
    "matching": dict(
        family="static-communication",
        label="Matching index",
        label_zh="匹配指数",
        description="Shared-neighbour ratio, (A^2)_ij / (k_i + k_j - 2 A_ij).",
        needs="matrix",
    ),
    "communicability": dict(
        family="static-communication",
        label="Communicability",
        label_zh="通信性",
        description="(expm(A))_ij, weighted walks of every length (Estrada–Hatano).",
        needs="matrix",
    ),
    "degree_product": dict(
        family="static-communication",
        label="Degree product",
        label_zh="度乘积",
        description="Product of binary degrees, k_i * k_j (homophily baseline).",
        needs="matrix",
    ),
    "mfpt": dict(
        family="static-communication",
        label="Mean first-passage time",
        label_zh="平均首次通过时间",
        description="Inverse-normalised random-walk hitting time.",
        needs="matrix",
    ),
    "path_ease": dict(
        family="static-communication",
        label="Path-length ease",
        label_zh="路径长度易度",
        description="1 / (1 + shortest-path length) on weighted shortest paths.",
        needs="matrix",
    ),
    "diffusion_efficiency": dict(
        family="static-communication",
        label="Diffusion efficiency",
        label_zh="扩散效率",
        description="1 / L_ij, the inverse weighted shortest-path length.",
        needs="matrix",
    ),
    "search_info": dict(
        family="static-communication",
        label="Search information",
        label_zh="搜索信息",
        description="Cost of the shortest path under the random-walk prior.",
        needs="matrix",
    ),
    "flow": dict(
        family="static-communication",
        label="Two-step flow",
        label_zh="两步流",
        description="(D^-1 A)^2, two-step random-walk reachability, symmetrised.",
        needs="matrix",
    ),
    # -- 4 reliance (regional, BOLD) --
    "reliance": dict(
        family="profile-regression",
        label="FC-on-SC profile slope",
        label_zh="FC对SC剖面斜率",
        description="Per-subject OLS slope of a region's FC profile on its SC profile, "
                    "re-expressed at edge level (static; not the BOLD-based reliance "
                    "of the main analysis).",
        needs="edges",
    ),
    "reliance_raw": dict(
        family="profile-regression",
        label="Reliance (raw)",
        label_zh="结构依赖（原始）",
        description="corr(ts_{t-1} W, ts_t) per region; structural reconstruction of activity.",
        needs="bold",
    ),
    "reliance_deriv": dict(
        family="profile-regression",
        label="Reliance (derivative)",
        label_zh="结构依赖（导数）",
        description="corr(-ts + ts @ W, dts/dt) per region; structural reconstruction of change.",
        needs="bold",
    ),
    "reliance_raw_nonlinear": dict(
        family="profile-regression",
        label="Reliance (nonlinear raw)",
        label_zh="结构依赖（非线性原始）",
        description="corr((1-a)ts + a*tanh(ts@W), ts_t) with leak + tanh activation.",
        needs="bold",
    ),
    "reliance_deriv_nonlinear": dict(
        family="profile-regression",
        label="Reliance (nonlinear derivative)",
        label_zh="结构依赖（非线性导数）",
        description="corr(a*(tanh(ts@W)-ts), dts/dt) with leak + tanh activation.",
        needs="bold",
    ),
    # -- 5 hybrid (regional, FC + SC) --
    "hybrid_global": dict(
        family="hybrid-model",
        label="Hybrid (global)",
        label_zh="混合模型（全局）",
        description="One set of communication weights across all regions, optimised per subject.",
        needs="matrix+fc",
    ),
    "hybrid_per_region": dict(
        family="hybrid-model",
        label="Hybrid (per-region)",
        label_zh="混合模型（逐区域）",
        description="One weight vector per region, optimised per subject.",
        needs="matrix+fc",
    ),
    "hybrid_cv": dict(
        family="hybrid-model",
        label="Hybrid (cross-validated)",
        label_zh="混合模型（交叉验证）",
        description="5-fold cross-validated global hybrid weights.",
        needs="matrix+fc",
    ),
    "hybrid_ridge": dict(
        family="hybrid-model",
        label="Hybrid (ridge)",
        label_zh="混合模型（岭回归）",
        description="Ridge regression of FC on communication profiles per region.",
        needs="matrix+fc",
    ),
    "hybrid_network": dict(
        family="hybrid-model",
        label="Hybrid (network)",
        label_zh="混合模型（网络约束）",
        description="Network-constrained hybrid with K-means clusters.",
        needs="matrix+fc",
    ),
    # -- 4 frontier (regional, BOLD + SC) --
    "information_geometry": dict(
        family="frontier",
        label="Information geometry",
        label_zh="信息几何",
        description="Log-Euclidean geodesic coupling on the SPD manifold.",
        needs="bold",
    ),
    "multi_scale_transfer_entropy": dict(
        family="frontier",
        label="Multi-scale transfer entropy",
        label_zh="多尺度转移熵",
        description="Conditional mutual information SC->FC across temporal scales.",
        needs="bold",
    ),
    "criticality": dict(
        family="frontier",
        label="Criticality index",
        label_zh="临界性指数",
        description="Distance of SC eigenvalue spectrum from lambda=1, modulating per-region coupling.",
        needs="matrix+fc",
    ),
    "granger": dict(
        family="frontier",
        label="Granger causality SC->FC",
        label_zh="格兰杰因果 SC->FC",
        description="F-statistic for SC-weighted past predicting BOLD future beyond autoregression.",
        needs="bold",
    ),
    # -- 2 innovative (regional, BOLD + SC) --
    "rnc": dict(
        family="innovative",
        label="Reliance-Null Coupling",
        label_zh="依赖-零模型耦合",
        description="Reliance minus weight-shuffle null; isolates topology-specific structural constraint.",
        needs="bold",
    ),
    "trr": dict(
        family="innovative",
        label="Temporal-Reach Reliance",
        label_zh="时间可达依赖",
        description="Weighted sum of multi-lag reliance (1/k weights); captures temporal persistence.",
        needs="bold",
    ),
}

PIPELINES = {
    "pearson": "Full correlation between regional BOLD series.",
    "precision": "Partial correlation after Ledoit-Wolf shrinkage.",
}

ESTIMATORS = {
    "regional": "SD across subjects of regional coupling c_i (one number per region).",
    "node": "Mean over the edges a region carries of edge-level SD (node projection).",
    "edge": "SD across subjects of a single edge (one number per edge).",
}

SCALES = {
    "regional": "Regional summary (n_regions units).",
    "network": "Network summary (n_networks units).",
    "edge": "Single edges (n_edges units).",
    "gap": "Within-minus-between gap (one scalar per subject).",
}

# --------------------------------------------------------------- presets ----

PRESETS = {
    "default": dict(
        metric="communicability", pipeline="pearson", estimator="regional", scale="regional",
        note="The paper's default: a static communication measure read at the regional scale.",
    ),
    "edge-signal": dict(
        metric="communicability", pipeline="pearson", estimator="edge", scale="edge",
        note="Where individual differences are readable: edge-level deviation-spectrum slope.",
    ),
    "node-contrast": dict(
        metric="communicability", pipeline="pearson", estimator="node", scale="regional",
        note="Same data, node projection: the axis is expected to flip sign.",
    ),
    "precision-fingerprint": dict(
        metric="direct", pipeline="precision", estimator="regional", scale="regional",
        note="Precision pipeline: individual identity rather than group prediction.",
    ),
}


class CoordinateError(ValueError):
    """Raised when a coordinate cannot be resolved."""


@dataclass(frozen=True)
class Coordinate:
    """One point in the four-layer readout coordinate system."""

    metric: str = "communicability"
    pipeline: str = "pearson"
    estimator: str = "regional"
    scale: str = "regional"

    def __post_init__(self):
        if self.metric not in METRICS:
            raise CoordinateError(
                "unknown metric %r; available: %s" % (self.metric, ", ".join(sorted(METRICS)))
            )
        if self.pipeline not in PIPELINES:
            raise CoordinateError(
                "unknown pipeline %r; available: %s" % (self.pipeline, ", ".join(sorted(PIPELINES)))
            )
        if self.estimator not in ESTIMATORS:
            raise CoordinateError(
                "unknown estimator %r; available: %s" % (self.estimator, ", ".join(sorted(ESTIMATORS)))
            )
        if self.scale not in SCALES:
            raise CoordinateError(
                "unknown scale %r; available: %s" % (self.scale, ", ".join(sorted(SCALES)))
            )

    @property
    def label(self) -> str:
        return "%s | %s | %s-level %s | %s" % (
            self.metric, self.pipeline, self.estimator, "SD", self.scale)

    @property
    def slug(self) -> str:
        return "%s__%s__%s__%s" % (self.metric, self.pipeline, self.estimator, self.scale)

    def as_dict(self) -> dict:
        d = asdict(self)
        d["label"] = self.label
        d["family"] = METRICS[self.metric]["family"]
        d["metric_description"] = METRICS[self.metric]["description"]
        return d

    @classmethod
    def from_preset(cls, name: str) -> "Coordinate":
        if name not in PRESETS:
            raise CoordinateError(
                "unknown preset %r; available: %s" % (name, ", ".join(sorted(PRESETS)))
            )
        p = dict(PRESETS[name])
        p.pop("note", None)
        return cls(**p)

    def with_(self, **kw) -> "Coordinate":
        """Return a copy with some layers replaced."""
        d = asdict(self)
        unknown = set(kw) - set(d)
        if unknown:
            raise CoordinateError("unknown layer(s): %s" % ", ".join(sorted(unknown)))
        d.update(kw)
        return Coordinate(**d)


def describe() -> str:
    """Human-readable description of the coordinate system."""
    out = ["SFC readout coordinate system", ""]
    out.append("Layer 1 - definition")
    for k, v in METRICS.items():
        out.append("  %-24s %-24s %s" % (k, v["family"], v["description"]))
    out.append("")
    out.append("Layer 2 - functional pipeline")
    for k, v in PIPELINES.items():
        out.append("  %-16s %s" % (k, v))
    out.append("")
    out.append("Layer 3 - dispersion estimator")
    for k, v in ESTIMATORS.items():
        out.append("  %-16s %s" % (k, v))
    out.append("")
    out.append("Layer 4 - observation scale")
    for k, v in SCALES.items():
        out.append("  %-16s %s" % (k, v))
    out.append("")
    out.append("Presets")
    for k, v in PRESETS.items():
        out.append("  %-22s %s" % (k, v["note"]))
    return "\n".join(out)


# ------------------------------------------------------- custom metrics ----

def register_metric(name, family="custom", label=None, label_zh=None, description="", needs="edges"):
    """Make a custom definition a first-class member of layer 1."""
    METRICS[name] = dict(family=family, label=label or name,
                         label_zh=label_zh or label or name,
                         description=description, needs=needs)
    return name


def unregister_metric(name):
    METRICS.pop(name, None)


# Metrics whose per-subject cost is high (SLSQP fits, matrix logarithms,
# time-series models).  They are fully supported, but they are excluded from
# the default sweep so that "run everything" never turns into an hour of
# compute on a large cohort.
SLOW_METRICS = {
    "reliance_raw", "reliance_deriv", "reliance_raw_nonlinear", "reliance_deriv_nonlinear",
    "hybrid_global", "hybrid_per_region", "hybrid_cv", "hybrid_ridge", "hybrid_network",
    "information_geometry", "multi_scale_transfer_entropy", "criticality", "granger",
    "rnc", "trr",
}


def fast_metrics():
    """Definitions cheap enough for a default coordinate sweep."""
    return [m for m in METRICS if m not in SLOW_METRICS]


def catalog() -> dict:
    """Everything a front-end needs to populate its controls."""
    return {
        "metrics": [dict(id=k, family=v["family"], label=v["label"],
                         label_zh=v.get("label_zh", v["label"]),
                         description=v["description"], needs=v["needs"],
                         cost=("slow" if k in SLOW_METRICS else "fast"))
                    for k, v in METRICS.items()],
        "pipelines": [dict(id=k, description=v) for k, v in PIPELINES.items()],
        "estimators": [dict(id=k, description=v) for k, v in ESTIMATORS.items()],
        "scales": [dict(id=k, description=v) for k, v in SCALES.items()],
        "presets": [dict(id=k, **{kk: vv for kk, vv in v.items()})
                    for k, v in PRESETS.items()],
        "nulls": [
            dict(id="auto", description="spin when coordinates are present, otherwise msr"),
            dict(id="spin", description="spherical rotation; keeps spatial autocorrelation"),
            dict(id="msr", description="Moran spectral randomisation; keeps autocorrelation"),
            dict(id="permutation", description="plain shuffle; smoke tests only"),
        ],
    }
