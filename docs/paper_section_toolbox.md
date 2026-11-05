# Manuscript-ready text for the SFCX toolbox

> 用途：给 D 线稿件补一节「Software / reusability」，把论文从"一个发现"包装成"一个可复用框架"。
> 以下英文可直接放进 Methods 末尾或作为独立小节。方括号处按需替换。

---

## 1. Software subsection (Methods, ~180 words)

**A coordinate system for SFC analyses.**
Because an SFC value is a function of the data *and* of the analysis coordinate, we implemented the
four coordinate layers as an explicit, executable object (`sfcx`; https://github.com/[repo]/sfcx).
A coordinate is the tuple (definition, pipeline, estimator, scale); the package refuses to run a
partial coordinate and records the full tuple with every result. Three analysis routes are provided
as first-class operations: (i) *topography*, the dispersion map at a chosen estimator, its alignment
with the cortical hierarchy, and its focal epicentres; (ii) *individual signal*, the
deviation-spectrum slope computed at regional, network, edge and within-minus-between scales, with a
three-tier partial-correlation ladder and a leakage-controlled incremental-prediction benchmark; and
(iii) *separability*, the partial association between SFC dispersion and the hierarchy after
controlling FC dispersion.

Spatial inference uses one of three null engines, and the engine is written into every reported
p value: spherical spin permutations when parcel coordinates are supplied, Moran spectral
randomisation otherwise, and plain permutation only for smoke tests. A single command
(`sfcx sweep`) recomputes the entire definition × estimator matrix, which is the analysis behind
Fig. [X]; a second command (`sfcx demo`) regenerates a synthetic cohort and its full report, so the
manuscript's central claim can be reproduced without access to the original data.

---

## 2. Availability statement

```
Code availability
SFCX is released under the MIT licence at https://github.com/[repo]/sfcx. The version used here is
archived at https://doi.org/[zenodo-doi]. The repository contains the full coordinate system, the
three analysis routes, three null engines (spin / Moran spectral randomisation / permutation), a
synthetic-cohort generator, a command-line interface and a self-contained HTML report. A single
command (`sfcx demo`) reproduces the coordinate sweep shown in Fig. [X] from synthetic data in
approximately 3 s; analysis of the study cohorts uses the same interface with the archived edge
matrices.
```

---

## 3. What to claim, and what not to

| 可以写 | 不要写 |
|---|---|
| "the coordinate is a first-class object; results are not interpretable without it" | "we solve the reproducibility problem of SFC" |
| "the same data yield axis values from −0.7 to +0.5 depending on the coordinate" | "the axis does not exist" |
| "the toolbox makes the estimator explicit and portable" | "the toolbox replaces existing pipelines" |
| "spin and MSR agree in direction but differ in p" | "p values are now correct" |

---

## 4. Suggested figure caption (if a toolbox panel is added)

> **Fig. X | The readout coordinate system is executable.** a, A coordinate is the tuple
> (definition, pipeline, estimator, scale); any partial specification is rejected. b, Axis strength
> for eight definitions × three estimators computed on one cohort with one null model; the sign of
> the axis is a property of the coordinate, not of the data. c, The three analysis routes
> (topography, individual signal, separability) and the outputs each returns. d, One-command
> reproduction: synthetic cohort, route, sweep and report.

---

## 5. Cover-letter sentence (for a journal that values tools)

> To make the coordinate argument usable rather than merely stated, we release SFCX: the four
> coordinate layers are an executable object, the definition × estimator sweep is one command, and a
> synthetic cohort reproduces the central figure without access to the study data. Reviewers can
> therefore check, on their own data, whether a reported connectional axis survives a change of
> estimator — which is the question the field has so far answered inconsistently.
