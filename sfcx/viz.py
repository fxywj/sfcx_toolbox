# -*- coding: utf-8 -*-
"""Figures.  Light background, thin spines, no chartjunk."""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

INK = "#1f3a5f"
TEAL = "#2a9d8f"
ORANGE = "#e76f51"
SAND = "#e9c46a"
GREY = "#8d99ae"

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#c9d1d9",
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})


def _save(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return path


def coordinate_matrix(mat, rows, cols, path, pmat=None, title="", cbar_label="Spearman r"):
    """Heatmap of an axis statistic across metric x estimator (or scale)."""
    m = np.asarray(mat, dtype=float)
    fig, ax = plt.subplots(figsize=(max(4.2, 1.15 * len(cols) + 1.6),
                                    max(3.0, 0.55 * len(rows) + 1.4)))
    vmax = np.nanmax(np.abs(m)) or 1.0
    im = ax.imshow(m, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
                   aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=20, ha="right")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if not np.isfinite(m[i, j]):
                continue
            star = ""
            if pmat is not None and np.isfinite(pmat[i, j]):
                star = "*" if pmat[i, j] < 0.05 else ""
            ax.text(j, i, "%.2f%s" % (m[i, j], star), ha="center", va="center",
                    fontsize=7.5, color="white" if abs(m[i, j]) > 0.55 * vmax else INK)
    ax.set_title(title, fontsize=10, pad=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(cbar_label)
    cb.outline.set_visible(False)
    return _save(fig, path)


def spectrum(v, path, slope=None, label="Regional dispersion", color=INK, title=""):
    v = np.asarray(v, dtype=float)
    order = np.argsort(-np.nan_to_num(v, nan=-np.inf))
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    x = np.arange(v.size)
    base = np.nanmin(v)
    ax.plot(x, v[order], color=color, lw=1.6)
    ax.fill_between(x, base, v[order], color=color, alpha=0.08)
    if slope is not None:
        ax.plot(x[[0, -1]], [v[order][0], v[order][-1]], ls="--", lw=1, color=ORANGE,
                label="slope = %.3g" % slope)
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("rank")
    ax.set_ylabel(label)
    ax.set_title("Descending spectrum", fontsize=10)
    ax2 = axes[1]
    ax2.hist(v[np.isfinite(v)], bins=28, color=TEAL, alpha=0.75, edgecolor="white")
    ax2.set_xlabel(label)
    ax2.set_ylabel("count")
    ax2.set_title("Distribution", fontsize=10)
    fig.suptitle(title, fontsize=10, y=1.03)
    return _save(fig, path)


def layers_bar(within, between, path, p=None, title=""):
    fig, ax = plt.subplots(figsize=(3.4, 2.9))
    vals = [within, between]
    bars = ax.bar(["within", "between"], vals, color=[ORANGE, GREY], width=0.55)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v, "%.3f" % v, ha="center", va="bottom",
                fontsize=8)
    ax.set_ylabel("mean edge dispersion")
    if p is not None:
        ax.set_title("delta = %.3f, p = %.3g" % (within - between, p), fontsize=9.5)
    else:
        ax.set_title(title, fontsize=9.5)
    return _save(fig, path)


def epicentre_strip(v, networks, focal, path, title=""):
    """Network-ordered dispersion profile plus a per-network bar.

    No cortical surface is implied: this is an ordering by network, not anatomy.
    """
    v = np.asarray(v, dtype=float)
    net = np.asarray(networks)
    order = np.lexsort((-v, net))
    x = np.arange(v.size)
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 4.0),
                             gridspec_kw={"height_ratios": [1.25, 1]})

    ax = axes[0]
    ax.plot(x, v[order], color=INK, lw=1.4)
    ax.fill_between(x, np.nanmin(v), v[order], color=INK, alpha=0.07)
    for b in np.flatnonzero(np.diff(net[order]) != 0):
        ax.axvline(b + 0.5, color=GREY, lw=0.8, ls=":")
    pos = {int(r): k for k, r in enumerate(order)}
    fx = [pos[int(r)] for r in focal if int(r) in pos]
    fy = [v[order][pos[int(r)]] for r in focal if int(r) in pos]
    ax.scatter(fx, fy, s=26, color=ORANGE, zorder=5, edgecolor="white", linewidth=0.6,
               label="focal epicentre")
    ax.set_xlim(-1, v.size)
    ax.set_ylabel("dispersion")
    ax.set_xticks([])
    ax.set_title("Dispersion profile, ordered by network then value", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax2 = axes[1]
    gs = np.unique(net)
    means = [np.nanmean(v[net == g]) if np.any(net == g) else np.nan for g in gs]
    ax2.bar([str(int(g)) for g in gs], means, color=TEAL, width=0.6)
    for i, m in enumerate(means):
        if np.isfinite(m):
            ax2.text(i, m, "%.3f" % m, ha="center", va="bottom", fontsize=7.5)
    ax2.set_ylabel("network mean")
    ax2.set_xlabel("network")
    ax2.set_ylim(0, np.nanmax(means) * 1.18)
    fig.suptitle(title, fontsize=10, y=1.0)
    return _save(fig, path)


def reliability_scatter(va, vb, path, r=None):
    fig, ax = plt.subplots(figsize=(3.2, 3.0))
    ax.scatter(va, vb, s=6, color=INK, alpha=0.35, linewidths=0)
    lo = np.nanmin([np.nanmin(va), np.nanmin(vb)])
    hi = np.nanmax([np.nanmax(va), np.nanmax(vb)])
    ax.plot([lo, hi], [lo, hi], ls="--", lw=1, color=ORANGE)
    ax.set_xlabel("half 1")
    ax.set_ylabel("half 2")
    ax.set_title("split-half r = %.3f" % r if r is not None else "split-half", fontsize=9.5)
    return _save(fig, path)


def individual_scatter(x, y, path, r=None, p=None, xlabel="slope", ylabel="behaviour"):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    ax.scatter(x[ok], y[ok], s=9, color=TEAL, alpha=0.45, linewidths=0)
    if ok.sum() > 2:
        b = np.polyfit(x[ok], y[ok], 1)
        xs = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 50)
        ax.plot(xs, np.polyval(b, xs), color=ORANGE, lw=1.4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if r is not None:
        ax.set_title("r = %.3f%s" % (r, "" if p is None else ", p = %.3g" % p), fontsize=9.5)
    return _save(fig, path)


def route_card(text, path, title="Active coordinate"):
    fig, ax = plt.subplots(figsize=(5.2, 2.2))
    ax.axis("off")
    ax.add_patch(plt.Rectangle((0.01, 0.05), 0.98, 0.9, facecolor="#f6f8fa",
                               edgecolor="#c9d1d9", transform=ax.transAxes))
    ax.text(0.05, 0.86, title, fontsize=10, fontweight="bold", color=INK,
            transform=ax.transAxes)
    ax.text(0.05, 0.62, text, fontsize=8.6, color=INK, va="top",
            transform=ax.transAxes, linespacing=1.5)
    return _save(fig, path)
