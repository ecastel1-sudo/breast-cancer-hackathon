"""
Plotting — generalised from the original lanthanide notebook.

Kept the structure you already had (corner pair-plots with threshold lines,
per-group KDE panels, a pass-rate summary), but:
  * `campaign`  -> `group`
  * metal-coordination physics panel -> interface physics panel
  * thresholds/labels now come from config, not hard-coded metal chemistry.
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

from . import config as C


# --------------------------------------------------------------------------- #
# palette / order helpers
# --------------------------------------------------------------------------- #
def group_order(df):
    return [g for g in df["group"].cat.categories if (df["group"] == g).any()]


def palette(df):
    order = group_order(df)
    base = sns.color_palette("colorblind", len(order))
    pal = {}
    for i, g in enumerate(order):
        pal[g] = C.GROUP_COLORS.get(g, base[i])
    return pal


def _label(col):
    if col in C.AXIS_LABELS:
        return C.AXIS_LABELS[col]
    if col in C.PHYSICS_METRICS:
        return C.PHYSICS_METRICS[col].label
    return col


def _thresholds(col):
    if col in C.CONF_THRESHOLDS:
        return [C.CONF_THRESHOLDS[col]]
    if col in C.PHYSICS_METRICS and C.PHYSICS_METRICS[col].threshold is not None:
        return [C.PHYSICS_METRICS[col].threshold]
    if col in C.RMSD_THRESHOLDS:
        return [C.RMSD_THRESHOLDS[col]]
    return []


# --------------------------------------------------------------------------- #
# 1. confidence pair-plot (one figure per group)
# --------------------------------------------------------------------------- #
def _add_threshold_lines(grid, color="crimson", lw=0.9, ls="--", alpha=0.85):
    for i, yvar in enumerate(grid.y_vars):
        for j, xvar in enumerate(grid.x_vars):
            ax = grid.axes[i][j]
            if ax is None:
                continue
            for t in _thresholds(xvar):
                ax.axvline(t, color=color, lw=lw, ls=ls, alpha=alpha, zorder=5)
            if i != j:
                for t in _thresholds(yvar):
                    ax.axhline(t, color=color, lw=lw, ls=ls, alpha=alpha, zorder=5)


def _relabel(grid, diag_ylabel=""):
    for i, yvar in enumerate(grid.y_vars):
        for j, xvar in enumerate(grid.x_vars):
            ax = grid.axes[i][j]
            if ax is None:
                continue
            if ax.get_xlabel():
                ax.set_xlabel(_label(xvar))
            if i == j:
                ax.set_ylabel(diag_ylabel)
            elif ax.get_ylabel():
                ax.set_ylabel(_label(yvar))


def pairplot_group(df, group, metrics=None, sample_n=C.SAMPLE_N, figdir=C.FIGDIR):
    metrics = metrics or [m for m in C.CONF_METRICS if m in df]
    g = df[df["group"] == group]
    g = g[metrics].dropna(how="all")
    if sample_n and len(g) > sample_n:
        g = g.sample(sample_n, random_state=0)
    if g.empty:
        return None
    grid = sns.PairGrid(g, vars=metrics, corner=True)
    grid.map_lower(sns.scatterplot, s=6, alpha=0.3, linewidth=0,
                   color=palette(df)[group], rasterized=True)
    grid.map_diag(sns.kdeplot, fill=True, color=palette(df)[group])
    _add_threshold_lines(grid)
    _relabel(grid)
    grid.figure.suptitle(f"Confidence metrics — {group}", y=1.02)
    Path(figdir).mkdir(parents=True, exist_ok=True)
    grid.figure.savefig(Path(figdir) / f"pairplot_{group}.png", bbox_inches="tight", dpi=150)
    return grid


# --------------------------------------------------------------------------- #
# 2. interface physics panel (replaces the metal-coordination panel)
# --------------------------------------------------------------------------- #
def physics_panel(df, metrics=None, figdir=C.FIGDIR):
    metrics = metrics or [m for m in C.PHYSICS_METRICS if m in df and df[m].notna().any()]
    pal = palette(df)
    order = group_order(df)
    ncol = 3
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.0 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, col in zip(axes, metrics):
        for grp in order:
            s = df.loc[df["group"] == grp, col].dropna()
            if len(s) < 2 or s.nunique() < 2:
                # too few / constant: draw a rug so it's not invisible
                if len(s):
                    ax.plot(s, np.zeros_like(s), "|", color=pal[grp], alpha=0.6)
                continue
            sns.kdeplot(x=s, ax=ax, color=pal[grp], fill=True, alpha=0.12, linewidth=1.2)
        for t in _thresholds(col):
            ax.axvline(t, color="crimson", ls="--", lw=0.9, alpha=0.85, zorder=5)
        ax.set_xlabel(_label(col))
        ax.set_ylabel("density")
    for ax in axes[len(metrics):]:
        ax.set_visible(False)

    handles = [Patch(facecolor=pal[g], edgecolor=pal[g], alpha=0.55, label=g) for g in order]
    fig.legend(handles=handles, loc="upper center", ncol=min(len(order), 6),
               frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.04))
    fig.suptitle("Interface physics metrics — by group", y=1.0)
    fig.tight_layout()
    Path(figdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(figdir) / "physics_panel.png", bbox_inches="tight", dpi=150)
    return fig


# --------------------------------------------------------------------------- #
# 3. hotspot recovery heatmap (which epitope residues get engaged)
# --------------------------------------------------------------------------- #
def hotspot_heatmap(df, figdir=C.FIGDIR):
    hs = df.attrs.get("hotspots")
    if hs is None or hs.empty:
        return None
    frac = (hs.groupby(["group", "hotspot"], observed=True)["engaged"]
              .mean().unstack("hotspot"))
    fig, ax = plt.subplots(figsize=(0.5 * frac.shape[1] + 3, 0.5 * frac.shape[0] + 1.5))
    sns.heatmap(frac, annot=True, fmt=".0%", cmap="mako", vmin=0, vmax=1,
                cbar_kws={"label": "fraction of designs engaging hotspot"}, ax=ax)
    ax.set_xlabel("antigen hotspot residue"); ax.set_ylabel("group")
    ax.set_title("Hotspot engagement by group")
    fig.tight_layout()
    Path(figdir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(figdir) / "hotspot_heatmap.png", bbox_inches="tight", dpi=150)
    return fig
