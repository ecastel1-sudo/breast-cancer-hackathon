"""
Statistics — carried over from the original notebook and generalised.

Same tools you were using (Kruskal–Wallis across groups, pairwise Mann–Whitney
with Cliff's delta effect size, bootstrap CIs, pass-rate tables) but driven by
the generic metric config rather than the metal-RMSD metrics.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations

from . import config as C


def boot_ci(x, stat=np.median, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    b = np.array([stat(rng.choice(x, len(x), replace=True)) for _ in range(n_boot)])
    return np.percentile(b, [2.5, 97.5])


def cliffs_delta(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return (gt - lt) / (len(x) * len(y))


def _mag(d):
    d = abs(d)
    return ("negligible" if d < 0.147 else "small" if d < 0.33
            else "medium" if d < 0.474 else "large")


def descriptive_table(df, metrics=None, min_n=5):
    metrics = metrics or _default_metrics(df)
    rows = []
    for m in metrics:
        for grp in df["group"].cat.categories:
            s = df.loc[df["group"] == grp, m].dropna()
            if len(s) < min_n:
                continue
            q1, q3 = s.quantile([.25, .75])
            lo, hi = boot_ci(s)
            rows.append({"metric": _label(m), "group": grp, "n": len(s),
                         "median": round(float(s.median()), 3),
                         "IQR": f"{q1:.3g}–{q3:.3g}",
                         "median 95% CI": f"{lo:.3g}–{hi:.3g}",
                         "range": f"{s.min():.3g}–{s.max():.3g}"})
    return pd.DataFrame(rows).set_index(["metric", "group"])


def compare_groups(df, metric, min_n=20, max_pairwise_n=4000, seed=0):
    """Kruskal–Wallis across groups + pairwise Mann–Whitney with Cliff's delta."""
    rng = np.random.default_rng(seed)
    groups = {g: df.loc[df["group"] == g, metric].dropna().values
              for g in df["group"].cat.categories}
    groups = {g: v for g, v in groups.items() if len(v) >= min_n}
    out = {"metric": _label(metric), "groups": list(groups)}
    if len(groups) < 2:
        out["note"] = "need ≥2 groups with enough data"
        return out, pd.DataFrame()

    H, p = kruskal(*groups.values())
    out["kruskal_H"] = round(float(H), 3)
    out["kruskal_p"] = float(p)

    rows = []
    for a, b in combinations(groups, 2):
        xa, xb = groups[a], groups[b]
        if len(xa) > max_pairwise_n:
            xa = rng.choice(xa, max_pairwise_n, replace=False)
        if len(xb) > max_pairwise_n:
            xb = rng.choice(xb, max_pairwise_n, replace=False)
        U, pp = mannwhitneyu(xa, xb, alternative="two-sided")
        d = cliffs_delta(xa, xb)
        rows.append({"A": a, "B": b, "n_A": len(groups[a]), "n_B": len(groups[b]),
                     "U": float(U), "p": float(pp),
                     "cliffs_delta": round(float(d), 3), "effect": _mag(d)})
    pw = pd.DataFrame(rows)
    if not pw.empty:
        # Holm correction across the pairwise family
        pw = pw.sort_values("p").reset_index(drop=True)
        mtot = len(pw)
        pw["p_holm"] = np.minimum.accumulate(
            (pw["p"].values * (mtot - np.arange(mtot)))[::-1])[::-1].clip(0, 1)
    return out, pw


def success_table(df, min_n=5, n_boot=2000, seed=0):
    """Pass-rate per metric per group, with a binomial bootstrap CI."""
    rng = np.random.default_rng(seed)
    rows = []
    specs = {m: (C.CONF_THRESHOLDS[m], True) for m in C.CONF_THRESHOLDS if m in df}
    for m, spec in C.PHYSICS_METRICS.items():
        if m in df and spec.threshold is not None:
            specs[m] = (spec.threshold, spec.higher_better)
    for m, t in C.RMSD_THRESHOLDS.items():
        if m in df:
            specs[m] = (t, False)

    for m, (t, higher) in specs.items():
        for grp in df["group"].cat.categories:
            s = df.loc[df["group"] == grp, m].dropna()
            if len(s) < min_n:
                continue
            ok = (s >= t) if higher else (s <= t)
            k = int(ok.sum()); p = k / len(s)
            b = rng.binomial(len(s), p, n_boot) / len(s)
            lo, hi = np.percentile(b, [2.5, 97.5])
            rows.append({"metric": _label(m),
                         "threshold": f"{'≥' if higher else '≤'} {t}",
                         "group": grp, "n": len(s), "n_pass": k,
                         "%_pass": round(100 * p, 1),
                         "95% CI": f"{100*lo:.1f}–{100*hi:.1f}"})
    return pd.DataFrame(rows).set_index(["metric", "group"])


def _default_metrics(df):
    cols = [m for m in C.CONF_METRICS if m in df]
    cols += [m for m in C.PHYSICS_METRICS if m in df and df[m].notna().any()]
    return cols


def _label(col):
    if col in C.AXIS_LABELS:
        return C.AXIS_LABELS[col]
    if col in C.PHYSICS_METRICS:
        return C.PHYSICS_METRICS[col].label
    return col
