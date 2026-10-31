"""
Scan a directory of designs and build the per-design score table.

Layout expected (flexible — see GROUP_FROM in config):

    designs/
      rfdiffusion/
        design_0001.pdb         (+ optional design_0001_confidence.json)
        design_0002.cif.gz
      bindcraft/
        ...

Each structure is one binder–antigen complex. The group label (used to compare
sets in the plots/stats) is derived per config.GROUP_FROM.
"""

from __future__ import annotations
from pathlib import Path
import re
import json
import numpy as np
import pandas as pd

from . import config as C
from .structure_io import load_structure, split_binder_antigen
from .confidence import load_confidence
from .interface_metrics import score_structure


STRUCTURE_GLOBS = ("*.pdb", "*.pdb.gz", "*.cif", "*.cif.gz", "*.ent", "*.mmcif")


def find_structures(root: str | Path) -> list[Path]:
    root = Path(root)
    files: list[Path] = []
    for pat in STRUCTURE_GLOBS:
        files += root.rglob(pat)
    # de-dup (a .cif and .cif.gz pair is unusual but guard anyway) & sort
    return sorted(set(files))


def group_label(path: Path, root: Path) -> str:
    if C.GROUP_FROM == "parent_dir":
        rel = path.relative_to(root)
        return rel.parts[0] if len(rel.parts) > 1 else "all"
    if C.GROUP_FROM == "filename_prefix":
        return path.name.split("_", 1)[0]
    if C.GROUP_FROM == "regex":
        m = re.search(C.GROUP_REGEX, path.name)
        return m.group(1) if m else "ungrouped"
    return "all"


def score_directory(root: str | Path, target: C.Target,
                    binder_chains=C.BINDER_CHAINS, verbose=True) -> pd.DataFrame:
    root = Path(root)
    structures = find_structures(root)
    if not structures:
        raise FileNotFoundError(f"No structures ({STRUCTURE_GLOBS}) under {root}")

    records = []
    hotspot_rows = []
    for i, sp in enumerate(structures, 1):
        rec = {"design": sp.stem, "path": str(sp),
               "group": group_label(sp, root), "target": target.name}
        try:
            arr = load_structure(sp)
            binder, antigen = split_binder_antigen(arr, target.antigen_chain, binder_chains)
            conf = load_confidence(sp)
            m = score_structure(binder, antigen, target, pae=conf.get("pae"))
            per_hotspot = m.pop("_hotspot_min_dist_A")
            rec.update(m)
            rec["ptm"] = conf["ptm"]
            rec["iptm"] = conf["iptm"]
            rec["error"] = ""
            for rid, d in per_hotspot.items():
                hotspot_rows.append({"design": sp.stem, "group": rec["group"],
                                     "hotspot": rid, "min_dist_A": d,
                                     "engaged": (d is not None and not np.isnan(d)
                                                 and d <= target.contact_cutoff_A)})
        except Exception as e:                      # one bad file never kills the run
            rec["error"] = f"{type(e).__name__}: {e}"
            if verbose:
                print(f"  ! {sp.name}: {rec['error']}")
        records.append(rec)
        if verbose and (i % 50 == 0 or i == len(structures)):
            print(f"  scored {i}/{len(structures)}")

    df = pd.DataFrame(records)
    _apply_group_order(df)
    df.attrs["hotspots"] = pd.DataFrame(hotspot_rows)
    df.attrs["target"] = target
    return df


def _apply_group_order(df: pd.DataFrame):
    groups = list(C.GROUP_ORDER)
    for g in df["group"].unique():
        if g not in groups:
            groups.append(g)
    df["group"] = pd.Categorical(df["group"], categories=groups, ordered=True)


def add_pass_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Add boolean pass columns + an overall `passes` flag from config thresholds."""
    masks = []
    for m, t in C.CONF_THRESHOLDS.items():
        if m in df:
            ok = df[m] >= t
            df[f"pass_{m}"] = ok
            masks.append(ok.fillna(False))
    for m, spec in C.PHYSICS_METRICS.items():
        if spec.threshold is None or m not in df:
            continue
        ok = df[m] >= spec.threshold if spec.higher_better else df[m] <= spec.threshold
        df[f"pass_{m}"] = ok
        masks.append(ok.fillna(False))
    df["passes"] = np.logical_and.reduce(masks) if masks else True
    return df
