"""
Parse folding-model confidence outputs.

Supports the common layouts written by AlphaFold3, Chai-1 and Boltz-1/2:
every one of them emits a per-model JSON with pTM / ipTM and a PAE matrix, and
writes per-residue pLDDT into the B-factor column of the structure file.

The functions here are deliberately permissive: they look for several key
spellings and fall back to NaN rather than crashing, so a missing metric shows
up as a blank in the scores table instead of killing the run.

What you get back
-----------------
    ptm, iptm              global confidences (float, NaN if absent)
    mean_plddt             mean over all binder residues
    interface_plddt        mean pLDDT over binder interface residues (set later,
                           once the interface is known — see interface_metrics)
    pae_interface          mean PAE over binder<->antigen residue pairs
    pae_hotspot            mean PAE over binder<->hotspot residue pairs
"""

from __future__ import annotations
from pathlib import Path
import json
import math
import numpy as np


_PTM_KEYS  = ("ptm", "pTM", "ptm_score")
_IPTM_KEYS = ("iptm", "ipTM", "iptm_score", "interface_ptm")
_PAE_KEYS  = ("pae", "predicted_aligned_error", "pae_matrix")


def _first(d: dict, keys) -> float:
    for k in keys:
        if k in d and d[k] is not None:
            v = d[k]
            try:
                return float(v)
            except (TypeError, ValueError):
                pass
    return math.nan


def find_confidence_json(structure_path: str | Path) -> Path | None:
    """Heuristically locate the confidence JSON that pairs with a structure.

    Looks next to the structure for files whose name shares the stem or that
    match the usual AF3/Chai/Boltz confidence-file patterns.
    """
    structure_path = Path(structure_path)
    d = structure_path.parent
    stem = structure_path.name
    for suf in (".cif.gz", ".pdb.gz", ".cif", ".pdb", ".ent", ".mmcif"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break

    candidates = []
    for pat in (
        f"{stem}*confidence*.json", f"{stem}*summary*.json", f"{stem}*.json",
        "*summary_confidences*.json", "*confidence*.json", "*scores*.json",
    ):
        candidates += sorted(d.glob(pat))
    for c in candidates:
        if c.suffix == ".json":
            return c
    return None


def load_confidence(structure_path: str | Path) -> dict:
    """Return {ptm, iptm, pae} for a structure, parsing its paired JSON if any.

    `pae` is a square np.ndarray (residue x residue) or None.
    """
    out = {"ptm": math.nan, "iptm": math.nan, "pae": None}
    jpath = find_confidence_json(structure_path)
    if jpath is None:
        return out
    try:
        data = json.loads(Path(jpath).read_text())
    except (json.JSONDecodeError, OSError):
        return out

    # Some tools nest everything under a top-level key; flatten one level.
    flat = dict(data)
    for v in list(data.values()):
        if isinstance(v, dict):
            for k2, v2 in v.items():
                flat.setdefault(k2, v2)

    out["ptm"] = _first(flat, _PTM_KEYS)
    out["iptm"] = _first(flat, _IPTM_KEYS)
    for k in _PAE_KEYS:
        if k in flat and flat[k] is not None:
            try:
                arr = np.asarray(flat[k], dtype=float)
                if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
                    out["pae"] = arr
                    break
            except (TypeError, ValueError):
                pass
    return out
