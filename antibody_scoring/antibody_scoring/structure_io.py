"""
Structure loading and chain splitting.

Reads PDB / mmCIF (optionally .gz) into a biotite AtomArray, and splits it into
"binder" and "antigen" selections according to config. Model-confidence B-factor
(pLDDT) is preserved in the `b_factor` annotation, which is how AF3/Chai/Boltz
write per-residue pLDDT into the coordinate file.
"""

from __future__ import annotations
from pathlib import Path
import gzip
import numpy as np

import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb


def load_structure(path: str | Path) -> struc.AtomArray:
    """Load the first model of a PDB/mmCIF/.gz file as an AtomArray."""
    path = Path(path)
    name = path.name.lower()
    opener = gzip.open if name.endswith(".gz") else open
    stem = name[:-3] if name.endswith(".gz") else name

    if stem.endswith((".cif", ".mmcif", ".bcif")):
        with opener(path, "rt") as fh:
            f = pdbx.CIFFile.read(fh)
        arr = pdbx.get_structure(f, model=1, extra_fields=["b_factor", "charge"])
    elif stem.endswith((".pdb", ".ent")):
        with opener(path, "rt") as fh:
            f = pdb.PDBFile.read(fh)
        arr = f.get_structure(model=1, extra_fields=["b_factor", "charge"])
    else:
        raise ValueError(f"Unrecognised structure extension: {path.name}")

    # Keep only the polymer; drop waters/ions/ligands for interface geometry.
    arr = arr[struc.filter_amino_acids(arr)]
    return arr


def chain_ids(arr: struc.AtomArray) -> list[str]:
    return list(dict.fromkeys(arr.chain_id.tolist()))


def split_binder_antigen(
    arr: struc.AtomArray,
    antigen_chain: str,
    binder_chains: tuple[str, ...] | None,
) -> tuple[struc.AtomArray, struc.AtomArray]:
    """Return (binder_atoms, antigen_atoms).

    If `binder_chains` is None, the binder is every chain that is not the
    antigen chain.
    """
    present = set(chain_ids(arr))
    if antigen_chain not in present:
        raise ValueError(
            f"Antigen chain '{antigen_chain}' not in structure "
            f"(chains present: {sorted(present)})"
        )
    antigen = arr[arr.chain_id == antigen_chain]

    if binder_chains is None:
        binder = arr[arr.chain_id != antigen_chain]
    else:
        missing = set(binder_chains) - present
        if missing:
            raise ValueError(f"Binder chain(s) {sorted(missing)} not in structure")
        binder = arr[np.isin(arr.chain_id, list(binder_chains))]

    if binder.array_length() == 0:
        raise ValueError("Binder selection is empty — check BINDER_CHAINS / antigen_chain")
    return binder, antigen
