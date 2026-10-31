"""
Physics- / geometry-based metrics for a binder–antigen interface.

These replace the metal-coordination physics of the original lanthanide
pipeline (bond length, coordination number, oxygen fraction, ...) with the
standard quantities used to judge a protein–protein / antibody–antigen
interface. None of them needs Rosetta or any external binary — everything is
computed from coordinates with biotite + numpy.

Metrics
-------
hotspot_recovery          fraction of configured hotspot residues contacted
n_hotspot_contacts        count of configured hotspots contacted
hotspot_min_dist_A        per-hotspot nearest binder heavy-atom distance (dict)
avoid_contacts            count of 'avoid' residues contacted (want 0)
interface_bsa_A2          buried surface area = SASA_binder + SASA_antigen
                          - SASA_complex  (full buried area, both sides)
n_interface_res           residues with any heavy atom within contact cutoff
n_hbonds                  cross-interface hydrogen bonds (biotite geometry)
n_salt_bridges            cross-interface charged pairs (<4 Å)
n_clashes                 cross-interface heavy-atom pairs closer than clash cut
sc_proxy                  shape-complementarity proxy in [0,1] (see note)
interface_hydrophobic_frac  hydrophobic share of the buried interface area
binder_radius_gyration_A  compactness sanity check on the binder

Note on sc_proxy
----------------
True Lawrence–Colman shape complementarity (Sc) needs the molecular surface and
is normally taken from CCP4 `sc` or Rosetta. We provide a fast, dependency-free
proxy: the fraction of binder interface atoms that have a complementary antigen
atom within a snug shell (good packing) rather than a gap. It correlates with Sc
and is fine for ranking, but is not the published Sc statistic — swap in real Sc
if you have it (the column name stays the same).
"""

from __future__ import annotations
import numpy as np
import biotite.structure as struc

from .config import Target, PHYSICS_METRICS


# Kyte–Doolittle "hydrophobic" set (side-chain character, coarse but standard).
_HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
_POS = {"LYS", "ARG", "HIS"}
_NEG = {"ASP", "GLU"}
# Atoms carrying the formal charge used for salt bridges.
_POS_ATOMS = {"LYS": {"NZ"}, "ARG": {"NH1", "NH2", "NE"}, "HIS": {"ND1", "NE2"}}
_NEG_ATOMS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}

_CLASH_CUT_A = 2.2      # heavy-atom pair closer than this across the interface
_SALT_CUT_A = 4.0
_SC_SHELL_LO, _SC_SHELL_HI = 3.0, 5.0   # "snug" packing shell for the proxy


def _heavy(arr: struc.AtomArray) -> struc.AtomArray:
    return arr[arr.element != "H"]


def _residue_contacts(binder, antigen, cutoff):
    """Boolean: for each antigen residue, is any heavy atom within `cutoff`
    of any binder heavy atom? Returns (res_ids, contacted_mask)."""
    b = _heavy(binder).coord
    res_ids = np.unique(antigen.res_id)
    contacted = np.zeros(len(res_ids), dtype=bool)
    ac = _heavy(antigen)
    for i, rid in enumerate(res_ids):
        rc = ac.coord[ac.res_id == rid]
        if rc.size == 0:
            continue
        # min pairwise distance res-atoms <-> binder-atoms
        d2 = ((rc[:, None, :] - b[None, :, :]) ** 2).sum(-1)
        if d2.size and np.sqrt(d2.min()) <= cutoff:
            contacted[i] = True
    return res_ids, contacted


def hotspot_metrics(binder, antigen, target: Target) -> dict:
    """Hotspot recovery + per-hotspot nearest-distance + avoid-contacts."""
    b = _heavy(binder).coord
    ac = _heavy(antigen)
    out_min = {}
    engaged = 0
    for rid in target.hotspots:
        rc = ac.coord[ac.res_id == rid]
        if rc.size == 0:
            out_min[int(rid)] = np.nan          # hotspot not present in structure
            continue
        d = np.sqrt(((rc[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min()
        out_min[int(rid)] = float(d)
        if d <= target.contact_cutoff_A:
            engaged += 1
    n_present = sum(1 for v in out_min.values() if not np.isnan(v))
    recovery = engaged / n_present if n_present else np.nan

    avoid_hit = 0
    for rid in target.avoid:
        rc = ac.coord[ac.res_id == rid]
        if rc.size == 0:
            continue
        d = np.sqrt(((rc[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min()
        if d <= target.contact_cutoff_A:
            avoid_hit += 1

    return {
        "hotspot_recovery": recovery,
        "n_hotspot_contacts": engaged,
        "hotspot_min_dist_A": out_min,
        "avoid_contacts": avoid_hit,
    }


def buried_surface_area(binder, antigen):
    """Full buried area and the hydrophobic fraction of it.

    Uses Shrake–Rupley SASA. Hydrophobic fraction is the share of buried area
    contributed by hydrophobic residues (both sides)."""
    complex_arr = binder + antigen
    try:
        sasa_c = struc.sasa(complex_arr, vdw_radii="Single")
        sasa_b = struc.sasa(binder, vdw_radii="Single")
        sasa_a = struc.sasa(antigen, vdw_radii="Single")
    except Exception:
        return {"interface_bsa_A2": np.nan, "interface_hydrophobic_frac": np.nan}

    sasa_c = np.nan_to_num(sasa_c)
    sasa_b = np.nan_to_num(sasa_b)
    sasa_a = np.nan_to_num(sasa_a)

    nb = binder.array_length()
    # Per-atom buried area: free (isolated) minus bound (in complex).
    buried_b = sasa_b - sasa_c[:nb]
    buried_a = sasa_a - sasa_c[nb:]
    buried_b[buried_b < 0] = 0
    buried_a[buried_a < 0] = 0
    bsa = float(buried_b.sum() + buried_a.sum())

    # Hydrophobic share of buried area.
    hb = np.isin(binder.res_name, list(_HYDROPHOBIC))
    ha = np.isin(antigen.res_name, list(_HYDROPHOBIC))
    hydro = float(buried_b[hb].sum() + buried_a[ha].sum())
    frac = hydro / bsa if bsa > 0 else np.nan
    return {"interface_bsa_A2": bsa, "interface_hydrophobic_frac": frac}


def interface_residues(binder, antigen, cutoff):
    """Indices/ids of binder and antigen residues at the interface."""
    bid, bmask = _residue_contacts(antigen, binder, cutoff)   # binder res contacted by antigen
    aid, amask = _residue_contacts(binder, antigen, cutoff)   # antigen res contacted by binder
    return {
        "binder_res": bid[bmask],
        "antigen_res": aid[amask],
        "n_interface_res": int(bmask.sum() + amask.sum()),
    }


def count_hbonds(binder, antigen):
    """Cross-interface hydrogen bonds via biotite's geometric detector."""
    complex_arr = binder + antigen
    nb = binder.array_length()
    try:
        triplets = struc.hbond(complex_arr)   # (donor, hydrogen, acceptor) atom idx
    except Exception:
        return np.nan
    if len(triplets) == 0:
        return 0
    donor, _, acceptor = triplets[:, 0], triplets[:, 1], triplets[:, 2]
    # cross-interface = donor and acceptor on opposite sides
    cross = ((donor < nb) & (acceptor >= nb)) | ((donor >= nb) & (acceptor < nb))
    return int(cross.sum())


def count_salt_bridges(binder, antigen, cutoff=_SALT_CUT_A):
    def charged_atoms(arr, table):
        sel = []
        for rn, atoms in table.items():
            m = (arr.res_name == rn) & np.isin(arr.atom_name, list(atoms))
            if m.any():
                sel.append(arr.coord[m])
        return np.vstack(sel) if sel else np.empty((0, 3))

    b_pos = charged_atoms(binder, _POS_ATOMS); b_neg = charged_atoms(binder, _NEG_ATOMS)
    a_pos = charged_atoms(antigen, _POS_ATOMS); a_neg = charged_atoms(antigen, _NEG_ATOMS)

    def pairs(x, y):
        if len(x) == 0 or len(y) == 0:
            return 0
        d = np.sqrt(((x[:, None, :] - y[None, :, :]) ** 2).sum(-1))
        return int((d <= cutoff).any(axis=1).sum())

    # binder+ <-> antigen-   and   binder- <-> antigen+
    return pairs(b_pos, a_neg) + pairs(b_neg, a_pos)


def count_clashes(binder, antigen, cutoff=_CLASH_CUT_A):
    b = _heavy(binder).coord
    a = _heavy(antigen).coord
    if len(b) == 0 or len(a) == 0:
        return 0
    # chunk to keep memory bounded for large complexes
    n = 0
    step = 2000
    for i in range(0, len(b), step):
        chunk = b[i:i + step]
        d2 = ((chunk[:, None, :] - a[None, :, :]) ** 2).sum(-1)
        n += int((d2 < cutoff ** 2).sum())
    return n


def sc_proxy(binder, antigen, lo=_SC_SHELL_LO, hi=_SC_SHELL_HI):
    """Packing proxy for shape complementarity in [0,1] (see module docstring)."""
    b = _heavy(binder).coord
    a = _heavy(antigen).coord
    if len(b) == 0 or len(a) == 0:
        return np.nan
    # interface binder atoms = those with an antigen atom within hi
    good = 0
    inter = 0
    step = 2000
    for i in range(0, len(b), step):
        chunk = b[i:i + step]
        dmin = np.sqrt(((chunk[:, None, :] - a[None, :, :]) ** 2).sum(-1)).min(axis=1)
        at_interface = dmin <= hi
        inter += int(at_interface.sum())
        good += int(((dmin >= lo) & (dmin <= hi)).sum())
    return good / inter if inter else np.nan


def binder_rg(binder):
    c = binder[binder.atom_name == "CA"].coord
    if len(c) == 0:
        c = _heavy(binder).coord
    if len(c) == 0:
        return np.nan
    com = c.mean(0)
    return float(np.sqrt(((c - com) ** 2).sum(1).mean()))


def score_structure(binder, antigen, target: Target, pae=None) -> dict:
    """Compute every physics metric for one binder–antigen complex.

    `pae` (optional) is the folding model's residue×residue PAE matrix; if given,
    interface/hotspot PAE is averaged over the relevant residue pairs.
    """
    cutoff = target.contact_cutoff_A
    m = {}
    m.update(hotspot_metrics(binder, antigen, target))
    m.update(buried_surface_area(binder, antigen))
    iface = interface_residues(binder, antigen, cutoff)
    m["n_interface_res"] = iface["n_interface_res"]
    m["n_hbonds"] = count_hbonds(binder, antigen)
    m["n_salt_bridges"] = count_salt_bridges(binder, antigen)
    m["n_clashes"] = count_clashes(binder, antigen)
    m["sc_proxy"] = sc_proxy(binder, antigen)
    m["binder_radius_gyration_A"] = binder_rg(binder)

    # Interface pLDDT: mean B-factor over binder interface residues.
    bi = iface["binder_res"]
    if len(bi):
        mask = np.isin(binder.res_id, bi)
        m["interface_plddt"] = float(np.nanmean(binder.b_factor[mask]))
    else:
        m["interface_plddt"] = np.nan
    m["mean_plddt"] = float(np.nanmean(binder.b_factor)) if binder.array_length() else np.nan

    # PAE summaries if a matrix is available and shapes are consistent.
    m["pae_interface"] = np.nan
    m["pae_hotspot"] = np.nan
    if pae is not None:
        m.update(_pae_summaries(binder, antigen, iface, target, pae))

    # keep only JSON-friendly scalars in the flat record; the per-hotspot dict
    # is handy but stored separately by the pipeline
    m["_hotspot_min_dist_A"] = m.pop("hotspot_min_dist_A")
    return m


def _pae_summaries(binder, antigen, iface, target, pae):
    """Mean PAE over binder<->antigen and binder<->hotspot residue pairs.

    Assumes the PAE matrix is indexed by a simple concatenation of binder then
    antigen residues in file order — the common AF3/Boltz convention. If your
    indexing differs, adjust here.
    """
    out = {"pae_interface": np.nan, "pae_hotspot": np.nan}
    nb = len(np.unique(binder.res_id))
    na = len(np.unique(antigen.res_id))
    if pae.shape[0] != nb + na:
        return out   # can't line up safely; leave NaN
    bi = np.arange(nb)
    ai = np.arange(nb, nb + na)
    block = pae[np.ix_(bi, ai)]
    out["pae_interface"] = float(np.nanmean(block)) if block.size else np.nan

    ag_res = np.unique(antigen.res_id)
    hs_idx = [nb + k for k, rid in enumerate(ag_res) if rid in set(target.hotspots)]
    if hs_idx:
        hblock = pae[np.ix_(bi, np.array(hs_idx))]
        out["pae_hotspot"] = float(np.nanmean(hblock)) if hblock.size else np.nan
    return out
