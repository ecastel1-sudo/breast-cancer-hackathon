"""
================================================================================
 EDIT THIS FILE.  Almost everything a team needs to change lives here.
================================================================================

This is the single place to configure a scoring run for your own target.
The two things you will touch most often are:

    1. TARGETS          -> which antigen residues are the "hotspots" you want
                           the designed binder/antibody to engage.
    2. CHAINS           -> which chain is the binder and which is the antigen.

Everything downstream (metrics, plots, statistics) reads from here, so you do
NOT need to edit the analysis code to retarget the pipeline.

Residue numbering note
----------------------
Hotspot residues are given as author/PDB residue numbers on the *target*
(antigen) chain, exactly as they appear in your structure files. If your
folding model renumbers chains, check one structure in a viewer first and use
the numbers you actually see there.
"""

from __future__ import annotations
from dataclasses import dataclass, field


# ------------------------------------------------------------------------------
# 1. HOTSPOT TARGETS  <-- the thing most teams change
# ------------------------------------------------------------------------------
# Define one Target per antigen you are designing against. A Target names the
# antigen chain and the list of epitope "hotspot" residues the binder should
# contact. The pipeline then reports, per design, how many of these hotspots
# are actually engaged (hotspot recovery) and how buried/close each one is.
#
# The example below is a placeholder (HER2-like numbering) -- REPLACE the
# residue numbers with your real epitope. You can keep several targets here and
# select one per run with --target on the command line.

@dataclass
class Target:
    name: str                      # free label, used in filenames/plots
    antigen_chain: str             # chain ID of the antigen in the structures
    hotspots: tuple[int, ...]      # author residue numbers of epitope hotspots
    # Optional: residues you specifically do NOT want engaged (e.g. to steer
    # away from a conserved/!safe patch). Reported but never required.
    avoid: tuple[int, ...] = ()
    # A hotspot counts as "engaged" if any binder heavy atom is within this many
    # angstrom of any heavy atom of the hotspot residue.
    contact_cutoff_A: float = 5.0


TARGETS: dict[str, Target] = {

    # ---- REPLACE ME -------------------------------------------------------
    "example_antigen": Target(
        name="example_antigen",
        antigen_chain="A",
        hotspots=(42, 44, 66, 68, 70, 99, 101),   # <-- your epitope residues
        avoid=(),
        contact_cutoff_A=5.0,
    ),
    # -----------------------------------------------------------------------

    # Add more targets as needed, e.g.:
    # "her2_domainIV": Target(name="her2_domainIV", antigen_chain="A",
    #                         hotspots=(557, 558, 579, 583, 592, 593, 615)),
}

# Which target to use when --target is not passed on the CLI.
DEFAULT_TARGET = "example_antigen"


# ------------------------------------------------------------------------------
# 2. CHAINS — who is the binder, who is the antigen
# ------------------------------------------------------------------------------
# The antigen chain is taken from the selected Target (above). The binder is
# everything else by default; if your complexes have extra chains (e.g. a
# light + heavy chain antibody, or a scaffold), list the binder chain IDs
# explicitly here. Leave as None to mean "all chains that are not the antigen".
BINDER_CHAINS: tuple[str, ...] | None = None        # e.g. ("H", "L")


# ------------------------------------------------------------------------------
# 3. GROUPING — how designs are bucketed for comparison
# ------------------------------------------------------------------------------
# The original pipeline grouped by "campaign". Here it is a generic "group":
# it can be the design method (rfdiffusion vs bindcraft), the target, the
# backbone family — whatever you want to compare in the plots/stats.
#
# GROUP_FROM controls how the group label is derived from each design's path:
#   "parent_dir"  -> immediate parent folder name (designs/<group>/<id>/...)
#   "filename_prefix" -> text before the first '_' in the file name
#   "regex"       -> use GROUP_REGEX (must contain one capture group)
GROUP_FROM = "parent_dir"
GROUP_REGEX = r"^([A-Za-z0-9]+)_"      # only used when GROUP_FROM == "regex"

# Fixed left-to-right order + colours for groups in plots. Any group not listed
# here is appended automatically with an auto colour.
GROUP_ORDER: list[str] = []            # e.g. ["rfdiffusion", "bindcraft", "native"]
GROUP_COLORS: dict[str, str] = {}      # e.g. {"rfdiffusion": "#4C72B0"}


# ------------------------------------------------------------------------------
# 4. METRICS — names, display labels, and pass/fail thresholds
# ------------------------------------------------------------------------------
# CONFIDENCE metrics come from the folding model (AF3 / Chai-1 / Boltz).
# Higher is better for all of these.
CONF_METRICS = [
    "iptm", "ptm", "mean_plddt", "interface_plddt",
]

# PHYSICS / geometry metrics computed from the structure itself.
# `higher_better` tells the threshold + stats code which direction is good.
@dataclass
class MetricSpec:
    label: str
    higher_better: bool
    threshold: float | None = None     # pass/fail line drawn on plots; None = none

PHYSICS_METRICS: dict[str, MetricSpec] = {
    "hotspot_recovery":    MetricSpec("Hotspot recovery (fraction)", True,  0.75),
    "n_hotspot_contacts":  MetricSpec("Hotspots engaged (count)",     True,  None),
    "interface_bsa_A2":    MetricSpec("Buried interface area (Å²)",   True,  800.0),
    "n_interface_res":     MetricSpec("Interface residues",           True,  None),
    "n_hbonds":            MetricSpec("Interface H-bonds",            True,  3.0),
    "n_salt_bridges":      MetricSpec("Interface salt bridges",       True,  None),
    "n_clashes":           MetricSpec("Interface clashes",           False, 0.0),
    "sc_proxy":            MetricSpec("Shape-complementarity proxy", True,  0.55),
    "interface_hydrophobic_frac":
                           MetricSpec("Hydrophobic fraction of interface", False, None),
    "binder_radius_gyration_A":
                           MetricSpec("Binder radius of gyration (Å)", True, None),
}

# RMSD-to-reference metrics (self-consistency / stability). Lower is better.
# Populated only if you provide reference structures (see README). Left here so
# the stats/plots code keeps working if you supply them.
RMSD_METRICS = ["interface_rmsd", "binder_ca_rmsd"]
RMSD_THRESHOLDS = {"interface_rmsd": 1.0, "binder_ca_rmsd": 2.0}   # Å, Watson 2023-style

# Confidence thresholds (higher = better); used for the pass-rate summary.
CONF_THRESHOLDS = {
    "iptm":            0.80,
    "ptm":             0.70,
    "mean_plddt":      70.0,
    "interface_plddt": 70.0,
}

# Display labels for confidence + RMSD metrics (physics labels live in the
# MetricSpec above).
AXIS_LABELS = {
    "iptm":            "ipTM",
    "ptm":             "pTM",
    "mean_plddt":      "Mean pLDDT",
    "interface_plddt": "Interface pLDDT",
    "interface_rmsd":  "Interface RMSD (Å)",
    "binder_ca_rmsd":  "Binder Cα RMSD (Å)",
}


# ------------------------------------------------------------------------------
# 5. OUTPUT
# ------------------------------------------------------------------------------
FIGDIR = "figures"          # where plots are written
SCORES_CSV = "scores.csv"   # the per-design metric table
SAMPLE_N = 19000            # points per group in scatter panels; None = all
