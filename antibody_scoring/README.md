# antibody_scoring

A general scoring pipeline for **predicted binder / antibody structures** against a
configurable antigen with defined **hotspot (epitope) residues**.

Point it at a folder of predicted complexes (PDB or mmCIF, optionally gzipped),
tell it which antigen residues are the hotspots you want engaged, and it writes
a per-design metric table, a pass/fail summary, cross-group statistics, and a
set of figures.

It is a generalisation of a structure-scoring notebook originally built for
metal-binding protein design: the grouping, the "site" definition, the metrics,
and the thresholds are now all **data in `config.py`**, so you retarget it
without touching the analysis code.

---

## 1. The one thing you'll change: hotspot targets

Open **`antibody_scoring/config.py`** and edit the `TARGETS` block:

```python
TARGETS = {
    "her2_domainIV": Target(
        name="her2_domainIV",
        antigen_chain="A",                       # chain ID of the antigen
        hotspots=(557, 558, 579, 583, 592, 615), # <-- your epitope residues
        avoid=(),                                # residues to steer away from
        contact_cutoff_A=5.0,                    # "contact" = heavy atom within 5 Å
    ),
}
DEFAULT_TARGET = "her2_domainIV"
```

- `hotspots` are **author/PDB residue numbers on the antigen chain**, exactly as
  they appear in your structure files. Open one structure in a viewer and use the
  numbers you see there.
- Keep several targets in the dict and pick one per run with `--target NAME`.
- A hotspot counts as *engaged* when any binder heavy atom is within
  `contact_cutoff_A` of any heavy atom of that residue. The pipeline reports, per
  design, how many hotspots are engaged (**hotspot recovery**) and the nearest
  binder distance to each one.

That's the whole retargeting step. If your complexes have named antibody chains
(e.g. heavy **H** + light **L**), also set `BINDER_CHAINS = ("H", "L")`;
otherwise the binder is "every chain that isn't the antigen".

---

## 2. Install

```bash
pip install -r requirements.txt
```

Only standard scientific packages plus **biotite** (structure IO, Shrake–Rupley
SASA, H-bond geometry). No Rosetta, PyRosetta, CCP4, or any external binary is
required.

---

## 3. Arrange your designs

One structure file = one binder–antigen complex. Group designs into subfolders;
the subfolder name becomes the **group** used to compare sets in the plots and
statistics (change this with `GROUP_FROM` in config):

```
designs/
├── rfdiffusion/
│   ├── design_0001.pdb
│   ├── design_0001_confidence.json      # optional, see below
│   └── design_0002.cif.gz
├── bindcraft/
│   └── ...
└── native/                              # e.g. a positive-control set
    └── ...
```

**Confidence files (optional but recommended).** If a folding-model JSON sits
next to a structure (AlphaFold3 `*_summary_confidences.json`, Chai-1, or Boltz
output), the pipeline reads `pTM`, `ipTM` and the PAE matrix from it and computes
interface/hotspot PAE. Per-residue **pLDDT is read from the B-factor column** of
the structure, which is where all three tools write it. No JSON → those columns
are left blank, everything else still runs.

---

## 4. Run

```bash
python -m antibody_scoring.run_scoring designs --target her2_domainIV --out scoring_output
```

Outputs under `scoring_output/`:

| file | contents |
|------|----------|
| `scores.csv` | one row per design, every metric + pass flags |
| `hotspot_contacts.csv` | per-design, per-hotspot nearest distance + engaged flag |
| `summary.txt` | pass-rate table + cross-group Kruskal–Wallis / Mann–Whitney |
| `figures/physics_panel.png` | interface-metric distributions by group |
| `figures/hotspot_heatmap.png` | which epitope residues each group engages |
| `figures/pairplot_<group>.png` | confidence-metric pair plots per group |

A bad structure file never kills the run — it gets an `error` note in `scores.csv`
and is skipped.

### From Python / a notebook

```python
from antibody_scoring import config, pipeline, plots, stats

tgt = config.TARGETS["her2_domainIV"]
df  = pipeline.score_directory("designs", tgt)
df  = pipeline.add_pass_flags(df)

plots.physics_panel(df)
plots.hotspot_heatmap(df)
print(stats.success_table(df))                 # pass rates, with bootstrap CIs
info, pairwise = stats.compare_groups(df, "interface_bsa_A2")
```

`df` is a plain pandas DataFrame — filter, sort, and rank however you like
(`df.sort_values("hotspot_recovery", ascending=False)`).

---

## 5. Metrics

### Confidence (from the folding model — higher is better)
| column | meaning |
|--------|---------|
| `iptm` | interface pTM — the single most useful binder-quality signal |
| `ptm` | global pTM |
| `mean_plddt` | mean pLDDT over the binder |
| `interface_plddt` | mean pLDDT over binder **interface** residues |
| `pae_interface` / `pae_hotspot` | mean PAE over binder↔antigen / binder↔hotspot pairs (needs PAE matrix) |

### Physics / geometry (computed from coordinates)
| column | meaning | direction |
|--------|---------|-----------|
| `hotspot_recovery` | fraction of configured hotspots engaged | higher |
| `n_hotspot_contacts` | count of hotspots engaged | higher |
| `avoid_contacts` | count of `avoid` residues touched | lower (want 0) |
| `interface_bsa_A2` | buried surface area, both sides (Shrake–Rupley ΔSASA) | higher |
| `n_interface_res` | residues with a heavy atom across the interface | — |
| `n_hbonds` | cross-interface hydrogen bonds (geometric) | higher |
| `n_salt_bridges` | cross-interface charged pairs < 4 Å | higher |
| `n_clashes` | cross-interface heavy-atom pairs < 2.2 Å | lower (want 0) |
| `sc_proxy` | shape-complementarity proxy in [0, 1] (see note) | higher |
| `interface_hydrophobic_frac` | hydrophobic share of buried area | context (very high → sticky) |
| `binder_radius_gyration_A` | binder compactness sanity check | — |

**Hydrogens:** `n_hbonds` needs hydrogens in the structure. AF3/Chai/Boltz
outputs are heavy-atom only, so either protonate first (`reduce`, `pdb2pqr`,
OpenMM) or treat `n_hbonds` as "0 = unknown" and lean on BSA + salt bridges +
`sc_proxy`.

**`sc_proxy` is a proxy, not true Sc.** Real Lawrence–Colman shape
complementarity needs the molecular surface (CCP4 `sc` or Rosetta). This proxy —
the fraction of interface binder atoms packed in a snug 3–5 Å shell rather than
leaving a gap — correlates with Sc and is fine for ranking. If you compute real
Sc elsewhere, drop it into the `sc_proxy` column and the plots/stats pick it up
unchanged.

### Ideas for stronger physics metrics (drop-in slots)
The architecture makes these easy to add in `interface_metrics.py` (add a key to
`PHYSICS_METRICS` in config and a line in `score_structure`):

- **ΔΔG / interface energy** via PyRosetta `InterfaceAnalyzerMover` or a
  FoldX/OpenMM MM-GBSA single point — the standard "is this interface actually
  favourable" number.
- **SAP score** (spatial aggregation propensity) and **net interface charge** to
  flag developability / non-specific binders.
- **Unsatisfied buried polar atoms** (buried donors/acceptors with no partner) —
  a classic false-positive filter for designed interfaces.
- **Rosetta `ddg` + `sc` + `interface_buried_sasa`** together, if you already run
  Rosetta — map each to a column and the rest of the pipeline is unchanged.
- **Self-consistency RMSD**: refold each design's binder sequence alone and
  measure Cα RMSD to the designed pose (`binder_ca_rmsd`), and interface RMSD to
  a reference (`interface_rmsd`). The `RMSD_METRICS` columns and their thresholds
  already exist — fill them in and the stats/plots use them.

---

## 6. Thresholds & pass/fail

Every threshold lives in `config.py` (`CONF_THRESHOLDS`, and the `threshold`
field of each `MetricSpec` in `PHYSICS_METRICS`). `add_pass_flags()` adds a
boolean `pass_<metric>` column per gated metric and an overall `passes` column
(AND of all gates). The defaults are deliberately generic starting points —
calibrate them against a positive control set (known binders / native complexes)
before trusting the pass rate.

---

## 7. Layout

```
antibody_scoring/
├── README.md
├── requirements.txt
└── antibody_scoring/
    ├── config.py            ← EDIT THIS (targets, chains, metrics, thresholds)
    ├── structure_io.py      load PDB/mmCIF, split binder/antigen chains
    ├── confidence.py        parse AF3 / Chai-1 / Boltz pTM/ipTM/PAE
    ├── interface_metrics.py physics / geometry metrics for the interface
    ├── pipeline.py          scan a folder → scores DataFrame
    ├── plots.py             pair plots, physics panel, hotspot heatmap
    ├── stats.py             Kruskal–Wallis, Mann–Whitney, Cliff's δ, bootstrap
    └── run_scoring.py       command-line entry point
```

## Assumptions worth checking
- **Residue numbering** of your hotspots matches the numbering in the structure
  files.
- **PAE indexing** is assumed to be binder-residues-then-antigen-residues in file
  order (the common AF3/Boltz convention). If your tool differs, adjust
  `_pae_summaries` in `interface_metrics.py`; interface PAE is left blank rather
  than guessed when shapes don't line up.
- For antibodies, hotspots are set on the **antigen**; the binder can be single-
  chain (nanobody/scFv) or multi-chain (Fv) via `BINDER_CHAINS`.
