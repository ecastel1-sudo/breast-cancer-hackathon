"""
Command-line entry point.

    python -m antibody_scoring.run_scoring DESIGNS_DIR [--target NAME] [--out DIR]

Writes:
    <out>/scores.csv            per-design metric table
    <out>/hotspot_contacts.csv  per-design per-hotspot distances
    <out>/figures/*.png         physics panel, hotspot heatmap, pair-plots
    <out>/summary.txt           pass-rate + cross-group stats
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

from . import config as C
from . import pipeline, plots, stats


def main(argv=None):
    ap = argparse.ArgumentParser(description="Score predicted binder/antibody structures.")
    ap.add_argument("designs", help="directory of structures (searched recursively)")
    ap.add_argument("--target", default=C.DEFAULT_TARGET,
                    help=f"target name from config.TARGETS (default: {C.DEFAULT_TARGET})")
    ap.add_argument("--out", default="scoring_output", help="output directory")
    ap.add_argument("--no-plots", action="store_true", help="skip figures")
    args = ap.parse_args(argv)

    if args.target not in C.TARGETS:
        sys.exit(f"Unknown target '{args.target}'. Available: {list(C.TARGETS)}")
    target = C.TARGETS[args.target]
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    figdir = out / C.FIGDIR

    print(f"Target '{target.name}': antigen chain {target.antigen_chain}, "
          f"{len(target.hotspots)} hotspots {target.hotspots}")
    print(f"Scoring structures under {args.designs} ...")
    df = pipeline.score_directory(args.designs, target)
    df = pipeline.add_pass_flags(df)

    df.drop(columns=["path"], errors="ignore").to_csv(out / C.SCORES_CSV, index=False)
    hs = df.attrs.get("hotspots")
    if hs is not None and not hs.empty:
        hs.to_csv(out / "hotspot_contacts.csv", index=False)
    print(f"Wrote {out / C.SCORES_CSV}  ({len(df)} designs, "
          f"{int(df.get('passes', False).sum()) if 'passes' in df else 0} passing all gates)")

    with open(out / "summary.txt", "w") as fh:
        fh.write(f"Target: {target.name}  (antigen {target.antigen_chain}, "
                 f"hotspots {target.hotspots})\n\n")
        fh.write("=== pass-rate by metric & group ===\n")
        fh.write(stats.success_table(df).to_string() + "\n\n")
        for m in [x for x in C.PHYSICS_METRICS if x in df and df[x].notna().any()]:
            info, pw = stats.compare_groups(df, m)
            fh.write(f"--- {info.get('metric', m)} ---\n")
            if "kruskal_p" in info:
                fh.write(f"Kruskal–Wallis H={info['kruskal_H']}  p={info['kruskal_p']:.3g}\n")
                if not pw.empty:
                    fh.write(pw.to_string(index=False) + "\n")
            fh.write("\n")
    print(f"Wrote {out / 'summary.txt'}")

    if not args.no_plots:
        for g in plots.group_order(df):
            plots.pairplot_group(df, g, figdir=figdir)
        plots.physics_panel(df, figdir=figdir)
        plots.hotspot_heatmap(df, figdir=figdir)
        print(f"Wrote figures to {figdir}")


if __name__ == "__main__":
    main()
