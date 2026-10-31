"""
antibody_scoring — a general scoring pipeline for predicted binder/antibody
structures against a configurable antigen with defined hotspot residues.

Quick start (Python):
    from antibody_scoring import config, pipeline, plots, stats
    tgt = config.TARGETS[config.DEFAULT_TARGET]
    df  = pipeline.score_directory("designs", tgt)
    df  = pipeline.add_pass_flags(df)
    plots.physics_panel(df); plots.hotspot_heatmap(df)
    print(stats.success_table(df))

Or from the shell:
    python -m antibody_scoring.run_scoring designs --target example_antigen
"""
from . import config, structure_io, confidence, interface_metrics, pipeline, plots, stats  # noqa

__all__ = ["config", "structure_io", "confidence", "interface_metrics",
           "pipeline", "plots", "stats"]
__version__ = "0.1.0"
