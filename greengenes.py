"""
greengenes.py  –  Stage 4: GreenGenes2 mapping + core metrics.

Downloads the GreenGenes2 backbone and phylogeny, maps the merged table to
GG2, runs core-metrics-phylogenetic, and tidies all outputs into output/.

Usage (standalone):
    python greengenes.py

or import and call greengenes().
"""

from __future__ import annotations

import shutil
from pathlib import Path

import config


def greengenes() -> None:
    backbone_name  = config.GG2_BACKBONE.split("/")[-1]
    phylogeny_name = config.GG2_PHYLOGENY.split("/")[-1]

    config.run(f"wget -c {config.GG2_BACKBONE}")
    config.run(f"wget -c {config.GG2_PHYLOGENY}")

    config.qiime(
        f"qiime greengenes2 non-v4-16s "
        f"--i-table merged_table_filtered.qza "
        f"--i-sequences merged_rep_seqs.qza "
        f"--p-threads 16 "
        f"--i-backbone {backbone_name} "
        f"--o-mapped-table merged_table_gg2.qza "
        f"--o-representatives merged_seqs_gg2.qza"
    )

    config.qiime(
        f"qiime diversity core-metrics-phylogenetic "
        f"--i-phylogeny {phylogeny_name} "
        f"--i-table merged_table_gg2.qza "
        f"--p-sampling-depth 2000 "
        f"--m-metadata-file merged_metadata.txt "
        f"--p-n-jobs-or-threads auto "
        f"--output-dir core-metrics-results"
    )

    Path(backbone_name).unlink(missing_ok=True)
    Path(phylogeny_name).unlink(missing_ok=True)

    Path("merged_files").mkdir(exist_ok=True)
    for f in sorted(Path(".").glob("merged*")):
        shutil.move(str(f), f"merged_files/{f.name}")

    Path("output").mkdir(exist_ok=True)
    for src, dst in [
        ("finishedtables",       "output/projects"),
        ("missingfiles",         "output/error_projects"),
        ("core-metrics-results", "output/visualizations"),
        ("merged_files",         "output/merged_files"),
    ]:
        if Path(src).exists():
            shutil.move(src, dst)


if __name__ == "__main__":
    greengenes()
