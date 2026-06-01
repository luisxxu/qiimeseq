"""
merge.py  –  Stage 3: merge all finished studies.

Collects every finished study's deblur table, rep-seqs, and combined
metadata, then merges them into a single feature table, rep-seqs artifact,
and metadata file (filtering the table to samples present in the metadata).

Usage (standalone):
    python merge.py

or import and call merge().
"""

from __future__ import annotations

import shutil
from pathlib import Path

import config


def merge() -> None:
    finished = Path("finishedtables")
    staging  = finished / "deblurtables"
    staging.mkdir(exist_ok=True)

    finished_studies = (finished / "finisheddeblurtables.txt").read_text().splitlines()
    for name in finished_studies:
        src = finished / name
        for pattern in ("*combined-metadata.txt", "*-table-deblur*", "*-rep-seqs*"):
            for f in src.rglob(pattern):
                shutil.copy(f, staging / f.name)

    # Merge metadata: first file keeps the header; subsequent files skip it
    merged_meta = staging / "merged_metadata.txt"
    first = True
    for f in sorted(staging.glob("*combined-metadata*")):
        lines = f.read_text().splitlines()
        if first:
            merged_meta.write_text("\n".join(lines) + "\n")
            first = False
        else:
            with merged_meta.open("a") as out:
                out.write("\n".join(lines[1:]) + "\n")

    tables = " ".join(
        f"--i-tables {f.name}" for f in sorted(staging.glob("*-table-deblur*"))
    )
    config.qiime(
        f"qiime feature-table merge {tables} --o-merged-table merged_table.qza",
        cwd=staging,
    )

    seqs = " ".join(
        f"--i-data {f.name}" for f in sorted(staging.glob("*-rep-seqs*"))
    )
    config.qiime(
        f"qiime feature-table merge-seqs {seqs} --o-merged-data merged_rep_seqs.qza",
        cwd=staging,
    )

    config.qiime(
        "qiime feature-table filter-samples "
        "--i-table merged_table.qza "
        "--m-metadata-file merged_metadata.txt "
        "--o-filtered-table merged_table_filtered.qza",
        cwd=staging,
    )

    for name in ("merged_metadata.txt", "merged_table_filtered.qza", "merged_rep_seqs.qza"):
        src = staging / name
        if src.exists():
            shutil.move(str(src), Path(".") / name)

    shutil.rmtree(staging)


if __name__ == "__main__":
    merge()
