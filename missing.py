"""
missing.py  –  Stage 2: audit study output directories.

Sorts every study folder into finishedtables/ (has a deblur table) or
missingfiles/ (missing seqs and/or table), and writes the corresponding
tracking lists.

Usage (standalone):
    python missing.py

or import and call check_missing().
"""

from __future__ import annotations

import shutil
from pathlib import Path


def check_missing() -> None:
    finished    = Path("finishedtables")
    missing_dir = Path("missingfiles")
    finished.mkdir(exist_ok=True)
    missing_dir.mkdir(exist_ok=True)

    skip = {finished.name, missing_dir.name, "data", "_sbatch_scripts"}
    no_seqs, no_tables, complete = [], [], []

    for d in Path(".").iterdir():
        if not d.is_dir() or d.name in skip:
            continue
        has_txt   = any(d.rglob("*.txt"))
        has_seqs  = any(d.rglob("*-seqs-deblur.qza"))
        has_table = any(d.rglob("*-table-deblur.qza"))

        if has_txt and not has_seqs and not has_table:
            no_seqs.append(d.name)
            shutil.move(str(d), missing_dir / d.name)
        elif has_txt and has_seqs and not has_table:
            no_tables.append(d.name)
            shutil.move(str(d), missing_dir / d.name)
        elif has_table:
            complete.append(d.name)
            shutil.move(str(d), finished / d.name)

    print(
        f"  {len(no_seqs)} missing seqs, "
        f"{len(no_tables)} missing tables, "
        f"{len(complete)} complete"
    )
    if no_seqs:
        Path("missingseqs.txt").write_text("\n".join(no_seqs))
    if no_tables:
        Path("missingdeblurtables.txt").write_text("\n".join(no_tables))
    if complete:
        (finished / "finisheddeblurtables.txt").write_text("\n".join(complete))


if __name__ == "__main__":
    check_missing()
