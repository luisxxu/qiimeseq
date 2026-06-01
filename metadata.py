"""
metadata.py  –  Combine manifest + custom metadata + ENA metadata.

Produces <name>-combined-metadata.txt by mapping each run_accession (parsed
from the manifest filepaths) to its SampleID and merging in the full ENA
metadata table.

Usage (standalone, from inside the study directory):
    python metadata.py <name>

or import and call combine(name, study_dir).
"""

from __future__ import annotations

import sys
from pathlib import Path


def combine(name: str, study_dir: Path | str = ".") -> Path:
    import pandas as pd  # imported lazily so the orchestrator works without pandas

    study_dir = Path(study_dir)

    manifest_file = study_dir / f"{name}-manifest.txt"
    custom_file   = study_dir / f"{name}-metadata.txt"
    ena_file      = study_dir / f"{name}-ena-metadata.txt"
    output_file   = study_dir / f"{name}-combined-metadata.txt"

    manifest    = pd.read_csv(manifest_file, sep="\t", header=None,
                              names=["SampleID", "filepath"])
    custom_meta = pd.read_csv(custom_file, sep="\t", keep_default_na=True)
    ena_meta    = pd.read_csv(ena_file, sep="\t")

    run_to_sample = {}
    for _, row in manifest.iterrows():
        # extract run_accession from the filepath (e.g. SRR8143822)
        run_accession = row["filepath"].split("/")[-1].split("_")[0]
        run_to_sample[run_accession] = row["SampleID"]

    ena_filtered = ena_meta[ena_meta["run_accession"].isin(run_to_sample)].copy()
    ena_filtered["SampleID"] = ena_filtered["run_accession"].map(run_to_sample)

    merged = pd.merge(
        custom_meta,
        ena_filtered.drop(columns=["filepath"], errors="ignore"),
        on="SampleID",
        how="left",
    )
    merged.to_csv(output_file, sep="\t", index=False, na_rep="NaN")
    return output_file


if __name__ == "__main__":
    combine(sys.argv[1])
