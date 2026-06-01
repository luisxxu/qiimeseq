"""
config.py  –  Shared configuration and helpers for the qiimeseq pipeline.

Imported by downloaddeblur.py, missing.py, merge.py, greengenes.py, and
qiimeseq.py so that every stage uses the same QIIME2 environment, ENA
endpoints, and subprocess wrappers.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path

# ── QIIME2 / environment ──────────────────────────────────────────────────────

CONDA_INIT = "source ~/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV  = "qiime2-2023.2"
TMPDIR     = "/ddn_scratch/lxxu/tmp"

# ── GreenGenes2 reference ─────────────────────────────────────────────────────

GG2_BACKBONE  = "https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza"
GG2_PHYLOGENY = "https://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.id.nwk.qza"

# ── ENA portal ────────────────────────────────────────────────────────────────

ENA_BASE         = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FASTQ_FIELDS = "library_strategy,fastq_ftp"
ENA_META_FIELDS  = (
    "sample_accession,library_name,secondary_sample_accession,run_accession,"
    "experiment_accession,fastq_bytes,fastq_ftp,fastq_md5,library_source,"
    "instrument_platform,submitted_format,library_strategy,library_layout,"
    "tax_id,scientific_name,instrument_model,library_selection,center_name,"
    "experiment_title,study_title,study_alias,experiment_alias,sample_alias,"
    "sample_title,study_accession"
)

# ── Download tuning ───────────────────────────────────────────────────────────

MAX_RETRIES        = 20
PARALLEL_DOWNLOADS = 8

# ── SLURM ─────────────────────────────────────────────────────────────────────

SBATCH_LOGS = "/home/lxxu/logs"
SBATCH_MAIL = "lxxu@ucsd.edu"


# ── Helpers ───────────────────────────────────────────────────────────────────

def ena_url(accession: str, fields: str) -> str:
    """Build an ENA filereport API URL for the given accession and fields."""
    return (
        f"{ENA_BASE}?accession={accession}&result=read_run"
        f"&fields={fields}&format=tsv&download=true&limit=0"
    )


def run(cmd: str, cwd: Path | None = None, check: bool = True) -> int:
    """Run a shell command, echoing a truncated form of it first."""
    print(f"  $ {cmd[:120]}")
    r = subprocess.run(cmd, shell=True, cwd=cwd)
    if check and r.returncode != 0:
        raise RuntimeError(f"Command failed (exit {r.returncode}): {cmd[:80]}")
    return r.returncode


def qiime(cmd: str, cwd: Path | None = None) -> int:
    """Run a QIIME2 command inside the configured conda environment."""
    return run(
        f"{CONDA_INIT} && conda activate {CONDA_ENV} && export TMPDIR={TMPDIR} && {cmd}",
        cwd=cwd,
    )


def load_samples(path: str) -> tuple[list[str], list[dict[str, str]]]:
    """Load a TSV sample sheet, returning (headers, list-of-row-dicts)."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = list(reader.fieldnames or [])
        rows    = list(reader)
    return headers, rows
