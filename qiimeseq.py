#!/usr/bin/env python3
"""
qiimeseq.py  –  Orchestrator for the QIIME2 16S amplicon pipeline.

Stages (each lives in its own module, runnable standalone too):
  1. downloaddeblur.py : ENA download → manifest/metadata → import/filter/deblur
  2. missing.py        : sort finished / incomplete study folders
  3. merge.py          : merge tables, rep-seqs, and metadata
  4. greengenes.py     : GreenGenes2 mapping + core-metrics-phylogenetic

Usage
-----
    # Run every stage locally, in order:
    python qiimeseq.py data/sample.tsv

    # Submit each stage as a chained SLURM job:
    python qiimeseq.py data/sample.tsv --sbatch
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import config
import downloaddeblur
import missing
import merge as merge_stage
import greengenes


# ── SBATCH mode ───────────────────────────────────────────────────────────────

def _sbatch_header(job_name: str, time: str, mem: str, cpus: int, array: str = "") -> str:
    array_line = f"#SBATCH --array={array}\n" if array else ""
    return (
        f"#!/bin/bash -l\n"
        f"#SBATCH --job-name={job_name}\n"
        f"#SBATCH --output={config.SBATCH_LOGS}/%x-%j.out\n"
        f"#SBATCH --error={config.SBATCH_LOGS}/%x-%j.err\n"
        f"#SBATCH --mail-type=ALL\n"
        f"#SBATCH --mail-user={config.SBATCH_MAIL}\n"
        f"#SBATCH --time={time}\n"
        f"#SBATCH --mem={mem}\n"
        f"#SBATCH --cpus-per-task={cpus}\n"
        f"{array_line}\n"
    )


def submit_sbatch(n_studies: int, sample_file: str) -> None:
    """Generate one .sbatch wrapper per stage and submit them chained."""
    here    = Path(__file__).resolve().parent
    sfile   = Path(sample_file).resolve()
    scripts = Path("_sbatch_scripts")
    scripts.mkdir(exist_ok=True)

    def write_script(name: str, header: str, body: str) -> Path:
        p = scripts / name
        p.write_text(header + body + "\n")
        return p

    # Each sbatch script just cd's to the repo and runs the stage module.
    j1 = write_script(
        "downloaddeblur.sbatch",
        _sbatch_header("qsq_download", "10:00:00", "32G", 8, array=f"2-{n_studies + 1}%4"),
        f"cd {here}\n"
        f"python3 downloaddeblur.py {sfile} $SLURM_ARRAY_TASK_ID",
    )
    j2 = write_script(
        "missing.sbatch",
        _sbatch_header("qsq_missing", "1:00:00", "16G", 4),
        f"cd {here}\npython3 missing.py",
    )
    j3 = write_script(
        "merge.sbatch",
        _sbatch_header("qsq_merge", "10:00:00", "16G", 4),
        f"cd {here}\npython3 merge.py",
    )
    j4 = write_script(
        "greengenes.sbatch",
        _sbatch_header("qsq_greengenes", "7-00", "32G", 16),
        f"cd {here}\npython3 greengenes.py",
    )

    def sbatch(*extra: str) -> str:
        r = subprocess.run(["sbatch", *extra], capture_output=True, text=True, check=True)
        return r.stdout.strip().split()[-1]

    id1 = sbatch(str(j1))
    print(f"Job 1 (download+deblur array): {id1}")
    id2 = sbatch(f"--dependency=afterany:{id1}", str(j2))
    print(f"Job 2 (audit missing):         {id2}")
    id3 = sbatch(f"--dependency=afterok:{id2}", str(j3))
    print(f"Job 3 (merge):                 {id3}")
    id4 = sbatch(f"--dependency=afterok:{id3}", str(j4))
    print(f"Job 4 (greengenes):            {id4}")
    print("All jobs submitted. Good luck.")


# ── Local mode ────────────────────────────────────────────────────────────────

def run_local(sample_file: str) -> None:
    headers, samples = config.load_samples(sample_file)

    print(f"Running {len(samples)} studies locally...\n")
    for row in samples:
        downloaddeblur.process_study(row, headers)

    print("\n-- Auditing outputs --")
    missing.check_missing()

    print("\n-- Merging --")
    merge_stage.merge()

    print("\n-- GreenGenes2 + core metrics --")
    greengenes.greengenes()

    print("\nPipeline complete.")


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="QIIME2 16S amplicon pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python qiimeseq.py data/sample.tsv\n"
            "  python qiimeseq.py data/sample.tsv --sbatch"
        ),
    )
    parser.add_argument("sample_file",
                        help="TSV with at least 'study' and 'accession' columns")
    parser.add_argument("--sbatch", action="store_true",
                        help="Submit pipeline stages as SLURM jobs instead of running locally")

    args = parser.parse_args()

    if args.sbatch:
        _, samples = config.load_samples(args.sample_file)
        submit_sbatch(len(samples), args.sample_file)
    else:
        run_local(args.sample_file)


if __name__ == "__main__":
    main()
