#!/usr/bin/env python3
"""
qiimeseq.py  –  Automated QIIME2 pipeline for 16S amplicon data

Stages
------
1. Per-study : download AMPLICON FASTQs from ENA, build manifest/metadata,
               run QIIME2 import → quality-filter → deblur.
2. Audit     : sort finished / incomplete study folders.
3. Merge     : combine all feature tables, rep-seqs, and metadata.
4. GreenGenes2: map to GG2 backbone, compute core-metrics-phylogenetic.

Usage
-----
    # Run everything locally:
    python qiimeseq.py data/sample.tsv

    # Submit each stage as a SLURM job (dependencies chained automatically):
    python qiimeseq.py data/sample.tsv --sbatch
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────

CONDA_INIT = "source ~/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV  = "qiime2-2023.2"
TMPDIR     = "/ddn_scratch/lxxu/tmp"

GG2_BACKBONE  = "https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza"
GG2_PHYLOGENY = "https://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.id.nwk.qza"

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

MAX_RETRIES        = 20
PARALLEL_DOWNLOADS = 8

SBATCH_LOGS = "/home/lxxu/logs"
SBATCH_MAIL = "lxxu@ucsd.edu"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _ena_url(accession: str, fields: str) -> str:
    return (
        f"{ENA_BASE}?accession={accession}&result=read_run"
        f"&fields={fields}&format=tsv&download=true&limit=0"
    )


def _run(cmd: str, cwd: Path | None = None, check: bool = True) -> int:
    print(f"  $ {cmd[:120]}")
    r = subprocess.run(cmd, shell=True, cwd=cwd)
    if check and r.returncode != 0:
        raise RuntimeError(f"Command failed (exit {r.returncode}): {cmd[:80]}")
    return r.returncode


def _qiime(cmd: str, cwd: Path | None = None) -> int:
    return _run(
        f"{CONDA_INIT} && conda activate {CONDA_ENV} && export TMPDIR={TMPDIR} && {cmd}",
        cwd=cwd,
    )


def _load_samples(path: str) -> tuple[list[str], list[dict[str, str]]]:
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = list(reader.fieldnames or [])
        rows    = list(reader)
    return headers, rows


# ── Stage 1: per-study download + deblur ──────────────────────────────────────

def _fetch_fastq_urls(accession: str, dest_dir: Path) -> list[str]:
    txt = dest_dir / f"{accession}.txt"
    urllib.request.urlretrieve(_ena_url(accession, ENA_FASTQ_FIELDS), txt)
    urls: list[str] = []
    with open(txt) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("library_strategy") == "AMPLICON":
                for u in row.get("fastq_ftp", "").split(";"):
                    if u.strip():
                        urls.append(u.strip())
    return urls


def _download_one(url: str, dest_dir: Path) -> bool:
    full_url = f"ftp://{url}" if not url.startswith("ftp://") else url
    r = subprocess.run(
        ["wget", "--tries=5", "--waitretry=5", "-c", full_url, "-P", str(dest_dir)],
        capture_output=True,
    )
    ok = r.returncode == 0 and (dest_dir / url.split("/")[-1]).exists()
    if not ok:
        print(f"    Warning: download failed for {url.split('/')[-1]}")
    return ok


def _download_all(urls: list[str], dest_dir: Path) -> None:
    url_by_name = {u.split("/")[-1]: u for u in urls}
    expected    = set(url_by_name)

    with ThreadPoolExecutor(max_workers=PARALLEL_DOWNLOADS) as ex:
        list(ex.map(lambda u: _download_one(u, dest_dir), urls))

    for attempt in range(1, MAX_RETRIES + 1):
        downloaded = {f.name for f in dest_dir.glob("*.fastq.gz")}
        missing    = expected - downloaded
        if not missing:
            print(f"  All {len(expected)} files present.")
            return
        print(f"  Retry {attempt}/{MAX_RETRIES}: {len(missing)} files still missing...")
        with ThreadPoolExecutor(max_workers=PARALLEL_DOWNLOADS) as ex:
            list(ex.map(lambda n: _download_one(url_by_name[n], dest_dir), missing))

    still = expected - {f.name for f in dest_dir.glob("*.fastq.gz")}
    if still:
        print(f"  Warning: {len(still)} files missing after {MAX_RETRIES} retries")
        (dest_dir / "missing_files.txt").write_text("\n".join(sorted(still)))


def _rename_single_end(acc_dir: Path) -> None:
    """Add _1 suffix to .fastq.gz files that lack a paired-end suffix."""
    for f in list(acc_dir.glob("*.fastq.gz")):
        stem = f.name[: -len(".fastq.gz")]
        if not (stem.endswith("_1") or stem.endswith("_2")):
            f.rename(acc_dir / f"{stem}_1.fastq.gz")


def _write_manifest_and_metadata(
    study: str,
    acc_dir: Path,
    row: dict[str, str],
    headers: list[str],
) -> None:
    study_dir = acc_dir.parent
    manifest_path = study_dir / f"{study}-manifest.txt"
    meta_path     = study_dir / f"{study}-metadata.txt"

    fastq_files = sorted(acc_dir.glob("*_1.*"))
    sample_ids  = [f"{study}-{i}" for i in range(1, len(fastq_files) + 1)]

    manifest_lines = ["sample-id\tabsolute-filepath"]
    for sid, f in zip(sample_ids, fastq_files):
        manifest_lines.append(f"{sid}\t{f.resolve()}")
    manifest_path.write_text("\n".join(manifest_lines) + "\n")

    # Metadata: first column is always SampleID
    meta_lines = ["SampleID\t" + "\t".join(headers[1:])]
    for sid in sample_ids:
        values = "\t".join(row.get(h, "") for h in headers[1:])
        meta_lines.append(f"{sid}\t{values}")
    meta_path.write_text("\n".join(meta_lines) + "\n")


def _merge_ena_metadata(study: str, study_dir: Path) -> None:
    """Merge manifest + custom metadata + ENA metadata → combined-metadata.txt."""
    try:
        import pandas as pd
    except ImportError:
        raise RuntimeError("pandas is required: pip install pandas")

    manifest    = pd.read_csv(study_dir / f"{study}-manifest.txt",
                              sep="\t", header=None, names=["SampleID", "filepath"])
    custom_meta = pd.read_csv(study_dir / f"{study}-metadata.txt",
                              sep="\t", keep_default_na=True)
    ena_meta    = pd.read_csv(study_dir / f"{study}-ena-metadata.txt", sep="\t")

    run_to_sample = {
        row["filepath"].split("/")[-1].split("_")[0]: row["SampleID"]
        for _, row in manifest.iterrows()
    }

    ena_filtered = ena_meta[ena_meta["run_accession"].isin(run_to_sample)].copy()
    ena_filtered["SampleID"] = ena_filtered["run_accession"].map(run_to_sample)

    merged = pd.merge(
        custom_meta,
        ena_filtered.drop(columns=["filepath"], errors="ignore"),
        on="SampleID",
        how="left",
    )
    merged.to_csv(study_dir / f"{study}-combined-metadata.txt",
                  sep="\t", index=False, na_rep="NaN")


def _qiime_import(study: str, study_dir: Path) -> None:
    out = study_dir / f"{study}-seqs.qza"
    cmd = (
        f"qiime tools import "
        f"--type 'SampleData[SequencesWithQuality]' "
        f"--input-path {study}-manifest.txt "
        f"--output-path {study}-seqs.qza "
        f"--input-format SingleEndFastqManifestPhred33V2"
    )
    for attempt in range(1, 4):
        _qiime(cmd, cwd=study_dir)
        if out.exists():
            return
        print(f"  Import attempt {attempt} failed, retrying...")
    raise RuntimeError(f"qiime import failed for {study}")


def _qiime_quality_filter(study: str, study_dir: Path) -> None:
    out = study_dir / f"{study}-demux-filtered.qza"
    cmd = (
        f"qiime quality-filter q-score "
        f"--i-demux {study}-seqs.qza "
        f"--o-filtered-sequences {study}-demux-filtered.qza "
        f"--o-filter-stats {study}-demux-filter-stats.qza"
    )
    for attempt in range(1, 4):
        _qiime(cmd, cwd=study_dir)
        if out.exists():
            return
        print(f"  Quality-filter attempt {attempt} failed, retrying...")
    raise RuntimeError(f"quality-filter failed for {study}")


def _qiime_deblur(study: str, study_dir: Path) -> None:
    rep_seqs = study_dir / f"{study}-rep-seqs-deblur.qza"
    table    = study_dir / f"{study}-table-deblur.qza"
    cmd = (
        f"qiime deblur denoise-16S "
        f"--i-demultiplexed-seqs {study}-demux-filtered.qza "
        f"--p-trim-length 150 "
        f"--p-left-trim-len 0 "
        f"--p-jobs-to-start 8 "
        f"--o-representative-sequences {study}-rep-seqs-deblur.qza "
        f"--o-table {study}-table-deblur.qza "
        f"--p-sample-stats "
        f"--o-stats {study}-deblur-stats.qza "
        f"--verbose"
    )
    for attempt in range(1, 4):
        _qiime(cmd, cwd=study_dir)
        if rep_seqs.exists() and table.exists():
            return
        print(f"  Deblur attempt {attempt} failed, retrying...")
    raise RuntimeError(f"deblur failed for {study}")


def process_study(row: dict[str, str], headers: list[str]) -> None:
    """Run the full per-study pipeline: download → manifest → QIIME2."""
    study     = row["study"].strip()
    accession = row["accession"].strip()

    print(f"\n=== {study}  ({accession}) ===")

    study_dir = Path(study)
    study_dir.mkdir(exist_ok=True)
    acc_dir = study_dir / accession
    acc_dir.mkdir(exist_ok=True)

    print("  Fetching FASTQ URL list from ENA...")
    urls = _fetch_fastq_urls(accession, acc_dir)
    if not urls:
        print(f"  No AMPLICON files found for {accession}; skipping.")
        return

    print(f"  Downloading {len(urls)} files...")
    _download_all(urls, acc_dir)

    _rename_single_end(acc_dir)
    _write_manifest_and_metadata(study, acc_dir, row, headers)

    print("  Fetching full ENA metadata...")
    urllib.request.urlretrieve(
        _ena_url(accession, ENA_META_FIELDS),
        study_dir / f"{study}-ena-metadata.txt",
    )
    _merge_ena_metadata(study, study_dir)

    print("  QIIME2: import...")
    _qiime_import(study, study_dir)

    print("  QIIME2: quality-filter...")
    _qiime_quality_filter(study, study_dir)

    print("  QIIME2: deblur...")
    _qiime_deblur(study, study_dir)

    print(f"  Done: {study}")


# ── Stage 2: audit ────────────────────────────────────────────────────────────

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


# ── Stage 3: merge ────────────────────────────────────────────────────────────

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
    _qiime(f"qiime feature-table merge {tables} --o-merged-table merged_table.qza",
           cwd=staging)

    seqs = " ".join(
        f"--i-data {f.name}" for f in sorted(staging.glob("*-rep-seqs*"))
    )
    _qiime(f"qiime feature-table merge-seqs {seqs} --o-merged-data merged_rep_seqs.qza",
           cwd=staging)

    _qiime(
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


# ── Stage 4: GreenGenes2 ──────────────────────────────────────────────────────

def greengenes() -> None:
    backbone_name  = GG2_BACKBONE.split("/")[-1]
    phylogeny_name = GG2_PHYLOGENY.split("/")[-1]

    _run(f"wget -c {GG2_BACKBONE}")
    _run(f"wget -c {GG2_PHYLOGENY}")

    _qiime(
        f"qiime greengenes2 non-v4-16s "
        f"--i-table merged_table_filtered.qza "
        f"--i-sequences merged_rep_seqs.qza "
        f"--p-threads 16 "
        f"--i-backbone {backbone_name} "
        f"--o-mapped-table merged_table_gg2.qza "
        f"--o-representatives merged_seqs_gg2.qza"
    )

    _qiime(
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


# ── SBATCH mode ───────────────────────────────────────────────────────────────

def _sbatch_header(job_name: str, time: str, mem: str, cpus: int, array: str = "") -> str:
    array_line = f"#SBATCH --array={array}\n" if array else ""
    return (
        f"#!/bin/bash -l\n"
        f"#SBATCH --job-name={job_name}\n"
        f"#SBATCH --output={SBATCH_LOGS}/%x-%j.out\n"
        f"#SBATCH --error={SBATCH_LOGS}/%x-%j.err\n"
        f"#SBATCH --mail-type=ALL\n"
        f"#SBATCH --mail-user={SBATCH_MAIL}\n"
        f"#SBATCH --time={time}\n"
        f"#SBATCH --mem={mem}\n"
        f"#SBATCH --cpus-per-task={cpus}\n"
        f"{array_line}\n"
    )


def submit_sbatch(n_studies: int, sample_file: str) -> None:
    script  = Path(__file__).resolve()
    sfile   = Path(sample_file).resolve()
    scripts = Path("_sbatch_scripts")
    scripts.mkdir(exist_ok=True)

    def write_script(name: str, content: str) -> Path:
        p = scripts / name
        p.write_text(content)
        return p

    # Job 1: array over all studies (4 at a time)
    j1_script = write_script("downloaddeblur.sbatch",
        _sbatch_header("qsq_download", "10:00:00", "32G", 8,
                       array=f"2-{n_studies + 1}%4") +
        f"python3 {script} {sfile} --_run-study $SLURM_ARRAY_TASK_ID\n"
    )
    # Job 2: audit missing
    j2_script = write_script("missing.sbatch",
        _sbatch_header("qsq_missing", "1:00:00", "16G", 4) +
        f"python3 {script} {sfile} --_run-missing\n"
    )
    # Job 3: merge
    j3_script = write_script("merge.sbatch",
        _sbatch_header("qsq_merge", "10:00:00", "16G", 4) +
        f"python3 {script} {sfile} --_run-merge\n"
    )
    # Job 4: greengenes
    j4_script = write_script("greengenes.sbatch",
        _sbatch_header("qsq_greengenes", "7-00", "32G", 16) +
        f"python3 {script} {sfile} --_run-greengenes\n"
    )

    def sbatch(*extra: str) -> str:
        r = subprocess.run(["sbatch", *extra], capture_output=True, text=True, check=True)
        return r.stdout.strip().split()[-1]

    j1 = sbatch(str(j1_script))
    print(f"Job 1 (download+deblur array): {j1}")
    j2 = sbatch(f"--dependency=afterany:{j1}", str(j2_script))
    print(f"Job 2 (audit missing):         {j2}")
    j3 = sbatch(f"--dependency=afterok:{j2}", str(j3_script))
    print(f"Job 3 (merge):                 {j3}")
    j4 = sbatch(f"--dependency=afterok:{j3}", str(j4_script))
    print(f"Job 4 (greengenes):            {j4}")
    print("All jobs submitted. Good luck.")


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

    # Internal flags used by the generated sbatch scripts — not shown in help
    parser.add_argument("--_run-study",      type=int, dest="run_study",
                        metavar="TASK_ID",   help=argparse.SUPPRESS)
    parser.add_argument("--_run-missing",    action="store_true", dest="run_missing",
                        help=argparse.SUPPRESS)
    parser.add_argument("--_run-merge",      action="store_true", dest="run_merge",
                        help=argparse.SUPPRESS)
    parser.add_argument("--_run-greengenes", action="store_true", dest="run_greengenes",
                        help=argparse.SUPPRESS)

    args = parser.parse_args()
    headers, samples = _load_samples(args.sample_file)

    # ── Worker dispatch (called by sbatch scripts) ───────────────────────────
    if args.run_study is not None:
        # SLURM task IDs start at 2 (row 1 is the header)
        process_study(samples[args.run_study - 2], headers)
        return
    if args.run_missing:
        check_missing()
        return
    if args.run_merge:
        merge()
        return
    if args.run_greengenes:
        greengenes()
        return

    # ── Orchestration ────────────────────────────────────────────────────────
    if args.sbatch:
        submit_sbatch(len(samples), args.sample_file)
        return

    print(f"Running {len(samples)} studies locally...\n")
    for row in samples:
        process_study(row, headers)

    print("\n-- Auditing outputs --")
    check_missing()

    print("\n-- Merging --")
    merge()

    print("\n-- GreenGenes2 + core metrics --")
    greengenes()

    print("\nPipeline complete.")


if __name__ == "__main__":
    main()
