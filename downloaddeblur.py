"""
downloaddeblur.py  –  Stage 1: per-study download + deblur.

For one study: download its AMPLICON FASTQs from ENA, build the
manifest/metadata files, then run QIIME2 import → quality-filter → deblur.

Usage (standalone, equivalent to one SLURM array task):
    python downloaddeblur.py <sample.tsv> <task_id>

where <task_id> follows the SLURM convention (2 = first study row, since
row 1 is the header). Or import and call process_study(row, headers).
"""

from __future__ import annotations

import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import config
import metadata


# ── Download ──────────────────────────────────────────────────────────────────

def fetch_fastq_urls(accession: str, dest_dir: Path) -> list[str]:
    """Query ENA and return the AMPLICON FASTQ FTP URLs for an accession."""
    import csv

    txt = dest_dir / f"{accession}.txt"
    urllib.request.urlretrieve(config.ena_url(accession, config.ENA_FASTQ_FIELDS), txt)

    urls: list[str] = []
    with open(txt) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("library_strategy") == "AMPLICON":
                for u in row.get("fastq_ftp", "").split(";"):
                    if u.strip():
                        urls.append(u.strip())
    return urls


def download_one(url: str, dest_dir: Path) -> bool:
    import subprocess

    full_url = f"ftp://{url}" if not url.startswith("ftp://") else url
    r = subprocess.run(
        ["wget", "--tries=5", "--waitretry=5", "-c", full_url, "-P", str(dest_dir)],
        capture_output=True,
    )
    ok = r.returncode == 0 and (dest_dir / url.split("/")[-1]).exists()
    if not ok:
        print(f"    Warning: download failed for {url.split('/')[-1]}")
    return ok


def download_all(urls: list[str], dest_dir: Path) -> None:
    """Download all URLs in parallel, retrying any that go missing."""
    url_by_name = {u.split("/")[-1]: u for u in urls}
    expected    = set(url_by_name)

    with ThreadPoolExecutor(max_workers=config.PARALLEL_DOWNLOADS) as ex:
        list(ex.map(lambda u: download_one(u, dest_dir), urls))

    for attempt in range(1, config.MAX_RETRIES + 1):
        downloaded = {f.name for f in dest_dir.glob("*.fastq.gz")}
        missing    = expected - downloaded
        if not missing:
            print(f"  All {len(expected)} files present.")
            return
        print(f"  Retry {attempt}/{config.MAX_RETRIES}: {len(missing)} files still missing...")
        with ThreadPoolExecutor(max_workers=config.PARALLEL_DOWNLOADS) as ex:
            list(ex.map(lambda n: download_one(url_by_name[n], dest_dir), missing))

    still = expected - {f.name for f in dest_dir.glob("*.fastq.gz")}
    if still:
        print(f"  Warning: {len(still)} files missing after {config.MAX_RETRIES} retries")
        (dest_dir / "missing_files.txt").write_text("\n".join(sorted(still)))


# ── Manifest / metadata ───────────────────────────────────────────────────────

def rename_single_end(acc_dir: Path) -> None:
    """Add _1 suffix to .fastq.gz files that lack a paired-end suffix."""
    for f in list(acc_dir.glob("*.fastq.gz")):
        stem = f.name[: -len(".fastq.gz")]
        if not (stem.endswith("_1") or stem.endswith("_2")):
            f.rename(acc_dir / f"{stem}_1.fastq.gz")


def write_manifest_and_metadata(
    study: str,
    acc_dir: Path,
    row: dict[str, str],
    headers: list[str],
) -> None:
    study_dir     = acc_dir.parent
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


# ── QIIME2 steps ──────────────────────────────────────────────────────────────

def qiime_import(study: str, study_dir: Path) -> None:
    out = study_dir / f"{study}-seqs.qza"
    cmd = (
        f"qiime tools import "
        f"--type 'SampleData[SequencesWithQuality]' "
        f"--input-path {study}-manifest.txt "
        f"--output-path {study}-seqs.qza "
        f"--input-format SingleEndFastqManifestPhred33V2"
    )
    for attempt in range(1, 4):
        config.qiime(cmd, cwd=study_dir)
        if out.exists():
            return
        print(f"  Import attempt {attempt} failed, retrying...")
    raise RuntimeError(f"qiime import failed for {study}")


def qiime_quality_filter(study: str, study_dir: Path) -> None:
    out = study_dir / f"{study}-demux-filtered.qza"
    cmd = (
        f"qiime quality-filter q-score "
        f"--i-demux {study}-seqs.qza "
        f"--o-filtered-sequences {study}-demux-filtered.qza "
        f"--o-filter-stats {study}-demux-filter-stats.qza"
    )
    for attempt in range(1, 4):
        config.qiime(cmd, cwd=study_dir)
        if out.exists():
            return
        print(f"  Quality-filter attempt {attempt} failed, retrying...")
    raise RuntimeError(f"quality-filter failed for {study}")


def qiime_deblur(study: str, study_dir: Path) -> None:
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
        config.qiime(cmd, cwd=study_dir)
        if rep_seqs.exists() and table.exists():
            return
        print(f"  Deblur attempt {attempt} failed, retrying...")
    raise RuntimeError(f"deblur failed for {study}")


# ── Orchestration for one study ───────────────────────────────────────────────

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
    urls = fetch_fastq_urls(accession, acc_dir)
    if not urls:
        print(f"  No AMPLICON files found for {accession}; skipping.")
        return

    print(f"  Downloading {len(urls)} files...")
    download_all(urls, acc_dir)

    rename_single_end(acc_dir)
    write_manifest_and_metadata(study, acc_dir, row, headers)

    print("  Fetching full ENA metadata...")
    urllib.request.urlretrieve(
        config.ena_url(accession, config.ENA_META_FIELDS),
        study_dir / f"{study}-ena-metadata.txt",
    )
    metadata.combine(study, study_dir)

    print("  QIIME2: import...")
    qiime_import(study, study_dir)

    print("  QIIME2: quality-filter...")
    qiime_quality_filter(study, study_dir)

    print("  QIIME2: deblur...")
    qiime_deblur(study, study_dir)

    print(f"  Done: {study}")


if __name__ == "__main__":
    sample_file = sys.argv[1]
    task_id     = int(sys.argv[2])
    headers, samples = config.load_samples(sample_file)
    # SLURM task IDs start at 2 (row 1 is the header)
    process_study(samples[task_id - 2], headers)
