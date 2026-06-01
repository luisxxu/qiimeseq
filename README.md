# qiimeseq

Downloading and sanitizing data for analysis is an incredibly arduous task that pleases only the most Sisyphean. To find joy in this task, one must automate the process as revolt against the human condition.

The pipeline takes a metadata sheet requiring columns named "accession" and "study", then downloads the fastq files and uses Qiime2's tools (import, quality-filter, deblur, greengenes2) to produce visualizations and other files necessary for analysis.

Special thanks to Sam, Noah, and Ananya for their contributions to this project.

Make sure to install Qiime2-2023.2 before running the pipeline.

# Qiime2-2023.2 Installation
[Qiime2 Installation Instructions](https://docs.qiime2.org/2023.2/install/)

# Layout

This is the Python rewrite. Each stage is its own module and can be run on its own, or all together via the orchestrator.

| File | Stage |
|------|-------|
| `qiimeseq.py` | Orchestrator — runs every stage in order, or submits them to SLURM with `--sbatch` |
| `downloaddeblur.py` | 1. ENA download → manifest/metadata → QIIME2 import / quality-filter / deblur |
| `missing.py` | 2. Sort finished vs. incomplete study folders |
| `merge.py` | 3. Merge feature tables, rep-seqs, and metadata |
| `greengenes.py` | 4. GreenGenes2 mapping + core-metrics-phylogenetic |
| `metadata.py` | Combine manifest + custom + ENA metadata (used by stage 1) |
| `config.py` | Shared configuration and helpers |

# Running the pipeline

Pass the location of your metadata sheet. Run everything locally:

```
python qiimeseq.py data/sample.tsv
```

Or submit each stage as a chained SLURM job (download is a job array, 4 studies at a time; later stages depend on earlier ones):

```
python qiimeseq.py data/sample.tsv --sbatch
```

The `--sbatch` flag generates one `.sbatch` wrapper per stage in `_sbatch_scripts/` and submits them with the right dependencies; without it, the same stages run directly in the current shell.

You can also run any single stage by hand, e.g.:

```
python downloaddeblur.py data/sample.tsv 2   # process the first study row
python missing.py
python merge.py
python greengenes.py
```

A small example sheet lives at `data/sample.tsv`.

Ensure that the metadata files don't have any punctuation, but especially **DO NOT** use any punctuation for the study/title (e.g. KannEtAlArgentina2021 is okay instead of Kann_et_al_Argentina_2021).

ok bye
