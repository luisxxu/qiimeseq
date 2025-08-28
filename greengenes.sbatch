#!/bin/bash -l
#SBATCH --job-name=qsqgreengenes
#SBATCH --output=/home/lxxu/logs/%x-%j.out
#SBATCH --error=/home/lxxu/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lxxu@ucsd.edu
#SBATCH --time=7-00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

wget https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
wget https://ftp.microbio.me/greengenes_release/2022.10/2022.10.phylogeny.id.nwk.qza

source ~/miniforge3/etc/profile.d/conda.sh
conda activate qiime2-2023.2

qiime greengenes2 non-v4-16s \
    --i-table merged_table_filtered.qza \
    --i-sequences merged_rep_seqs.qza \
    --p-threads 16 \
    --i-backbone 2022.10.backbone.full-length.fna.qza \
    --o-mapped-table merged_table_gg2.qza \
    --o-representatives merged_seqs_gg2.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 2022.10.phylogeny.id.nwk.qza \
  --i-table merged_table_gg2.qza \
  --p-sampling-depth 2000 \
  --m-metadata-file merged_metadata.txt \
  --p-n-jobs-or-threads auto \
  --output-dir core-metrics-results

rm 2022.10.backbone.full-length.fna.qza
rm 2022.10.phylogeny.id.nwk.qza