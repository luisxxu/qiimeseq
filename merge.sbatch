#!/bin/bash
#SBATCH --job-name=qsqmerge
#SBATCH --output=/home/lxxu/logs/%x-%j.out
#SBATCH --error=/home/lxxu/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lxxu@ucsd.edu
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

cd finishedtables
mkdir deblurtables

while IFS= read -r line; do
    line=$(echo $line | tr -d './')
    cp ${line}/*combined-metadata.txt ${line}/*-table-deblur* ${line}/*-rep-seqs* deblurtables
done < "finisheddeblurtables.txt"

cd "deblurtables"

for file in *combined-metadata*; do
    if [[ -f "merged_metadata.txt" ]]; then
        tail -n +2 $file >> "merged_metadata.txt" 
    else
        cat $file > "merged_metadata.txt"
    fi
done

source ~/miniforge3/etc/profile.d/conda.sh
conda activate qiime2-2023.2

qiime feature-table merge \
    --i-tables *-table* \
    --o-merged-table merged_table.qza

qiime feature-table merge-seqs \
  --i-data *-rep-seqs* \
  --o-merged-data merged_rep_seqs.qza

qiime feature-table filter-samples \
  --i-table merged_table.qza \
  --m-metadata-file merged_metadata.txt \
  --o-filtered-table merged_table_filtered.qza

mv merged_metadata.txt merged_table_filtered.qza merged_rep_seqs.qza ..
cd ..
rm -rf deblurtables
mv merged_metadata.txt merged_table_filtered.qza merged_rep_seqs.qza ..
cd ..
