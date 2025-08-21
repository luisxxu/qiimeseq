#!/bin/bash
#SBATCH --job-name=qsqmissing
#SBATCH --output=/home/lxxu/logs/%x-%j.out
#SBATCH --error=/home/lxxu/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lxxu@ucsd.edu
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

mkdir finishedtables
mkdir missingfiles

find . -mindepth 1 -maxdepth 1 -type d | while read dir; do
    has_txt=$(find "$dir" -type f -name "*.txt" | head -n 1)
    has_seqs=$(find "$dir" -type f -name "*-seqs-deblur.qza" | head -n 1)
    has_table=$(find "$dir" -type f -name "*-table-deblur.qza" | head -n 1)

    if [[ -n "$has_txt" && -z "$has_seqs" && -z "$has_table" ]]; then
        echo "$dir" >> missingseqs.txt
        mv $dir missingfiles
    fi
    if [[ -n "$has_txt" && -n "$has_seqs" && -z "$has_table" ]]; then
        echo "$dir" >> missingdeblurtables.txt
        mv $dir missingfiles
    fi
    if [[ -n "$has_table" ]]; then
        echo "$dir" >> finisheddeblurtables.txt
        mv $dir finishedtables
    fi

done

# Count them
count=$(wc -l < missingseqs.txt)
echo "$count folders have .txt files but are missing seqs"

count=$(wc -l < missingdeblurtables.txt)
echo "$count folders have seqs but are missing tables"

count=$(wc -l < finisheddeblurtables.txt)
echo "$count folders have tables"

mv finisheddeblurtables.txt finishedtables