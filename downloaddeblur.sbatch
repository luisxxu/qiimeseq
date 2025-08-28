#!/bin/bash -l
#SBATCH --job-name=qsqdownloaddeblur
#SBATCH --output=/home/lxxu/logs/%x-%j.out
#SBATCH --error=/home/lxxu/logs/%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lxxu@ucsd.edu
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8

sheet=$1

# read header line into array
IFS=$'\t' read -r -a headers < $sheet

# get the raw line with header and study metadata from the slurm array task id
header_line=$(sed -n "1p" "$sheet" | tr -d '\n')
data_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sheet")

# read metadata line into array
IFS=$'\t' read -r -a values <<< "$data_line"

# export metadata as variables
for i in "${!headers[@]}"; do
  # clean variable name
  varname=$(echo "${headers[$i]}" | tr -d '[:punct:]' | tr ' ' '_' | sed 's/^_//;s/_$//')

  # assign metadata value to variable
  declare "$varname=${values[$i]}"
done

accession=$(echo "$accession" | tr -d '[:space:]')

echo "Processing study: $study with accession: $accession"

echo "No previous study directory found. Creating new directory."

mkdir "$study"
cd "$study"
mkdir "$accession"

# download fastq files from study txt using ftp links
count=0
{
    curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=library_strategy,fastq_ftp&format=tsv&download=true&limit=0" >  "${accession}/${accession}.txt"
    tail -n +2 "${accession}/${accession}.txt" | while read -r run ftp strategy; do
        if [ "$strategy" = "AMPLICON" ]; then
            echo "$ftp" | tr ';' '\n' | while read -r url; do
                if [ -n "$url" ]; then
                    wget --tries=5 --waitretry=5 -c "ftp://$url" -P "${accession}" &
                    ((count++))
                    if ((count % 8 == 0)); then
                        wait
                    fi
                fi
            done
        fi
    done
} > /dev/null 2>&1

echo "Initial downloads complete. Retrying downloads on missing files now."

awk -F'\t' 'NR>1 && $3=="AMPLICON" {n=split($2, urls, ";"); for(i=1;i<=n;i++) print urls[i]}' "${accession}/${accession}.txt" > "${accession}/urls.txt"

URL_FILE="${accession}/urls.txt"
DOWNLOAD_DIR="${accession}"
MISSING_NAMES="${accession}/missing_files.txt"
MISSING_URLS="${accession}/missing_urls.txt"
FINAL_NAMES="${accession}/final_missing_files.txt"

retries=0

# Loop until there are no missing files
while true; do
    echo "Checking for missing files..."

    # Extract expected filenames
    basename -a $(cat "$URL_FILE") | sort > "${accession}/expected_files.txt"

    # Get downloaded filenames
    ls "$DOWNLOAD_DIR" | sort > "${accession}/downloaded_files.txt"

    # Find missing files by comparing expected vs downloaded
    comm -23 "${accession}/expected_files.txt" "${accession}/downloaded_files.txt" > "$MISSING_NAMES"

    # If no missing files, exit loop
    if [[ ! -s "$MISSING_NAMES" ]]; then
        echo "All files downloaded."
        break
    fi

    if (( retries >= 20 )); then
        echo "Maximum retries reached. Exiting."
        echo "Missing files: $(cat "$MISSING_NAMES")"
        break
    fi

    # Extract missing file URLs
    grep -Ff "$MISSING_NAMES" "$URL_FILE" > "$MISSING_URLS"

    echo "Retrying download of missing files..."
    count=0
    while read -r url; do
        wget --tries=5 --waitretry=5 -c "ftp://$url" -P "$DOWNLOAD_DIR" &
        ((count++))
        if ((count % 8 == 0)); then
            wait
        fi
    done < "$MISSING_URLS"

    wait  # Ensure all background downloads finish before next loop
    retries=$((retries + 1))
done

echo "Done with download step. Generating metadata and manifest files."

# generate metadata file name
metadata_file="$study""-""metadata.txt"

# generate manifest file name
manifest_file="$study""-""manifest.txt"

# headers for manifest file
echo -e "sample-id\tabsolute-filepath" > "$manifest_file"
# headers for metadata file
echo -e "${header_line}" > "$metadata_file"

name=$study
max_files=500

for file in "$accession"/*.fastq.gz; do
    # Skip files that already have _1 or _2 before .fastq.gz
    if [[ "$file" =~ _[12]\.fastq\.gz$ ]]; then
        continue
    fi

    # Extract the base name without the extension
    base="${file%.fastq.gz}"
    
    # Rename the file to add _1 before the extension
    mv "$file" "${base}_1.fastq.gz"
done

# write sample-ids/filepaths to manifest file
# write sample-ids and other info to metadata file
counter=1
for file in "$accession"/*_1.*; do
    filepath="$PWD/$file"
    filename="$name-$counter"
    sample_id="${filename%.*}"
    echo -e "$sample_id\t$filepath" >> "$manifest_file"
    echo -e "$data_line" >> "$metadata_file"
    ((counter+=1))
    if [ $counter -ge $max_files ]; then
        echo "Reached limit of $max_files files. Stopping."
        break
    fi
done

# number metadata rows and change header
linecount=0
while IFS=$'\t' read -r first rest; do
    if [[ $linecount -eq 0 ]]; then
        echo -e "SampleID\t${rest}" >> ${name}-metadata_counts.txt
        linecount=1
        continue
    fi
    echo -e "${first}-${linecount}\t${rest}" >> ${name}-metadata_counts.txt
    linecount=$((linecount + 1))
done < $metadata_file

mv "${name}-metadata_counts.txt" "$metadata_file"

echo "Manifest file '$manifest_file' created."
echo "Metadata file '$metadata_file' created."

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=sample_accession,library_name,secondary_sample_accession,run_accession,experiment_accession,fastq_bytes,fastq_ftp,fastq_md5,library_source,instrument_platform,submitted_format,library_strategy,library_layout,tax_id,scientific_name,instrument_model,library_selection,center_name,experiment_title,study_title,study_alias,experiment_alias,sample_alias,sample_title,study_accession&format=tsv&download=true&limit=0" > "${name}-ena-metadata.txt"

cp ../metadata.py .
python metadata.py $name
rm metadata.py

# activate qiime
source ~/miniforge3/etc/profile.d/conda.sh
conda activate qiime2-2023.2   

export TMPDIR=/ddn_scratch/lxxu/tmp

# generate demux
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$name"-manifest.txt \
    --output-path "$name"-seqs.qza \
    --input-format SingleEndFastqManifestPhred33V2

echo "Denoising with deblur..."
# generate quality scores for deblur
qiime quality-filter q-score \
    --i-demux "$name"-seqs.qza \
    --o-filtered-sequences "$name"-demux-filtered.qza \
    --o-filter-stats "$name"-demux-filter-stats.qza

# run deblur at 150bp
qiime deblur denoise-16S \
    --i-demultiplexed-seqs "$name"-demux-filtered.qza \
    --p-trim-length 150 \
    --p-left-trim-len 0 \
    --p-jobs-to-start 8 \
    --o-representative-sequences "$name"-rep-seqs-deblur.qza \
    --o-table "$name"-table-deblur.qza \
    --p-sample-stats \
    --o-stats "$name"-deblur-stats.qza \
    --verbose