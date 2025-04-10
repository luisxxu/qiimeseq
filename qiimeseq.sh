#!/bin/bash

toolkitpath=$1

export PATH=${HOME}/edirect:${PATH}
export PATH=$toolkitpath:$PATH

# prompt user for info
read -p "Study name (e.g. ColstonJackson-et-al-2021) (NO UNDERSCORES): " name
read -p "Country: " location
read -p "Continent: " continent
read -p "Accession Number: " accession
read -p "Notes: " notes
read -p "Forward primer: " forward
read -p "Reverse primer: " reverse
read -p "Primer (V4 or V3-V4): " primer
read -p "Age (Adults/Children/Both): " age
read -p "Sex (Male/Female/Both): " sex
read -p "Disease: " disease
read -p "Disease Group: " disease_group
read -p "Initials: " initials

# generate temp directory for fastq files
mkdir "$name"-fastq-files
mkdir "$accession"
cd "$accession"

# grab 30 sample accessions from project accession
esearch -db sra -query $accession | efetch -format runinfo | cut -d "," -f 1 > $accession".csv"
tail -n +2 $accession".csv" > temp.csv && mv temp.csv $accession".csv"
head -n 30 $accession".csv" > temp.csv && mv temp.csv $accession".csv"

# download sample accessions
while IFS= read -r sample
do
  prefetch -q --max-size 1g $sample
  fasterq-dump -q $sample
  if ! ls *_1.fastq 1> /dev/null 2>&1; then
    echo "WARNING: No paired end reads found."
  else
    rm "$sample""_2.fastq"
  fi
  mv *.fastq ../"$name"-fastq-files
done < $accession".csv"

# cleanup
cd ..
rm -rf "$accession"

cd "$name"-fastq-files

# generate metadata file name
metadata_file="$name""-""metadata.txt"

# generate manifest file name
manifest_file="$name""-""manifest.txt"
counter=1

# headers for manifest file
echo -e "sample-id\tabsolute-filepath" > "$manifest_file"

# write sample-ids/filepaths to manifest file
# write sample-ids and other info to metadata file
for file in *fastq; do
    filepath="$PWD""/""$file"
    filename="$name""-""$counter"
    sample_id="${filename%.*}"
    echo -e "$sample_id\t$filepath" >> "$manifest_file"
    echo -e "$sample_id\t$notes\t$primer\t$forward\t$reverse\t$location\t$initials\t$age\t$sex\t$disease\t$accession" >> "$metadata_file"

    ((counter+=1))
done

mv *.txt ..
cd ..

echo "Manifest file '$manifest_file' created."
echo "Metadata file '$metadata_file' created."

# activate qiime
source ~/.zshrc
conda activate qiime2-2023.2

# generate demux
qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$name"-manifest.txt \
    --output-path "$name"-demux.qza \
    --input-format SingleEndFastqManifestPhred33V2

# generate demux visualization
qiime demux summarize \
    --i-data "$name"-demux.qza \
    --o-visualization "$name"-demux.qzv

# generate quality scores for deblur
qiime quality-filter q-score \
    --i-demux "$name"-demux.qza \
    --o-filtered-sequences "$name"-demux-filtered.qza \
    --o-filter-stats "$name"-demux-filter-stats.qza

# run deblur at 150bp
qiime deblur denoise-16S \
  --i-demultiplexed-seqs "$name"-demux-filtered.qza \
  --p-trim-length 150 \
  --p-left-trim-len 0 \
  --p-jobs-to-start 4 \
  --o-representative-sequences "$name"-rep-seqs-deblur.qza \
  --o-table "$name"-table-deblur.qza \
  --p-sample-stats \
  --o-stats "$name"-deblur-stats.qza

# move all files to single directory
mkdir $name
mv "$name"-fastq-files $name
mv deblur.log $name
mv $name-* $name