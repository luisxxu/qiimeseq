import pandas as pd
import sys

name = sys.argv[1]

manifest_file = f"{name}-manifest.txt"
custom_metadata_file = f"{name}-metadata.txt"
ena_metadata_file = f"{name}-ena-metadata.txt"

output_file = f"{name}-combined-metadata.txt"

manifest = pd.read_csv(manifest_file, sep="\t", header=None, names=["SampleID","filepath"])

custom_meta = pd.read_csv(custom_metadata_file, sep="\t", keep_default_na=True)

ena_meta = pd.read_csv(ena_metadata_file, sep="\t")

run_to_sample = {}
for _, row in manifest.iterrows():
    # extract run_accession from the filepath (filename, last part)
    run_accession = row['filepath'].split("/")[-1].split("_")[0]  # e.g., SRR8143822
    run_to_sample[run_accession] = row['SampleID']

ena_meta_filtered = ena_meta[ena_meta['run_accession'].isin(run_to_sample.keys())].copy()

# Map run_accession to sample_id
ena_meta_filtered['SampleID'] = ena_meta_filtered['run_accession'].map(run_to_sample)

# Merge with custom metadata on sample_id
merged_meta = pd.merge(custom_meta, ena_meta_filtered.drop(columns=['filepath'], errors='ignore'), 
                       on='SampleID', how='left')

# Save merged metadata
merged_meta.to_csv(output_file, sep="\t", index=False, na_rep="NaN")