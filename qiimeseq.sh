wget -O gmtol.tsv "https://docs.google.com/spreadsheets/d/1uxOu99l7wGm_hghIer-HcNrSp1DF7qE2b6INBigrwx8/export?format=tsv" > /dev/null 2>&1

num_jobs=$(wc -l < gmtol.tsv)
num_jobs=$((num_jobs + 1)) # add one for header line
echo $num_jobs

#echo "Starting download of deblur data for $((num_jobs - 1)) studies."
JOB_1=$(sbatch --array=2-$num_jobs%2 downloaddeblur.sbatch gmtol.tsv | cut -d' ' -f4)
echo "Job 1: $JOB_1"

#echo "Download of deblur data completed. Checking missing files."
JOB_2=$(sbatch --dependency=afterok:${JOB_1} missing.sbatch | cut -d' ' -f4)
echo "Job 2: $JOB_2"

#echo "Merging data."
JOB_3=$(sbatch --dependency=afterok:${JOB_2} merge.sbatch | cut -d' ' -f4)
echo "Job 3: $JOB_3"

#echo "Running greengenes2 mapping and generating core metrics."
JOB_4=$(sbatch --dependency=afterok:${JOB_3} greengenes.sbatch | cut -d' ' -f4)
echo "Job 4: $JOB_4"

echo "All jobs submitted. Good luck."