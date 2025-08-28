# qiimeseq

Downloading and sanitizing data for analysis is an incredibly arduous task that pleases only the most Sisyphean. To find joy in this task, one must automate the process as revolt against the human condition. 

The script takes a metadata sheet with, requiring a column named "accession" and "study", then downloads the fastq files and uses Qiime2's tools (demux, quality-filter, deblur) to produce visualizations and other files necessary for analysis.

Special thanks to Sam, Noah, and Ananya for their contributions to this project.

Make sure to install Qiime2-2023.2 before running the script.

# Qiime2-2023.2 Installation
[Qiime2 Installation Instructions](https://docs.qiime2.org/2023.2/install/)

# Running the script
The script only takes one parameter: the file location of the metadata sheet. For example, this is how the script can be run:
```
bash qiimeseq.sh ~/gmtol/gmtol.tsv
```
The name/location of your metadata file may be different, so change accordingly.

Ensure that the metadata files don't have any punctuation, but especially **DO NOT** use any punctuation for the study/title (e.g. KannEtAlArgentina2021 is okay instead of Kann_et_al_Argentina_2021).

ok bye
