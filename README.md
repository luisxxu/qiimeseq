# qiimeseq

Downloading and sanitizing data for analysis is an incredibly arduous task that pleases only the most Sisyphean. To find joy in this task, one must automate the process as revolt against the human condition. 

The script takes metadata fields and a project accession number, then downloads the fastq files and uses Qiime2's tools (demux, quality-filter, deblur) to produce visualizations and other files necessary for analysis.

Make sure to install EDirect, SRAToolkit, and Qiime2-2023.2 before running the script.

# EDirect Installation
```
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

# SRAToolkit Installation
[SRAToolkit Installation Instructions](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Make sure to keep track of where you unzip SRAToolkit, this will be important later!


# Qiime2-2023.2 Installation
[Qiime2 Installation Instructions](https://docs.qiime2.org/2023.2/install/)

# Running the script
The script only takes one parameter: the file location of sratoolkit/bin. For example, this is how my code is run:
```
bash qiimeseq.sh /Users/luisxu/Documents/HMTOL/sratoolkit.3.2.1-mac-arm64/bin
```
The name/location of your sratoolkit file may be different, so change accordingly. You can also run "pwd" once you're inside the sratoolkit/bin folder and copy/paste the folder location.

After running the script, you'll be prompted to enter in metadata fields. **DO NOT** use any punctuation for the authors/name (e.g. KannEtAlArgentina2021 instead of Kann_et_al_Argentina_2021).

Make sure you placed your bash script in its own separate folder. After the script completes, it'll create a new folder that contains all the output files. Each output is expected to be around 1-2gb, so ensure you have enough storage space on your computer.

The script should run anywhere from 20-40 minutes, but you can also run multiple terminal instances of this script processing multiple studies at the same time.

ok bye
