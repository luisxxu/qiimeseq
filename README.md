# qiimeseq

Make sure to install EDirect, SRAToolkit, and Qiime2-2023.7 before running the script.

# EDirect Installation
```
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```

# SRAToolkit Installation
[SRAToolkit Installation Instructions](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)

Make sure to keep track of where you unzip SRAToolkit, this will be important later!


# Qiime2-2023.7 Installation
[Qiime2 Installation Instructions](https://docs.qiime2.org/2023.7/install/)

# Running the script
The script only takes one parameter: the file location of sratoolkit/bin. For example, this is how my code is run:
```
bash qiime2-pipe.sh /Users/luisxu/Documents/HMTOL/sratoolkit.3.2.1-mac-arm64/bin
```
The name/location of your sratoolkit file may be different, so change accordingly. You can also run "pwd" once you're inside the sratoolkit/bin folder.

After running the script, you'll be prompted to enter in metadata fields. **DO NOT** use any underscores for the authors/name (e.g. Kann-et-al-2021 instead of Kann_et_al_2021).

The script should run anywhere from 20-40 minutes, but you can also run multiple terminal instances of this script at the same time.

ok bye
