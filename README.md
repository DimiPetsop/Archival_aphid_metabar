# Archival_aphid_metabar

### This repository includes all scripts required to reproduce the analyses, statistics and figures included in the article: [Identifying archived insect bulk samples using DNA metabarcoding: A case study using the long-term Rothamsted Insect Survey](https://doi.org/10.1002/edn3.542)

You can download the raw sequence files from [Zenodo](https://zenodo.org/records/10995475)

There are multiple steps required which are summarised bellow:
- To demultiplex we used MetaBeat: https://github.com/HullUni-bioinformatics/metaBEAT.
After installation you can use the tool on the command line with the pattern file provided in Raw directory. For example you can run:
```metabeat
/usr/bin/metaBEAT_global.py -Q ../metadata/pattern_metabeat --trim_qual 1 -n 6
```

  
- After demultiplexing run: [Trimming_dada.R](/Scripts/01_Trimming_dada.R)
- Download the MIDORI2 database from here: [MIDORI2_259](https://www.reference-midori.info/download/Databases/GenBank259_2023-12-17/BLAST/uniq/MIDORI2_UNIQ_NUC_GB259_CO1_BLAST.zip)
- Then run blastn on the [command line](/Scripts/blast_assignment.txt) (Please follow the installation steps according to your OS [here](https://www.ncbi.nlm.nih.gov/books/NBK569861/))
- From this point onwards the scripts are numbered and can be run in that particular order although flexible.
