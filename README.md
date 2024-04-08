# Archival_aphid_metabar

### This repository includes all scripts required to reproduce the analyses, statistics and figures included in the article: ###

You can download the raw sequence files here: 

There are multiple steps required which can be followed bellow:
- To demultiplex we used MetaBeat: https://github.com/HullUni-bioinformatics/metaBEAT
- After demultiplexing run: [Trimming_dada.R](/Scripts/Trimming_dada.R)
- Download the MIDORI2 database from here: [MIDORI2_259)(https://www.reference-midori.info/download/Databases/GenBank259_2023-12-17/BLAST/uniq/MIDORI2_UNIQ_NUC_GB259_CO1_BLAST.zip)
- Then run blastn on the [command line](/Scripts/blast_assignment.txt)
