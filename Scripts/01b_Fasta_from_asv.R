library(ShortRead)
library(Biostrings)
library(taxreturn)
rm(list=ls())

asv_table_22 <- read.table("data/processed/seqtab_nochim_22_seed_aphidsb.tsv")


fasta <- DNAStringSet(colnames(asv_table_22))

names_fasta<- paste0("ASV",seq_along(colnames(asv_table_22)))

names(fasta) <- names_fasta               

write_fasta(fasta,file="data/processed/aphids_fasta.fasta")
