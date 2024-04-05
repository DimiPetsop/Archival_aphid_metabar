#### Script to check how many sequences were thrown out at each step ####
## Check and install required packages 
rm(list=ls())
list.of.packages <- c("dplyr", "stringr","taxreturn","Biostrings")

lapply(list.of.packages, require, character.only = TRUE)

result <- read.csv("Data/Processed/output_aphids_259",sep="\t",header=FALSE)

seqtab_final <- read.table("Data/Processed/seqtab_nochim_22_seed_aphidsb.tsv")

# Copy names 
blast_names <-c("qseqid","sseqid","stitle","pident","length","mismatch","gapopen", 
                "qstart","qend","qlen","sstart","send","evalue","bitscore","qcovs")


colnames(result) <-blast_names

assigned <- read.csv("Data/Processed/taxa_asv_259.csv")


seqtab_finalb <- seqtab_final
seqtab_finalb <- as.data.frame(seqtab_finalb)
names(seqtab_finalb) <-paste0("ASV",seq_along(seqtab_finalb))


seqtab_finalb <- as.data.frame(t(seqtab_finalb))
seqtab_finalb$ASV <- rownames(seqtab_finalb)

# How many have lower pident and qcovs
thrown_filt <-result %>% dplyr::filter(!is.na(sseqid)) %>%
  dplyr::filter(pident < 97 | qcovs < 95)

thrown_filt <- thrown_filt %>% filter(!(qseqid %in% assigned$OTU))





# How many have smaller length than 305 or bigger than 315
throw_length <- result %>% dplyr::filter(qlen < 305 | qlen >315)

throw_length <- throw_length %>% filter(!(qseqid %in% assigned$OTU))

#sum(colSums(seqtab_finalb[seqtab_finalb$ASV %in% unique(throw_length$qseqid),-188]))
# How many did not get assigned in total from blast? 

filtered_threshold_names <- unique(thrown_filt$qseqid)
filtered_length <- unique(throw_length$qseqid)

only_length <- setdiff(filtered_length,filtered_threshold_names)
only_threshold <- setdiff(filtered_threshold_names,filtered_length)

both_length_threhold <- intersect(filtered_length,filtered_threshold_names)

sum(colSums(seqtab_finalb[seqtab_finalb$ASV %in% only_threshold,-188]))
sum(colSums(seqtab_finalb[seqtab_finalb$ASV %in% only_length,-188]))
sum(colSums(seqtab_finalb[seqtab_finalb$ASV %in% both_length_threhold,-188]))


# How many reads in the samples found
sum(colSums(seqtab_finalb[seqtab_finalb$ASV %in% unique(result$qseqid),-188]))

# # Check how many reads were thrown out: 
# Not assigned at all
sum(colSums(seqtab_finalb[!(seqtab_finalb$ASV %in% result$qseqid),-188]))
not_assigned <- seqtab_finalb[!(seqtab_finalb$ASV %in% result$qseqid),-188]







asv_throw <- c(unique(throw_length$qseqid),unique(thrown_filt$qseqid),unique(rownames(not_assigned)))
dif_asv_throw <-setdiff(seqtab_finalb$ASV,asv_throw)
dif_ass_throw <- setdiff(dif_asv_throw,assigned$OTU)

sum(colSums(seqtab_finalb[seqtab_finalb$ASV %in% dif_ass_throw,-188]))
