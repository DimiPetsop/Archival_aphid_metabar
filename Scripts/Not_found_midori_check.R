##### Find how many aphids genera and species are actually within Midori database #####
rm(list=ls())

library(taxreturn)
library(Biostrings)
library(tidyverse)
## Get ncbi database
db <- get_ncbi_taxonomy()

## Get MIDORI database fasta and convert to CSV

midori_db <- readDNAStringSet("Data/Raw/MIDORI2_UNIQ_NUC_GB259_CO1_DADA2.fasta")

seq_name = names(midori_db)
sequence = paste(midori_db)
midori_df <- data.frame(seq_name, sequence)

tax_id <- sub(".*_(\\d+);.*","\\1",midori_df$seq_name)
tax_id <- as.data.frame(tax_id)
tax_id$tax_id <- as.integer(tax_id$tax_id)

# Join ncbi with tax_id 
midori_taxa_info <- inner_join(tax_id,db,by="tax_id")

## Get aphid data 
mota_data <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")

#Number of genera not found in database
setdiff(mota_data$genus,midori_taxa_info$genus)

mota_data <- mota_data[!grepl("spp", mota_data$binomial),]
mota_data <- mota_data %>% mutate(binomial = case_when(binomial =="Aphis fabae gp." ~ "Aphis fabae",
                                                         binomial !="Aphis fabae" ~ binomial))
# Species not represented in database
setdiff(mota_data$binomial,midori_taxa_info$tax_name)

not_found_species_remove <- mota_data %>% filter(!(binomial %in% setdiff(mota_data$binomial,midori_taxa_info$tax_name)))

write.csv(not_found_species_remove,"Data/Processed/removed_species_long.csv",row.names = F)



## Get blast output data

result <- read.csv("Data/Processed/output_aphids_259",sep="\t",header=FALSE) # Change this to required database output

blast_names <-c("qseqid","sseqid","stitle","pident","length","mismatch","gapopen", 
                "qstart","qend","qlen","sstart","send","evalue","bitscore","qcovs")


colnames(result) <-blast_names

result$stitle <-as.character(result$stitle)

# Parse blast assignment
result$stitle <-str_extract(result$stitle,"\\w+_\\w+_\\d+")

result <- result %>% mutate(genus = stitle %>% str_remove("_.*$"),
       spp = stitle %>% str_match("_\\w+_") %>% str_match("[^_].+[^_]")) 
colnames(result)[17] <- "spp"


not_found_taxa <- read.csv("Data/Processed/not_found_taxa.csv")

not_found_blast <- result[result$genus %in% not_found_taxa$genus,]

# What's the mean pident for every genus
not_found_blast %>% group_by(genus) %>% summarise(mean=mean(pident))

# What's the top pid for every genus
not_found_blast %>% group_by(genus) %>% summarise(max=max(pident))

