#### Create long aphid dataframe for species and genera
rm(list=ls())
library(tidyverse)


# Load wide format
joined_tablesb <- read.csv("data/processed/taxa_asv_259.csv")

# Fix name strings
colnames(joined_tablesb) <- sub("_forward.paired.fastq.gz", "", colnames(joined_tablesb))
colnames(joined_tablesb) <- sub("X", "", colnames(joined_tablesb))

#### Get species long dataframe

# Get ASVs with genus level identification
genera <- joined_tablesb %>% filter(taxon =="species" |
                                      taxon == "genus")

# Get all non-aphid families
genera <- genera %>% filter(family != "Adelgidae", 
                            family != "Anoeciidae",
                            family != "Aphididae",
                            family != "Pemphigidae",
                            family != "Thelaxidae")

#Remove non aphids
genera <- genera %>% filter(genus!="class")


genera_long <- gather(genera,date,measure,"00110R":"0EE127R",factor_key = TRUE)
genera_long <- genera_long[genera_long$measure != 0,]
genera_long <- genera_long %>% aggregate(measure~genus+date,FUN=sum)


taxa_only <- unique(joined_tablesb[c(2,4:9)],)

genera_long <- inner_join(genera_long,taxa_only,by="genus")


samples_per_order <- genera_long %>% aggregate(measure~order+date,FUN=sum)
samples_per_order %>% group_by(order) %>% summarise(n())


# Merge long dataframes with metadata
meta_data <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")
meta_data <- meta_data %>% select(unique_id,CalDate,year,month)
meta_data <- unique(meta_data)
colnames(meta_data)[1] <-"date"

# Merge species df
species_long_meta <- inner_join(species_long,meta_data,by="date")

# Merge genera df
genera_long_meta <- inner_join(genera_long,meta_data,by="date")



# Write both to csv

write.csv(genera_long_meta,"Data/Processed/non_target_long.csv", 
          row.names = FALSE)

arthropoda <- genera_long[genera_long$phylum == "Arthropoda",]
sum(arthropoda$measure)


non_arthropod <- genera_long[genera_long$phylum != "Arthropoda",]
sum(non_arthropod$measure)
