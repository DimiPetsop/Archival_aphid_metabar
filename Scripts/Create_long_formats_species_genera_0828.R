#### Create long aphid dataframe for species and genera
rm(list=ls())

# Install and load packages required
install.packages("tidyverse")
library(tidyverse)
# Load wide format
joined_tablesb <- read.csv("Data/Processed/taxa_asv_259.csv")

# Fix name strings
colnames(joined_tablesb) <- sub("_forward.paired.fastq.gz", "", colnames(joined_tablesb))
colnames(joined_tablesb) <- sub("X", "", colnames(joined_tablesb))

#### Get species long dataframe

# Get ASVs with genus, species identification
species <- joined_tablesb %>% filter(taxon =="species")

# Get only Hemiptera


# Filter Hemiptera that are not within the major Aphid families
species <- species %>% filter(order == "Hemiptera")
species <- species %>% filter(family != "Miridae", 
                              family != "Lygaeidae")


#Remove non aphids
species <- species %>% filter(genus!="class")


species_long <- gather(species,date,measure,colnames(species)[12]:colnames(species)[198],factor_key = TRUE)
species_long <- species_long[species_long$measure != 0,]
species_long <- species_long %>% aggregate(measure~species+date,FUN=sum)


#### Get genera long dataframe

# Get ASVs with genus level identification
genera <- joined_tablesb %>% filter(taxon =="species" |
                                       taxon == "genus")

# Get only Hemiptera
genera <- genera %>% filter(order == "Hemiptera")
genera <- genera %>% filter(family != "Miridae", 
                            family != "Lygaeidae")

#Remove non aphids
genera <- genera %>% filter(genus!="class")


genera_long <- gather(genera,date,measure,"00110R":"0EE127R",factor_key = TRUE)
genera_long <- genera_long[genera_long$measure != 0,]
genera_long <- genera_long %>% aggregate(measure~genus+date,FUN=sum)


# Merge long dataframes with metadata
meta_data <- read.csv("data/raw/Joined.aphids.cleaned_23.csv")
meta_data <- meta_data %>% select(unique_id,CalDate,year,month)
meta_data <- unique(meta_data)
colnames(meta_data)[1] <-"date"

# Merge species df
species_long_meta <- inner_join(species_long,meta_data,by="date")

# Merge genera df
genera_long_meta <- inner_join(genera_long,meta_data,by="date")






# Write both to csvs

write.csv(species_long_meta, 
          "Data/Processed/species_long_23.csv",
          row.names = FALSE)

write.csv(genera_long_meta,"Data/Processed/genera_long_23.csv", 
          row.names = FALSE)
