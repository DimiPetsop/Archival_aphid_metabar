#### Calculate abundance of non-identified taxa from metabarcoding ####
library(tidyverse)

mota_long <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")

meta_long <- read.csv("Data/Processed/genera_long_23.csv")


genus_not_found <- setdiff(mota_long$genus,meta_long$genus)

genus_abundance <- mota_long[mota_long$genus %in% genus_not_found,]

genus_abundance <- genus_abundance %>% aggregate(DailyCount~genus,FUN=sum)

write.csv(genus_abundance,"Data/Processed/not_found_taxa.csv",row.names = F)
