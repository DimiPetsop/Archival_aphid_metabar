#### Create venn diagram ####

#library(ggVennDiagram)


list.of.packages <- c("ggvenn", "ggVennDiagram","VennDiagram","eulerr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)

meta_long <- read.csv("Data/Processed/genera_long_23.csv")
meta_genus <- unique(meta_long$genus)

mota_long <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")



mota_genus <- unique(mota_long$genus)



x <- list(Metabarcoding= meta_genus,
          Morphology= mota_genus)

svg(file="Figures/Venn_genera.svg", width=6.5,height=5.1)

plot(euler(x),fill=c("orange","lightblue"),lty=0,labels=list(3),quantities = T)
dev.off()

#### Do same for species ####
  meta_long_sp <- read.csv("Data/Processed/species_long_23.csv")


mota_long <- mota_long[!grepl("spp", mota_long$binomial),]
Joined_agg_sp <- mota_long %>% mutate(binomial = case_when(binomial =="Aphis fabae gp." ~ "Aphis fabae",
                                                         binomial !="Aphis fabae" ~ binomial))


y <- list(Metabarcoding= unique(meta_long_sp$species),
          Morphology= unique(Joined_agg_sp$binomial))
svg(file="Figures/Venn_species.svg",width=6.5,height=5.1)

plot(euler(y),fill=c("orange","lightblue"),lty=0,labels=list(3),quantities = T)

dev.off()
