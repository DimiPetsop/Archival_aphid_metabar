### NMDS analysis 

rm(list=ls())

library(tidyverse)
library(viridis)

library(vegan)
library(ggrepel)

# set seed for NMDS
set.seed(31301992)

# Get dataset asv dataset
asv_tab <- read.csv("Data/Processed/taxa_asv_259.csv")

# Parse column namee
colnames(asv_tab) <- sub("_forward.paired.fastq.gz", "", colnames(asv_tab))
colnames(asv_tab) <- sub("X", "", colnames(asv_tab))


#Get ASVs with genus level identification
genera <- asv_tab %>% filter(taxon =="species" | taxon == "genus")
genera <- genera %>% filter(genus !="class")

# Get only Hemiptera
genera <- genera %>% filter(order == "Hemiptera")
genera <- genera %>% filter(family != "Miridae" | family != "Lygaeidae" )
#Remove non aphids


#Long metabarcoding dataset
metabar_long <- read.csv("Data/Processed/genera_long_23.csv")


##Create matrix for morphological and metabarcoding data

morpho_data <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")
meta_data <- morpho_data[,c(1,2,6,7)]
meta_data <-unique(meta_data)
colnames(meta_data)[2]<-"CalDate"
colnames(meta_data)[1]<-"date"

metabar_by_year <- aggregate(measure~genus+year,data=metabar_long,FUN=sum)

morpho_by_year <- aggregate(DailyCount~genus+year,data=morpho_data,FUN=sum)


matrix_1 <- spread(metabar_by_year,genus,measure)
matrix_1 <- data.frame(matrix_1[,-1],row.names=matrix_1[,1])
matrix_1[is.na(matrix_1)] <- 0
matrix_1[matrix_1>=1] <- 1

matrix_2 <- spread(morpho_by_year,genus,DailyCount)
matrix_2 <- data.frame(matrix_2[,-1],row.names=matrix_2[,1])
matrix_2[is.na(matrix_2)] <- 0
matrix_2[matrix_2>=1] <- 1

insect.dist_meta <- decostand(matrix_1,method="total")
dist_mat_meta <- vegdist(insect.dist_meta,method="jaccard")
insec_nmds_meta <- metaMDS(dist_mat_meta,distance="jaccard")

insect.dist_morpho <- decostand(matrix_2,method="total")
dist_mat_morpho <- vegdist(insect.dist_morpho,method="jaccard")
insec_nmds_morpho <- metaMDS(dist_mat_morpho,distance="jaccard")


nmds_meta <- as.data.frame(scores(insec_nmds_meta))
nmds_morpho <-as.data.frame(scores(insec_nmds_morpho))

nmds_meta$Dataset <- "META"
nmds_morpho$Dataset <- "MOTA"
nmds_meta$year <-rownames(nmds_meta)
nmds_morpho$year <- rownames(nmds_morpho)

nmds_joined <- rbind(nmds_meta,nmds_morpho)


pcoa_meta <- cmdscale(dist_mat_meta,eig=TRUE)
pca_meta<- as.data.frame(scores(pcoa_meta))

pcoa_morpho <-cmdscale(dist_mat_morpho,eig=TRUE)
pca_morpho <- as.data.frame(scores(pcoa_morpho))

pca_meta$group <- "meta"
pca_morpho$group <- "morpho"
pca_meta$year <-rownames(pca_meta)
pca_morpho$year <- rownames(pca_morpho)

pca_joined <- rbind(pca_meta,pca_morpho)

ggplot(data=pca_joined,aes(x=Dim1,y=Dim2,colour=group))+geom_point()+geom_text(aes(label=year))

# Plot NMDS for both
svg(file="Figures/Figure_4.svg", width=8.3,height=6.1)

ggplot(data=nmds_joined,aes(x=NMDS1,y=NMDS2,colour=Dataset)) +
  geom_point(size=4,aes(shape=Dataset))+geom_text_repel(aes(label=year),size=4) +
  geom_line(aes(group=year),colour="black",linetype="dotted",linewidth=0.55)+theme_bw() +
  theme(axis.text = element_text(size=14),text = element_text(size=14))
dev.off()



#p + labs(fill="Dataset")
#ggplot(data=pca_joined,aes(x=Dim1,y=Dim2,colour=group)) +
 # geom_point(size=4,aes(shape=group))+geom_text_repel(aes(label=year),size=3) +
#  geom_line(aes(group=year),colour="black",linetype="dotted",linewidth=0.40)+theme_bw() +
#  theme(axis.text = element_text(size=14),text = element_text(size=14))



