### Filter datasets at different levels ####
rm(list=ls())


# Load and install required packages
list.of.packages <- c("tidyverse", "viridis","ggdist")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)
source("Scripts/Functions_archival_aphids.R")





# Load long format dataset plus morpho data long format
long_metabar <- read.csv("Data/Processed/species_long_23.csv")

Joined_agg <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")
Joined_agg <- Joined_agg[!grepl("spp", Joined_agg$binomial),]
Joined_agg <- Joined_agg %>% mutate(binomial = case_when(binomial =="Aphis fabae gp." ~ "Aphis fabae",
                                                         binomial !="Aphis fabae" ~ binomial))



# Calculate total reads per sample 
reads_per_sample <- read.table("Data/Processed/seqtab_nochim_22_seed_aphidsb.tsv")

rownames(reads_per_sample) <- sub("_forward.paired.fastq.gz", "",
                                  rownames(reads_per_sample))
rownames(reads_per_sample) <- sub("X", "", rownames(reads_per_sample))
reads_per_sample <- as.data.frame(rowSums(reads_per_sample))
colnames(reads_per_sample) <-"Total_reads"
reads_per_sample$date <-rownames(reads_per_sample)


#Join datasets with total reads and long-meta
reads_join <- inner_join(long_metabar,reads_per_sample,by="date")

#Remove "X" column
reads_join <- reads_join[,c(1,2,3,7)]
long_metabar <- long_metabar[,c(1,2,3)]

Joined_agg <- Joined_agg[,c(1,3,4)]
colnames(Joined_agg)[1] <-"date"
colnames(Joined_agg)[2] <- "species"

#### Calculate all stats for no filter ####
percentage_common <- congruence_find_species(long_metabar,Joined_agg)


#### Calculate all stats for 0.1% filter ####
sample_01_perc_spe <- reads_join[reads_join$measure >=(0.01*reads_join$Total_reads),]

percentage_common_1 <- congruence_find_species(sample_01_perc_spe,Joined_agg)
#### Calculate all stats for 0.05% filter ####

sample_05_perc_spe <- reads_join[reads_join$measure >=(0.005*reads_join$Total_reads),]

percentage_common_05 <- congruence_find_species(sample_05_perc_spe,Joined_agg)

sample_02_perc_spe <- reads_join[reads_join$measure >=(0.002*reads_join$Total_reads),]

#### Calculate all stats for 0.01% filter ####
percentage_common_02 <- congruence_find_species(sample_02_perc_spe,Joined_agg)


#### Create boxplots for figures ####
congruence_nofilter <- percentage_common[,c(1,6)]
congruence_01filter <- percentage_common_1[,c(1,6)]
congruence_05filter <- percentage_common_05[,c(1,6)]
congruence_02filter <- percentage_common_02[,c(1,6)]

jaccard_nofilter <- percentage_common[,c(1,7)]
jaccard_01filter <- percentage_common_1[,c(1,7)]
jaccard_05filter <- percentage_common_05[,c(1,7)]
jaccard_02filter <- percentage_common_02[,c(1,7)]

jaccard_nofilter$group <- "No filtering"
jaccard_01filter$group <-"1% filter"
jaccard_05filter$group <- "0.5% filter"
jaccard_02filter$group <- "0.2% filter"

colnames(jaccard_nofilter)[2] <-"percent"
colnames(jaccard_01filter)[2] <- "percent"
colnames(jaccard_05filter)[2] <- "percent"
colnames(jaccard_02filter)[2] <- "percent"

congruence_nofilter$group <- "No filtering"
congruence_01filter$group <-"1% filter"
congruence_05filter$group <- "0.5% filter"
congruence_02filter$group <- "0.2% filter"

congruence_nofilter$percent <- congruence_nofilter$percent/100
congruence_01filter$percent <- congruence_01filter$percent/100
congruence_05filter$percent <- congruence_05filter$percent/100
congruence_02filter$percent <- congruence_02filter$percent/100

jaccard_nofilter$percent <- round(jaccard_nofilter$percent,digits=2)
jaccard_01filter$percent <-round(jaccard_01filter$percent,digits=2)
jaccard_05filter$percent <- round(jaccard_05filter$percent,digits=2)
jaccard_02filter$percent <- round(jaccard_02filter$percent,digits=2)

jaccard_all <- rbind(jaccard_01filter,jaccard_05filter,jaccard_nofilter,jaccard_02filter)
jaccard_all$method="jaccard"

congruence_all <- rbind(congruence_nofilter,congruence_05filter,congruence_01filter,congruence_02filter)
congruence_all$method="congruence"
congruence <- rbind(jaccard_all,congruence_all)

## Plot supplementary Figure congruence: 

svg(file="Figures/Congruence_species_detect.svg", width=6.7,height=5.2)

na.omit(congruence[congruence$method=="congruence",])%>% arrange(percent) %>% mutate(name=factor(group,levels=c("1% filter","0.5% filter", "0.2% filter","No filtering"))) %>%
  ggplot(aes(x=name,y=percent*100,fill=name,group=interaction(method,name)))+
  stat_halfeye(adjust=.5,width=.6,justification=-.3,alpha=0.8)+
  geom_boxplot(width=0.25,alpha=0.7,position=position_dodge(1))+
  theme_bw()+geom_point(size=1.1,alpha=.2,position=position_jitter(width=.1))+
  scale_fill_viridis(discrete=TRUE,option="H")+
  xlab("Filter thresholds")+
  ylab("Congruence as detectability")+coord_flip()+
  guides(fill="none")+
  theme(axis.text = element_text(size=14),text = element_text(size=14))

dev.off()
svg(file="Figures/Congruence_species_jaccard.svg", width=6.7,height=5.2)

na.omit(congruence[congruence$method=="jaccard",])%>% arrange(percent) %>% mutate(name=factor(group,levels=c("1% filter","0.5% filter", "0.2% filter","No filtering"))) %>%
  ggplot(aes(x=name,y=percent*100,fill=name,group=interaction(method,name)))+
  stat_halfeye(adjust=.5,width=.6,justification=-.3,alpha=0.8)+
  geom_boxplot(width=0.25,alpha=0.7,position=position_dodge(1))+
  theme_bw()+geom_point(size=1.1,alpha=.2,position=position_jitter(width=.1))+
  scale_fill_viridis(discrete=TRUE,option="H")+
  xlab("Filter thresholds")+
  ylab("Congruence as similarity")+coord_flip()+
  guides(fill="none")+
  theme(axis.text = element_text(size=14),text = element_text(size=14))

dev.off()
#### Calculate false positives and negatives at species level


false_negative <- percentage_common[,c(1,4)]
false_negative01 <- percentage_common_1[,c(1,4)]
false_negative02 <- percentage_common_05[,c(1,4)]
false_negative03 <- percentage_common_02[,c(1,4)]

colnames(false_negative)[2] <- "n"
colnames(false_negative01)[2] <- "n"
colnames(false_negative02)[2] <- "n"
colnames(false_negative03)[2] <- "n"

false_negative$group <- "False negative"
false_negative01$group <-"False negative"
false_negative02$group <-"False negative"
false_negative03$group <-"False negative"

false_negative$filter <- "No filtering"
false_negative01$filter <-"1% filter"
false_negative02$filter <-"0.5% filter"
false_negative03$filter <-"0.2% filter"



#congruence_02filter <- percentage_common_02[,c(1,6)]

false_positive <- percentage_common[,c(1,3)]
false_positive01 <- percentage_common_1[,c(1,3)]
false_positive02 <- percentage_common_05[,c(1,3)]
false_positive03 <- percentage_common_02[,c(1,3)]


colnames(false_positive)[2] <- "n"
colnames(false_positive01)[2] <- "n"
colnames(false_positive02)[2] <- "n"
colnames(false_positive03)[2] <- "n"

false_positive$group <- "False positive"
false_positive01$group <-"False positive"
false_positive02$group <-"False positive"
false_positive03$group <-"False positive"


false_positive$filter <- "No filtering"
false_positive01$filter <-"1% filter"
false_positive02$filter <-"0.5% filter"
false_positive03$filter <-"0.2% filter"

false_positives <- rbind(false_positive,false_positive01,false_positive02,false_positive03)
false_negatives <- rbind(false_negative,false_negative01,false_negative02,false_negative03)

all_controls <- rbind(false_positives,false_negatives)



### Supplementary plot errors 
svg(file="Figures/Errors_species.svg", width=6.5,height=5.3)
na.omit(all_controls)%>% arrange(n) %>% mutate(name=factor(filter,levels=c("No filtering","0.2% filter", "0.5% filter","1% filter"))) %>%
  ggplot(aes(x=name,y=n,fill=group,group=interaction(group,name)))+
  geom_boxplot(width=0.25,alpha=0.7)+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE,option="H")+
  xlab("Filter thresholds")+
  ylab("Number of taxa identified as false positives/negatives")+
  guides(fill=guide_legend(title="Errors"))+
  theme(axis.text = element_text(size=14),text = element_text(size=14))
dev.off()

#### Summary stats for document ####

# Mean of false positives
mean(false_positive$n)

# Mean of false negatives
mean(false_negative$n)


#### Compare false positives and negatives before and after the exclusion of taxa ####
all_controls_removed <- read.csv("Data/Processed/all_controls_removed_species.csv")
all_controls_removed$method = "remove"
all_controls$method = "remain"

all_controls_combined <- rbind(all_controls,all_controls_removed)



drop_na(all_controls_combined[all_controls_combined$group =="False positive",]) %>% group_by(filter,method)%>%
  summarise(mean(n))

