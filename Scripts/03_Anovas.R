#### Anova comparing non destructive vs destructive on 1. Total reads, 2. Congruence
rm(list=ls())
library(tidyverse)
library(ggrepel)
#Load meta data info
group_non_destr <- read.csv("Data/Raw/Subset_1_F.csv")
group_destr <- read.csv("Data/Raw/Subset_2_F.csv")

#Load asv_table
asv_table <- read.table("Data/Processed/seqtab_nochim_22_seed_aphidsb.tsv")


reads_per_sample <- as.data.frame(colSums(t(asv_table)))
colnames(reads_per_sample) <- "total_reads" 
rownames(reads_per_sample) <- sub("_forward.paired.fastq.gz", "", rownames(reads_per_sample))
rownames(reads_per_sample) <- sub("X", "", rownames(reads_per_sample))

reads_per_sample$sample <- rownames(reads_per_sample)

group_non_destr$method <- "Non_destructive"
group_destr$method <- "Destructive"

group_non_destr <- group_non_destr %>% select(CalDate,Unique.id,Grouping)

group_destr <- group_destr %>% select(CalDate,Unique.id,Grouping)


#Anova for differences in total reads #

grouped_info_df <- rbind(group_destr,group_non_destr)
colnames(grouped_info_df)[2] <- "sample"

joined_df <- inner_join(reads_per_sample,grouped_info_df,by="sample")


anova_total_reads <- aov(total_reads ~ as.factor(Grouping),data=joined_df)
summary(anova_total_reads)



#Anova for differences in congruence (as similarity & congruence measured by jaccar's index)
congruence <- read.csv("Data/Processed/congruence_table.csv")

colnames(congruence)[1] <- "sample"


congruence_joined <- inner_join(congruence,joined_df,by="sample")

# Anova for congruence based on "detectability@
anova_congruence_similarity <- aov(percent~Grouping,data=congruence_joined)
summary(anova_congruence_similarity)

# Anova for congruence based on jaccard index
anova_congruence_jaccard <- aov(jaccard*100~Grouping,data=congruence_joined)
summary(anova_congruence_jaccard)


congruence_joined$CalDate <- as.Date(congruence_joined$CalDate , format="%d/%m/%Y")

congruence_joined <- congruence_joined %>%
  dplyr::mutate(year = lubridate::year(CalDate))

# GLM for congruence + year + reads + year:reads

m1 <- glm(percent/100 ~year+log(total_reads)+year:log(total_reads),family="binomial", data=congruence_joined)

m2 <- glm(jaccard ~year+log(total_reads)+ year:log(total_reads),family="binomial", data=congruence_joined)


# Get long format meta data and long format mota data and get linear model plus plot
meta_long <- read.csv("Data/Processed/genera_long_23.csv")
mota_long <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")

meta_genus_agg <- meta_long %>% aggregate(measure~genus,FUN=sum)
mota_genus_agg <- mota_long %>% aggregate(DailyCount~genus,FUN=sum)

meta_genus_agg$relative_per <- round(meta_genus_agg$measure/sum(meta_genus_agg$measure),digits=4)
mota_genus_agg$relative_per <- round(mota_genus_agg$DailyCount/sum(mota_genus_agg$DailyCount),digits=4)

joined_meta_mota <- inner_join(mota_genus_agg,meta_genus_agg,by="genus")

# Regression for reads/abundance
regres_reads <- lm(log(measure)~log(DailyCount),data=joined_meta_mota)
summary(regres_reads)


# Supplementary plot 
svg(file="Figures/GLM_supplementary.svg", width=11,height=5.5)

ggplot(joined_meta_mota,aes(y=log(measure),x=log(DailyCount))) + 
  geom_point()+geom_text_repel(aes(label=genus),size = 3.5)+
  geom_smooth(method="lm",se=T,color="blue")+ theme_bw()+
  labs(x="Logged transformed reads from MOTA",y="Logged transformed counts from META") +
  theme(axis.text = element_text(size=14),text = element_text(size=14))
dev.off()
# Table Supplementary A3
meta_genus_agg_02 <- meta_genus_agg[meta_genus_agg$relative_per >=0.02,]

# Table Supplementary A4
mota_genus_agg_02 <- mota_genus_agg[mota_genus_agg$relative_per >=0.02,]


# Table A7 column 2

mean_congruence_per_year <- congruence_joined %>% aggregate(percent~year,FUN=mean)
sd_congruence_per_year <-congruence_joined %>% aggregate(percent~year,FUN=sd)

table_23 <- inner_join(mean_congruence_per_year,sd_congruence_per_year,by="year")


## Table A7 column 3
aggregated_counts <- mota_long %>% group_by(genus,year) %>%
  summarise(total=mean(DailyCount))

aggregated_counts %>% group_by(year) %>%
  summarise(n())



## Table A7 column 4
aggregated_reads <- meta_long %>% group_by(genus,year) %>%
  summarise(total=mean(measure))

aggregated_reads %>% group_by(year) %>%
  summarise(n())

