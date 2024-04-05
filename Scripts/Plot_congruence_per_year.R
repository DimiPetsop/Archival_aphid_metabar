#### Plot congruence per year #####

rm(list=ls())


# Load and install required packages
list.of.packages <- c("tidyverse", "viridis","viridisLite","lme4")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


congruence <- read.csv("Data/Processed/congruence_table.csv")

meta_data <- read.csv("Data/Raw/Joined.aphids.cleaned_23.csv")

meta_data <- meta_data %>% select(unique_id,year)

meta_data <- unique(meta_data)

colnames(meta_data)[1] <-colnames(congruence)[1]

congruence_joined <- inner_join(congruence,meta_data,by="Sample")

# Plot congrunece
svg(file="Figures/Congruence_per_year.svg", width=12,height=5.3)
congruence_joined %>% ggplot(aes(x=as.factor(year), y=percent, fill=as.factor(year))) +
  geom_boxplot() +
  scale_fill_viridis(discrete=TRUE,alpha=0.6,option="magma") +
  geom_jitter(color="black", size=0.8, alpha=0.9) +xlab("Year")+
  ylab("Congruence as %") + theme_bw()+ 
  theme(axis.text=element_text(size=14), text=element_text(size=14)) +guides(fill="none")
dev.off()

# Plot jaccards similarity
congruence_joined %>% ggplot(aes(x=as.factor(year), y=jaccard*100, fill=as.factor(year))) +
  geom_boxplot() +
  scale_fill_viridis(discrete=TRUE,alpha=0.6,option="magma") +
  geom_jitter(color="black", size=0.4, alpha=0.9) +xlab("Year")+
  ylab("Congruence as %") + theme_bw()+ 
  theme(axis.text=element_text(size=14), text=element_text(size=14)) +guides(fill="none")



## Congruence GLM ##
congruence_mut <- congruence_joined %>% mutate(percent = percent/100)


m1 <- glm(percent ~year,family="binomial", data=congruence_mut)

m2 <- glmer(percent ~ year +(1|year),data=congruence_mut, family=binomial)



