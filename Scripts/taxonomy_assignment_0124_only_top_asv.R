#### Taxonomic Assignment script, Parses files from blast and db and matches ASV table to taxonomy #####


## Check and install required packages 
list.of.packages <- c("dplyr", "stringr","taxreturn","Biostrings")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
devtools::install_github("alexpiper/taxreturn")

lapply(list.of.packages, require, character.only = TRUE)

#Load blast output
result <- read.csv("Data/Processed/output_aphids_259",sep="\t",header=FALSE) # Change this to required database output

#Load DADA2 asv table
seqtab_final <- read.table("Data/Processed/seqtab_nochim_22_seed_aphidsb.tsv")

# Copy names 
blast_names <-c("qseqid","sseqid","stitle","pident","length","mismatch","gapopen", 
                "qstart","qend","qlen","sstart","send","evalue","bitscore","qcovs")


colnames(result) <-blast_names


#### Process only asv's with stringent criteria ####
# Filter at 97% percent identity and 95 query cover%

results <-result %>% dplyr::filter(!is.na(sseqid)) %>%
  dplyr::filter(pident >= 97, qcovs >= 95)

# Filter reads between 305 and 315 inclusive (10bp window)
#rm(result)
results <- results %>% dplyr::filter(qlen >= 305, qlen <=315)


results$stitle <-as.character(results$stitle)

# Parse blast assignment
results$stitle <-str_extract(results$stitle,"\\w+_\\w+_\\d+")

#Get only top hits
top_hit <- results %>%
  dplyr::group_by(qseqid) %>%
  dplyr::top_n(1, bitscore) %>%
  mutate(genus = stitle %>% str_remove("_.*$"),
         spp = stitle %>% str_match("_\\w+_") %>% str_match("[^_].+[^_]"))  %>%   
  dplyr::summarise(spp = paste(sort(unique(spp)), collapse = "/"), genus, pident, qcovs, evalue, bitscore, .groups = "keep") %>% 
  dplyr::mutate(binomial = paste(genus, spp)) %>% 
  dplyr::distinct() %>% 
  dplyr::add_tally() %>% 
  dplyr::mutate(binomial = dplyr::case_when(
    n > 1 ~ as.character(NA), 
    n == 1 ~ binomial
  )) %>% 
  dplyr::select(OTU = qseqid, genus, species = binomial)


### Remove top hits where species assignment = NA

#Download taxonomy
# options(timeout=15000), enable only if there is a mistake when downloading the database
db <- get_ncbi_taxonomy()


# Collapse top_hit & heirarchy from database
heirarchy <- top_hit %>% 
  dplyr::select(genus) %>%
  distinct() %>% 
  left_join(db %>% dplyr::select(genus, family, order, class, phylum, kingdom, superkingdom), by="genus") %>%
  distinct()  %>% mutate(family = case_when(is.na(family) ~ "NA", !is.na(family)  ~ family)) %>%
  dplyr::select(OTU, genus, family) %>%
  left_join(db %>% dplyr::select(family, order, class, phylum, kingdom, superkingdom), by="family") %>%
  distinct() 


# Note how heirarchy has 10923 rows, (mistake from database)   
## Count  asv's numbers (rows)
heirarchy <- heirarchy %>%  group_by(OTU) %>% mutate(counts=n())

# Find asv with max counts
heirarchy[max(heirarchy$counts),]$OTU


heirarchy <-  heirarchy %>% dplyr::filter(OTU !=heirarchy[max(heirarchy$counts),]$OTU) 

#Remove plant misidentifications Viridiplantae, mismatch due to database. Insect genera get Viridiplantae due to synonyms
heirarchy <- heirarchy %>% dplyr::filter(kingdom !="Viridiplantae" )


# Count again 
heirarchy <- heirarchy %>%  group_by(OTU) %>% mutate(counts=n())
# We will accept only unique ASV's 
heirarchy <- heirarchy %>% dplyr::filter(counts == 1)


# Get ASV names where counts is not 1

not_assigned <- heirarchy %>% dplyr::filter(counts!=1)

seqtab_finalb <- seqtab_final
seqtab_finalb <- as.data.frame(seqtab_finalb)
names(seqtab_finalb) <-paste0("ASV",seq_along(seqtab_finalb))

#### Collapse ASV's with taxonomy ####
tax <- tibble::enframe(colnames(seqtab_finalb)) %>%
  dplyr::select(-name, OTU = value) %>%
  left_join(top_hit)%>%
  left_join(heirarchy) %>%
  distinct()


# Remove species assignments with "/"
tax[grep("/",tax$species),]$species <- NA

 

# Get those that don't have NA in all columns but OTU
taxb <- tax[rowSums(is.na(tax[,-1])) != ncol(tax[,-1]),]


taxb <- tax[rowSums(is.na(tax[,-c(1,2)])) != ncol(tax[,-c(1,2)]),]

# Remove rows where there is a genus but no other heirarchy (Gillateridae-synonym for Adelges)
taxb <- taxb %>%  group_by(OTU) %>% mutate(counts=n())

taxb <- taxb %>% dplyr::filter(counts == 1)
taxb <- taxb %>% dplyr::filter(genus != "class")
taxb <- taxb %>% mutate(taxon = case_when(is.na(species) ~ "genus",
                                          !is.na(genus)  ~ "species"))



# Replace taxonomy on asv table

seqtab_finalb <- t(seqtab_finalb)
seqtab_finalb <- as.data.frame(seqtab_finalb)
seqtab_finalb$OTU <- rownames(seqtab_finalb)

# How many reads did the ASV's that had multiple hits?
sum(colSums(seqtab_finalb[seqtab_finalb$OTU %in% unique(not_assigned$OTU),-188]))
#118657

joined_tablesb <- inner_join(taxb,seqtab_finalb,by="OTU")


write.csv(joined_tablesb,"Data/Processed/taxa_asv_259.csv",row.names = F)


##### Create fasta from ASV table

fasta <- DNAStringSet(colnames(seqtab_final))

names_fasta<- paste0("ASV",seq_along(colnames(seqtab_final)))

names(fasta) <- names_fasta               

write_fasta(fasta,file="Data/Processed/aphids_fasta.fasta")




##### Check if controls are within the samples ####


### Check for cross contamination of PCR negative or PCR positive ### 
a <-seqtab_finalb["00ExPos1_forward.paired.fastq.gz"]

df2=a[which(apply(a,1,function(x) x != 0)),]
ex_pos_1 <- as.data.frame(seqtab_finalb[c("ASV85","ASV1171","ASV1257"),])
colSums(ex_pos_1[,1:187])


b <- seqtab_finalb["00ExPos2_forward.paired.fastq.gz"]
df3=b[which(apply(b,1,function(x) x != 0)),]
ex_pos_2 <- as.data.frame(seqtab_finalb[c("ASV44","ASV1091"),])
colSums(ex_pos_2[,1:187])

c <-seqtab_finalb["00PCRpos1_forward.paired.fastq.gz"]

df4=a[which(apply(a,1,function(x) x != 0)),]
pcr_pos_1 <- as.data.frame(seqtab_finalb["ASV66",])
colSums(ex_pos_1[,1:187])

d<-seqtab_finalb["00PCRpos2_forward.paired.fastq.gz"]

df4=a[which(apply(a,1,function(x) x != 0)),]
pcr_pos_2 <- as.data.frame(seqtab_finalb["ASV66"])
colSums(ex_pos_1[,1:187])

















