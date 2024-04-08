##########
## Script for DADA2 trimming and cleaning of Rothamsted samples 
##########


library(dada2)
packageVersion("dada2")
library(tidyverse)
packageVersion("tidyverse")
library(Biostrings)
packageVersion("Biostrings")
library(ShortRead)
packageVersion("ShortRead")

path <- "/data/Raw/extracted_fastqs/"

set.seed(01109291)

head(list.files(path))

fnFs <- sort(list.files(path, pattern="_forward.paired.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_reverse.paired.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "\\_"), `[`, 1)
sample.names

#Primer Sequences
FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"
        
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

setwd(path)


fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

fnFs.filtN <- fnFs.filtN[file.exists(fnFs.filtN)] # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- fnRs.filtN[file.exists(fnRs.filtN)] 


#### add new names after the ones deleted
path <- paste0(path,"/filtN")

fnFs <- sort(list.files(path, pattern="_forward.paired.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_reverse.paired.fastq.gz", full.names = TRUE))


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Check that cutadapt is findable
#cutadapt <- "~/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2("cutadapt", args = "--version") # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt

for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--cores=40", # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_forward.paired.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_reverse.paired.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
sample.names <- sapply(strsplit(basename(cutFs), "\\_"), `[`, 1)
head(sample.names)

#pdf(paste(home,"QualityProfile.pdf", sep=""))
#plotQualityProfile(fnFs[3])
#plotQualityProfile(fnRs[3])
#dev.off()


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_forward.paired.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_reverse.paired.fastq.gz"))


#filtFs <-filtFs[file.exists(filtFs)]
#filtRs <-filtRs[file.exists(filtRs)]
#cutFs <- cutFs[file.exists(cutFs)]
#cutRs <- cutRs[file.exists(cutRs)]


sample.names
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 200, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

out 

filtFs <-filtFs[file.exists(filtFs)]
filtRs <-filtRs[file.exists(filtRs)]



errF <- learnErrors(filtFs,multithread = T)
errR <- learnErrors(filtRs,multithread = T)

dadaFs <- dada(filtFs,err=errF,multithread = T)
dadaRs <- dada(filtRs,err = errR,multithread = T)

mergers <- mergePairs(dadaFs,filtFs,dadaRs,filtRs,verbose=F)

seqtab <- makeSequenceTable(mergers)

table(nchar(getSequences(seqtab)))


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))

track_filter <- out[file.exists(filtFs),]

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")

head(track)
track <- as.data.frame(track)
track_filter <- as.data.frame(track_filter)

track_combined <- cbind(out[file.exists(filtFs),],sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#Write files in intermediate

write.table(seqtab.nochim,"./seqtab_nochim_22_seed_aphids_24.tsv")
write.table(track,"./seqtab_stats_22_seed_aphids_24.tsv")
write.table(track_filter,"~/seqtab_stats_filter_24.tsv")
write.table(track_combined,"~/track_combined_24.tsv")
list.files(getwd())
write.table(seqtab.nochim,"~/seqtab_nochim_22_seed_aphids_24.tsv")
write.table(track,"~/seqtab_stats_22_seed_aphids_24.tsv")
getwd()
list.files(getwd())

