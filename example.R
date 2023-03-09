#### MOCK EXAMPLE####
# Illumina Test - Test with two samples 
#
# written by: José A. Castillo
#
# Last modified Jan,
# First written March,2023
# The analysis is based in the tutorial (http://rstudio-pubs-static.s3.amazonaws.com/273862_8258982fdfac49a8ba9d6ed651b9b539.html)

####  Introduction  ####


##### Establecer el directorio del Proyecto ####
getwd()

##### Configurar sesión de R ####
#### 1. Enlistar los paquetes necesarios en diferentes vectores dependindo de su reporitorio de origen ####

## Repositorio CRAN
cran_packages <- c("bookdown", "knitr", "tidyverse","plyr","grid","gridExtra", "kableExtra","xtable","ggpubr")

## Repositorio Bioconductor
bio_packages <- c("phyloseq","dada2","DECIPHER","phangorn","ggpubr","BiocManager","DESeq2","microbiome","philr","ShortRead", "Biostrings")

# Repositorio GitHub
git_source <- c('twbattaglia/btools',"gmteunisse/fantaxtic","MadsAlbertsen/ampvis2","opisthokonta/tsnemicrobiota")

# fuente/nombre del paquete

git_packages <- c("btools", "fantaxtic", "ampvis2","tsnemicrobiota" )



#### 2. Instalar los paquetes definidos arriba usando la funcion corrspondiente a cada repositorio ####
remotes::install_github("kasperskytte/ampvis2")
install.packages('devtools')

if(!"devtools" %in% installed.packages()){
  install.packages("devtools")
}
devtools::install_github("gmteunisse/fantaxtic")  
devtools::install_github('twbattaglia/btools')
library(btools)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead",force=TRUE)


install_github("opisthokonta/tsnemicrobiota")


# Instalar paquetes CRAN
.inst <- cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(cran_packages[!.inst])
}

# Instalar paquetes BioConductor
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BioManager")
.inst <- bio_packages %in% installed.packages()
if(any(!.inst)) {
  BiocManager::install(bio_packages[!.inst])
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genefilter")
BiocManager::install("phyloseq")


install.packages('devtools')
devtools::install_github('twbattaglia/btools')
library(btools)

#### 3. Cargar los paquetes requeridos a la sesion actual de R ####
sapply(c(cran_packages, bioc_packages, git_packages), require, character.only = TRUE)

# NOTE: Si los paquetes ya estan instalados en tu computadora, solo necesitas enlistar (1er script) y cargar (3er script) los paquetes

#### Functions for Primer Identification ####

# If your dada2 version shouldbe 1.4 or higher.
# The data we will work with are the same as those in the Mothur Miseq SOP 
# walkthrough. Download the example data and unzip. These files represent 
# longitudinal samples from a mouse post-weaning and one mock community control. 
# For now just consider them paired-end fastq files to be processed. Define the 
# following path variable so that it points to the extracted directory on your 
# machine:

#This example uses data from


path <- "/home/umbreon/Documents/Estancia_Inves_UTPL/NGS/Example/data/MiSeq_SOP"
list.files


##### FILTER AND TRIM ####

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)



plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


##### Perform filtering and trimming ####

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))



out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


head(out)


# Start the clock
ptm <- proc.time()
# Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
# Stop the clock
proc.time() - ptm



# Start the clock
ptm <- proc.time()
# Learn error rates
errR <- learnErrors(filtRs, multithread=TRUE)
# Stop the clock
proc.time() - ptm



plotErrors(errF, nominalQ=TRUE)



##### Dereplication ####
derepFs <- derepFastq(filtFs, verbose=TRUE)

derepRs <- derepFastq(filtRs, verbose=TRUE)


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#### Sample Inference ####
# Start the clock
ptm <- proc.time()
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Stop the clock
proc.time() - ptm


dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]


##### Merge paired reads #####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])



seqtab <- makeSequenceTable(mergers)
dim(seqtab)


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))



hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")


##### Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


##### Track reads through the pipeline ####

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

#### Assign taxonomy  #####
### Formatting custom databases ###
taxtrain <- "/home/umbreon/Documents/Estancia_Inves_UTPL/NGS/Example/data/silva_nr_v123_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, taxtrain, multithread=TRUE)
unname(head(taxa))


#### Handoff to phyloseq ####
#### Import into phyloseq: ####

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")


# Make a data.frame holding the sample data
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps

#### Plot the species richness ####

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When") + theme_bw()



##### Create ordination plots ####
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

plot_ordination(ps, ord.nmds.bray, color="When", title="Bray NMDS")



#### Bar plot ####
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")







