# PNNL R Skills Assessment Script
# This script is going to map Proteomic data
# Visualilzation will be (optimisitcally) implemanted in a shiny app
#   as the data includes phosphosite positions over 7000 genes


library("rio")
library("dplyr")
library("R.utils")
library("BiocManager")
library("readr")
library("stringr")
library("seqinr")
library("biomaRt")
library("SparkR")
library("vegan")

# data frame of refseq and peptide id
phosphopeptides <- read.delim("./phosphopeptides.txt", header = TRUE, sep = "\t", dec = ".")

# data maping peptide id to gene id
# getting the biomaRt marts to access the refseq to gene mapping
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# Peptide RefSeq to gene mappipng, called inline later to save memory
# getBM(attributes = c("refseq_peptide", "external_gene_name"),
#       filters = "refseq_peptide", values = unique(phosphopeptides$RefSeq)
#       mart = ensembl)

# adding genes each refseq is in to our data frame
phosphopeptides <- dplyr::left_join(phosphopeptides, getBM(attributes = c("refseq_peptide", "external_gene_name"),
                                                           filters = "refseq_peptide", values = unique(phosphopeptides$RefSeq),
                                                           mart = ensembl),
                                    by = c("RefSeq" = "refseq_peptide"))

####################################################################################

# removes all RefSeq ids that do not have a corresponding gene in the database loaded
phosphopeptides <- phosphopeptides[!is.na(phosphopeptides$external_gene_name),]

####################################################################################

# Finding locations of the peptides

# loading the fasta file which is the amino acid sequences associated to each refseq
# and putting those into a data frame
RefSeq <- read.fasta("./H_sapiens_RefSeq.fasta.gz", seqtype = "AA", seqonly = FALSE,
                     strip.desc = TRUE, as.string = TRUE)
RefSeqDF <- data.frame("RefSeq" = names(RefSeq), "Sequence" = paste(RefSeq), "Seq_length" = str_length(paste(RefSeq)))
# creating a master data frame by RefSeq id
masterDF <- dplyr::left_join(phosphopeptides, RefSeqDF, by = c("RefSeq" = "RefSeq"))

# creating PepString in the masterDF:
# the Peptide sequence without the punctuation to search the AA sequence later
masterDF$PepString = str_replace_all(masterDF$Peptide, "[:punct:]", "")

# creating LocationSeq in masterDF:
# The Peptide sequence without any punctuation except for the "*" which is the
# phosphorylation sites for each respective peptide
masterDF$LocationSeq = str_replace_all(phosphopeptides$Peptide, "\\.", "")
#masterDF$LocationSeq = str_replace_all(phosphopeptides$LocationSeq, "[\-]", "")

#####################################################################################

# mapping phosphorylation locations

##############################################################################

# First site in each peptide

# position of the first sites for each peptide
masterDF$SeqStart <-(as.vector(str_locate(masterDF$Sequence, masterDF$PepString)[,1]) +
                       as.vector(str_locate(masterDF$LocationSeq, "[*]")[,1]))
# creating the phosphorylation site to the masterDF
masterDF$temp <- paste0(str_sub(masterDF$Sequence, masterDF$SeqStart-1, masterDF$SeqStart-1),
                        masterDF$SeqStart)
masterDF$Site <- paste0(masterDF$external_gene_name, "-", masterDF$temp)

# adding the phosphorylation site to the plotting data frame
# keeping track of all the sites in a long string of the all specific sites
# separated by commas in order to process them in plotting later using str_split
gene_plot_DF <- dplyr::select(masterDF, external_gene_name, temp, Seq_length) %>%
                dplyr::group_by(external_gene_name) %>%
                dplyr::mutate(Sites = paste(temp, collapse = ",")) %>%
                dplyr::select(external_gene_name, Seq_length, Sites) %>%
                dplyr::distinct()

############################################################################

# Second phosphorylation site mapping for all applicable peptides

# LocationSeq substring excluding everything before the first *
# in order to have only the second peptide bond location
# or the sequence does not contain "*"
shortLocationSeq <- str_sub(masterDF$LocationSeq, 
                            as.vector(str_locate(masterDF$LocationSeq, "[*]")[,1])+1)


# if the peptide has multiple sites:
# find the location of the second phosphorylation site
# add the second site to the masterDF

masterDF$SeqStart <-  ifelse(str_detect(shortLocationSeq, "[*]"), 
                             masterDF$SeqStart + as.vector(str_locate(shortLocationSeq, "[*]")[,1]), 0)

masterDF$temp <- ifelse(masterDF$SeqStart !=0, 
                        paste0(str_sub(masterDF$Sequence, masterDF$SeqStart-1, 
                                       masterDF$SeqStart-1), masterDF$SeqStart),
                        "")

# completes the site mapping on masterDF
masterDF$Site <- paste0(masterDF$Site, masterDF$temp)




# second site data frame
gene_plot_DF2 <- dplyr::select(masterDF, external_gene_name, temp) %>%
  dplyr::filter(temp != "") %>%
  dplyr::group_by(external_gene_name) %>%
  dplyr::mutate(Sites2 = paste(temp, collapse = ",")) %>%
  dplyr::select(external_gene_name, Sites2) %>%
  dplyr::distinct()

gene_plot_DF3 <- dplyr::left_join(gene_plot_DF, gene_plot_DF2, by = "external_gene_name") 

# combining all the sites if there are more to add
gene_plot_DF3$allSites <- ifelse(!is.na(gene_plot_DF3$Sites2), 
                                 paste(gene_plot_DF3$Sites, gene_plot_DF3$Sites2, sep = ","),
                                 gene_plot_DF3$Sites)

names(masterDF)[names(masterDF) == "external_gene_name"] <- "Gene"

dplyr::select(masterDF, Gene, RefSeq, Peptide, Site) %>%
  write.csv(file = "Mapping.csv")

dplyr::select(gene_plot_DF3, external_gene_name, Seq_length, Sites) %>%
  write.csv(file = "plottingdata.csv")



