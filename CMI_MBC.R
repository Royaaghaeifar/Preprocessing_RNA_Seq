# Project: CMI MBC
rm(list = ls())

#Install packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.16")
# 
# BiocManager::install(c("TCGAbiolinks"))

#loading libraries:
library(TCGAbiolinks)
library(limma)
library(tidyverse)
library(biomaRt)


#Downloading the data with TCGABiolink
query <- GDCquery(project = "CMI-MBC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  access = "open"
)


GDCdownload(query, method = "api")
data <- GDCprepare(query = query, save = TRUE, save.filename = "exp.rda")
rna <- as.data.frame(SummarizedExperiment::assay(data))
clinical <- data.frame(data@colData)

#Preprocessing step
clinical$primary_site <- vapply(clinical$primary_site, paste, collapse = " ", character(1L))
clinical$disease_type <- vapply(clinical$disease_type, paste, collapse = " ", character(1L))


# Check if data is normalized
ggplot(rna, aes(x= MBCProject_3946_T4_RNA)) + geom_density(color="darkblue")+
  labs(title="Raw Data Expression Value Distribution",x="Intensity Value", y = "Density")+
  theme(plot.title = element_text(hjust = 0.5))


# perform log2
rna_norm <- sapply(rna, function(x) log2(as.numeric(x)))
rna_norm <- as.data.frame(rna_norm)
rna_norm[rna_norm == -Inf] <- 0
  


ggplot(rna_norm, aes(x= MBCProject_3896_T1_RNA)) + geom_density(color="darkblue")+
  labs(title="Log2 Data Expression Vlue Distribution",x="Intensity Value", y = "Density")+
  theme(plot.title = element_text(hjust = 0.5))

#10595 genes have no expression in the rna matrix
table(rowSums(rna_norm[,0:203]) == 0)

# add gene names
rna_norm$gene <- rownames(rna)

#remove genes with no-expression in more 50% of samples:
no_exp <- data.frame(count = apply(rna, 1, function(x) length(which(x== 0))))
no_exp$gene <- row.names(no_exp)
no_exp <- no_exp %>% filter(count > dim(rna)[2]/2)
rna <- rna %>% filter(!(row.names(rna) %in% no_exp$gene)) 









# find genes associated with ensemble codes
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "ensembl_gene_id",
  values = rownames(rna),
  uniqueRows=TRUE)



library("illuminaHumanv4.db")
gene_id <- data.frame(Gene=unlist(mget(x = rownames(rna), envir = illuminaHumanv4SYMBOL)))



