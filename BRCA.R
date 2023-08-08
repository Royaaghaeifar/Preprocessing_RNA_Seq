
# loading librariies:
suppressMessages({
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(survival)
library(survminer)
library(tidyverse)
library(DESeq2)
  })

setwd("/Users/royaaghaeifar/Downloads/Disertation/gene/Clean Data/BRCA/")

# We can see list of projects in GDC 
projects <- getGDCprojects()


# Download the data
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  access = "open")

GDCdownload(query, method = "api")
data <- GDCprepare(query = query, save = TRUE, save.filename = "exp.rda")
rna <- as.data.frame(SummarizedExperiment::assay(data))
clinic <- data.frame(data@colData)


# Save the data in a folder
write.csv(rna, file =  "Raw/gene_brca.csv", row.names = T)
write.csv(clinic, file = "Raw/brca_clinic.csv", row.names = FALSE)


### import data
# I saved the data to save time 
clinic <- read_csv("Raw/brca_clinic.csv")


# I created two variables that can be used for the design variable in DEseq 
clinic <- clinic %>% mutate(normal = ifelse(substr(barcode,14,14) == '1', "N","T"),
                            normal = as.factor(normal))

clinic <- clinic %>% mutate(type = case_when(sample_type == "Metastatic" ~ "A",
                                            sample_type == "Primary Tumor" ~ "B",
                                            TRUE ~ "C"),
                            type = as.factor(type))


# define survival time and event in clinical dataset -------------------------
#Check data contains days in the clinical dataset
colnames(clinic)[grep("days", colnames(clinic))]


# Define time to event 
clinic <- clinic %>% 
  mutate(time_to_event = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death)) %>%
  mutate(time_to_event= ifelse(time_to_event == 0, NA, time_to_event))

# Data for event can be found under column vital_status 
table(clinic$vital_status)

# Remove a sample with missing vital status
clinic$barcode[is.na(clinic$vital_status)] 
clinic <- clinic %>% filter(barcode != "TCGA-BH-A0B2-01A-11R-A10J-07")

#recoding vital_status
clinic <- clinic %>%
  mutate(event = ifelse(vital_status == "Alive", 0, 1))
         

# create a subset from original clinical data
clinic <-clinic[, c("barcode", "sample_type", "age_at_diagnosis", "type", 
                    "age_at_index", "initial_weight", "race", "gender", "normal",
                     "ethnicity", "time_to_event", "event")]


# Analyze rna data
rna <- read.csv("Raw/brca_genes.csv", row.names = 1)

# save row names in a vector
genename <- rownames(rna)


# remove genes with no vital status
rna$`TCGA-BH-A0B2-01A-11R-A10J-07` <- NULL

# visualize the raw counts
plotDensities(rna, col = "darkblue", main="Raw Data Expression Distribution", xlab= "Intensity Values", legend=F)

# remove low expressed genes using the DESeq2 method
dds <- DESeqDataSetFromMatrix(countData = rna, colData = clinic, design = ~ normal)


des <- DESeq(dds)


res <- results(des, contrast = c("normal", "N", "T"))


# visualizing the counts of reads for a single gene across the groups
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="sample_type", 
                returnData=TRUE)

ggplot(d, aes(x= sample_type, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0), color = "darkblue") + 
  scale_y_log10(breaks=c(25,100,400))+
  xlab("Sample Type")+
  ylab("Log10(count)")+
  ggtitle("The Counts of reads across Sample Type")+
  theme(plot.title = element_text(hjust = 0.5, face="bold"))


# keep genes based on Dseq threshold
keep <- res$baseMean > metadata(res)$filterThreshold
                                           
# remove genes with low expression
dds <- dds[keep,]

# Perform VST normalization
dds_vst <- vst(dds)


# to check for batch effects and the like.
pcaData <- plotPCA(dds_vst, intgroup=c("sample_type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=sample_type, shape=sample_type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


#Access the normalized expression values
norm_data <- data.frame(assay(dds_vst))
hist(assay(dds_vst), xlab='Intensity Values', main= "VST Normalization")

plotDensities(assay(dds_vst), col = "darkblue", main="VST Normalization Visualization", xlab= "Intensity Values", legend=F)


# add gene name
gene_brca <- gene_brca[keep,]
norm_data$gene_id <- gene_brca$gene_id

library(vsn)
meanSdPlot(as.matrix(norm_data), ylim =c(0, 2.5))

write.csv(norm_data, file = "brca_genes.csv", row.names=FALSE)


