
# loading librariies:
suppressMessages({
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(tidyverse)
library(DESeq2)
  })

# Set work directory
dir <- "C:/Users/aghaer01/Downloads/WISE/Gene/BRCA"

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
write.csv(rna, file = paste0(dir, "Raw/gene_brca.csv"), row.names = T)
write.csv(clinic, file = paste0(dir, "Raw/brca_clinic.csv"), row.names = FALSE)


### import data
# I saved the data to save time 
# Analyze rna data
rna <- read.csv(paste0(dir, "Raw/brca_genes.csv"), row.names = 1)

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

write.csv(norm_data, file = paste0(dir, "brca_genes.csv"), row.names=FALSE)


#Analyze clinical data
clinic <- read_csv(paste0(dir, "Raw/brca_clinic.csv"))


# I created two variables that can be used for the design variable in DEseq 
clinic <- clinic %>% mutate(normal = ifelse(substr(barcode,14,14) == '1', "N","T"),
                            normal = as.factor(normal))

clinic <- clinic %>% mutate(type = case_when(sample_type == "Metastatic" ~ "M",
                                             sample_type == "Primary Tumor" ~ "T",
                                             TRUE ~ "N"),
                            type = as.factor(type))


# define survival time and event in clinical dataset -------------------------
#Check data contains days in the clinical dataset
colnames(clinic)[grep("days", colnames(clinic))]


# Define time to event 
clinic <- clinic %>% 
  mutate(time_to_event = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
         time_to_event= ifelse(time_to_event == 0, NA, time_to_event),
         time_to_event = ifelse(is.na(time_to_event), paper_days_to_last_followup, time_to_event))

# Data for event can be found under column vital_status 
table(clinic$vital_status)

# Remove a sample with missing vital status
clinic$barcode[is.na(clinic$vital_status)] 
clinic <- clinic %>% filter(barcode != "TCGA-BH-A0B2-01A-11R-A10J-07")

#recoding vital_status
clinic <- clinic %>%
  mutate(event = ifelse(vital_status == "Alive", 0, 1))

clinic <- clinic %>%
  mutate(age_at_diagnosis= round(age_at_diagnosis/365, 2))


# create a subset from original clinical data
# clinic <-clinic[, c("barcode", "sample_type", "age_at_diagnosis", "type", 
#                     "age_at_index", "initial_weight", "race", "gender", "normal",
#                      "ethnicity", "time_to_event", "event")]


# Survival Analysis
# creates a survival object 

# install.packages(c("ggsurvfit", "gtsummary", "tidycmprsk"))
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

# devtools::install_github("zabore/condsurv")
library(condsurv)


Surv(clinic$time_to_event, clinic$event)[1:10]

# survival curves using the Kaplan-Meier method. Time and surve are important output
#time: the timepoints at which the curve has a step, i.e. at least one event occurred
#surv: the estimate of survival at the corresponding time

s1 <- survfit(Surv(time_to_event, event) ~ 1, data = clinic)
str(s1)


survfit2(Surv(time_to_event, event) ~ 1, data = clinic) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  )+
  add_confidence_interval()+
  add_risktable() # the numbers at risk

# estimate the probability of surviving beyond a certain time
summary(survfit(Surv(time_to_event, event) ~ 1, data = clinic), times = 3500.25)


# Put the results on a table
survfit(Surv(time_to_event, event) ~ 1, data = clinic) %>% 
  tbl_survfit(
    times = 3500.25,
    label_header = "**1-year survival (95% CI)**"
  )


# Median survival time estimates
survfit(Surv(time_to_event, event) ~ 1, data = clinic) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )


# Compare survival between two groups
survfit2(Surv(time_to_event, event) ~ normal, data = clinic) %>% 
  ggsurvfit() +
  labs(
    x = "Days to Event",
    y = "Overall survival probability"
  ) +
  add_risktable()


# Check if there is a significant difference in overall survival according to sample type 
survdiff(Surv(time_to_event, event) ~ sample_type, data = clinic)


# CPH regression model; 
#exp = T gives the hazard ratio instead of log hazard ratio
#The HR shows the instantaneous rate of occurrence of the event of interest 
#in those who are still at risk for the event. It is regression parameter not a risk)
# HR > 1 indicates an increased hazard of death.
coxph(Surv(time_to_event, event) ~ normal + age_at_diagnosis + initial_weight,
        data = clinic) %>%
         tbl_regression(exp = TRUE)


# Assessing CPH
cph <- coxph(Surv(time_to_event, event) ~ sample_type + age_at_diagnosis + initial_weight,
             data = clinic)
cz <- cox.zph(cph)

#A significant p-value indicates that the proportional hazards assumption is violated
#If p-values >0.05, we do not reject H0
#We conclude that the PH assumption is satisfied for each covariate and for the model overall
print(cz)


plot(cz)


# Estimate the CPH for a specific group
coxph(Surv(time_to_event, event) ~  age_at_diagnosis+ initial_weight, 
           subset = sample_type == "Solid Tissue Normal", data = clinic) %>% 
             tbl_regression(exp = TRUE)






# Shows the cumulative incidence of death from BC
cuminc(Surv(time_to_event, as.factor(event)) ~ normal, data = clinic) %>% 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()

# test for a difference between groups over the entire follow-up period
cuminc(Surv(time_to_event, as.factor(event)) ~ normal, data = clinic) %>% 
  tbl_cuminc(
    times = 1826.25, 
    label_header = "**{time/365.25}-year cuminc**") %>% 
  add_p()


# estimate the subdistribution hazards: estimates CPH for those who have not yet experienced the event
crr(Surv(time_to_event, as.factor(event)) ~  age_at_diagnosis, data = clinic) %>% 
  tbl_regression(exp = TRUE)
