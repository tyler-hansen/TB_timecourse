library(limma)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(dendextend)
library(reshape2)
library(edgeR)
library(qvalue)
library(cowplot)
library(stats)
library(ngram)

setwd("/Users/Haley/Desktop/lab/code/time_course/QC/")

meta_data <- read.table("meta_data_time_course.txt", header = T, sep = ",", check.names = FALSE)
reads <- read.table("GE_matrix_time_course_uncorrected_raw_counts.txt", header = T, sep = ",")

length(which(colnames(reads)!=rownames(meta_data)))

meta_data$African_admixture[is.na(meta_data$African_admixture)] <- 0.000010
meta_data$flow_cell <- as.factor(meta_data$flow_cell)
meta_data$lane <- as.factor(meta_data$lane)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
meta_data$infection = factor(meta_data$infection, levels=c("NI","Mtb_MOI_5"))
meta_data$time_infection <- as.factor(paste0(meta_data$time_point_hr, meta_data$infection))

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")
# remove columns that are 0hr and MTB_moi_5 because this doesn't make any sense
# these were samples that were collected immediately after mtb infection, so should be functionally the same as the NI time 0 controls
remove <- c("0Mtb_MOI_5")
samples_to_remove <- meta_data[meta_data$time_infection %in% remove,]
names_to_remove <- as.character(samples_to_remove$sample_ID_GE)
reads <- reads[, -which(names(reads) %in% c(names_to_remove))]
meta_data <- meta_data[meta_data$sample_ID_GE %in% colnames(reads), ]
length(which(colnames(reads)!=rownames(meta_data)))

## remove all AF145 samples (10) because all clusters incorrectly on a UMAP
individual <- c("AF145")
ind_remove <- meta_data[meta_data$individual %in% individual,]
names_ind_remove <- as.character(ind_remove$sample_ID_GE)
reads <- reads[, -which(names(reads) %in% c(names_ind_remove))]
meta_data <- meta_data[meta_data$sample_ID_GE %in% colnames(reads), ]
length(which(colnames(reads)!=rownames(meta_data)))


## remove AF95 18hr ctl because this sample clusters with infected on PCA
ctl <- c("5_AF95_T18_NI")
ctl_remove <- meta_data[meta_data$sample_ID %in% ctl,]
names_ctl_remove <- as.character(ctl_remove$sample_ID_GE)
reads <- reads[, -which(names(reads) %in% c(names_ctl_remove))]
meta_data <- meta_data[meta_data$sample_ID_GE %in% colnames(reads), ]
length(which(colnames(reads)!=rownames(meta_data)))

## 361 samples
write.table(meta_data, "meta_data_GOOD_SAMPLES.txt", sep = ",", quote = FALSE)
write.table(reads, "GE_uncorrected_raw_counts_GOOD_SAMPLES.txt", sep = ",", quote = FALSE)


## calculate and store corrected expression matrix (log2)
dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)
design = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale, data = meta_data)

v <- voom(dge, design, plot = TRUE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

corr_expression <- v$E - vfit$coefficients[,"flow_cell2"]%*%t(design[,"flow_cell2"]) - vfit$coefficients[,"flow_cell3"]%*%t(design[,"flow_cell3"]) - vfit$coefficients[,"perc_Aligned_scale"]%*%t(design[,"perc_Aligned_scale"]) - vfit$coefficients[,"perc_GC_scale"]%*%t(design[,"perc_GC_scale"]) - vfit$coefficients[,"perc_Dups_scale"]%*%t(design[,"perc_Dups_scale"]) 

write.table(corr_expression, "corrected_expression_log2CPM_voom.txt", sep = ",", quote = FALSE)









