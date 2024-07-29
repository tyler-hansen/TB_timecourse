library(limma)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(edgeR)
library(qvalue)
library(cowplot)
library(stats)
library(ngram)
library(RColorBrewer)
library(dendextend)
library(reshape2)
library(rlist)
library(gridExtra)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")

meta_data = read.table("meta_data_GOOD_SAMPLES_with_ccScores.txt", header = T, sep = ",", check.names = FALSE)
reads = read.table("GE_uncorrected_raw_counts_GOOD_SAMPLES.txt", header = T, sep = ",")

length(which(colnames(reads)!=rownames(meta_data)))

meta_data$internal_ID <- as.factor(meta_data$internal_ID)
meta_data$sample_number <- as.factor(meta_data$sample_number)
meta_data$flow_cell <- as.factor(meta_data$flow_cell)
meta_data$lane <- as.factor(meta_data$lane)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
meta_data$infection = factor(meta_data$infection, levels=c("NI","Mtb_MOI_5"))

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_dup_ctls")
## duplicate non-stimulated samples for those in which we don't have control data
# first, let's duplicate the GE data 
rownames(meta_data[meta_data$infection == "NI",])

dup_reads <- reads

dup_reads$X16_EU238_T1_NI <- reads$X16_EU238_T0_NI
dup_reads$X16_EU238_T2_NI <- reads$X16_EU238_T0_NI
dup_reads$X16_EU238_T3_NI <- reads$X16_EU238_T4_NI
dup_reads$X16_EU238_T5_NI <- reads$X16_EU238_T4_NI
dup_reads$X16_EU238_T6_NI <- reads$X16_EU238_T4_NI
dup_reads$X16_EU238_T7_NI <- reads$X16_EU238_T8_NI
dup_reads$X16_EU238_T9_NI <- reads$X16_EU238_T8_NI
dup_reads$X16_EU238_T10_NI <- reads$X16_EU238_T12_NI
dup_reads$X16_EU238_T14_NI <- reads$X16_EU238_T12_NI
dup_reads$X16_EU238_T16_NI <- reads$X16_EU238_T18_NI
dup_reads$X16_EU238_T42_NI <- reads$X16_EU238_T36_NI

dup_reads$X15_EU148_T1_NI <- reads$X15_EU148_T0_NI
dup_reads$X15_EU148_T2_NI <- reads$X15_EU148_T0_NI
dup_reads$X15_EU148_T3_NI <- reads$X15_EU148_T4_NI
dup_reads$X15_EU148_T5_NI <- reads$X15_EU148_T4_NI
dup_reads$X15_EU148_T6_NI <- reads$X15_EU148_T4_NI
dup_reads$X15_EU148_T7_NI <- reads$X15_EU148_T8_NI
dup_reads$X15_EU148_T9_NI <- reads$X15_EU148_T8_NI
dup_reads$X15_EU148_T10_NI <- reads$X15_EU148_T12_NI
dup_reads$X15_EU148_T14_NI <- reads$X15_EU148_T12_NI
dup_reads$X15_EU148_T16_NI <- reads$X15_EU148_T18_NI
dup_reads$X15_EU148_T42_NI <- reads$X15_EU148_T36_NI

dup_reads$X14_EU144_T1_NI <- reads$X14_EU144_T4_NI
dup_reads$X14_EU144_T2_NI <- reads$X14_EU144_T4_NI
dup_reads$X14_EU144_T3_NI <- reads$X14_EU144_T4_NI
dup_reads$X14_EU144_T5_NI <- reads$X14_EU144_T4_NI
dup_reads$X14_EU144_T6_NI <- reads$X14_EU144_T4_NI
dup_reads$X14_EU144_T7_NI <- reads$X14_EU144_T8_NI
dup_reads$X14_EU144_T9_NI <- reads$X14_EU144_T8_NI
dup_reads$X14_EU144_T10_NI <- reads$X14_EU144_T12_NI
dup_reads$X14_EU144_T14_NI <- reads$X14_EU144_T12_NI
dup_reads$X14_EU144_T16_NI <- reads$X14_EU144_T18_NI
dup_reads$X14_EU144_T42_NI <- reads$X14_EU144_T36_NI

dup_reads$X13_EU140_T1_NI <- reads$X13_EU140_T0_NI
dup_reads$X13_EU140_T2_NI <- reads$X13_EU140_T0_NI
dup_reads$X13_EU140_T3_NI <- reads$X13_EU140_T4_NI
dup_reads$X13_EU140_T5_NI <- reads$X13_EU140_T4_NI
dup_reads$X13_EU140_T6_NI <- reads$X13_EU140_T4_NI
dup_reads$X13_EU140_T7_NI <- reads$X13_EU140_T8_NI
dup_reads$X13_EU140_T9_NI <- reads$X13_EU140_T8_NI
dup_reads$X13_EU140_T10_NI <- reads$X13_EU140_T12_NI
dup_reads$X13_EU140_T14_NI <- reads$X13_EU140_T12_NI
dup_reads$X13_EU140_T16_NI <- reads$X13_EU140_T18_NI
dup_reads$X13_EU140_T42_NI <- reads$X13_EU140_T36_NI

dup_reads$X12_EU122_T1_NI <- reads$X12_EU122_T0_NI
dup_reads$X12_EU122_T2_NI <- reads$X12_EU122_T0_NI
dup_reads$X12_EU122_T3_NI <- reads$X12_EU122_T4_NI
dup_reads$X12_EU122_T5_NI <- reads$X12_EU122_T4_NI
dup_reads$X12_EU122_T6_NI <- reads$X12_EU122_T4_NI
dup_reads$X12_EU122_T7_NI <- reads$X12_EU122_T8_NI
dup_reads$X12_EU122_T9_NI <- reads$X12_EU122_T8_NI
dup_reads$X12_EU122_T10_NI <- reads$X12_EU122_T12_NI
dup_reads$X12_EU122_T14_NI <- reads$X12_EU122_T12_NI
dup_reads$X12_EU122_T16_NI <- reads$X12_EU122_T18_NI
dup_reads$X12_EU122_T42_NI <- reads$X12_EU122_T36_NI

dup_reads$X10_EU118_T1_NI <- reads$X10_EU118_T0_NI
dup_reads$X10_EU118_T2_NI <- reads$X10_EU118_T0_NI
dup_reads$X10_EU118_T3_NI <- reads$X10_EU118_T4_NI
dup_reads$X10_EU118_T5_NI <- reads$X10_EU118_T4_NI
dup_reads$X10_EU118_T6_NI <- reads$X10_EU118_T4_NI
dup_reads$X10_EU118_T7_NI <- reads$X10_EU118_T8_NI
dup_reads$X10_EU118_T9_NI <- reads$X10_EU118_T8_NI
dup_reads$X10_EU118_T10_NI <- reads$X10_EU118_T12_NI
dup_reads$X10_EU118_T14_NI <- reads$X10_EU118_T12_NI
dup_reads$X10_EU118_T16_NI <- reads$X10_EU118_T18_NI
dup_reads$X10_EU118_T42_NI <- reads$X10_EU118_T36_NI

dup_reads$X18_EU03_T1_NI <- reads$X18_EU03_T0_NI
dup_reads$X18_EU03_T2_NI <- reads$X18_EU03_T0_NI
dup_reads$X18_EU03_T3_NI <- reads$X18_EU03_T4_NI
dup_reads$X18_EU03_T5_NI <- reads$X18_EU03_T4_NI
dup_reads$X18_EU03_T6_NI <- reads$X18_EU03_T4_NI
dup_reads$X18_EU03_T7_NI <- reads$X18_EU03_T8_NI
dup_reads$X18_EU03_T9_NI <- reads$X18_EU03_T8_NI
dup_reads$X18_EU03_T10_NI <- reads$X18_EU03_T12_NI
dup_reads$X18_EU03_T14_NI <- reads$X18_EU03_T12_NI
dup_reads$X18_EU03_T16_NI <- reads$X18_EU03_T18_NI
dup_reads$X18_EU03_T42_NI <- reads$X18_EU03_T36_NI

dup_reads$X17_EU02_T1_NI <- reads$X17_EU02_T0_NI
dup_reads$X17_EU02_T2_NI <- reads$X17_EU02_T0_NI
dup_reads$X17_EU02_T3_NI <- reads$X17_EU02_T4_NI
dup_reads$X17_EU02_T5_NI <- reads$X17_EU02_T4_NI
dup_reads$X17_EU02_T6_NI <- reads$X17_EU02_T4_NI
dup_reads$X17_EU02_T7_NI <- reads$X17_EU02_T8_NI
dup_reads$X17_EU02_T9_NI <- reads$X17_EU02_T8_NI
dup_reads$X17_EU02_T10_NI <- reads$X17_EU02_T12_NI
dup_reads$X17_EU02_T14_NI <- reads$X17_EU02_T12_NI
dup_reads$X17_EU02_T16_NI <- reads$X17_EU02_T18_NI
dup_reads$X17_EU02_T42_NI <- reads$X17_EU02_T36_NI

dup_reads$X5_AF95_T1_NI <- reads$X5_AF95_T0_NI
dup_reads$X5_AF95_T2_NI <- reads$X5_AF95_T0_NI
dup_reads$X5_AF95_T3_NI <- reads$X5_AF95_T4_NI
dup_reads$X5_AF95_T5_NI <- reads$X5_AF95_T4_NI
dup_reads$X5_AF95_T6_NI <- reads$X5_AF95_T4_NI
dup_reads$X5_AF95_T7_NI <- reads$X5_AF95_T8_NI
dup_reads$X5_AF95_T9_NI <- reads$X5_AF95_T8_NI
dup_reads$X5_AF95_T10_NI <- reads$X5_AF95_T12_NI
dup_reads$X5_AF95_T14_NI <- reads$X5_AF95_T12_NI
dup_reads$X5_AF95_T16_NI <- reads$X5_AF95_T12_NI
dup_reads$X5_AF95_T42_NI <- reads$X5_AF95_T36_NI

dup_reads$X4_AF69_T1_NI <- reads$X4_AF69_T0_NI
dup_reads$X4_AF69_T2_NI <- reads$X4_AF69_T0_NI
dup_reads$X4_AF69_T3_NI <- reads$X4_AF69_T4_NI
dup_reads$X4_AF69_T5_NI <- reads$X4_AF69_T4_NI
dup_reads$X4_AF69_T6_NI <- reads$X4_AF69_T4_NI
dup_reads$X4_AF69_T7_NI <- reads$X4_AF69_T8_NI
dup_reads$X4_AF69_T9_NI <- reads$X4_AF69_T8_NI
dup_reads$X4_AF69_T10_NI <- reads$X4_AF69_T12_NI
dup_reads$X4_AF69_T14_NI <- reads$X4_AF69_T12_NI
dup_reads$X4_AF69_T16_NI <- reads$X4_AF69_T18_NI
dup_reads$X4_AF69_T42_NI <- reads$X4_AF69_T36_NI

dup_reads$X3_AF55_T1_NI <- reads$X3_AF55_T0_NI
dup_reads$X3_AF55_T2_NI <- reads$X3_AF55_T0_NI
dup_reads$X3_AF55_T3_NI <- reads$X3_AF55_T4_NI
dup_reads$X3_AF55_T5_NI <- reads$X3_AF55_T4_NI
dup_reads$X3_AF55_T6_NI <- reads$X3_AF55_T4_NI
dup_reads$X3_AF55_T7_NI <- reads$X3_AF55_T8_NI
dup_reads$X3_AF55_T9_NI <- reads$X3_AF55_T8_NI
dup_reads$X3_AF55_T10_NI <- reads$X3_AF55_T12_NI
dup_reads$X3_AF55_T14_NI <- reads$X3_AF55_T12_NI
dup_reads$X3_AF55_T16_NI <- reads$X3_AF55_T18_NI
dup_reads$X3_AF55_T42_NI <- reads$X3_AF55_T36_NI

dup_reads$X9_AF193_T1_NI <- reads$X9_AF193_T0_NI
dup_reads$X9_AF193_T2_NI <- reads$X9_AF193_T0_NI
dup_reads$X9_AF193_T3_NI <- reads$X9_AF193_T4_NI
dup_reads$X9_AF193_T5_NI <- reads$X9_AF193_T4_NI
dup_reads$X9_AF193_T6_NI <- reads$X9_AF193_T4_NI
dup_reads$X9_AF193_T7_NI <- reads$X9_AF193_T8_NI
dup_reads$X9_AF193_T9_NI <- reads$X9_AF193_T8_NI
dup_reads$X9_AF193_T10_NI <- reads$X9_AF193_T12_NI
dup_reads$X9_AF193_T16_NI <- reads$X9_AF193_T18_NI
dup_reads$X9_AF193_T42_NI <- reads$X9_AF193_T36_NI

dup_reads$X8_AF183_T1_NI <- reads$X8_AF183_T0_NI
dup_reads$X8_AF183_T2_NI <- reads$X8_AF183_T0_NI
dup_reads$X8_AF183_T3_NI <- reads$X8_AF183_T4_NI
dup_reads$X8_AF183_T5_NI <- reads$X8_AF183_T4_NI
dup_reads$X8_AF183_T6_NI <- reads$X8_AF183_T4_NI
dup_reads$X8_AF183_T7_NI <- reads$X8_AF183_T8_NI
dup_reads$X8_AF183_T9_NI <- reads$X8_AF183_T8_NI
dup_reads$X8_AF183_T10_NI <- reads$X8_AF183_T12_NI
dup_reads$X8_AF183_T14_NI <- reads$X8_AF183_T12_NI
dup_reads$X8_AF183_T16_NI <- reads$X8_AF183_T18_NI
dup_reads$X8_AF183_T42_NI <- reads$X8_AF183_T36_NI



# next, duplicate the meta data 
dup_meta_data <- subset(meta_data, select = -c(sample_ID, sample_ID_GE, sample_number, internal_ID))

## X16_EU238
X16_EU238_T1_NI <- dup_meta_data[c("X16_EU238_T0_NI"),]
X16_EU238_T1_NI$time_point_hr[1] <- 1
rownames(X16_EU238_T1_NI) <- c("X16_EU238_T1_NI")

X16_EU238_T2_NI <- dup_meta_data[c("X16_EU238_T0_NI"),]
X16_EU238_T2_NI$time_point_hr[1] <- 2
rownames(X16_EU238_T2_NI) <- c("X16_EU238_T2_NI")

X16_EU238_T3_NI <- dup_meta_data[c("X16_EU238_T4_NI"),]
X16_EU238_T3_NI$time_point_hr[1] <- 3
rownames(X16_EU238_T3_NI) <- c("X16_EU238_T3_NI")

X16_EU238_T5_NI <- dup_meta_data[c("X16_EU238_T4_NI"),]
X16_EU238_T5_NI$time_point_hr[1] <- 5
rownames(X16_EU238_T5_NI) <- c("X16_EU238_T5_NI")

X16_EU238_T6_NI <- dup_meta_data[c("X16_EU238_T4_NI"),]
X16_EU238_T6_NI$time_point_hr[1] <- 6
rownames(X16_EU238_T6_NI) <- c("X16_EU238_T6_NI")

X16_EU238_T7_NI <- dup_meta_data[c("X16_EU238_T8_NI"),]
X16_EU238_T7_NI$time_point_hr[1] <- 7
rownames(X16_EU238_T7_NI) <- c("X16_EU238_T7_NI")

X16_EU238_T9_NI <- dup_meta_data[c("X16_EU238_T8_NI"),]
X16_EU238_T9_NI$time_point_hr[1] <- 9
rownames(X16_EU238_T9_NI) <- c("X16_EU238_T9_NI")

X16_EU238_T10_NI <- dup_meta_data[c("X16_EU238_T12_NI"),]
X16_EU238_T10_NI$time_point_hr[1] <- 10
rownames(X16_EU238_T10_NI) <- c("X16_EU238_T10_NI")

X16_EU238_T14_NI <- dup_meta_data[c("X16_EU238_T12_NI"),]
X16_EU238_T14_NI$time_point_hr[1] <- 14
rownames(X16_EU238_T14_NI) <- c("X16_EU238_T14_NI")

X16_EU238_T16_NI <- dup_meta_data[c("X16_EU238_T18_NI"),]
X16_EU238_T16_NI$time_point_hr[1] <- 16
rownames(X16_EU238_T16_NI) <- c("X16_EU238_T16_NI")

X16_EU238_T42_NI <- dup_meta_data[c("X16_EU238_T36_NI"),]
X16_EU238_T42_NI$time_point_hr[1] <- 42
rownames(X16_EU238_T42_NI) <- c("X16_EU238_T42_NI")

X16_EU238 <- rbind(X16_EU238_T1_NI, X16_EU238_T2_NI, X16_EU238_T3_NI, X16_EU238_T5_NI, X16_EU238_T6_NI, X16_EU238_T7_NI, X16_EU238_T9_NI, X16_EU238_T10_NI, X16_EU238_T14_NI, X16_EU238_T16_NI, X16_EU238_T42_NI)




## X15_EU148
X15_EU148_T1_NI <- dup_meta_data[c("X15_EU148_T0_NI"),]
X15_EU148_T1_NI$time_point_hr[1] <- 1
rownames(X15_EU148_T1_NI) <- c("X15_EU148_T1_NI")

X15_EU148_T2_NI <- dup_meta_data[c("X15_EU148_T0_NI"),]
X15_EU148_T2_NI$time_point_hr[1] <- 2
rownames(X15_EU148_T2_NI) <- c("X15_EU148_T2_NI")

X15_EU148_T3_NI <- dup_meta_data[c("X15_EU148_T4_NI"),]
X15_EU148_T3_NI$time_point_hr[1] <- 3
rownames(X15_EU148_T3_NI) <- c("X15_EU148_T3_NI")

X15_EU148_T5_NI <- dup_meta_data[c("X15_EU148_T4_NI"),]
X15_EU148_T5_NI$time_point_hr[1] <- 5
rownames(X15_EU148_T5_NI) <- c("X15_EU148_T5_NI")

X15_EU148_T6_NI <- dup_meta_data[c("X15_EU148_T4_NI"),]
X15_EU148_T6_NI$time_point_hr[1] <- 6
rownames(X15_EU148_T6_NI) <- c("X15_EU148_T6_NI")

X15_EU148_T7_NI <- dup_meta_data[c("X15_EU148_T8_NI"),]
X15_EU148_T7_NI$time_point_hr[1] <- 7
rownames(X15_EU148_T7_NI) <- c("X15_EU148_T7_NI")

X15_EU148_T9_NI <- dup_meta_data[c("X15_EU148_T8_NI"),]
X15_EU148_T9_NI$time_point_hr[1] <- 9
rownames(X15_EU148_T9_NI) <- c("X15_EU148_T9_NI")

X15_EU148_T10_NI <- dup_meta_data[c("X15_EU148_T12_NI"),]
X15_EU148_T10_NI$time_point_hr[1] <- 10
rownames(X15_EU148_T10_NI) <- c("X15_EU148_T10_NI")

X15_EU148_T14_NI <- dup_meta_data[c("X15_EU148_T12_NI"),]
X15_EU148_T14_NI$time_point_hr[1] <- 14
rownames(X15_EU148_T14_NI) <- c("X15_EU148_T14_NI")

X15_EU148_T16_NI <- dup_meta_data[c("X15_EU148_T18_NI"),]
X15_EU148_T16_NI$time_point_hr[1] <- 16
rownames(X15_EU148_T16_NI) <- c("X15_EU148_T16_NI")

X15_EU148_T42_NI <- dup_meta_data[c("X15_EU148_T36_NI"),]
X15_EU148_T42_NI$time_point_hr[1] <- 42
rownames(X15_EU148_T42_NI) <- c("X15_EU148_T42_NI")

X15_EU148 <- rbind(X15_EU148_T1_NI, X15_EU148_T2_NI, X15_EU148_T3_NI, X15_EU148_T5_NI, X15_EU148_T6_NI, X15_EU148_T7_NI, X15_EU148_T9_NI, X15_EU148_T10_NI, X15_EU148_T14_NI, X15_EU148_T16_NI, X15_EU148_T42_NI)




## X14_EU144
X14_EU144_T1_NI <- dup_meta_data[c("X14_EU144_T4_NI"),]
X14_EU144_T1_NI$time_point_hr[1] <- 1
rownames(X14_EU144_T1_NI) <- c("X14_EU144_T1_NI")

X14_EU144_T2_NI <- dup_meta_data[c("X14_EU144_T4_NI"),]
X14_EU144_T2_NI$time_point_hr[1] <- 2
rownames(X14_EU144_T2_NI) <- c("X14_EU144_T2_NI")

X14_EU144_T3_NI <- dup_meta_data[c("X14_EU144_T4_NI"),]
X14_EU144_T3_NI$time_point_hr[1] <- 3
rownames(X14_EU144_T3_NI) <- c("X14_EU144_T3_NI")

X14_EU144_T5_NI <- dup_meta_data[c("X14_EU144_T4_NI"),]
X14_EU144_T5_NI$time_point_hr[1] <- 5
rownames(X14_EU144_T5_NI) <- c("X14_EU144_T5_NI")

X14_EU144_T6_NI <- dup_meta_data[c("X14_EU144_T4_NI"),]
X14_EU144_T6_NI$time_point_hr[1] <- 6
rownames(X14_EU144_T6_NI) <- c("X14_EU144_T6_NI")

X14_EU144_T7_NI <- dup_meta_data[c("X14_EU144_T8_NI"),]
X14_EU144_T7_NI$time_point_hr[1] <- 7
rownames(X14_EU144_T7_NI) <- c("X14_EU144_T7_NI")

X14_EU144_T9_NI <- dup_meta_data[c("X14_EU144_T8_NI"),]
X14_EU144_T9_NI$time_point_hr[1] <- 9
rownames(X14_EU144_T9_NI) <- c("X14_EU144_T9_NI")

X14_EU144_T10_NI <- dup_meta_data[c("X14_EU144_T12_NI"),]
X14_EU144_T10_NI$time_point_hr[1] <- 10
rownames(X14_EU144_T10_NI) <- c("X14_EU144_T10_NI")

X14_EU144_T14_NI <- dup_meta_data[c("X14_EU144_T12_NI"),]
X14_EU144_T14_NI$time_point_hr[1] <- 14
rownames(X14_EU144_T14_NI) <- c("X14_EU144_T14_NI")

X14_EU144_T16_NI <- dup_meta_data[c("X14_EU144_T18_NI"),]
X14_EU144_T16_NI$time_point_hr[1] <- 16
rownames(X14_EU144_T16_NI) <- c("X14_EU144_T16_NI")

X14_EU144_T42_NI <- dup_meta_data[c("X14_EU144_T36_NI"),]
X14_EU144_T42_NI$time_point_hr[1] <- 42
rownames(X14_EU144_T42_NI) <- c("X14_EU144_T42_NI")

X14_EU144 <- rbind(X14_EU144_T1_NI, X14_EU144_T2_NI, X14_EU144_T3_NI, X14_EU144_T5_NI, X14_EU144_T6_NI, X14_EU144_T7_NI, X14_EU144_T9_NI, X14_EU144_T10_NI, X14_EU144_T14_NI, X14_EU144_T16_NI, X14_EU144_T42_NI)




## X13_EU140
X13_EU140_T1_NI <- dup_meta_data[c("X13_EU140_T0_NI"),]
X13_EU140_T1_NI$time_point_hr[1] <- 1
rownames(X13_EU140_T1_NI) <- c("X13_EU140_T1_NI")

X13_EU140_T2_NI <- dup_meta_data[c("X13_EU140_T0_NI"),]
X13_EU140_T2_NI$time_point_hr[1] <- 2
rownames(X13_EU140_T2_NI) <- c("X13_EU140_T2_NI")

X13_EU140_T3_NI <- dup_meta_data[c("X13_EU140_T4_NI"),]
X13_EU140_T3_NI$time_point_hr[1] <- 3
rownames(X13_EU140_T3_NI) <- c("X13_EU140_T3_NI")

X13_EU140_T5_NI <- dup_meta_data[c("X13_EU140_T4_NI"),]
X13_EU140_T5_NI$time_point_hr[1] <- 5
rownames(X13_EU140_T5_NI) <- c("X13_EU140_T5_NI")

X13_EU140_T6_NI <- dup_meta_data[c("X13_EU140_T4_NI"),]
X13_EU140_T6_NI$time_point_hr[1] <- 6
rownames(X13_EU140_T6_NI) <- c("X13_EU140_T6_NI")

X13_EU140_T7_NI <- dup_meta_data[c("X13_EU140_T8_NI"),]
X13_EU140_T7_NI$time_point_hr[1] <- 7
rownames(X13_EU140_T7_NI) <- c("X13_EU140_T7_NI")

X13_EU140_T9_NI <- dup_meta_data[c("X13_EU140_T8_NI"),]
X13_EU140_T9_NI$time_point_hr[1] <- 9
rownames(X13_EU140_T9_NI) <- c("X13_EU140_T9_NI")

X13_EU140_T10_NI <- dup_meta_data[c("X13_EU140_T12_NI"),]
X13_EU140_T10_NI$time_point_hr[1] <- 10
rownames(X13_EU140_T10_NI) <- c("X13_EU140_T10_NI")

X13_EU140_T14_NI <- dup_meta_data[c("X13_EU140_T12_NI"),]
X13_EU140_T14_NI$time_point_hr[1] <- 14
rownames(X13_EU140_T14_NI) <- c("X13_EU140_T14_NI")

X13_EU140_T16_NI <- dup_meta_data[c("X13_EU140_T18_NI"),]
X13_EU140_T16_NI$time_point_hr[1] <- 16
rownames(X13_EU140_T16_NI) <- c("X13_EU140_T16_NI")

X13_EU140_T42_NI <- dup_meta_data[c("X13_EU140_T36_NI"),]
X13_EU140_T42_NI$time_point_hr[1] <- 42
rownames(X13_EU140_T42_NI) <- c("X13_EU140_T42_NI")

X13_EU140 <- rbind(X13_EU140_T1_NI, X13_EU140_T2_NI, X13_EU140_T3_NI, X13_EU140_T5_NI, X13_EU140_T6_NI, X13_EU140_T7_NI, X13_EU140_T9_NI, X13_EU140_T10_NI, X13_EU140_T14_NI, X13_EU140_T16_NI, X13_EU140_T42_NI)




## X12_EU122
X12_EU122_T1_NI <- dup_meta_data[c("X12_EU122_T0_NI"),]
X12_EU122_T1_NI$time_point_hr[1] <- 1
rownames(X12_EU122_T1_NI) <- c("X12_EU122_T1_NI")

X12_EU122_T2_NI <- dup_meta_data[c("X12_EU122_T0_NI"),]
X12_EU122_T2_NI$time_point_hr[1] <- 2
rownames(X12_EU122_T2_NI) <- c("X12_EU122_T2_NI")

X12_EU122_T3_NI <- dup_meta_data[c("X12_EU122_T4_NI"),]
X12_EU122_T3_NI$time_point_hr[1] <- 3
rownames(X12_EU122_T3_NI) <- c("X12_EU122_T3_NI")

X12_EU122_T5_NI <- dup_meta_data[c("X12_EU122_T4_NI"),]
X12_EU122_T5_NI$time_point_hr[1] <- 5
rownames(X12_EU122_T5_NI) <- c("X12_EU122_T5_NI")

X12_EU122_T6_NI <- dup_meta_data[c("X12_EU122_T4_NI"),]
X12_EU122_T6_NI$time_point_hr[1] <- 6
rownames(X12_EU122_T6_NI) <- c("X12_EU122_T6_NI")

X12_EU122_T7_NI <- dup_meta_data[c("X12_EU122_T8_NI"),]
X12_EU122_T7_NI$time_point_hr[1] <- 7
rownames(X12_EU122_T7_NI) <- c("X12_EU122_T7_NI")

X12_EU122_T9_NI <- dup_meta_data[c("X12_EU122_T8_NI"),]
X12_EU122_T9_NI$time_point_hr[1] <- 9
rownames(X12_EU122_T9_NI) <- c("X12_EU122_T9_NI")

X12_EU122_T10_NI <- dup_meta_data[c("X12_EU122_T12_NI"),]
X12_EU122_T10_NI$time_point_hr[1] <- 10
rownames(X12_EU122_T10_NI) <- c("X12_EU122_T10_NI")

X12_EU122_T14_NI <- dup_meta_data[c("X12_EU122_T12_NI"),]
X12_EU122_T14_NI$time_point_hr[1] <- 14
rownames(X12_EU122_T14_NI) <- c("X12_EU122_T14_NI")

X12_EU122_T16_NI <- dup_meta_data[c("X12_EU122_T18_NI"),]
X12_EU122_T16_NI$time_point_hr[1] <- 16
rownames(X12_EU122_T16_NI) <- c("X12_EU122_T16_NI")

X12_EU122_T42_NI <- dup_meta_data[c("X12_EU122_T36_NI"),]
X12_EU122_T42_NI$time_point_hr[1] <- 42
rownames(X12_EU122_T42_NI) <- c("X12_EU122_T42_NI")

X12_EU122 <- rbind(X12_EU122_T1_NI, X12_EU122_T2_NI, X12_EU122_T3_NI, X12_EU122_T5_NI, X12_EU122_T6_NI, X12_EU122_T7_NI, X12_EU122_T9_NI, X12_EU122_T10_NI, X12_EU122_T14_NI, X12_EU122_T16_NI, X12_EU122_T42_NI)




## X10_EU118
X10_EU118_T1_NI <- dup_meta_data[c("X10_EU118_T0_NI"),]
X10_EU118_T1_NI$time_point_hr[1] <- 1
rownames(X10_EU118_T1_NI) <- c("X10_EU118_T1_NI")

X10_EU118_T2_NI <- dup_meta_data[c("X10_EU118_T0_NI"),]
X10_EU118_T2_NI$time_point_hr[1] <- 2
rownames(X10_EU118_T2_NI) <- c("X10_EU118_T2_NI")

X10_EU118_T3_NI <- dup_meta_data[c("X10_EU118_T4_NI"),]
X10_EU118_T3_NI$time_point_hr[1] <- 3
rownames(X10_EU118_T3_NI) <- c("X10_EU118_T3_NI")

X10_EU118_T5_NI <- dup_meta_data[c("X10_EU118_T4_NI"),]
X10_EU118_T5_NI$time_point_hr[1] <- 5
rownames(X10_EU118_T5_NI) <- c("X10_EU118_T5_NI")

X10_EU118_T6_NI <- dup_meta_data[c("X10_EU118_T4_NI"),]
X10_EU118_T6_NI$time_point_hr[1] <- 6
rownames(X10_EU118_T6_NI) <- c("X10_EU118_T6_NI")

X10_EU118_T7_NI <- dup_meta_data[c("X10_EU118_T8_NI"),]
X10_EU118_T7_NI$time_point_hr[1] <- 7
rownames(X10_EU118_T7_NI) <- c("X10_EU118_T7_NI")

X10_EU118_T9_NI <- dup_meta_data[c("X10_EU118_T8_NI"),]
X10_EU118_T9_NI$time_point_hr[1] <- 9
rownames(X10_EU118_T9_NI) <- c("X10_EU118_T9_NI")

X10_EU118_T10_NI <- dup_meta_data[c("X10_EU118_T12_NI"),]
X10_EU118_T10_NI$time_point_hr[1] <- 10
rownames(X10_EU118_T10_NI) <- c("X10_EU118_T10_NI")

X10_EU118_T14_NI <- dup_meta_data[c("X10_EU118_T12_NI"),]
X10_EU118_T14_NI$time_point_hr[1] <- 14
rownames(X10_EU118_T14_NI) <- c("X10_EU118_T14_NI")

X10_EU118_T16_NI <- dup_meta_data[c("X10_EU118_T18_NI"),]
X10_EU118_T16_NI$time_point_hr[1] <- 16
rownames(X10_EU118_T16_NI) <- c("X10_EU118_T16_NI")

X10_EU118_T42_NI <- dup_meta_data[c("X10_EU118_T36_NI"),]
X10_EU118_T42_NI$time_point_hr[1] <- 42
rownames(X10_EU118_T42_NI) <- c("X10_EU118_T42_NI")

X10_EU118 <- rbind(X10_EU118_T1_NI, X10_EU118_T2_NI, X10_EU118_T3_NI, X10_EU118_T5_NI, X10_EU118_T6_NI, X10_EU118_T7_NI, X10_EU118_T9_NI, X10_EU118_T10_NI, X10_EU118_T14_NI, X10_EU118_T16_NI, X10_EU118_T42_NI)




## X18_EU03
X18_EU03_T1_NI <- dup_meta_data[c("X18_EU03_T0_NI"),]
X18_EU03_T1_NI$time_point_hr[1] <- 1
rownames(X18_EU03_T1_NI) <- c("X18_EU03_T1_NI")

X18_EU03_T2_NI <- dup_meta_data[c("X18_EU03_T0_NI"),]
X18_EU03_T2_NI$time_point_hr[1] <- 2
rownames(X18_EU03_T2_NI) <- c("X18_EU03_T2_NI")

X18_EU03_T3_NI <- dup_meta_data[c("X18_EU03_T4_NI"),]
X18_EU03_T3_NI$time_point_hr[1] <- 3
rownames(X18_EU03_T3_NI) <- c("X18_EU03_T3_NI")

X18_EU03_T5_NI <- dup_meta_data[c("X18_EU03_T4_NI"),]
X18_EU03_T5_NI$time_point_hr[1] <- 5
rownames(X18_EU03_T5_NI) <- c("X18_EU03_T5_NI")

X18_EU03_T6_NI <- dup_meta_data[c("X18_EU03_T4_NI"),]
X18_EU03_T6_NI$time_point_hr[1] <- 6
rownames(X18_EU03_T6_NI) <- c("X18_EU03_T6_NI")

X18_EU03_T7_NI <- dup_meta_data[c("X18_EU03_T8_NI"),]
X18_EU03_T7_NI$time_point_hr[1] <- 7
rownames(X18_EU03_T7_NI) <- c("X18_EU03_T7_NI")

X18_EU03_T9_NI <- dup_meta_data[c("X18_EU03_T8_NI"),]
X18_EU03_T9_NI$time_point_hr[1] <- 9
rownames(X18_EU03_T9_NI) <- c("X18_EU03_T9_NI")

X18_EU03_T10_NI <- dup_meta_data[c("X18_EU03_T12_NI"),]
X18_EU03_T10_NI$time_point_hr[1] <- 10
rownames(X18_EU03_T10_NI) <- c("X18_EU03_T10_NI")

X18_EU03_T14_NI <- dup_meta_data[c("X18_EU03_T12_NI"),]
X18_EU03_T14_NI$time_point_hr[1] <- 14
rownames(X18_EU03_T14_NI) <- c("X18_EU03_T14_NI")

X18_EU03_T16_NI <- dup_meta_data[c("X18_EU03_T18_NI"),]
X18_EU03_T16_NI$time_point_hr[1] <- 16
rownames(X18_EU03_T16_NI) <- c("X18_EU03_T16_NI")

X18_EU03_T42_NI <- dup_meta_data[c("X18_EU03_T36_NI"),]
X18_EU03_T42_NI$time_point_hr[1] <- 42
rownames(X18_EU03_T42_NI) <- c("X18_EU03_T42_NI")

X18_EU03 <- rbind(X18_EU03_T1_NI, X18_EU03_T2_NI, X18_EU03_T3_NI, X18_EU03_T5_NI, X18_EU03_T6_NI, X18_EU03_T7_NI, X18_EU03_T9_NI, X18_EU03_T10_NI, X18_EU03_T14_NI, X18_EU03_T16_NI, X18_EU03_T42_NI)




## X17_EU02
X17_EU02_T1_NI <- dup_meta_data[c("X17_EU02_T0_NI"),]
X17_EU02_T1_NI$time_point_hr[1] <- 1
rownames(X17_EU02_T1_NI) <- c("X17_EU02_T1_NI")

X17_EU02_T2_NI <- dup_meta_data[c("X17_EU02_T0_NI"),]
X17_EU02_T2_NI$time_point_hr[1] <- 2
rownames(X17_EU02_T2_NI) <- c("X17_EU02_T2_NI")

X17_EU02_T3_NI <- dup_meta_data[c("X17_EU02_T4_NI"),]
X17_EU02_T3_NI$time_point_hr[1] <- 3
rownames(X17_EU02_T3_NI) <- c("X17_EU02_T3_NI")

X17_EU02_T5_NI <- dup_meta_data[c("X17_EU02_T4_NI"),]
X17_EU02_T5_NI$time_point_hr[1] <- 5
rownames(X17_EU02_T5_NI) <- c("X17_EU02_T5_NI")

X17_EU02_T6_NI <- dup_meta_data[c("X17_EU02_T4_NI"),]
X17_EU02_T6_NI$time_point_hr[1] <- 6
rownames(X17_EU02_T6_NI) <- c("X17_EU02_T6_NI")

X17_EU02_T7_NI <- dup_meta_data[c("X17_EU02_T8_NI"),]
X17_EU02_T7_NI$time_point_hr[1] <- 7
rownames(X17_EU02_T7_NI) <- c("X17_EU02_T7_NI")

X17_EU02_T9_NI <- dup_meta_data[c("X17_EU02_T8_NI"),]
X17_EU02_T9_NI$time_point_hr[1] <- 9
rownames(X17_EU02_T9_NI) <- c("X17_EU02_T9_NI")

X17_EU02_T10_NI <- dup_meta_data[c("X17_EU02_T12_NI"),]
X17_EU02_T10_NI$time_point_hr[1] <- 10
rownames(X17_EU02_T10_NI) <- c("X17_EU02_T10_NI")

X17_EU02_T14_NI <- dup_meta_data[c("X17_EU02_T12_NI"),]
X17_EU02_T14_NI$time_point_hr[1] <- 14
rownames(X17_EU02_T14_NI) <- c("X17_EU02_T14_NI")

X17_EU02_T16_NI <- dup_meta_data[c("X17_EU02_T18_NI"),]
X17_EU02_T16_NI$time_point_hr[1] <- 16
rownames(X17_EU02_T16_NI) <- c("X17_EU02_T16_NI")

X17_EU02_T42_NI <- dup_meta_data[c("X17_EU02_T36_NI"),]
X17_EU02_T42_NI$time_point_hr[1] <- 42
rownames(X17_EU02_T42_NI) <- c("X17_EU02_T42_NI")

X17_EU02 <- rbind(X17_EU02_T1_NI, X17_EU02_T2_NI, X17_EU02_T3_NI, X17_EU02_T5_NI, X17_EU02_T6_NI, X17_EU02_T7_NI, X17_EU02_T9_NI, X17_EU02_T10_NI, X17_EU02_T14_NI, X17_EU02_T16_NI, X17_EU02_T42_NI)




## X5_AF95
X5_AF95_T1_NI <- dup_meta_data[c("X5_AF95_T0_NI"),]
X5_AF95_T1_NI$time_point_hr[1] <- 1
rownames(X5_AF95_T1_NI) <- c("X5_AF95_T1_NI")

X5_AF95_T2_NI <- dup_meta_data[c("X5_AF95_T0_NI"),]
X5_AF95_T2_NI$time_point_hr[1] <- 2
rownames(X5_AF95_T2_NI) <- c("X5_AF95_T2_NI")

X5_AF95_T3_NI <- dup_meta_data[c("X5_AF95_T4_NI"),]
X5_AF95_T3_NI$time_point_hr[1] <- 3
rownames(X5_AF95_T3_NI) <- c("X5_AF95_T3_NI")

X5_AF95_T5_NI <- dup_meta_data[c("X5_AF95_T4_NI"),]
X5_AF95_T5_NI$time_point_hr[1] <- 5
rownames(X5_AF95_T5_NI) <- c("X5_AF95_T5_NI")

X5_AF95_T6_NI <- dup_meta_data[c("X5_AF95_T4_NI"),]
X5_AF95_T6_NI$time_point_hr[1] <- 6
rownames(X5_AF95_T6_NI) <- c("X5_AF95_T6_NI")

X5_AF95_T7_NI <- dup_meta_data[c("X5_AF95_T8_NI"),]
X5_AF95_T7_NI$time_point_hr[1] <- 7
rownames(X5_AF95_T7_NI) <- c("X5_AF95_T7_NI")

X5_AF95_T9_NI <- dup_meta_data[c("X5_AF95_T8_NI"),]
X5_AF95_T9_NI$time_point_hr[1] <- 9
rownames(X5_AF95_T9_NI) <- c("X5_AF95_T9_NI")

X5_AF95_T10_NI <- dup_meta_data[c("X5_AF95_T12_NI"),]
X5_AF95_T10_NI$time_point_hr[1] <- 10
rownames(X5_AF95_T10_NI) <- c("X5_AF95_T10_NI")

X5_AF95_T14_NI <- dup_meta_data[c("X5_AF95_T12_NI"),]
X5_AF95_T14_NI$time_point_hr[1] <- 14
rownames(X5_AF95_T14_NI) <- c("X5_AF95_T14_NI")

X5_AF95_T16_NI <- dup_meta_data[c("X5_AF95_T12_NI"),]
X5_AF95_T16_NI$time_point_hr[1] <- 16
rownames(X5_AF95_T16_NI) <- c("X5_AF95_T16_NI")

X5_AF95_T42_NI <- dup_meta_data[c("X5_AF95_T36_NI"),]
X5_AF95_T42_NI$time_point_hr[1] <- 42
rownames(X5_AF95_T42_NI) <- c("X5_AF95_T42_NI")

X5_AF95 <- rbind(X5_AF95_T1_NI, X5_AF95_T2_NI, X5_AF95_T3_NI, X5_AF95_T5_NI, X5_AF95_T6_NI, X5_AF95_T7_NI, X5_AF95_T9_NI, X5_AF95_T10_NI, X5_AF95_T14_NI, X5_AF95_T16_NI, X5_AF95_T42_NI)




## X4_AF69
X4_AF69_T1_NI <- dup_meta_data[c("X4_AF69_T0_NI"),]
X4_AF69_T1_NI$time_point_hr[1] <- 1
rownames(X4_AF69_T1_NI) <- c("X4_AF69_T1_NI")

X4_AF69_T2_NI <- dup_meta_data[c("X4_AF69_T0_NI"),]
X4_AF69_T2_NI$time_point_hr[1] <- 2
rownames(X4_AF69_T2_NI) <- c("X4_AF69_T2_NI")

X4_AF69_T3_NI <- dup_meta_data[c("X4_AF69_T4_NI"),]
X4_AF69_T3_NI$time_point_hr[1] <- 3
rownames(X4_AF69_T3_NI) <- c("X4_AF69_T3_NI")

X4_AF69_T5_NI <- dup_meta_data[c("X4_AF69_T4_NI"),]
X4_AF69_T5_NI$time_point_hr[1] <- 5
rownames(X4_AF69_T5_NI) <- c("X4_AF69_T5_NI")

X4_AF69_T6_NI <- dup_meta_data[c("X4_AF69_T4_NI"),]
X4_AF69_T6_NI$time_point_hr[1] <- 6
rownames(X4_AF69_T6_NI) <- c("X4_AF69_T6_NI")

X4_AF69_T7_NI <- dup_meta_data[c("X4_AF69_T8_NI"),]
X4_AF69_T7_NI$time_point_hr[1] <- 7
rownames(X4_AF69_T7_NI) <- c("X4_AF69_T7_NI")

X4_AF69_T9_NI <- dup_meta_data[c("X4_AF69_T8_NI"),]
X4_AF69_T9_NI$time_point_hr[1] <- 9
rownames(X4_AF69_T9_NI) <- c("X4_AF69_T9_NI")

X4_AF69_T10_NI <- dup_meta_data[c("X4_AF69_T12_NI"),]
X4_AF69_T10_NI$time_point_hr[1] <- 10
rownames(X4_AF69_T10_NI) <- c("X4_AF69_T10_NI")

X4_AF69_T14_NI <- dup_meta_data[c("X4_AF69_T12_NI"),]
X4_AF69_T14_NI$time_point_hr[1] <- 14
rownames(X4_AF69_T14_NI) <- c("X4_AF69_T14_NI")

X4_AF69_T16_NI <- dup_meta_data[c("X4_AF69_T18_NI"),]
X4_AF69_T16_NI$time_point_hr[1] <- 16
rownames(X4_AF69_T16_NI) <- c("X4_AF69_T16_NI")

X4_AF69_T42_NI <- dup_meta_data[c("X4_AF69_T36_NI"),]
X4_AF69_T42_NI$time_point_hr[1] <- 42
rownames(X4_AF69_T42_NI) <- c("X4_AF69_T42_NI")

X4_AF69 <- rbind(X4_AF69_T1_NI, X4_AF69_T2_NI, X4_AF69_T3_NI, X4_AF69_T5_NI, X4_AF69_T6_NI, X4_AF69_T7_NI, X4_AF69_T9_NI, X4_AF69_T10_NI, X4_AF69_T14_NI, X4_AF69_T16_NI, X4_AF69_T42_NI)




## X3_AF55
X3_AF55_T1_NI <- dup_meta_data[c("X3_AF55_T0_NI"),]
X3_AF55_T1_NI$time_point_hr[1] <- 1
rownames(X3_AF55_T1_NI) <- c("X3_AF55_T1_NI")

X3_AF55_T2_NI <- dup_meta_data[c("X3_AF55_T0_NI"),]
X3_AF55_T2_NI$time_point_hr[1] <- 2
rownames(X3_AF55_T2_NI) <- c("X3_AF55_T2_NI")

X3_AF55_T3_NI <- dup_meta_data[c("X3_AF55_T4_NI"),]
X3_AF55_T3_NI$time_point_hr[1] <- 3
rownames(X3_AF55_T3_NI) <- c("X3_AF55_T3_NI")

X3_AF55_T5_NI <- dup_meta_data[c("X3_AF55_T4_NI"),]
X3_AF55_T5_NI$time_point_hr[1] <- 5
rownames(X3_AF55_T5_NI) <- c("X3_AF55_T5_NI")

X3_AF55_T6_NI <- dup_meta_data[c("X3_AF55_T4_NI"),]
X3_AF55_T6_NI$time_point_hr[1] <- 6
rownames(X3_AF55_T6_NI) <- c("X3_AF55_T6_NI")

X3_AF55_T7_NI <- dup_meta_data[c("X3_AF55_T8_NI"),]
X3_AF55_T7_NI$time_point_hr[1] <- 7
rownames(X3_AF55_T7_NI) <- c("X3_AF55_T7_NI")

X3_AF55_T9_NI <- dup_meta_data[c("X3_AF55_T8_NI"),]
X3_AF55_T9_NI$time_point_hr[1] <- 9
rownames(X3_AF55_T9_NI) <- c("X3_AF55_T9_NI")

X3_AF55_T10_NI <- dup_meta_data[c("X3_AF55_T12_NI"),]
X3_AF55_T10_NI$time_point_hr[1] <- 10
rownames(X3_AF55_T10_NI) <- c("X3_AF55_T10_NI")

X3_AF55_T14_NI <- dup_meta_data[c("X3_AF55_T12_NI"),]
X3_AF55_T14_NI$time_point_hr[1] <- 14
rownames(X3_AF55_T14_NI) <- c("X3_AF55_T14_NI")

X3_AF55_T16_NI <- dup_meta_data[c("X3_AF55_T18_NI"),]
X3_AF55_T16_NI$time_point_hr[1] <- 16
rownames(X3_AF55_T16_NI) <- c("X3_AF55_T16_NI")

X3_AF55_T42_NI <- dup_meta_data[c("X3_AF55_T36_NI"),]
X3_AF55_T42_NI$time_point_hr[1] <- 42
rownames(X3_AF55_T42_NI) <- c("X3_AF55_T42_NI")

X3_AF55 <- rbind(X3_AF55_T1_NI, X3_AF55_T2_NI, X3_AF55_T3_NI, X3_AF55_T5_NI, X3_AF55_T6_NI, X3_AF55_T7_NI, X3_AF55_T9_NI, X3_AF55_T10_NI, X3_AF55_T14_NI, X3_AF55_T16_NI, X3_AF55_T42_NI)




## X9_AF193
X9_AF193_T1_NI <- dup_meta_data[c("X9_AF193_T0_NI"),]
X9_AF193_T1_NI$time_point_hr[1] <- 1
rownames(X9_AF193_T1_NI) <- c("X9_AF193_T1_NI")

X9_AF193_T2_NI <- dup_meta_data[c("X9_AF193_T0_NI"),]
X9_AF193_T2_NI$time_point_hr[1] <- 2
rownames(X9_AF193_T2_NI) <- c("X9_AF193_T2_NI")

X9_AF193_T3_NI <- dup_meta_data[c("X9_AF193_T4_NI"),]
X9_AF193_T3_NI$time_point_hr[1] <- 3
rownames(X9_AF193_T3_NI) <- c("X9_AF193_T3_NI")

X9_AF193_T5_NI <- dup_meta_data[c("X9_AF193_T4_NI"),]
X9_AF193_T5_NI$time_point_hr[1] <- 5
rownames(X9_AF193_T5_NI) <- c("X9_AF193_T5_NI")

X9_AF193_T6_NI <- dup_meta_data[c("X9_AF193_T4_NI"),]
X9_AF193_T6_NI$time_point_hr[1] <- 6
rownames(X9_AF193_T6_NI) <- c("X9_AF193_T6_NI")

X9_AF193_T7_NI <- dup_meta_data[c("X9_AF193_T8_NI"),]
X9_AF193_T7_NI$time_point_hr[1] <- 7
rownames(X9_AF193_T7_NI) <- c("X9_AF193_T7_NI")

X9_AF193_T9_NI <- dup_meta_data[c("X9_AF193_T8_NI"),]
X9_AF193_T9_NI$time_point_hr[1] <- 9
rownames(X9_AF193_T9_NI) <- c("X9_AF193_T9_NI")

X9_AF193_T10_NI <- dup_meta_data[c("X9_AF193_T12_NI"),]
X9_AF193_T10_NI$time_point_hr[1] <- 10
rownames(X9_AF193_T10_NI) <- c("X9_AF193_T10_NI")

X9_AF193_T16_NI <- dup_meta_data[c("X9_AF193_T18_NI"),]
X9_AF193_T16_NI$time_point_hr[1] <- 16
rownames(X9_AF193_T16_NI) <- c("X9_AF193_T16_NI")

X9_AF193_T42_NI <- dup_meta_data[c("X9_AF193_T36_NI"),]
X9_AF193_T42_NI$time_point_hr[1] <- 42
rownames(X9_AF193_T42_NI) <- c("X9_AF193_T42_NI")

X9_AF193 <- rbind(X9_AF193_T1_NI, X9_AF193_T2_NI, X9_AF193_T3_NI, X9_AF193_T5_NI, X9_AF193_T6_NI, X9_AF193_T7_NI, X9_AF193_T9_NI, X9_AF193_T10_NI, X9_AF193_T16_NI, X9_AF193_T42_NI)




## X8_AF183
X8_AF183_T1_NI <- dup_meta_data[c("X8_AF183_T0_NI"),]
X8_AF183_T1_NI$time_point_hr[1] <- 1
rownames(X8_AF183_T1_NI) <- c("X8_AF183_T1_NI")

X8_AF183_T2_NI <- dup_meta_data[c("X8_AF183_T0_NI"),]
X8_AF183_T2_NI$time_point_hr[1] <- 2
rownames(X8_AF183_T2_NI) <- c("X8_AF183_T2_NI")

X8_AF183_T3_NI <- dup_meta_data[c("X8_AF183_T4_NI"),]
X8_AF183_T3_NI$time_point_hr[1] <- 3
rownames(X8_AF183_T3_NI) <- c("X8_AF183_T3_NI")

X8_AF183_T5_NI <- dup_meta_data[c("X8_AF183_T4_NI"),]
X8_AF183_T5_NI$time_point_hr[1] <- 5
rownames(X8_AF183_T5_NI) <- c("X8_AF183_T5_NI")

X8_AF183_T6_NI <- dup_meta_data[c("X8_AF183_T4_NI"),]
X8_AF183_T6_NI$time_point_hr[1] <- 6
rownames(X8_AF183_T6_NI) <- c("X8_AF183_T6_NI")

X8_AF183_T7_NI <- dup_meta_data[c("X8_AF183_T8_NI"),]
X8_AF183_T7_NI$time_point_hr[1] <- 7
rownames(X8_AF183_T7_NI) <- c("X8_AF183_T7_NI")

X8_AF183_T9_NI <- dup_meta_data[c("X8_AF183_T8_NI"),]
X8_AF183_T9_NI$time_point_hr[1] <- 9
rownames(X8_AF183_T9_NI) <- c("X8_AF183_T9_NI")

X8_AF183_T10_NI <- dup_meta_data[c("X8_AF183_T12_NI"),]
X8_AF183_T10_NI$time_point_hr[1] <- 10
rownames(X8_AF183_T10_NI) <- c("X8_AF183_T10_NI")

X8_AF183_T14_NI <- dup_meta_data[c("X8_AF183_T12_NI"),]
X8_AF183_T14_NI$time_point_hr[1] <- 14
rownames(X8_AF183_T14_NI) <- c("X8_AF183_T14_NI")

X8_AF183_T16_NI <- dup_meta_data[c("X8_AF183_T18_NI"),]
X8_AF183_T16_NI$time_point_hr[1] <- 16
rownames(X8_AF183_T16_NI) <- c("X8_AF183_T16_NI")

X8_AF183_T42_NI <- dup_meta_data[c("X8_AF183_T36_NI"),]
X8_AF183_T42_NI$time_point_hr[1] <- 42
rownames(X8_AF183_T42_NI) <- c("X8_AF183_T42_NI")

X8_AF183 <- rbind(X8_AF183_T1_NI, X8_AF183_T2_NI, X8_AF183_T3_NI, X8_AF183_T5_NI, X8_AF183_T6_NI, X8_AF183_T7_NI, X8_AF183_T9_NI, X8_AF183_T10_NI, X8_AF183_T14_NI, X8_AF183_T16_NI, X8_AF183_T42_NI)


dup_meta_data <- rbind(dup_meta_data, X16_EU238, X15_EU148, X14_EU144, X13_EU140, X12_EU122, X10_EU118, X18_EU03, X17_EU02, X5_AF95, X4_AF69, X3_AF55, X9_AF193, X8_AF183)


length(which(colnames(dup_reads)!=rownames(dup_meta_data)))

write.table(dup_reads, "uncorrected_raw_counts_duplicated_controls_ABS_ctl.txt", quote = FALSE)
write.table(dup_meta_data, "meta_data_duplicated_contols_ABS_ctl_with_ccScores.txt", quote = FALSE)

## get corrected expression matrix
dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)
design = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale, data = dup_meta_data)

v <- voom(dge, design, plot = TRUE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

corr_expression <- v$E - vfit$coefficients[,"flow_cell2"]%*%t(design[,"flow_cell2"]) - vfit$coefficients[,"flow_cell3"]%*%t(design[,"flow_cell3"]) - vfit$coefficients[,"perc_Aligned_scale"]%*%t(design[,"perc_Aligned_scale"]) - vfit$coefficients[,"perc_GC_scale"]%*%t(design[,"perc_GC_scale"]) - vfit$coefficients[,"perc_Dups_scale"]%*%t(design[,"perc_Dups_scale"]) 

write.table(corr_expression, "corrected_expression_duplicated_controls_ABS_ctl.txt", sep = ",", quote = FALSE)



#################################
## begin actual analysis here ###
#################################
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_dup_ctls")
dup_reads <- read.table("uncorrected_raw_counts_duplicated_controls_ABS_ctl.txt")
dup_meta_data <- read.table("meta_data_duplicated_contols_ABS_ctl.txt")

length(which(colnames(dup_reads)!=rownames(dup_meta_data)))

dup_meta_data$flow_cell <- as.factor(dup_meta_data$flow_cell)
dup_meta_data$lane <- as.factor(dup_meta_data$lane)
dup_meta_data$time_point_hr <- as.factor(dup_meta_data$time_point_hr)
dup_meta_data$infection = factor(dup_meta_data$infection, levels=c("NI","Mtb_MOI_5"))
dup_meta_data$time_infection <- as.factor(paste0(dup_meta_data$time_point_hr, dup_meta_data$infection))

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/dupCtl_DE_response_Mtb")
dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)

design <- model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + individual + time_point_hr + infection:time_point_hr, data = dup_meta_data)

design <- design[, -c(38)]

v <- voom(dge, design, plot = FALSE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

betas = vfit$coefficients[, 19:56]
p_values = vfit$p.value[, 19:56]
fdrs = p_values
    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }


## let's subset on the MTB infection fdrs because these are what we care about right now
fdrs_inf = data.frame(fdrs[ , c(20:38)])
    
number_tp_sig <- function(x){
  sum(x < 0.05)
}

# what about genes that are significant in half of the time points? (10)
fdrs_inf$half_sig_fdr = apply(fdrs_inf, 1, number_tp_sig)
hits = rownames(fdrs_inf)[which(fdrs_inf$half_sig_fdr >= 10)]
# 8339 for 10 (50%) with a BH = 0.05

write.table(hits, "genes_DE_at_least_10_tps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# number of genes significant at each time point with an FDR (BH) of 5%
time_points <- as.data.frame(c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48))
colnames(time_points) <- c("time_point")
rownames(time_points) <- colnames(fdrs_inf[,1:19])

DE_genes <- c()

for(i in 1:19){
	DE_genes[i] <- sum(fdrs_inf[, i] < 0.05)
}

time_points <- cbind(time_points, DE_genes)

pdf("DE_genes_FDR_BH_0.05.pdf", height = 4, width = 8)
ggplot(time_points, aes(x = time_point, y = DE_genes, label = DE_genes)) +
		geom_line(alpha = 0.4) + 
		geom_point() +
		geom_text_repel(direction = "y", nudge_y = 1, nudge_x = 1) +
		xlab("Mtb infection (hr)") +
		ylab("number of DE genes") + 
		scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48)) +
		theme_bw() +
		theme(#panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
		ylim(0,11000) +
        geom_hline(yintercept = 0, linetype="dashed", color = "grey70")
dev.off()








######################################################
###### CHANGE ALL OF THIS DEPENDING ON NEW DATA ######
######################################################

## normalize data
tab = betas[which(rownames(betas) %in% hits), c(20:38)]
ttab = t(tab)
norm_ttab = sapply(1:ncol(ttab), function(x){ttab[,x]/max(abs(ttab[,x]))})
# normalizing logFC per gene, divide logFC by maximal ones so that max response time point is = 1 (genes that respond with same qualitative shape without considering intensity of response/no matter how strongly)
tab_norm = data.frame(t(norm_ttab))
rownames(tab_norm) = rownames(tab)


## dendrogram and clustering
clusters <- hclust(dist(tab_norm))
dend1 <- as.dendrogram(clusters)

## let's cut it at the height that leads to 10 clusters (this is arbitrary)
dend1 <- color_branches(dend1, k = 10)
dend1 <- color_labels(dend1, k = 10)

pdf("tree_labeled_k10.pdf", width = 40, height = 10)
plot(dend1)
dev.off()

ggd1 = as.ggdend(dend1)
arbol1 = ggplot(ggd1, labels = FALSE)

pdf("tree_k10.pdf", height = 30, width = 10)
print(arbol1)
dev.off()

pdf("radial_tree_k10.pdf", height = 30, width = 10)
arbol_radial = ggplot(ggd1, labels = FALSE, horiz = TRUE) + coord_polar(theta = "x")
dev.off()

# more parameters of dendrogram
ordered_tab <- tab_norm
ordered_genes <- labels(dend1)

ordered_tab = ordered_tab[ordered_genes, ]
ordered_tab$ordered_genes = factor(rownames(ordered_tab), levels = rownames(ordered_tab))
ordered_tab = ordered_tab[order(ordered_tab$ordered_genes), ]

# how many genes per cluster?
set_1 = which(ordered_tab$clusters_k_10==1)
# 82
set_2 = which(ordered_tab$clusters_k_10==2)
# 802
set_3 = which(ordered_tab$clusters_k_10==3)
# 862
set_4 = which(ordered_tab$clusters_k_10==4)
# 557
set_5 = which(ordered_tab$clusters_k_10==5)
# 441
set_6 = which(ordered_tab$clusters_k_10==6)
# 22
set_7 = which(ordered_tab$clusters_k_10==7)
# 149
set_8 = which(ordered_tab$clusters_k_10==8)
# 18
set_9 = which(ordered_tab$clusters_k_10==9)
# 14
set_10 = which(ordered_tab$clusters_k_10==10)
# 20


## cut to 10 clusters because this looks good
tab_norm$clusters_k_10 = data.frame(cutree(dend1, k = 10))[ ,1]

write.table(tab_norm, "normalized_betas_and_clusters_k_10_all_covars_regressed.txt")
results = list(fit = fit, betas = betas, p_values = p_values, fdrs = fdrs)
save(results, file = "DE_results_all_covars_regressed.Rdata")


## obtain colors for each cluster for downstream plotting
colors <- get_leaves_branches_col(dend1)
clusters <- ordered_tab$clusters_k_10
colors_df <- as.data.frame(cbind(clusters, colors))
colors_df <- distinct(colors_df)

# fill = colors_df$colors


## mean expression plots per cluster
values = unique(ordered_tab$clusters_k_10)
values = values[order(values)]
ordered_tab$cluster_def = 0

for(i in 1:length(values))
{
    set = which(ordered_tab$clusters_k_10 == values[i])
    ordered_tab$cluster_def[set] = i
}

clusters_as_they_appear = unique(ordered_tab$cluster_def)

plot_template = data.frame(timepoint = c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h"), time = c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48))
plot_template$mean = 0
plot_template$sd = 0
plot_template$sem = 0
plot_template$CI_low = 0
plot_template$CI_high = 0
plot_template$timepoint = factor(plot_template$timepoint, levels = unique(plot_template$timepoint))

get_summary=function(i){
    out = plot_template
    chunk = ordered_tab[which(ordered_tab$cluster_def == i), 1:20]

    title = paste0("cluster ",i,", ",nrow(chunk)," genes")
    out$mean = as.vector(apply(chunk, 2, mean))
    out$sd = as.vector(apply(chunk, 2, sd))
    out$sem = as.vector(apply(chunk, 2, sd)/sqrt(nrow(chunk)))
    out$CI_low = as.vector(apply(chunk, 2, function(x){quantile(x,0.025)}))
    out$CI_high = as.vector(apply(chunk, 2, function(x){quantile(x,0.975)}))

    col <- as.character(as.matrix(colors_df[colors_df[, "clusters"] == i,][2]))
    
    pl = ggplot(out) + geom_point(aes(x = time, y = mean)) + geom_line(aes(x = time, y = mean)) + geom_errorbar(aes(x = time, ymin = mean -sd, ymax = mean + sd), width = 0.2) + ylim(-1.2, 1.2) + geom_hline(yintercept = 0) + geom_ribbon(aes(x = time, ymin = mean - sd,ymax = mean + sd), fill = col, alpha = "0.5") + theme_classic() + xlab("infection (hrs)") + ylab("mean normalized logFC across genes") +
    theme(legend.position="top",
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    legend.title=element_blank(),
    legend.text=element_text(size=16),
    legend.key.size = unit(1, 'lines')
    ) + ggtitle(title)
    
    return(pl)
}


## plot trajectory for all clusters
plots_list <- list()

for(i in 1:length(clusters_as_they_appear))
{
    plots_list[[i]] <- get_summary(i)
}

trajectories <- grid.arrange(plots_list[[1]], plots_list[[6]], plots_list[[5]], plots_list[[3]], plots_list[[4]], plots_list[[9]], plots_list[[7]], plots_list[[2]], plots_list[[10]], plots_list[[8]], ncol = 5)

pdf(paste0("cluster_trajectories_k_10_correct.pdf"), width = 16, height = 10)
plot_grid(arbol1, trajectories, rel_heights = c(1/4, 3/4), nrow = 2)
dev.off()


## let's look at enrichments of clusters
# background is all genes in the dataset

background <- rownames(ordered_tab)
write.table(background, "gene_lists/background.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_1 <- rownames(ordered_tab[set_1, ])
write.table(genes_cluster_1, "gene_lists/genes_cluster_1.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_2 <- rownames(ordered_tab[set_2, ])
write.table(genes_cluster_2, "gene_lists/genes_cluster_2.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_3 <- rownames(ordered_tab[set_3, ])
write.table(genes_cluster_3, "gene_lists/genes_cluster_3.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_4 <- rownames(ordered_tab[set_4, ])
write.table(genes_cluster_4, "gene_lists/genes_cluster_4.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_5 <- rownames(ordered_tab[set_5, ])
write.table(genes_cluster_5, "gene_lists/genes_cluster_5.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_6 <- rownames(ordered_tab[set_6, ])
write.table(genes_cluster_6, "gene_lists/genes_cluster_6.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_7 <- rownames(ordered_tab[set_7, ])
write.table(genes_cluster_7, "gene_lists/genes_cluster_7.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_8 <- rownames(ordered_tab[set_8, ])
write.table(genes_cluster_8, "gene_lists/genes_cluster_8.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_9 <- rownames(ordered_tab[set_9, ])
write.table(genes_cluster_9, "gene_lists/genes_cluster_9.csv", row.names = FALSE, col.names = FALSE, sep = ",")

genes_cluster_10 <- rownames(ordered_tab[set_10, ])
write.table(genes_cluster_10, "gene_lists/genes_cluster_10.csv", row.names = FALSE, col.names = FALSE, sep = ",")


## heatmap of cluster 2 genes
library(gplots)
library(viridis)
clust2 <- as.matrix(ordered_tab[rownames(ordered_tab) %in% genes_cluster_2, 1:20])

colnames(clust2) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster2_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()




######################################################################################
# redo the cluster trajectory plots but do it for EUR and AFR individuals separately #
######################################################################################
eur_list <- c("european")
eur_meta_data <- dup_meta_data[dup_meta_data$ethnicity %in% eur_list,]

eur_reads <- dup_reads[, colnames(dup_reads) %in% rownames(eur_meta_data)]

dge_eur <- DGEList(counts = eur_reads)
dge_eur <- calcNormFactors(dge_eur)

design_eur <- model.matrix(~ perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + time_point_hr + infection:time_point_hr, data = eur_meta_data)

v <- voom(dge_eur, design_eur, plot = FALSE)
vfit <-lmFit(v, design_eur)
vfit <- eBayes(vfit)

betas_eur = vfit$coefficients[, 5:43]
p_values_eur = vfit$p.value[, 5:43]
fdrs_eur = p_values_eur
    
for(i in 1:ncol(fdrs_eur))
    {
        fdrs_eur[,i] = p.adjust(p_values_eur[, i], method = "bonferroni")
    }

tab = betas_eur[which(rownames(betas_eur) %in% hits), c(20:39)]
ttab = t(tab)
norm_ttab = sapply(1:ncol(ttab), function(x){ttab[,x]/max(abs(ttab[,x]))})
# normalizing logFC per gene, divide logFC by maximal ones so that max response time point is = 1 (genes that respond with same qualitative shape without considering intensity of response/no matter how strongly)
tab_norm_eur = data.frame(t(norm_ttab))
rownames(tab_norm_eur) = rownames(tab)

ordered_genes <- labels(dend1)

ordered_tab_EUR = tab_norm_eur[ordered_genes, ]
ordered_tab_EUR$clusters_k_10 <- ordered_tab$clusters_k_10
ordered_tab_EUR$ordered_genes = factor(rownames(ordered_tab_EUR), levels = rownames(ordered_tab_EUR))



# now, AFR data frame with betas
afr_list <- c("african")
afr_meta_data <- dup_meta_data[dup_meta_data$ethnicity %in% afr_list,]

afr_reads <- dup_reads[, colnames(dup_reads) %in% rownames(afr_meta_data)]

dge_afr <- DGEList(counts = afr_reads)
dge_afr <- calcNormFactors(dge_afr)

design_afr <- model.matrix(~ perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + time_point_hr + infection:time_point_hr, data = afr_meta_data)

v <- voom(dge_afr, design_afr, plot = FALSE)
vfit <-lmFit(v, design_afr)
vfit <- eBayes(vfit)

betas_afr = vfit$coefficients[, 5:43]
p_values_afr = vfit$p.value[, 5:43]
fdrs_afr = p_values_afr
    
for(i in 1:ncol(fdrs_afr))
    {
        fdrs_afr[,i] = p.adjust(p_values_afr[, i], method = "bonferroni")
    }

tab = betas_afr[which(rownames(betas_afr) %in% hits), c(20:39)]
ttab = t(tab)
norm_ttab = sapply(1:ncol(ttab), function(x){ttab[,x]/max(abs(ttab[,x]))})
# normalizing logFC per gene, divide logFC by maximal ones so that max response time point is = 1 (genes that respond with same qualitative shape without considering intensity of response/no matter how strongly)
tab_norm_afr = data.frame(t(norm_ttab))
rownames(tab_norm_afr) = rownames(tab)

ordered_tab_AFR = tab_norm_afr[ordered_genes, ]
ordered_tab_AFR$clusters_k_10 <- ordered_tab$clusters_k_10
ordered_tab_AFR$ordered_genes = factor(rownames(ordered_tab_AFR), levels = rownames(ordered_tab_AFR))

clusters_as_they_appear = unique(ordered_tab_AFR$clusters_k_10)


## ok let's plot to see if there are obvious differences
plot_template = data.frame(timepoint = c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h"), time = c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48))
plot_template$mean_EUR = 0
plot_template$mean_AFR = 0
plot_template$sd_EUR = 0
plot_template$sd_AFR = 0
plot_template$timepoint = factor(plot_template$timepoint, levels = unique(plot_template$timepoint))

EUR_dataset <- ordered_tab_EUR
AFR_dataset <- ordered_tab_AFR

get_summary_ethnicity=function(i){
    out = plot_template
    chunk_EUR = EUR_dataset[which(EUR_dataset$clusters_k_10 == i), 1:20]
    chunk_AFR = AFR_dataset[which(AFR_dataset$clusters_k_10 == i), 1:20]

    title = paste0("cluster ",i,", ",nrow(chunk_EUR)," genes")
    out$mean_EUR = as.vector(apply(chunk_EUR, 2, mean))
    out$sd_EUR = as.vector(apply(chunk_EUR, 2, sd))

    out$mean_AFR = as.vector(apply(chunk_AFR, 2, mean))
    out$sd_AFR = as.vector(apply(chunk_AFR, 2, sd))

    col <- as.character(as.matrix(colors_df[colors_df[, "clusters"] == i,][2]))
    
    pl = ggplot(out) + 
    geom_point(aes(x = time, y = mean_EUR)) + geom_line(aes(x = time, y = mean_EUR, linetype = "solid")) + geom_errorbar(aes(x = time, ymin = mean_EUR - sd_EUR, ymax = mean_EUR + sd_EUR), width = 0.2) + geom_ribbon(aes(x = time, ymin = mean_EUR - sd_EUR, ymax = mean_EUR + sd_EUR), fill = col, alpha = "0.3") + 
    ## now, AFR
	geom_point(aes(x = time, y = mean_AFR)) + geom_line(aes(x = time, y = mean_AFR, linetype = "dashed")) + geom_errorbar(aes(x = time, ymin = mean_AFR - sd_AFR, ymax = mean_AFR + sd_AFR), width = 0.2) + geom_ribbon(aes(x = time, ymin = mean_AFR - sd_AFR, ymax = mean_AFR + sd_AFR), fill = col, alpha = "0.3") + 
    ylim(-1.2, 1.2) + geom_hline(yintercept = 0) +
    scale_linetype_manual(labels = c("AFR", "EUR"), values=c("solid" = 1, "dashed" = 2)) +
    theme_classic() + xlab("infection (hrs)") + ylab("mean normalized logFC") +
    theme(legend.position="top",
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    legend.title=element_blank(),
    legend.text=element_text(size=12),
    legend.key.size = unit(1, 'lines')
    ) + ggtitle(title)
    
    return(pl)
}

## plot trajectory for all clusters
plots_list_ethnicity <- list()

for(i in 1:length(clusters_as_they_appear))
{
    plots_list_ethnicity[[i]] <- get_summary_ethnicity(i)
}

trajectories_ethnicity <- grid.arrange(plots_list_ethnicity[[6]], plots_list_ethnicity[[1]], plots_list_ethnicity[[4]], plots_list_ethnicity[[3]], plots_list_ethnicity[[7]], plots_list_ethnicity[[10]], plots_list_ethnicity[[9]], plots_list_ethnicity[[8]], plots_list_ethnicity[[2]], plots_list_ethnicity[[5]], ncol = 5)

pdf(paste0("cluster_trajectories_k_10_AFR_vs_EUR_with_SD.pdf"), width = 16, height = 10)
plot_grid(arbol1, trajectories_ethnicity, rel_heights = c(1/4, 3/4), nrow = 2)
dev.off()




############################################################
## look into cluster 2/5 more i.e. cluster within cluster ##
############################################################
setwd("/Users/Haley/Desktop/lab/code/time_course/DE/closest_time_point_no_interpolation/cluster_2_plots")
clust2 <- ordered_tab[rownames(ordered_tab) %in% genes_cluster_2, 1:20]

clust2_clusters <- hclust(dist(clust2))
dend2_clusts <- as.dendrogram(clust2_clusters)

## let's cut it at the height that leads to 10 clusters (this is arbitrary)
dend2_clusts <- color_branches(dend2_clusts, k = 3)
dend2_clusts <- color_labels(dend2_clusts, k = 3)

clust2$clusters_k_3 = data.frame(cutree(dend2_clusts, k = 3))[ ,1]

pdf("cluster_2_clusts.pdf", width = 40, height = 10)
plot(dend2_clusts)
dev.off()

ggd2 = as.ggdend(dend2_clusts)
arbol2 = ggplot(dend2_clusts, labels = FALSE)

colors_2 <- get_leaves_branches_col(dend2_clusts)
clusters <- clust2$clusters_k_3
colors_df <- as.data.frame(cbind(clusters, colors_2))
colors_df <- distinct(colors_df)

ordered_clust2_norm <- clust2

ordered_genes <- labels(dend2_clusts)

ordered_clust2_norm = ordered_clust2_norm[ordered_genes, ]
ordered_clust2_norm$ordered_genes = factor(rownames(ordered_clust2_norm), levels = rownames(ordered_clust2_norm))

## gene sets
set_1 = which(clust2_norm$clusters_k_3==1)
genes_clust_1 <- rownames(clust2_norm[set_1, ])
write.table(genes_clust_1, "genes_cluster_1.csv", row.names = FALSE, col.names = FALSE, sep = ",")

set_2 = which(clust2_norm$clusters_k_3==2)
genes_clust_2 <- rownames(clust2_norm[set_2, ])
write.table(genes_clust_2, "genes_cluster_2.csv", row.names = FALSE, col.names = FALSE, sep = ",")

set_3 = which(clust2_norm$clusters_k_3==3)
genes_clust_3 <- rownames(clust2_norm[set_3, ])
write.table(genes_clust_3, "genes_cluster_3.csv", row.names = FALSE, col.names = FALSE, sep = ",")


## perform GO enrichments on these sets vs the background gene list 
# still nothing?
clusters_as_they_appear = unique(ordered_clust2_norm$clusters_k_3)

plot_template = data.frame(timepoint = c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h"), time = c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48))
plot_template$mean = 0
plot_template$sd = 0
plot_template$sem = 0
plot_template$CI_low = 0
plot_template$CI_high = 0
plot_template$timepoint = factor(plot_template$timepoint, levels = unique(plot_template$timepoint))

cluster_df = ordered_clust2_norm

get_summary_clust=function(i){
    out = plot_template
    chunk = cluster_df[which(cluster_df$clusters_k_3 == i), 1:20]

    title = paste0("cluster ",i,", ",nrow(chunk)," genes")
    out$mean = as.vector(apply(chunk, 2, mean))
    out$sd = as.vector(apply(chunk, 2, sd))
    out$sem = as.vector(apply(chunk, 2, sd)/sqrt(nrow(chunk)))
    out$CI_low = as.vector(apply(chunk, 2, function(x){quantile(x,0.025)}))
    out$CI_high = as.vector(apply(chunk, 2, function(x){quantile(x,0.975)}))

    col <- as.character(as.matrix(colors_df[colors_df[, "clusters"] == i,][2]))
    
    pl = ggplot(out) + geom_point(aes(x = time, y = mean)) + geom_line(aes(x = time, y = mean)) + geom_errorbar(aes(x = time, ymin = mean -sd, ymax = mean + sd), width = 0.2) + ylim(-1.2, 1.2) + geom_hline(yintercept = 0) + geom_ribbon(aes(x = time, ymin = mean - sd,ymax = mean + sd), fill = col, alpha = "0.5") + theme_classic() + xlab("infection (hrs)") + ylab("mean normalized logFC across genes") +
    theme(legend.position="top",
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    legend.title=element_blank(),
    legend.text=element_text(size=16),
    legend.key.size = unit(1, 'lines')
    ) + ggtitle(title)
    
    return(pl)
}


## plot trajectory for all clusters
plots_list <- list()

for(i in 1:length(clusters_as_they_appear))
{
    plots_list[[i]] <- get_summary_clust2(i)
}

trajectories <- grid.arrange(plots_list[[1]], plots_list[[2]], plots_list[[3]], ncol = 3)

pdf(paste0("cluster2_trajectories_reclustering.pdf"), width = 16, height = 10)
plot_grid(arbol2, trajectories, rel_heights = c(1/4, 3/4), nrow = 2)
dev.off()


## look at a heatmap for cluster 2.2
genes_cluster_2.2 <- rownames(ordered_clust2_norm[set_2, ])
clust2_2 <- as.matrix(ordered_clust2_norm[rownames(ordered_clust2_norm) %in% genes_cluster_2.2, 1:20])

colnames(clust2_2) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster2.2_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_2,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.2 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


# let's look at EARLIEST (within one hour) response genes in this cluster
# genes in which normalized logFC at 1hr is already > 0.4
clust2_2_1hr_early_response <- clust2_2[clust2_2[, "1h"] > 0.60, ]
clust2_2_1hr_and_2hr <- clust2_2[clust2_2[, "2h"] > 0.75, ]

clust2_2["IL1B",]
clust2_2["IL6",]

clust2_2_1hr_early_response_list <- rownames(clust2_2_1hr_early_response)
write.table(clust2_2_1hr_early_response_list, "clust2_2_1hr_early_response_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

clust2_2_1hr_and_2hr_list <- rownames(clust2_2_1hr_and_2hr)
write.table(clust2_2_1hr_and_2hr_list, "clust2_2_1hr_and_2hr_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

pdf("cluster2.2_earliest_response_1hr.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_2_1hr_early_response,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.2 earliest response genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()

pdf("cluster2.2_high_response_2hr.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_2_1hr_and_2hr,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.2 high response 2hr, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


# do these VERY early response genes differ among AFR/EUR?
subset_ordered_tab_AFR <- ordered_tab_AFR[rownames(ordered_tab_AFR) %in% clust2_2_1hr_early_response_list, ]
subset_ordered_tab_EUR <- ordered_tab_EUR[rownames(ordered_tab_AFR) %in% clust2_2_1hr_early_response_list, ]

plot_template = data.frame(timepoint = c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h"), time = c(0,1,2,3,4,5,6,7,8,9,10))
plot_template$mean_EUR = 0
plot_template$mean_AFR = 0
plot_template$sd_EUR = 0
plot_template$sd_AFR = 0
plot_template$timepoint = factor(plot_template$timepoint, levels = unique(plot_template$timepoint))

EUR_dataset <- subset_ordered_tab_AFR
AFR_dataset <- subset_ordered_tab_EUR

get_summary_ethnicity=function(i){
    out = plot_template
    chunk_EUR = EUR_dataset[which(EUR_dataset$clusters_k_10 == i), 1:11]
    chunk_AFR = AFR_dataset[which(AFR_dataset$clusters_k_10 == i), 1:11]

    title = paste0("cluster ",i,", ",nrow(chunk_EUR)," genes")
    out$mean_EUR = as.vector(apply(chunk_EUR, 2, mean))
    out$sd_EUR = as.vector(apply(chunk_EUR, 2, sd))

    out$mean_AFR = as.vector(apply(chunk_AFR, 2, mean))
    out$sd_AFR = as.vector(apply(chunk_AFR, 2, sd))

    col <- as.character(as.matrix(colors_df[colors_df[, "clusters"] == i,][2]))
    
    pl = ggplot(out) + 
    geom_point(aes(x = time, y = mean_EUR)) + geom_line(aes(x = time, y = mean_EUR, linetype = "solid")) + geom_errorbar(aes(x = time, ymin = mean_EUR - sd_EUR, ymax = mean_EUR + sd_EUR), width = 0.2) + geom_ribbon(aes(x = time, ymin = mean_EUR - sd_EUR, ymax = mean_EUR + sd_EUR), fill = col, alpha = "0.3") + 
    ## now, AFR
	geom_point(aes(x = time, y = mean_AFR)) + geom_line(aes(x = time, y = mean_AFR, linetype = "dashed")) + geom_errorbar(aes(x = time, ymin = mean_AFR - sd_AFR, ymax = mean_AFR + sd_AFR), width = 0.2) + geom_ribbon(aes(x = time, ymin = mean_AFR - sd_AFR, ymax = mean_AFR + sd_AFR), fill = col, alpha = "0.3") + 
    ylim(-1.2, 1.2) + geom_hline(yintercept = 0) +
    scale_linetype_manual(labels = c("AFR", "EUR"), values=c("solid" = 1, "dashed" = 2)) +
    theme_classic() + xlab("infection (hrs)") + ylab("mean normalized logFC") +
    theme(legend.position="top",
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
    legend.title=element_blank(),
    legend.text=element_text(size=12),
    legend.key.size = unit(1, 'lines')
    ) + ggtitle(title)
    
    return(pl)
}


plots_list_ethnicity <- list()
plots_list_ethnicity[[1]] <- get_summary_ethnicity(2)

pdf(paste0("earliest_response_genes_AFR_vs_EUR.pdf"), width = 8, height = 8)
plot_grid(plots_list_ethnicity[[1]])
dev.off()



## look at a heatmap for cluster 2.3
genes_cluster_2.3 <- rownames(ordered_clust2_norm[set_3, ])
clust2_3 <- as.matrix(ordered_clust2_norm[rownames(ordered_clust2_norm) %in% genes_cluster_2.3, 1:20])

colnames(clust2_3) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster_2.3/cluster2.3_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_3,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.3 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()

# let's look at EARLIEST (within one hour) response genes in this cluster
# genes in which normalized logFC at 1hr is already > 0.4
clust2_3_4hr <- clust2_3[clust2_3[, "4h"] > 0.75, ]

clust2_3_4hr_list <- rownames(clust2_3_4hr)
write.table(clust2_3_4hr_list , "clust2_3_4hr_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

pdf("cluster_2.3/clust2.3_4hr_response_genes.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_3_4hr,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.3 4hr response genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


# later response
clust2_3_8hr <- clust2_3[clust2_3[, "8h"] > 0.75 & clust2_3[, "6h"] < 0.75, ]

clust2_3_8hr_list <- rownames(clust2_3_8hr)
write.table(clust2_3_8hr_list , "clust2_3_8hr_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

pdf("cluster_2.3/clust2.3_8hr_response_genes.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_3_8hr,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.3 8hr response genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


## look at a heatmap for cluster 2.1
genes_cluster_2.1 <- rownames(ordered_clust2_norm[set_1, ])
clust2_1 <- as.matrix(ordered_clust2_norm[rownames(ordered_clust2_norm) %in% genes_cluster_2.1, 1:20])

colnames(clust2_1) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster_2.1/cluster2.3_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust2_1,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2.1 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()

# let's look at EARLIEST (within one hour) response genes in this cluster
clust2_1_list <- rownames(clust2_1)
write.table(clust2_1_list , "clust2_1_list_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


## what about very early response genes across all clusters?
subset_early <- ordered_clust2_norm[, 1:20]
colnames(subset_early) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

clust2_early <- subset_early[subset_early[, "1h"] > 0.7, ]

clust2_early_list <- rownames(clust2_early)
write.table(clust2_early_list, "clust2_early_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

pdf("cluster_2_plots/cluster2_earliest_response_1hr.pdf", width = 7, height = 8, useDingbats=FALSE)
hmap <- heatmap.2(as.matrix(clust2_early),
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 2 earliest response genes, \nnormalized logFC\n logFC > 0.5 at 1hr")
	#col = inferno(30, direction = 1))
dev.off()





## NEEDS TO BE CHANGED TO REFLECT FLOW CELL -- SHOULD BE CLUSTER 4 ###



###########################
## what about cluster 5? ##
###########################
clust5 <- ordered_tab[rownames(ordered_tab) %in% genes_cluster_5, 1:20]

colnames(clust5) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster_5_plots/cluster5_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(as.matrix(clust5),
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 5 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


## let's do the subclustering for cluster 5
clust5_clusters <- hclust(dist(clust5))
dend5_clusts <- as.dendrogram(clust5_clusters)

## let's cut it at the height that leads to 10 clusters (this is arbitrary)
dend5_clusts <- color_branches(dend5_clusts, k = 3)
dend5_clusts <- color_labels(dend5_clusts, k = 3)

pdf("cluster_5_clusts.pdf", width = 40, height = 10)
plot(dend5_clusts)
dev.off()

clust5$clusters_k_3 = data.frame(cutree(dend5_clusts, k = 3))[ ,1]

ggd5 = as.ggdend(dend5_clusts)
arbol5 = ggplot(dend5_clusts, labels = FALSE)

colors_5 <- get_leaves_branches_col(dend5_clusts)
clusters <- clust5$clusters_k_3
colors_df <- as.data.frame(cbind(clusters, colors_5))
colors_df <- distinct(colors_df)

ordered_clust5_norm <- clust5

ordered_genes <- labels(dend5_clusts)

ordered_clust5_norm = ordered_clust5_norm[ordered_genes, ]
ordered_clust5_norm$ordered_genes = factor(rownames(ordered_clust5_norm), levels = rownames(ordered_clust5_norm))

## gene sets
set_1 = which(ordered_clust5_norm$clusters_k_3==1)
genes_clust_1 <- rownames(ordered_clust5_norm[set_1, ])
write.table(genes_clust_1, "cluster_5_plots/genes_cluster_1.csv", row.names = FALSE, col.names = FALSE, sep = ",")

set_2 = which(ordered_clust5_norm$clusters_k_3==2)
genes_clust_2 <- rownames(ordered_clust5_norm[set_5, ])
write.table(genes_clust_2, "cluster_5_plots/genes_cluster_2.csv", row.names = FALSE, col.names = FALSE, sep = ",")

set_3 = which(ordered_clust5_norm$clusters_k_3==3)
genes_clust_3 <- rownames(ordered_clust5_norm[set_3, ])
write.table(genes_clust_3, "cluster_5_plots/genes_cluster_3.csv", row.names = FALSE, col.names = FALSE, sep = ",")


## perform GO enrichments on these sets vs the background gene list 
clusters_as_they_appear = unique(ordered_clust5_norm$clusters_k_3)

plot_template = data.frame(timepoint = c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h"), time = c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48))
plot_template$mean = 0
plot_template$sd = 0
plot_template$sem = 0
plot_template$CI_low = 0
plot_template$CI_high = 0
plot_template$timepoint = factor(plot_template$timepoint, levels = unique(plot_template$timepoint))

cluster_df = ordered_clust5_norm

## plot trajectory for all clusters
plots_list <- list()

for(i in 1:length(clusters_as_they_appear))
{
    plots_list[[i]] <- get_summary_clust(i)
}

trajectories <- grid.arrange(plots_list[[1]], plots_list[[2]], plots_list[[3]], ncol = 3)

pdf(paste0("cluster_5_plots/cluster5_trajectories_reclustering.pdf"), width = 16, height = 10)
plot_grid(arbol5, trajectories, rel_heights = c(1/4, 3/4), nrow = 2)
dev.off()


## look at specific clusters
## cluster 5.2
genes_cluster_5.2 <- rownames(ordered_clust5_norm[set_2, ])
clust5_2 <- as.matrix(ordered_clust5_norm[rownames(ordered_clust5_norm) %in% genes_cluster_5.2, 1:20])

colnames(clust5_2) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster_5_plots/cluster5.2_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust5_2,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 5.2 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


# cluster 5.1
genes_cluster_5.1 <- rownames(ordered_clust5_norm[set_1, ])
clust5_1 <- as.matrix(ordered_clust5_norm[rownames(ordered_clust5_norm) %in% genes_cluster_5.1, 1:20])

colnames(clust5_1) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster_5_plots/cluster5.1_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust5_1,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 5.1 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()

# cluster 5.3
genes_cluster_5.3 <- rownames(ordered_clust5_norm[set_3, ])
clust5_3 <- as.matrix(ordered_clust5_norm[rownames(ordered_clust5_norm) %in% genes_cluster_5.3, 1:20])

colnames(clust5_3) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

pdf("cluster_5_plots/cluster5.3_normalized_logFC_heatmap.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(clust5_3,
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 5.3 genes, \nnormalized logFC")
	#col = inferno(30, direction = 1))
dev.off()


## what about very early response genes across all clusters?
subset_early <- ordered_clust5_norm[, 1:20]
colnames(subset_early) <- c("0h","1h","2h","3h","4h","5h","6h","7h","8h","9h","10h","12h","14h","16h","18h","24h","30h","36h","42h","48h")

clust5_early <- subset_early[subset_early[, "1h"] > 0.70, ]

clust5_early_list <- rownames(clust5_early)
write.table(clust5_early_list, "clust5_early_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

pdf("cluster_5_plots/cluster5_earliest_response_1hr.pdf", width = 6, height = 6, useDingbats=FALSE)
hmap <- heatmap.2(as.matrix(clust5_early),
	Colv = NULL,
	trace = "none",
	dendrogram = "row",
	density.info = "none", 
	col = inferno(30, direction = 1),
	main = "cluster 5 earliest response genes, \nnormalized logFC\n logFC > 0.7 at 1hr")
	#col = inferno(30, direction = 1))
dev.off()











