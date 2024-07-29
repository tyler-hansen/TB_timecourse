## PCA
library(ggfortify)
library(scales)
library(RColorBrewer)
library(devtools)
library(ggplot2)
library(reshape2)
library(limma)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")

meta_data = read.table("meta_data_GOOD_SAMPLES.txt", header = T, sep = ",", check.names = FALSE)
reads = read.table("GE_uncorrected_raw_counts_GOOD_SAMPLES.txt", header = T, sep = ",")

meta_data$internal_ID <- as.factor(meta_data$internal_ID)
meta_data$sample_number <- as.factor(meta_data$sample_number)
meta_data$flow_cell <- as.factor(meta_data$flow_cell)
meta_data$lane <- as.factor(meta_data$lane)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
meta_data$infection = factor(meta_data$infection, levels=c("NI","Mtb_MOI_5"))

dge <- DGEList(counts = reads)
dge <- calcNormFactors(dge)

design <- model.matrix(~ flow_cell, data = meta_data)

v <- voom(dge, plot = TRUE)

sampleTable.order = meta_data
voomed_reads.order = v$E

length((which(rownames(sampleTable.order)!=colnames(voomed_reads.order))))
## 0

## perform PCA on uncorrected expression values
pca = prcomp(t(voomed_reads.order))
loadings <- pca$rotation
scores <- pca$x

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/")
pdf("PCA/PCA_uncorrected_PC1_PC2.pdf")
autoplot(pca, data=sampleTable.order, colour='flow_cell') 
autoplot(pca,data=sampleTable.order,colour='infection') 
autoplot(pca,data=sampleTable.order,colour='ethnicity') 
autoplot(pca,data=sampleTable.order,colour='time_point_hr') 
autoplot(pca,data=sampleTable.order,colour='individual') 
dev.off()

pdf("PCA/PCA_uncorrected_PC2_PC3.pdf")
autoplot(pca,data=sampleTable.order,colour='flow_cell', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='infection', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='ethnicity', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='time_point_hr', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='individual', x = 2, y = 3) 
dev.off()


## corrected expression values
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")
corrected_expression <- read.table("corrected_expression_log2CPM_voom.txt", header = TRUE, sep = ",")

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/")
sampleTable.order = meta_data
voomed_reads.order = corrected_expression

pca = prcomp(t(voomed_reads.order))
loadings <- pca$rotation
scores <- pca$x

pdf("PCA/PCA_corrected_rela_impo_PC1_PC2.pdf")
autoplot(pca, data=sampleTable.order, colour='flow_cell') 
autoplot(pca,data=sampleTable.order,colour='infection') 
autoplot(pca,data=sampleTable.order,colour='ethnicity') 
autoplot(pca,data=sampleTable.order,colour='time_point_hr') 
autoplot(pca,data=sampleTable.order,colour='individual') 
dev.off()

pdf("PCA/PCA_corrected_rela_impo_PC2_PC3.pdf")
autoplot(pca,data=sampleTable.order,colour='flow_cell', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='infection', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='ethnicity', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='time_point_hr', x = 2, y = 3) 
autoplot(pca,data=sampleTable.order,colour='individual', x = 2, y = 3) 
dev.off()

pdf("PCA/PCA_corrected_rela_impo_PC3_PC4.pdf")
autoplot(pca,data=sampleTable.order,colour='flow_cell', x = 3, y = 4) 
autoplot(pca,data=sampleTable.order,colour='infection', x = 3, y = 4) 
autoplot(pca,data=sampleTable.order,colour='ethnicity', x = 3, y = 4) 
autoplot(pca,data=sampleTable.order,colour='time_point_hr', x = 3, y = 4) 
autoplot(pca,data=sampleTable.order,colour='individual', x = 3, y = 4) 
dev.off()



## ggplots
######## PC1 vs PC2 #############
pdf("PCA/PC1_PC2_by_flow_cell.pdf", width = 10, height = 8)
ggplot(scores, aes(PC1, PC2)) +
		geom_point(aes(colour = sampleTable.order$flow_cell), size = 2) +
		theme_classic() + 
		scale_color_manual(values = c("white", "grey70", "black")) +
		geom_point(colour = "black", size = 2.4, shape = 1) + 
		xlab(paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained")) + 
		ylab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(colour = "flow cell")
dev.off()

pdf("PCA/PC1_PC2_by_infection_and_ethnicity.pdf", width = 8, height = 6)
ggplot(scores, aes(PC1, PC2)) +
		geom_point(aes(colour = sampleTable.order$infection, shape = sampleTable.order$ethnicity), size = 2) +
		theme_classic() + 
		#geom_point(colour = "grey40", shape = sampleTable.order$ethnicity, size = 2.4) + 
		scale_shape_manual(values=c(16, 3)) + 
		xlab(paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained")) + 
		ylab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(shape = "ethnicity", colour = "infection status")
dev.off()

pdf("PCA/PC1_PC2_color_by_ethnicity.pdf", width = 8, height = 6)
ggplot(scores, aes(PC1, PC2)) +
		geom_point(aes(colour = sampleTable.order$ethnicity, shape = sampleTable.order$infection), size = 2) +
		theme_classic() + 
		#geom_point(colour = "grey40", shape = sampleTable.order$ethnicity, size = 2.4) + 
		scale_shape_manual(values=c(16, 3)) + 
		scale_color_manual(values = c("cyan4", "chocolate3")) +
		xlab(paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained")) + 
		ylab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(shape = "infection status", colour = "ethnicity")
dev.off()

sampleTable.order$condition <- with(sampleTable.order, ifelse(sampleTable.order$time_point_hr == "0", 0, " "))
sampleTable.order$condition <- as.factor(sampleTable.order$condition)

pdf("PCA/PC1_PC2_by_infection_and_timepoint_ctl_outlined.pdf", width = 8, height = 6)
ggplot(scores, aes(PC1, PC2)) +
		geom_point(aes(colour = sampleTable.order$time_point_hr, shape = sampleTable.order$infection), size = 2) +
		#geom_point(colour = "red", shape = as.factor(sampleTable.order$infection), size = 2.4) +
		geom_text(aes(label = sampleTable.order$condition), color = "red") +
		theme_classic() + 
		#geom_point(aes(fill = sampleTable.order$condition)) +
		#geom_point(colour = "grey70", shape = sampleTable.order$infection, size = 2) + 
		scale_shape_manual(values=c(16, 3)) +
		xlab(paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained")) + 
		ylab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) +
		scale_colour_viridis_d(end = 0.97) +
		#scale_colour_manual(values = c("red","black")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(shape = "infection status", color = "time point") 
dev.off()


## plot by cline in european ancestry
pdf("PCA/PC1_PC2_color_by_EUR_ancestry_only_continuous.pdf", width = 8, height = 6)
ggplot(scores, aes(PC1, PC2)) +
		geom_point(aes(colour = sampleTable.order$only_EUR_admix, shape = sampleTable.order$infection), size = 2) +
		theme_classic() + 
		#geom_point(colour = "grey40", shape = sampleTable.order$ethnicity, size = 2.4) + 
		scale_shape_manual(values=c(16, 3)) + 
		#scale_color_manual(values = c("cyan4", "chocolate3")) +
		xlab(paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained")) + 
		ylab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) +
		scale_colour_viridis_c(end = 0.8, option = "magma", na.value = "grey70") +
		#scale_colour_gradient(low = "brown1", high = "firebrick4") +
		ggtitle("only AFR individuals colored by EUR ancestry") + 
		labs(shape = "infection status", colour = "% european admixture")
dev.off()


pdf("PCA/PC1_PC2_color_by_individual.pdf", width = 8, height = 6)
ggplot(scores, aes(PC1, PC2)) +
		geom_point(aes(colour = sampleTable.order$individual, shape = sampleTable.order$infection), size = 2) +
		theme_classic() + 
		#geom_point(colour = "grey40", shape = sampleTable.order$ethnicity, size = 2.4) + 
		scale_shape_manual(values=c(16, 3)) + 
		#scale_color_manual(values = c("cyan4", "chocolate3")) +
		xlab(paste0("PC1: ", round(summary(pca)$importance[2,1],3)*100, "% variance explained")) + 
		ylab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) +
		#scale_colour_viridis_c(end = 0.8, option = "magma", na.value = "grey70") +
		#scale_colour_gradient(low = "brown1", high = "firebrick4") +
		ggtitle("individual") + 
		labs(shape = "infection status", colour = "individual")
dev.off()


######## PC2 vs PC3 #############
pdf("PCA/PC2_PC3_by_flow_cell.pdf", width = 8, height = 6)
ggplot(scores, aes(PC2, PC3)) +
		geom_point(aes(colour = sampleTable.order$flow_cell), size = 2) +
		theme_classic() + 
		scale_color_manual(values = c("white", "grey70", "black")) +
		geom_point(colour = "black", size = 2.4, shape = 1) + 
		xlab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) + 
		ylab(paste0("PC3: ", round(summary(pca)$importance[2,3],3)*100, "% variance explained")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(colour = "flow cell")
dev.off()

pdf("PCA/PC2_PC3_by_infection_and_ethnicity.pdf", width = 8, height = 6)
ggplot(scores, aes(PC2, PC3)) +
		geom_point(aes(colour = sampleTable.order$infection, shape = sampleTable.order$ethnicity), size = 2) +
		theme_classic() + 
		#geom_point(colour = "grey40", shape = sampleTable.order$ethnicity, size = 2.4) + 
		scale_shape_manual(values=c(16, 3)) + 
		xlab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) + 
		ylab(paste0("PC3: ", round(summary(pca)$importance[2,3],3)*100, "% variance explained")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(shape = "ethnicity", colour = "infection status")
dev.off()

pdf("PCA/PC2_PC3_color_by_ethnicity.pdf", width = 8, height = 6)
ggplot(scores, aes(PC2, PC3)) +
		geom_point(aes(colour = sampleTable.order$ethnicity, shape = sampleTable.order$infection), size = 2) +
		theme_classic() + 
		#geom_point(colour = "grey40", shape = sampleTable.order$ethnicity, size = 2.4) + 
		scale_shape_manual(values=c(16, 3)) + 
		scale_color_manual(values = c("cyan4", "chocolate3")) +
		xlab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) + 
		ylab(paste0("PC3: ", round(summary(pca)$importance[2,3],3)*100, "% variance explained")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(shape = "infection status", colour = "ethnicity")
dev.off()

pdf("PCA/PC2_PC3_by_infection_and_timepoint_ctl_outlined.pdf", width = 8, height = 6)
ggplot(scores, aes(PC2, PC3)) +
		geom_point(aes(colour = sampleTable.order$time_point_hr, shape = sampleTable.order$infection), size = 2) +
		#geom_point(colour = "red", shape = as.factor(sampleTable.order$infection), size = 2.4) +
		geom_text(aes(label = sampleTable.order$condition), color = "red") +
		theme_classic() + 
		#geom_point(aes(fill = sampleTable.order$condition)) +
		#geom_point(colour = "grey70", shape = sampleTable.order$infection, size = 2) + 
		scale_shape_manual(values=c(16, 3)) +
		xlab(paste0("PC2: ", round(summary(pca)$importance[2,2],3)*100, "% variance explained")) + 
		ylab(paste0("PC3: ", round(summary(pca)$importance[2,3],3)*100, "% variance explained")) +
		scale_colour_viridis_d(end = 0.97) +
		#scale_colour_manual(values = c("red","black")) +
		#scale_colour_gradient(low = "brown1", high = "firebrick4")
		labs(shape = "infection status", color = "time point") 
dev.off()





############
### UMAP ###
############
library(umap)
scores_subset <- scores[,1:35]

set.seed(2020)
exp_umap <- umap(scores_subset)
umap_loadings <- exp_umap$layout
colnames(umap_loadings) <- c("UMAP1","UMAP2")

pdf("UMAP/UMAP_time_point_and_infection.pdf", width = 8, height = 7)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$time_point_hr, shape = meta_data$infection)) +
	scale_colour_viridis_d(end = 0.9, option = "inferno") + 
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	labs(shape = "infection status", color = "time point")
dev.off()

pdf("UMAP/UMAP_ethnicity_and_infection.pdf", width = 7, height = 6)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$ethnicity, shape = meta_data$infection)) +
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	#scale_color_manual(values = c("steelblue3", "firebrick3")) +
	scale_color_manual(values = c("#37AFA9", "#00284b")) +
	labs(shape = "infection status", colour = "ethnicity")
dev.off()


pdf("UMAP/UMAP_individual_and_ethnicity.pdf", width = 7, height = 6)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$individual, shape = meta_data$infection)) +
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	labs(shape = "infection status", colour = "individual")
dev.off()

