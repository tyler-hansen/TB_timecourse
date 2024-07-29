## PCA
library(ggfortify)
library(scales)
library(RColorBrewer)
library(devtools)
library(ggplot2)
library(reshape2)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")

corrected_expression <- read.table("corrected_expression_ABS_ctl.txt", header = TRUE, sep = ",")
meta_data <- read.table("meta_data_GOOD_SAMPLES.txt", header = T, sep = ",", check.names = FALSE)

length(which(colnames(corrected_expression)!=rownames(meta_data)))

meta_data$flow_cell <- as.factor(meta_data$flow_cell)
meta_data$lane <- as.factor(meta_data$lane)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
meta_data$infection = factor(meta_data$infection, levels=c("NI","Mtb_MOI_5"))
meta_data$time_infection <- as.factor(paste0(meta_data$time_point_hr, meta_data$infection))

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/")
## subset by time point
infection_status <- c("NI","Mtb_MOI_5")

for (i in 1:length(infection_status)){
	infection <- infection_status[i]

	cols <- meta_data[meta_data$infection %in% infection,]
	reads <- corrected_expression[colnames(corrected_expression) %in% cols$sample_ID_GE]

	sampleTable.order = cols
	voomed_reads.order = reads

	pca = prcomp(t(voomed_reads.order))
	loadings <- pca$rotation
	scores <- pca$x

	pdf(paste0("PCA/by_infection/",infection,".pdf"), width = 10, height = 8)
	p1 <- autoplot(pca, data = sampleTable.order, colour = 'time_point_hr') + theme_bw() + ggtitle(paste0(infection)) + guides(color = guide_legend(ncol = 2))
	p2 <- autoplot(pca, data = sampleTable.order, colour = 'ethnicity') + theme_bw() + ggtitle(paste0(infection))
	p3 <- autoplot(pca, data = sampleTable.order, colour = 'individual') + theme_bw() + ggtitle(paste0(infection))
	p4 <- autoplot(pca, data = sampleTable.order, colour = 'flow_cell') + theme_bw() + ggtitle(paste0(infection))
	grid.arrange(p1, p2, p3, p4)
	dev.off()

	print(infection)
}





