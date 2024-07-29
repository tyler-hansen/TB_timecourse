library(ggfortify)
library(scales)
library(RColorBrewer)
library(devtools)
library(ggplot2)
library(reshape2)
library(limma)
library(Seurat)
library(umap)
library(tidyr)
library(gridExtra)
library(plyr)

## try UMAP first, on just cell cycle genes
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")
corrected_expression <- read.table("corrected_expression_log2CPM_voom.txt", header = TRUE, sep = ",")
meta_data <- read.table("meta_data_GOOD_SAMPLES.txt", header = T, sep = ",", check.names = FALSE)

meta_data$flow_cell <- as.factor(meta_data$flow_cell)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
meta_data$infection = factor(meta_data$infection, levels=c("NI","Mtb_MOI_5"))


## subset on cell cycle genes
# 97 total
# 94 in our data
S_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
cell_cycle_genes <- c(S_genes, g2m_genes)

cell_cycle_expression <- corrected_expression[which(rownames(corrected_expression) %in% cell_cycle_genes), ]

## run UMAP
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/cell_cycle_tests")
sampleTable.order = meta_data
voomed_reads.order = cell_cycle_expression
length((which(rownames(sampleTable.order)!=colnames(voomed_reads.order))))

pca = prcomp(t(voomed_reads.order))
loadings <- pca$rotation
scores <- as.data.frame(pca$x)

set.seed(2020)
exp_umap <- umap(scores)
umap_loadings <- as.data.frame(exp_umap$layout)
colnames(umap_loadings) <- c("UMAP1","UMAP2")

pdf("UMAP_time_infection_ccGenes_allCovars.pdf", width = 8, height = 7)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$time_point_hr, shape = meta_data$infection)) +
	scale_colour_viridis_d(end = 0.9, option = "inferno") + 
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	labs(shape = "infection status", color = "time point")
dev.off()

pdf("UMAP_ancestry_infection_ccGenes_allCovars.pdf", width = 7, height = 6)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$ethnicity, shape = meta_data$infection)) +
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	#scale_color_manual(values = c("steelblue3", "firebrick3")) +
	scale_color_manual(values = c("#37AFA9", "#00284b")) +
	labs(shape = "infection status", colour = "ethnicity")
dev.off()

pdf("PC1_PC2_infection_ancestry.pdf", width = 8, height = 6)
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

pdf("PC1_PC2_color_ancestry.pdf", width = 8, height = 6)
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


## plot cell cycle scores over time for both afr/eur in mtb and non-infected cells
cdt_list <- c("G2M_NI", "G2M_Mtb", "S_NI", "S_Mtb")
label_list <- c("G2M", "G2M", "S", "S")
genes_list <- list(g2m_genes, g2m_genes, S_genes, S_genes)
infection_list <- c("NI", "Mtb_MOI_5", "NI", "Mtb_MOI_5")
plots_list <- list()
plots_list_by_indiv <- list()

for(i in 1:length(cdt_list)){
	cdt <- cdt_list[i]
	cycle_genes <- genes_list[[i]]
	infection <- infection_list[i]
	label <- label_list[i]

	samples <- meta_data[(meta_data$infection %in% infection), ]
	expression <- corrected_expression[, colnames(corrected_expression) %in% samples$sample_ID_GE]
	cycle_expression <- expression[which(rownames(expression) %in% cycle_genes), ]

	## scale the expression per gene (by row)
	exp_mean_center <- t(scale(t(cycle_expression)))

	## do an average across time points across for each individual
	## convert to long format of data so perform averages across time points
	exp_mean_center_long <- melt(exp_mean_center)
	colnames(exp_mean_center_long) <- c("gene","individual","scaled_expression")

	exp_mean_center_long <- separate(exp_mean_center_long, individual, c("indiv_ID", "individual", "time_point", "infection"))
	exp_mean_center_long <- separate(exp_mean_center_long, individual, c("ancestry", "ID"), sep = "(?<=[A-Za-z])(?=[0-9])")

	exp_mean_center_long_avg <- exp_mean_center_long %>%
	                            group_by(indiv_ID, time_point) %>%
	                            summarise_at(vars(scaled_expression), funs(mean(.))) 

	#test <- subset(exp_mean_center_long, indiv_ID == "X10" & time_point == "T1")
	#test_means <- mean(test$scaled_expression)
	##store this as a score
	if(i == 1){
		G2M_NI_scores <- as.data.frame(exp_mean_center_long_avg)
	} else if(i == 2){
		G2M_Mtb_scores <- as.data.frame(exp_mean_center_long_avg)
	} else if(i == 3){
		S_NI_scores <- as.data.frame(exp_mean_center_long_avg)
	} else {
		S_Mtb_scores <- as.data.frame(exp_mean_center_long_avg)
	}


	exp_mean_center_avg <- as.data.frame(spread(exp_mean_center_long_avg, time_point, scaled_expression))

	## readd ancestry call
	exp_mean_center_long$individual <- paste0(exp_mean_center_long$ancestry, exp_mean_center_long$ID)
	ancestry <- unique(data.frame(exp_mean_center_long$indiv_ID, exp_mean_center_long$ancestry, exp_mean_center_long$individual))
	colnames(ancestry) <- c("indiv_ID", "ancestry", "individual")
	exp_mean_with_ancestry <- merge(exp_mean_center_avg, ancestry, by = "indiv_ID")

	rownames(exp_mean_with_ancestry) <- exp_mean_with_ancestry$indiv_ID
	exp_mean_with_ancestry <- exp_mean_with_ancestry[-c(1)]

	if(i == 1 || i == 3){
			exp_mean_with_ancestry <- exp_mean_with_ancestry[c("ancestry","individual","T0","T4","T8","T12","T18","T24","T30","T36","T48")]
		}else{
			exp_mean_with_ancestry <- exp_mean_with_ancestry[c("ancestry","individual","T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T12","T14","T16","T18","T24","T30","T36","T42","T48")]
		}

	exp_mean_to_plot <- melt(exp_mean_with_ancestry)

	plots_list[[i]] <- ggplot(exp_mean_to_plot, aes(x = ancestry, y = value)) + 
							geom_boxplot(aes(color = ancestry)) +
							facet_grid(. ~ variable) + 
							theme_bw() +
							theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5), legend.position = "none", plot.title = element_text(size = 10)) +
							ggtitle(cdt) + 
							ylab(paste0(label," score")) + 
							xlab("")

	plots_list_by_indiv[[i]] <- ggplot(exp_mean_to_plot, aes(x = ancestry, y = value)) + 
							facet_grid(. ~ variable) + 
							geom_point(aes(color = ifelse(individual == "EU144", "red", ifelse(individual == "AF193", "red", "black")))) +
							scale_color_identity() +
							theme_bw() +
							theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5), plot.title = element_text(size = 10)) +
							ggtitle(cdt) + 
							ylab(paste0(label," score")) + 
							xlab("")

}

pdf(paste0("cc_avg_scores.pdf"), width = 15, height = 6) 
grid.arrange(grobs = plots_list, nrow = 2, widths = c(1,2))
dev.off()

pdf(paste0("cc_avg_scores_by_individual.pdf"), width = 15, height = 6) 
grid.arrange(grobs = plots_list_by_indiv, nrow = 2, widths = c(1,2))
dev.off()


## if we regress these scores out, what happens?
# rearranage so everything can be merged into meta data
G2M_NI_scores$timeID <- paste0(G2M_NI_scores$time_point, G2M_NI_scores$indiv_ID)
G2M_NI_scores$infection <- "NI"
G2M_NI_scores$timeIDinf <- paste0(G2M_NI_scores$timeID, G2M_NI_scores$infection)
G2M_NI_scores <- select(G2M_NI_scores, timeIDinf, scaled_expression)

G2M_Mtb_scores$timeID <- paste0(G2M_Mtb_scores$time_point, G2M_Mtb_scores$indiv_ID)
G2M_Mtb_scores$infection <- "Mtb_MOI_5"
G2M_Mtb_scores$timeIDinf <- paste0(G2M_Mtb_scores$timeID, G2M_Mtb_scores$infection)
G2M_Mtb_scores <- select(G2M_Mtb_scores, timeIDinf, scaled_expression)

S_NI_scores$timeID <- paste0(S_NI_scores$time_point, S_NI_scores$indiv_ID)
S_NI_scores$infection <- "NI"
S_NI_scores$timeIDinf <- paste0(S_NI_scores$timeID, S_NI_scores$infection)
S_NI_scores <- select(S_NI_scores, timeIDinf, scaled_expression)

S_Mtb_scores$timeID <- paste0(S_Mtb_scores$time_point, S_Mtb_scores$indiv_ID)
S_Mtb_scores$infection <- "Mtb_MOI_5"
S_Mtb_scores$timeIDinf <- paste0(S_Mtb_scores$timeID, S_Mtb_scores$infection)
S_Mtb_scores <- select(S_Mtb_scores, timeIDinf, scaled_expression)

G2M_scores <- rbind(G2M_NI_scores, G2M_Mtb_scores)
colnames(G2M_scores) <- c("timeIDinf","G2M_scores")
S_scores <- rbind(S_NI_scores, S_Mtb_scores)
colnames(S_scores) <- c("timeIDinf","S_scores")

## rearrange meta data
meta_data <- separate(meta_data, sample_ID_GE, c("indiv_ID", "individual", "time_point", "infection_label"))
meta_data$timeID <- paste0(meta_data$time_point, meta_data$indiv_ID)
meta_data$timeIDinf <- paste0(meta_data$timeID, meta_data$infection)

## merge scores with the meta data
meta_data$sample_ID_GE <- paste0("X",meta_data$sample_ID)
meta_data <- join(meta_data, G2M_scores, by = "timeIDinf")
meta_data <- join(meta_data, S_scores, by = "timeIDinf")
rownames(meta_data) <- meta_data$sample_ID_GE
sample_ID_ccScores <- data.frame(meta_data$sample_ID_GE, meta_data$G2M_scores, meta_data$S_scores)
colnames(sample_ID_ccScores) <- c("sample_ID_GE", "G2M_scres", "S_scores")

pdf("G2M_S_score_correlation.pdf")
plot(meta_data$G2M_scores, meta_data$S_scores)
dev.off()

write.table(sample_ID_ccScores, "ccScores.txt", sep = ",", quote = FALSE, row.names = FALSE)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")
write.table(meta_data, "meta_data_GOOD_SAMPLES_with_ccScores.txt", sep = ",", quote = FALSE)



###################
## pop-DE genes ###
###################
# need to add sample_ID_G2Mscores to add into DUPLICATED meta data
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_dup_ctls")
dup_reads <- read.table("uncorrected_raw_counts_duplicated_controls_ABS_ctl.txt", header = TRUE)
dup_meta_data <- read.table("meta_data_duplicated_contols_ABS_ctl_with_ccScores.txt", header = T, check.names = FALSE)

length(which(colnames(dup_reads)!=rownames(dup_meta_data)))

dup_meta_data$flow_cell <- as.factor(dup_meta_data$flow_cell)
dup_meta_data$lane <- as.factor(dup_meta_data$lane)
dup_meta_data$time_point_hr <- as.factor(dup_meta_data$time_point_hr)
dup_meta_data$infection = factor(dup_meta_data$infection, levels=c("NI","Mtb_MOI_5"))
dup_meta_data$sample_ID_GE <- rownames(dup_meta_data)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/cell_cycle_tests/popDE_regressed_out_G2M")
## nested in time admixture model
## this model allows us to estimate global effects of the covariates INSTEAD of estimating them independently at each time point 
dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)

## for admixture effects within condition (by time point)
## "ancestry-related transcriptional differences in non-infected and infected macrophages"
design = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + G2M_scores + S_scores + (infection + African_admixture:infection):time_point_hr, data = dup_meta_data)

length(which(colnames(dup_reads)!=rownames(design)))
## error with Mtb time point 48 sample -- don't know why -- fix this later
design <- design[, -c(10,50)]
design <- design[, -c(47)]

## ex. infectionNI:African_admixture:time_point_hr0 and infectionMtb_MOI_5:African_admixture:time_point_hr0 
v <- voom(dge, design, plot = FALSE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

## results
betas = vfit$coefficients[, 47:85]
p_values = vfit$p.value[, 47:85]
fdrs = p_values
	    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }

#fdrs <- as.data.frame(fdrs)

## how many below threshold?
names <- c("NI_t0", "NI_t1", "Mtb_t1", "NI_t2", "Mtb_t2", "NI_t3", "Mtb_t3", "NI_t4", "Mtb_t4", "NI_t5", "Mtb_t5", "NI_t6", "Mtb_t6", "NI_t7", "Mtb_t7", "NI_t8", "Mtb_t8", "NI_t9", "Mtb_t9", "NI_t10", "Mtb_t10", "NI_t12", "Mtb_t12", "NI_t14", "Mtb_t14", "NI_t16", "Mtb_t16", "NI_t18", "Mtb_t18", "NI_t24", "Mtb_t24", "NI_t30", "Mtb_t30", "NI_t36", "Mtb_t36", "NI_t42", "Mtb_t42", "NI_t48", "Mtb_t48")
time_points <- c(0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,12,12,14,14,16,16,18,18,24,24,30,30,36,36,42,42,48,48)
time_points_sing <- c(0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48)

nested_pop_DE_numbers <- apply(fdrs, 2, function(x) sum(x < 0.05))
nested_pop_DE_numbers <- cbind.data.frame(time_points, nested_pop_DE_numbers)
nested_pop_DE_numbers$infection <- rownames(nested_pop_DE_numbers)

nested_pop_DE_numbers_to_plot <- separate(nested_pop_DE_numbers, "infection", "infection_status", sep = ":")
rownames(nested_pop_DE_numbers_to_plot) <- names
nested_pop_DE_numbers_to_plot$infection <- revalue(nested_pop_DE_numbers_to_plot$infection, c("infectionNI"="NI", "infectionMtb_MOI_5"="Mtb"))


## plot
pdf("popDE_nested_in_time_ABS_CTLS_ccRegressed.pdf", width = 9, height = 5.5)
ggplot(nested_pop_DE_numbers_to_plot, aes(x = time_points, y = nested_pop_DE_numbers, color = infection, label = nested_pop_DE_numbers)) +
	geom_line() +
	#geom_smooth() +
	geom_point() +
	geom_text_repel(direction = "y") +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
		  legend.title = element_blank()) + 
	ylab("popDE genes") +
	scale_x_continuous(breaks = time_points_sing) +
	xlab("infection time point")
dev.off()

popDE_list <- list()
for (i in 1:ncol(fdrs)){
	## subset on just upreg and significant genes
	list <- names(which(fdrs[, i] < 0.05))

	popDE_list[[i]] <- list
	print(i)
}
names(popDE_list) <- colnames(fdrs)




################################################################################################################################
## what's up with 30 and 36 hrs for the nested model? model this separately (do NOT nest time point, with and without scores) ##
################################################################################################################################
## nested in time admixture model
## this model allows us to estimate global effects of the covariates INSTEAD of estimating them independently at each time point
#### 30 HOURS ####
meta_data_30hr <- dup_meta_data[dup_meta_data$time_point_hr == 30, ]
meta_data_30hr$time_point_hr <- factor(meta_data_30hr$time_point_hr)
meta_data_30hr$flow_cell <- factor(meta_data_30hr$flow_cell)

reads_30hr <- reads[colnames(reads) %in% meta_data_30hr$sample_ID_GE]
length(which(colnames(reads_36hr)!=rownames(meta_data_30hr)))

dge <- DGEList(counts = reads_30hr)
dge <- calcNormFactors(dge)

design_tp = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + G2M_scores + S_scores + (infection + African_admixture:infection), data = meta_data_30hr)

v <- voom(dge, design_tp, plot = FALSE)
vfit <-lmFit(v, design_tp)
vfit <- eBayes(vfit)

## results
betas = vfit$coefficients[, 9:10]
p_values = vfit$p.value[, 9:10]
fdrs = p_values
	    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }

pop_DE_numbers_30hrs <- apply(fdrs, 2, function(x) sum(x < 0.05))
## infectionNI:African_admixture infectionMtb_MOI_5:African_admixture 
#                           2702                                    2

### WITHOUT CC SCORES
design_tp = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + (infection + African_admixture:infection), data = meta_data_30hr)

v <- voom(dge, design_tp, plot = FALSE)
vfit <-lmFit(v, design_tp)
vfit <- eBayes(vfit)

## results
betas = vfit$coefficients[, 7:8]
p_values = vfit$p.value[, 7:8]
fdrs = p_values
	    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }

pop_DE_numbers_30hrs <- apply(fdrs, 2, function(x) sum(x < 0.05))
## infectionNI:African_admixture infectionMtb_MOI_5:African_admixture 
#                              5                                    1



#### 36 HOURS ####
meta_data_36hr <- dup_meta_data[dup_meta_data$time_point_hr == 36, ]
meta_data_36hr$time_point_hr <- factor(meta_data_36hr$time_point_hr)
meta_data_36hr$flow_cell <- factor(meta_data_36hr$flow_cell)

reads_36hr <- reads[colnames(reads) %in% meta_data_36hr$sample_ID_GE]
length(which(colnames(reads_36hr)!=rownames(meta_data_36hr)))

dge <- DGEList(counts = reads_36hr)
dge <- calcNormFactors(dge)

design_tp = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + G2M_scores + S_scores + (infection + African_admixture:infection), data = meta_data_36hr)

v <- voom(dge, design_tp, plot = FALSE)
vfit <-lmFit(v, design_tp)
vfit <- eBayes(vfit)

## results
betas_36 = vfit$coefficients[, 9:10]
p_values_36 = vfit$p.value[, 9:10]
fdrs_36 = p_values
	    
for(i in 1:ncol(fdrs_36))
    {
        fdrs_36[,i] = p.adjust(p_values[, i], method = "BH")
    }

pop_DE_numbers_36hrs <- apply(fdrs_36, 2, function(x) sum(x < 0.05))
## infectionNI:African_admixture infectionMtb_MOI_5:African_admixture 
#                             44                                    0



## correlate nested betas for 36hr with non-nested-in-time betas for 36 hrs
pdf("betas_36hrNested_vs_36hrSingle.pdf")
plot(betas[,34], betas_36[,1], xlab = "nested model 36hr betas", ylab = "single tp model 36hr betas")
dev.off()




### WITHOUT CC SCORES
design_tp = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + (infection + African_admixture:infection), data = meta_data_36hr)

v <- voom(dge, design_tp, plot = FALSE)
vfit <-lmFit(v, design_tp)
vfit <- eBayes(vfit)

## results
betas = vfit$coefficients[, 7:8]
p_values = vfit$p.value[, 7:8]
fdrs = p_values
	    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }

pop_DE_numbers_36hrs <- apply(fdrs, 2, function(x) sum(x < 0.05))
## infectionNI:African_admixture infectionMtb_MOI_5:African_admixture 
#                              5                                    1



###################
## pop-DR genes ###
###################
dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)

design <- model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + G2M_scores + S_scores + individual + time_point_hr*African_admixture + (infection*African_admixture):time_point_hr, data = dup_meta_data)

## remove the global admixture effect which is confounded with the individual labels and the MTB_t0 samples (time_point_hr0:infectionMtb_MOI_5, time_point_hr0:African_admixture:infectionMtb_MOI_5)

design = design[,-c(40,60,80)]
v <- voom(dge, design, plot = FALSE)
fit <-lmFit(v, design)
fit <- eBayes(fit)

betas = fit$coefficients[, 78:96]
p_values = fit$p.value[, 78:96]
fdrs = p_values
	    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }

## how many below threshold?
names <- c("Mtb_t1", "Mtb_t2", "Mtb_t3", "Mtb_t4", "Mtb_t5", "Mtb_t6", "Mtb_t7", "Mtb_t8", "Mtb_t9", "Mtb_t10", "Mtb_t12", "Mtb_t14", "Mtb_t16", "Mtb_t18", "Mtb_t24", "Mtb_t30", "Mtb_t36", "Mtb_t42", "Mtb_t48")
time_points <- c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48)

DR_numbers <- apply(fdrs, 2, function(x) sum(x < 0.05))
DR_numbers <- cbind.data.frame(time_points, DR_numbers)
DR_numbers$infection <- rownames(DR_numbers)

rownames(DR_numbers) <- names

## plot
pdf("popDR_nested_in_time_ABS_CTLS_ccRegressed.pdf", width = 7, height = 5.5)
ggplot(DR_numbers, aes(x = time_points, y = DR_numbers)) +
	geom_line() +
	#geom_smooth() +
	geom_point() +
	geom_text_repel(aes(label = DR_numbers), direction = "y", nudge_y = 1) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		  panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) + 
	ylab("popDR genes") +
	scale_x_continuous(breaks = time_points) +
	xlab("infection time point")
dev.off()


popDR_list <- list()
for (i in 1:ncol(fdrs)){
	## subset on just upreg and significant genes
	list <- names(which(fdrs[, i] < 0.05))

	popDR_list[[i]] <- list
	print(i)
}
names(popDR_list) <- colnames(fdrs)



########################################
## redo umap with these regressed out ##
########################################
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")
uncorrected_expression <- read.table("GE_uncorrected_raw_counts_GOOD_SAMPLES.txt", header = TRUE, sep = ",")

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/cell_cycle_tests/")
# meta_data

length((which(rownames(meta_data)!=colnames(uncorrected_expression))))

dge <- DGEList(counts = uncorrected_expression)
dge <- calcNormFactors(dge)
design = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + G2M_scores + S_scores, data = meta_data)

v <- voom(dge, design, plot = TRUE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

expression_ALL_corrected <- v$E - vfit$coefficients[,"flow_cell2"]%*%t(design[,"flow_cell2"]) - vfit$coefficients[,"flow_cell3"]%*%t(design[,"flow_cell3"]) - vfit$coefficients[,"perc_Aligned_scale"]%*%t(design[,"perc_Aligned_scale"]) - vfit$coefficients[,"perc_GC_scale"]%*%t(design[,"perc_GC_scale"]) - vfit$coefficients[,"perc_Dups_scale"]%*%t(design[,"perc_Dups_scale"]) - vfit$coefficients[,"G2M_scores"]%*%t(design[,"G2M_scores"]) - vfit$coefficients[,"S_scores"]%*%t(design[,"S_scores"]) 

cell_cycle_CORRECTED_expression <- expression_ALL_corrected[which(rownames(expression_ALL_corrected) %in% cell_cycle_genes), ]

sampleTable.order = meta_data
voomed_reads.order = cell_cycle_CORRECTED_expression
length((which(rownames(sampleTable.order)!=colnames(voomed_reads.order))))

pca = prcomp(t(voomed_reads.order))
loadings <- pca$rotation
scores <- as.data.frame(pca$x)

set.seed(2020)
exp_umap <- umap(scores)
umap_loadings <- as.data.frame(exp_umap$layout)
colnames(umap_loadings) <- c("UMAP1","UMAP2")

pdf("UMAP_time_infection_ccGenes_ccRegressed.pdf", width = 8, height = 7)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$time_point_hr, shape = meta_data$infection)) +
	scale_colour_viridis_d(end = 0.9, option = "inferno") + 
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	labs(shape = "infection status", color = "time point")
dev.off()

pdf("UMAP_ancestry_infection_ccGenes_ccRegressed.pdf", width = 7, height = 6)
ggplot(umap_loadings, aes(UMAP1, UMAP2)) +
	geom_point(aes(color = meta_data$ethnicity, shape = meta_data$infection)) +
	scale_shape_manual(values=c(3, 16)) +
	theme_classic() + 
	#scale_color_manual(values = c("steelblue3", "firebrick3")) +
	scale_color_manual(values = c("#37AFA9", "#00284b")) +
	labs(shape = "infection status", colour = "ethnicity")
dev.off()

pdf("PC1_PC2_color_ancestry_ccRegressed.pdf", width = 8, height = 6)
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



## plot correlations of CC genes pairwise in sets to see if anti-correlation is true or not (or is the score sign being flipped somehow)
S_expression <- corrected_expression[which(rownames(corrected_expression) %in% S_genes), ]
G2M_expression <- corrected_expression[which(rownames(corrected_expression) %in% g2m_genes), ]
length((which(colnames(S_expression)!=colnames(G2M_expression))))

corr <- c()
correlations_all <- c()

for(i in 1:nrow(S_expression)){
	gene_i_S_gene <- t(as.vector(S_expression[i, ]))

	for(j in 1:nrow(G2M_expression)){
		gene_j_G2M_gene <- t(as.vector(G2M_expression[j,]))

		# calculate correlation
		corr <- cor(gene_i_S_gene, gene_j_G2M_gene)
		if(j == 1){
			corr_gene_i <- corr
		}else{
			corr_gene_i <- rbind(corr_gene_i, corr)
		}
	}

	if(i == 1){
		correlations_all <- corr_gene_i
	}else{
		correlations_all <- rbind(correlations_all, corr_gene_i)
		# 52 x 42 dimensions = 2184 rows
	}

}

colnames(correlations_all) <- "corrs"
correlations_all <- as.data.frame(correlations_all)

#S_test <- t(as.vector(S_expression[2,]))
#G2M_test <- t(as.vector(G2M_expression[3,]))
#cor(S_test, G2M_test)

pdf("distribution_pairwise_correlations_G2M_v_S_genes.pdf", width = 5, height = 5)
ggplot(correlations_all, aes(x = corrs)) +
	geom_density() +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
	geom_vline(xintercept = 0, color = c("grey50"), linetype = "dashed") +
	xlab("expression correlations")
dev.off()
















