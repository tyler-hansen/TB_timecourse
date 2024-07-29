library(ggplot2)
library(limma)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_dup_ctls")
dup_reads <- read.table("uncorrected_raw_counts_duplicated_controls_ABS_ctl.txt", header = TRUE)
dup_meta_data <- read.table("meta_data_duplicated_contols_ABS_ctl.txt", header = T, check.names = FALSE)

length(which(colnames(dup_reads)!=rownames(dup_meta_data)))

dup_meta_data$flow_cell <- as.factor(dup_meta_data$flow_cell)
dup_meta_data$lane <- as.factor(dup_meta_data$lane)
dup_meta_data$time_point_hr <- as.factor(dup_meta_data$time_point_hr)
dup_meta_data$infection = factor(dup_meta_data$infection, levels=c("NI","Mtb_MOI_5"))

###################
## pop-DR genes ###
###################
## let´s go after the genes for which the response, at each timepoint is different between EU and AF fellows
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/dupCtl_popDR")

design <- model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + individual + time_point_hr*African_admixture + (infection*African_admixture):time_point_hr, data = dup_meta_data)

## remove the global admixture effect which is confounded with the individual labels and the MTB_t0 samples (time_point_hr0:infectionMtb_MOI_5, time_point_hr0:African_admixture:infectionMtb_MOI_5)

design = design[,-c(38,58,78)]

dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = FALSE)
fit <-lmFit(v, design)
fit <- eBayes(fit)


####################
### FROM JOAQUIN ###
####################
## flip these considering admixture because it was initially modeled in europeans but now it is modeled in africans to be congruent with the popDE genes

## CONTRAST A: time_point_hr4: difference in expression between t=4 and t=0 FOR AFRICANS
## CONTRAST B: time_point_hr4:European_admixture: additional chunk to that difference that is present on europenas (i.e. effect of european admixture on the time effect t=4 minus t=0)
## CONTRAST C: time_point_hr4:infectionMtb_MOI_5: black infection effect at t=4 FOR AFRICANS
## CONTRAST D: time_point_hr4:European_admixture:infectionMtb_MOI_5: additional chunk that needs to be added to the previous contrast to get to the black infection effect at t=4 FOR EUROPEANS (i.e. effect of europen admixture on the black infection arrow)

## 77+ is popDR directly (contrast D)
## effect of african admixture in infected samples at x time point 
# betaNI + beta_excess = three way coefficient (NOT beta_NI and beta_TB separately)
# 95 columns 

betas = fit$coefficients[, 76:94]
p_values = fit$p.value[, 76:94]
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
pdf("popDR_nested_in_time_ABS_CTLS.pdf", width = 7, height = 5.5)
ggplot(DR_numbers, aes(x = time_points, y = DR_numbers)) +
	geom_line() +
	geom_smooth() +
	geom_point() +
	geom_text_repel(aes(label = DR_numbers), direction = "y", nudge_y = 1) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		  panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) + 
	ylab("popDR genes") +
	scale_x_continuous(breaks = time_points) +
	xlab("infection time point")
dev.off()

## merge betas and fdrs so can subset on both together
names_betas <- c("t1_betas", "t2_betas", "t3_betas", "t4_betas", "t5_betas", "t6_betas", "t7_betas", "t8_betas", "t9_betas", "t10_betas", "t12_betas", "t14_betas", "t16_betas", "t18_betas", "t24_betas", "t30_betas", "t36_betas", "t42_betas", "t48_betas")
names_fdrs <- c("t1_fdrs", "t2_fdrs", "t3_fdrs", "t4_fdrs", "t5_fdrs", "t6_fdrs", "t7_fdrs", "t8_fdrs", "t9_fdrs", "t10_fdrs", "t12_fdrs", "t14_fdrs", "t16_fdrs", "t18_fdrs", "t24_fdrs", "t30_fdrs", "t36_fdrs", "t42_fdrs", "t48_fdrs")

colnames(betas) <- names_betas
colnames(fdrs) <- names_fdrs

length(which(rownames(betas)!=rownames(fdrs)))


## logFC plots to look at distribution of effect sizes and how this changes over time
##################
## UPREGULATED ###
##################
## first, subset on genes that are ONLY upregulated and significant for each time point
infection_betas <- fit$coefficients[, 58:76]
popDR_fdrs <- fdrs
length(which(rownames(infection_betas)!=rownames(popDR_fdrs)))

popDR_lists_sep_tp_upreg <- list()

for (i in 1:ncol(popDR_fdrs)){
	## subset on just upreg and significant genes
	upreg_i <- names(which(popDR_fdrs[, i] < 0.05 & infection_betas[, i] > 0))

	popDR_lists_sep_tp_upreg[[i]] <- upreg_i
	print(i)
}

upreg_significant <- lengths(popDR_lists_sep_tp_upreg)
#[1]   0  13  25   3   1   1   4   2   3   1 342  58  58 102 325 413 476 636 905

## union of these sets
popDR_union_upreg <- Reduce(union, popDR_lists_sep_tp_upreg)
# 1764

## take the betas of these genes and plot per time point 
names_betas <- c(paste0("t1, n = ", upreg_significant[1]), paste0("t2, n = ", upreg_significant[2]), paste0("t3, n = ", upreg_significant[3]), paste0("t4, n = ", upreg_significant[4]), paste0("t5, n = ", upreg_significant[5]), paste0("t6, n = ", upreg_significant[6]), paste0("t7, n = ", upreg_significant[7]), paste0("t8, n = ", upreg_significant[8]), paste0("t9, n = ", upreg_significant[9]), paste0("t10, n = ", upreg_significant[10]), paste0("t12, n = ", upreg_significant[11]), paste0("t14, n = ", upreg_significant[12]), paste0("t16, n = ", upreg_significant[13]), paste0("t18, n = ", upreg_significant[14]), paste0("t24, n = ", upreg_significant[15]), paste0("t30, n = ", upreg_significant[16]), paste0("t36, n = ", upreg_significant[17]), paste0("t42, n = ", upreg_significant[18]), paste0("t48, n = ", upreg_significant[19]), "genes")

pop_DR_betas <- as.data.frame(betas[rownames(betas) %in% popDR_union_upreg,])
pop_DR_betas$genes <- rownames(pop_DR_betas)
colnames(pop_DR_betas) <- names_betas

pop_DR_betas_melt <- melt(pop_DR_betas)

pdf("beta_dist_popDR_upreg_fdr0.05_by_time_point.pdf", width = 14, height = 7)
ggplot(pop_DR_betas_melt, aes(x = value)) +
	geom_density(alpha = 0.3) +
	geom_vline(xintercept = 0) +
	facet_wrap(facets = vars(variable), ncol = 5) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
	ggtitle("subset on upregulated, popDR fdr < 0.05") +
	xlim(-3, 3)
dev.off()


#####################
### DOWNREGULATED ###
#####################
## first, subset on genes that are ONLY downregulated and significant for each time point
popDR_lists_sep_tp_downreg <- list()

for (i in 1:ncol(popDR_fdrs)){
	## subset on just upreg and significant genes
	downreg_i <- names(which(popDR_fdrs[, i] < 0.05 & infection_betas[, i] < 0))

	popDR_lists_sep_tp_downreg[[i]] <- downreg_i
	print(i)
}

downreg_significant <- lengths(popDR_lists_sep_tp_downreg)
#[1]    0    3   10    1    1    0    3    2    2    0  511  240  305  343  588 1020  911 1080 1043

## union of these sets
popDR_union_downreg <- Reduce(union, popDR_lists_sep_tp_downreg)
# 2485

## take the betas of these genes and plot per time point 
names_betas <- c(paste0("t1, n = ", downreg_significant[1]), paste0("t2, n = ", downreg_significant[2]), paste0("t3, n = ", downreg_significant[3]), paste0("t4, n = ", downreg_significant[4]), paste0("t5, n = ", downreg_significant[5]), paste0("t6, n = ", downreg_significant[6]), paste0("t7, n = ", downreg_significant[7]), paste0("t8, n = ", downreg_significant[8]), paste0("t9, n = ", downreg_significant[9]), paste0("t10, n = ", downreg_significant[10]), paste0("t12, n = ", downreg_significant[11]), paste0("t14, n = ", downreg_significant[12]), paste0("t16, n = ", downreg_significant[13]), paste0("t18, n = ", downreg_significant[14]), paste0("t24, n = ", downreg_significant[15]), paste0("t30, n = ", downreg_significant[16]), paste0("t36, n = ", downreg_significant[17]), paste0("t42, n = ", downreg_significant[18]), paste0("t48, n = ", downreg_significant[19]), "genes")

pop_DR_betas <- as.data.frame(betas[rownames(betas) %in% popDR_union_downreg,])
pop_DR_betas$genes <- rownames(pop_DR_betas)
colnames(pop_DR_betas) <- names_betas

pop_DR_betas_melt <- melt(pop_DR_betas)

pdf("beta_dist_popDR_downreg_fdr0.05_by_time_point.pdf", width = 14, height = 7)
ggplot(pop_DR_betas_melt, aes(x = value)) +
	geom_density(alpha = 0.3) +
	geom_vline(xintercept = 0) +
	facet_wrap(facets = vars(variable), ncol = 5) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
	ggtitle("subset on downregulated, fdr < 0.05") +
	xlim(-3.5, 3.5)
dev.off()


## ok so, now let's plot the logFC plot the way it is plotted in the cell paper (combination of upregulated and downregulated popDR genes in one plot)
## check this code off of joaquin's answer
## don't need to be perfect addition from the above vectors because some could be upreg/downreg depending on the time 
## there are 22 genes that seem to fall into this category
popDR_genes_ALL <- union(popDR_union_upreg, popDR_union_downreg)
# 4034
popDR_genes_opp <- intersect(popDR_union_downreg, popDR_union_upreg)
# 215


## global infection at time point X (time_point_hrX:infectionMtb_MOI_5) + (0 or 1)(time_point_hrX:African_admixture:infectionMtb_MOI_5)
## 1 would be coded as AFR, 0 coded as EUR
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)

## pull out betas of interest for a single time point 
## use subset dfs to iterate bc this is much easier 
global_infection_betas <- fit$coefficients[, 58:76]
global_infection_effects_names <- c("t1_infection", "t2_infection", "t3_infection", "t4_infection", "t5_infection", "t6_infection", "t7_infection", "t8_infection", "t9_infection", "t10_infection", "t12_infection", "t14_infection", "t16_infection", "t18_infection", "t24_infection", "t30_infection", "t36_infection", "t42_infection", "t48_infection")
colnames(global_infection_betas) <- global_infection_effects_names

ancestry_interaction_betas <- fit$coefficients[, 77:95]
ancestry_interaction_effects_names <- c("t1_anc", "t2_anc", "t3_anc", "t4_anc", "t5_anc", "t6_anc", "t7_anc", "t8_anc", "t9_anc", "t10_anc", "t12_anc", "t14_anc", "t16_anc", "t18_anc", "t24_anc", "t30_anc", "t36_anc", "t42_anc", "t48_anc")
colnames(ancestry_interaction_betas) <- ancestry_interaction_effects_names

## subset on popDR genes
global_infection_betas_popDR <- global_infection_betas[rownames(global_infection_betas) %in% popDR_genes_ALL,]
ancestry_interaction_betas_popDR <- ancestry_interaction_betas[rownames(ancestry_interaction_betas) %in% popDR_genes_ALL,]
length(which(rownames(global_infection_betas_popDR)!=rownames(ancestry_interaction_betas_popDR)))

## initialize matrices to store
logFC_pure_EUR <- matrix(nrow = length(popDR_genes_ALL), ncol = length(iteration))
rownames(logFC_pure_EUR) <- rownames(global_infection_betas_popDR)
colnames(logFC_pure_EUR) <- names

logFC_pure_AFR <- matrix(nrow = length(popDR_genes_ALL), ncol = length(iteration))
rownames(logFC_pure_AFR) <- rownames(global_infection_betas_popDR)
colnames(logFC_pure_AFR) <- names

differences_logFC <- matrix(nrow = length(popDR_genes_ALL), ncol = length(iteration))
rownames(differences_logFC) <- rownames(global_infection_betas_popDR)
colnames(differences_logFC) <- names

for (i in 1:length(iteration)){
	global_infec_beta_i <- global_infection_betas_popDR[,i]
	ancestry_interaction_beta_i <- ancestry_interaction_betas_popDR[,i]

	logFC_EUR_i <- global_infec_beta_i + (ancestry_interaction_beta_i * 0)
	logFC_AFR_i <- global_infec_beta_i + (ancestry_interaction_beta_i * 1)

	## store these in data frames to use separately
	logFC_pure_EUR[,i] <- logFC_EUR_i
	logFC_pure_AFR[,i] <- logFC_AFR_i

	## begin formula in cell paper
	abs_logFC_EUR_i <- abs(logFC_EUR_i)
	abs_logFC_AFR_i <- abs(logFC_AFR_i)
	difference <- abs_logFC_AFR_i - abs_logFC_EUR_i

	## store
	differences_logFC[,i] <- difference
	print(i)
}

differences_logFC <- as.data.frame(differences_logFC)
differences_logFC$genes <- rownames(differences_logFC)
differences_logFC_melt <- melt(differences_logFC)

pdf("logFC_differences_popDR_genes.pdf", width = 14, height = 7)
ggplot(differences_logFC_melt, aes(x = value)) +
	geom_density(alpha = 0.3) +
	geom_vline(xintercept = 0) +
	facet_wrap(facets = vars(variable), ncol = 5) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
	ggtitle("|logFC in AA| - |logFC in EA|") +
	xlim(-3, 3)
dev.off()



## try doing this with just popDR upreg list OR popDR downreg list
## popDR_union_upreg and popDR_union_downreg
global_infection_betas_upreg <- global_infection_betas[rownames(global_infection_betas) %in% popDR_union_upreg,]
hist(global_infection_betas_upreg[,2], breaks = 100)
hist(global_infection_betas_upreg[,18], breaks = 100)

global_infection_betas_downreg <- global_infection_betas[rownames(global_infection_betas) %in% popDR_union_downreg,]
hist(global_infection_betas_downreg[,2], breaks = 100)
hist(global_infection_betas_downreg[,18], breaks = 100)

##
### try with TNF example (upreg)
	global_infec <- global_infection_betas_popDR["TNF",]
	ancestry_interaction_beta <- ancestry_interaction_betas_popDR["TNF",]

	logFC_EUR <- global_infec + (ancestry_interaction_beta * 0)
	logFC_AFR <- global_infec + (ancestry_interaction_beta * 1)

	## begin formula in cell paper
	abs_logFC_EUR <- abs(logFC_EUR)
	abs_logFC_AFR <- abs(logFC_AFR)
	difference <- abs_logFC_AFR - abs_logFC_EUR


## plot for only upregulated
global_infection_betas_popDR <- global_infection_betas[rownames(global_infection_betas) %in% popDR_union_upreg,]
ancestry_interaction_betas_popDR <- ancestry_interaction_betas[rownames(ancestry_interaction_betas) %in% popDR_union_upreg,]
length(which(rownames(global_infection_betas_popDR)!=rownames(ancestry_interaction_betas_popDR)))

## initialize matrices to store
logFC_pure_EUR <- matrix(nrow = length(popDR_union_upreg), ncol = length(iteration))
rownames(logFC_pure_EUR) <- rownames(global_infection_betas_popDR)
colnames(logFC_pure_EUR) <- names

logFC_pure_AFR <- matrix(nrow = length(popDR_union_upreg), ncol = length(iteration))
rownames(logFC_pure_AFR) <- rownames(global_infection_betas_popDR)
colnames(logFC_pure_AFR) <- names

differences_logFC <- matrix(nrow = length(popDR_union_upreg), ncol = length(iteration))
rownames(differences_logFC) <- rownames(global_infection_betas_popDR)
colnames(differences_logFC) <- names

for (i in 1:length(iteration)){
	global_infec_beta_i <- global_infection_betas_popDR[,i]
	ancestry_interaction_beta_i <- ancestry_interaction_betas_popDR[,i]

	logFC_EUR_i <- global_infec_beta_i + (ancestry_interaction_beta_i * 0)
	logFC_AFR_i <- global_infec_beta_i + (ancestry_interaction_beta_i * 1)

	## begin formula in cell paper
	abs_logFC_EUR_i <- abs(logFC_EUR_i)
	abs_logFC_AFR_i <- abs(logFC_AFR_i)
	difference <- abs_logFC_AFR_i - abs_logFC_EUR_i

	## store
	differences_logFC[,i] <- difference
	print(i)
}

differences_logFC <- as.data.frame(differences_logFC)
differences_logFC$genes <- rownames(differences_logFC)
differences_logFC_melt <- melt(differences_logFC)

ggplot(differences_logFC_melt, aes(x = value)) +
	geom_density(alpha = 0.3) +
	geom_vline(xintercept = 0) +
	facet_wrap(facets = vars(variable), ncol = 5) +
	theme_bw() +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
	ggtitle("|logFC in AA| - |logFC in EA|") +
	xlim(-3, 3)



## plot upreg gene with very negative popDR effect
infection_betas <- fit$coefficients[, 58:76]
popDR_fdrs <- fdrs
 
upreg_strong_down_popDR_list <- list()

for (i in 1:ncol(popDR_fdrs)){
	## subset on just upreg and significant genes
	upreg_strong_down_popDR <- names(which(popDR_fdrs[, i] < 0.05 & infection_betas[, i] > 0 & betas[, i] < 0))

	upreg_strong_down_popDR_list[[i]] <- upreg_strong_down_popDR
	print(i)
}

## union of these sets
upreg_strong_down_popDR <- Reduce(union, upreg_strong_down_popDR_list)

#APOE 
#CCR5 ok
betas["TNF",]
infection_betas["C3",]

## plot the Mtb samples - NI samples for corrected expression (logFC)
dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)
design = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale, data = dup_meta_data)

v <- voom(dge, design, plot = FALSE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

corrected_expression <- v$E - vfit$coefficients[,"flow_cell2"]%*%t(design[,"flow_cell2"]) - vfit$coefficients[,"flow_cell3"]%*%t(design[,"flow_cell3"]) - vfit$coefficients[,"perc_Aligned_scale"]%*%t(design[,"perc_Aligned_scale"]) - vfit$coefficients[,"perc_GC_scale"]%*%t(design[,"perc_GC_scale"]) - vfit$coefficients[,"perc_Dups_scale"]%*%t(design[,"perc_Dups_scale"]) 

plot_logFC <- function(gene){
	exp_to_plot <- cbind.data.frame(corrected_expression[gene,], dup_meta_data$ethnicity, dup_meta_data$time_point_hr, dup_meta_data$infection)
	colnames(exp_to_plot) <- c("exp","genetic_ancestry","time_point","infection")
	exp_to_plot$genetic_ancestry <- revalue(exp_to_plot$genetic_ancestry, c("african"="AF", "european"="EU"))

	exp_to_plot_NI <- subset(exp_to_plot, infection == "NI" & !(time_point == 0))
	exp_to_plot_Mtb <- subset(exp_to_plot, infection == "Mtb_MOI_5")

	order_vec <- rownames(exp_to_plot_NI)
	order_vec <- gsub("NI", "Mtb_MOI_5", order_vec)

	exp_to_plot_Mtb <- exp_to_plot_Mtb[match(order_vec, rownames(exp_to_plot_Mtb)),]
	logFC <- exp_to_plot_Mtb$exp - exp_to_plot_NI$exp

	logFC_to_plot <- cbind.data.frame(logFC, exp_to_plot_NI$genetic_ancestry, exp_to_plot_NI$time_point)
	colnames(logFC_to_plot) <- c("logFC","genetic_ancestry","time_point")

	plot <- ggplot(logFC_to_plot, aes(genetic_ancestry, logFC)) +
		geom_boxplot() +
		facet_grid(cols = vars(time_point)) +
		ylab(gene)

	return(plot)
}

gene <- "TNF"

pdf(paste0(gene,"_logFC.pdf"), width = 10, height = 8)
plot_logFC(gene)
dev.off()


###############################################################################################################################
#### now do plots above (popDR up or down dist) but instead of taking the union of all time points, plot it per time point ####
###############################################################################################################################
infection_betas <- fit$coefficients[, 58:76]
popDR_fdrs <- fdrs
length(which(rownames(infection_betas)!=rownames(popDR_fdrs)))

popDR_lists_sep_tp_upreg <- list()

for (i in 1:ncol(popDR_fdrs)){
	## subset on just upreg and significant genes
	upreg_i <- names(which(popDR_fdrs[, i] < 0.05 & infection_betas[, i] > 0))

	popDR_lists_sep_tp_upreg[[i]] <- upreg_i
	print(i)
}

upreg_significant <- lengths(popDR_lists_sep_tp_upreg)
#[1]   0  13  25   3   1   1   4   2   3   1 342  58  58 102 325 413 476 636 905
## 11 to 19 


## take the betas of these genes and plot per time point 
time_points <- c(12, 14, 16, 18, 24, 30, 36, 42, 48)
plots_list <- list()

for (i in 1:length(time_points)){

	pop_DR_betas <- as.data.frame(betas[rownames(betas) %in% popDR_lists_sep_tp_upreg[[i + 10]],])
	pop_DR_betas <- as.data.frame(pop_DR_betas[, i + 10])
	colnames(pop_DR_betas) <- "tp_betas"

	plots_list[[i]] <- ggplot(pop_DR_betas, aes(x = tp_betas)) +
		geom_density(alpha = 0.3) +
		geom_vline(xintercept = 0) +
		theme_bw() +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
		xlim(-4, 4) + 
		ggtitle(paste0(time_points[i], " hr, genes = ", upreg_significant[i + 10]))

	print(time_points[i])

}

pdf(paste0("per_time_point_upregulated_popDR_significant.pdf"), width = 10, height = 10)
grid.arrange(grobs = plots_list, ncol = 3, top = textGrob("subset on upregulated, fdr < 0.05", gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()




##################
#### DOWNREG ####
##################
infection_betas <- fit$coefficients[, 58:76]
popDR_fdrs <- fdrs
length(which(rownames(infection_betas)!=rownames(popDR_fdrs)))

popDR_lists_sep_tp_downreg <- list()

for (i in 1:ncol(popDR_fdrs)){
	## subset on just downreg and significant genes
	downreg_i <- names(which(popDR_fdrs[, i] < 0.05 & infection_betas[, i] < 0))

	popDR_lists_sep_tp_downreg[[i]] <- downreg_i
	print(i)
}

downreg_significant <- lengths(popDR_lists_sep_tp_downreg)
#[1]    0    3   10    1    1    0    3    2    2    0  511  240  305  343  588 1020  911 1080 1043
## 11 to 19 


## take the betas of these genes and plot per time point 
time_points <- c(12, 14, 16, 18, 24, 30, 36, 42, 48)
plots_list <- list()

for (i in 1:length(time_points)){

	pop_DR_betas <- as.data.frame(betas[rownames(betas) %in% popDR_lists_sep_tp_downreg[[i + 10]],])
	pop_DR_betas <- as.data.frame(pop_DR_betas[, i + 10])
	colnames(pop_DR_betas) <- "tp_betas"

	plots_list[[i]] <- ggplot(pop_DR_betas, aes(x = tp_betas)) +
		geom_density(alpha = 0.3) +
		geom_vline(xintercept = 0) +
		theme_bw() +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5)) +
		xlim(-4, 4) + 
		ggtitle(paste0(time_points[i], " hr, genes = ", downreg_significant[i + 10]))

	print(time_points[i])

}

pdf(paste0("per_time_point_downregulated_popDR_significant.pdf"), width = 10, height = 10)
grid.arrange(grobs = plots_list, ncol = 3, top = textGrob("subset on downregulated, fdr < 0.05", gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()


##################################################################################################################################
## log2FC for all the popDR genes for the 18 hours time point, upregulated n=102 (print the beta for pop-dr on top of the plot) ##
##################################################################################################################################
popDR_betas <- betas
infection_betas <- fit$coefficients[, 58:76]
popDR_fdrs <- fdrs

plots_list <- list()

time_point_18_upreg <- as.vector(unlist(popDR_lists_sep_tp_upreg[14]))

for (i in 1:length(time_point_18_upreg)){

	gene <- time_point_18_upreg[i]
	popDR_beta <- popDR_betas[gene, 14]

	exp_to_plot <- cbind.data.frame(corrected_expression[gene,], dup_meta_data$ethnicity, dup_meta_data$time_point_hr, dup_meta_data$infection)
	colnames(exp_to_plot) <- c("exp","genetic_ancestry","time_point","infection")
	exp_to_plot$genetic_ancestry <- revalue(exp_to_plot$genetic_ancestry, c("african"="AF", "european"="EU"))

	exp_to_plot_NI <- subset(exp_to_plot, infection == "NI" & !(time_point == 0))
	exp_to_plot_Mtb <- subset(exp_to_plot, infection == "Mtb_MOI_5")

	order_vec <- rownames(exp_to_plot_NI)
	order_vec <- gsub("NI", "Mtb_MOI_5", order_vec)

	exp_to_plot_Mtb <- exp_to_plot_Mtb[match(order_vec, rownames(exp_to_plot_Mtb)),]
	logFC <- exp_to_plot_Mtb$exp - exp_to_plot_NI$exp

	logFC_to_plot <- cbind.data.frame(logFC, exp_to_plot_NI$genetic_ancestry, exp_to_plot_NI$time_point)
	colnames(logFC_to_plot) <- c("logFC","genetic_ancestry","time_point")
	logFC_to_plot <- subset(logFC_to_plot, time_point == 18)

	plots_list[[i]] <- ggplot(logFC_to_plot, aes(genetic_ancestry, logFC)) +
		geom_boxplot() +
		ylab(gene) +
		ggtitle(paste0("popDR beta = ", round(popDR_beta, 3))) +
		ylim(-4,4) +
		theme(plot.title = element_text(size = 12))

	if (i%%10 == 0) {
		print(i)
	}
}

pdf(paste0("timepoint_specific_logFC_popDR_betas/18hr_up.pdf"), width = 10, height = 10)
marrangeGrob(grobs = plots_list, ncol = 4, nrow = 4, top = textGrob("subset on 18hr upregulated, fdr < 0.05", gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()










##################################################################################################################################
### log2FC all the popDR genes for down-regulated genes that overlap time points 24h and 30 hours + those that are 30h specific (basically for the union of 24 + 30 hours just exclude those that are 24 hours-only) ###
##################################################################################################################################
## 24hr is 15, 30hr is 16
## what's in 24hr that's not in 30 hrs? (24hr-specific)
## union then above genes these out
popDR_betas <- betas
infection_betas <- fit$coefficients[, 58:76]
popDR_fdrs <- fdrs

plots_list <- list()

time_point_24_downreg <- as.vector(unlist(popDR_lists_sep_tp_downreg[15]))
time_point_30_downreg <- as.vector(unlist(popDR_lists_sep_tp_downreg[16]))

specific_24hrs <- setdiff(time_point_24_downreg, time_point_30_downreg)
union <- union(time_point_24_downreg, time_point_30_downreg)

union_and_30hr_specific <- union[!union %in% specific_24hrs]
union_and_30hr_specific <- sort(union_and_30hr_specific)

for (i in 1:length(union_and_30hr_specific)){

	gene <- union_and_30hr_specific[i]
	popDR_beta <- popDR_betas[gene, c(15,16)]

	exp_to_plot <- cbind.data.frame(corrected_expression[gene,], dup_meta_data$ethnicity, dup_meta_data$time_point_hr, dup_meta_data$infection)
	colnames(exp_to_plot) <- c("exp","genetic_ancestry","time_point","infection")
	exp_to_plot$genetic_ancestry <- revalue(exp_to_plot$genetic_ancestry, c("african"="AF", "european"="EU"))

	exp_to_plot_NI <- subset(exp_to_plot, infection == "NI" & !(time_point == 0))
	exp_to_plot_Mtb <- subset(exp_to_plot, infection == "Mtb_MOI_5")

	order_vec <- rownames(exp_to_plot_NI)
	order_vec <- gsub("NI", "Mtb_MOI_5", order_vec)

	exp_to_plot_Mtb <- exp_to_plot_Mtb[match(order_vec, rownames(exp_to_plot_Mtb)),]
	logFC <- exp_to_plot_Mtb$exp - exp_to_plot_NI$exp

	logFC_to_plot <- cbind.data.frame(logFC, exp_to_plot_NI$genetic_ancestry, exp_to_plot_NI$time_point)
	colnames(logFC_to_plot) <- c("logFC","genetic_ancestry","time_point")
	logFC_to_plot <- subset(logFC_to_plot, time_point == 24 | time_point == 30)

	plots_list[[i]] <- ggplot(logFC_to_plot, aes(genetic_ancestry, logFC)) +
		geom_boxplot() +
		ylab(gene) +
		facet_grid(cols = vars(time_point)) +
		ggtitle(paste0("24hr = ", round(popDR_beta[1], 3),", 30hr = ", round(popDR_beta[2], 3))) +
		ylim(-4,4) +
		theme(plot.title = element_text(size = 9), axis.title.x=element_blank())

	if (i%%100 == 0) {
		print(i)
	}
}


pdf(paste0("timepoint_specific_logFC_popDR_betas/30hr_down.pdf"), width = 10, height = 10)
marrangeGrob(grobs = plots_list, ncol = 4, nrow = 4, top = textGrob("subset on 30hr downregulated, fdr < 0.05", gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()









