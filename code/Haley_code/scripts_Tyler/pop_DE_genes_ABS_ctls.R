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
## pop-DE genes ###
###################
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES//dupCtl_popDE")

## nested in time admixture model
## this model allows us to estimate global effects of the covariates INSTEAD of estimating them independently at each time point 
dge <- DGEList(counts = dup_reads)
dge <- calcNormFactors(dge)

## for admixture effects within condition (by time point)
## "ancestry-related transcriptional differences in non-infected and infected macrophages"
design = model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + (infection + African_admixture:infection):time_point_hr, data = dup_meta_data)

## error with Mtb time point 48 sample -- don't know why -- fix this later
design <- design[, -c(8,48)]
design <- design[, -c(45)]

## ex. infectionNI:African_admixture:time_point_hr0 and infectionMtb_MOI_5:African_admixture:time_point_hr0 
v <- voom(dge, design, plot = FALSE)
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

## results
betas = vfit$coefficients[, 45:83]
p_values = vfit$p.value[, 45:83]
fdrs = p_values
	    
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }

fdrs <- as.data.frame(fdrs)

setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/dupCtl_popDE")
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
pdf("popDE_nested_in_time_ABS_CTLS_good.pdf", width = 7, height = 3.5, useDingbats=FALSE)
ggplot(nested_pop_DE_numbers_to_plot, aes(x = time_points, y = nested_pop_DE_numbers, color = infection, label = nested_pop_DE_numbers)) +
	geom_line() +
	geom_smooth() +
	geom_point() +
	geom_text_repel(direction = "y") +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		  panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
		  legend.title = element_blank()) + 
	ylab("popDE genes") +
	scale_x_continuous(breaks = time_points_sing) +
	xlab("time point") +
	scale_color_manual(values = c("#D95F02","#2171b5"))
dev.off()


cb181d

## which genes?
## get corrected expression matrix
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_dup_ctls")
corrected_expression <- read.table("corrected_expression_duplicated_controls_ABS_ctl_log2CPM.txt", sep = ",", header = TRUE)
length(which(colnames(corrected_expression)!=rownames(dup_meta_data)))
dup_meta_data <- as.data.frame(dup_meta_data)
corrected_expression <- as.data.frame(corrected_expression)

plot_ancestry_effect_by_time_point <- function(gene){
	exp_to_plot <- cbind.data.frame(t(corrected_expression[gene,]), dup_meta_data$ethnicity, dup_meta_data$time_point_hr, dup_meta_data$infection)
	colnames(exp_to_plot) <- c("exp","genetic_ancestry","time_point","infection")
	exp_to_plot$genetic_ancestry <- revalue(exp_to_plot$genetic_ancestry, c("african"="AF", "european"="EU"))

	plot <- ggplot(exp_to_plot, aes(genetic_ancestry, exp)) +
		geom_boxplot(aes(color = genetic_ancestry)) +
		facet_grid(infection ~ time_point) +
		ylab(gene) + 
		theme_bw() +
		theme(panel.grid.minor = element_blank(),
		  panel.border = element_rect(colour = "grey80", fill = NA, size = 1.5),
		  legend.title = element_blank(),
		  axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
		scale_color_manual(values = c("#88419d", "#238b45"))

	return(plot)
}


setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES//dupCtl_popDE")
## early genes that pop up as popDE at 2hrs
## CSF2, IL6, IL1A, CCL20, JUN

## from 0 to 6 hrs -- reverses at late time point (48hrs)
pdf("CLEC7A.pdf", height = 4, width = 9)
plot_ancestry_effect_by_time_point("CLEC7A")
dev.off()



## union of set of DE genes
## for NI first
NI_df <- fdrs[, grepl("NI", names(fdrs))]
Mtb_df <- fdrs[, grepl("Mtb", names(fdrs))]

get_popDE_gene_list_UNION <- function(df_type){
	DE_gene_list <- list()
	for(i in 1:ncol(df_type)){
		DE_gene_list[[i]] <- rownames(subset(fdrs, df_type[i] < 0.05))
	}

	union_list <- Reduce(union, DE_gene_list) 

	return(union_list)
}

NI_union_popDE <- get_popDE_gene_list_UNION(NI_df)
# 163

Mtb_union_popDE <- get_popDE_gene_list_UNION(Mtb_df)
# 5578

union_total_popDE <- union(NI_union_popDE, Mtb_union_popDE)





