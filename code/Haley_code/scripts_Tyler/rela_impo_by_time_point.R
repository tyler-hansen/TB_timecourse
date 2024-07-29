library(limma)
library(edgeR)
library(relaimpo)
library(robustbase)
library(reshape2)
library(ggplot2)
library(dplyr)

#setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_dup_ctls")
raw_counts <- read.table("uncorrected_raw_counts_duplicated_controls_ABS_ctl.txt", header = TRUE)
meta_data <- read.table("meta_data_duplicated_contols_ABS_ctl.txt", header = T, check.names = FALSE)

length(which(colnames(raw_counts)!=rownames(meta_data)))

cols <- meta_data

cols$flow_cell <- as.factor(cols$flow_cell)
cols$time_point_hr <- as.factor(cols$time_point_hr)
cols$infection <- as.factor(cols$infection)
cols$infection = factor(cols$infection, levels=c("NI","Mtb_MOI_5"))

## subset on raw reads
reads <- raw_counts[colnames(raw_counts) %in% rownames(cols)]

## initialize dfs
time_points <- c("1","2","3","4","5","6","7","8","9","10","12","14","16","18","24","30","36","42","48")
rela_impo_outs_ALL_time_points <- data.frame(matrix(nrow = 5, ncol = length(time_points)))
NORMALIZED_rela_impo_outs_ALL_time_points <- data.frame(matrix(nrow = 5, ncol = length(time_points)))
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/relaimpo_time_point")

for(j in 1:length(time_points)){

    ## subset both meta data and reads per time point
    ## make sure time point is a FACTOR

    tp <- time_points[j]
    cols_j <- cols[cols$time_point_hr %in% tp,]

    write.table(cols_j, paste0("tp_subsets_meta_data/",tp,"_meta_data.txt"), quote = FALSE)
    cols_j <- read.table(paste0("tp_subsets_meta_data/",tp,"_meta_data.txt"), header = TRUE)

    cols_j$time_point_hr <- as.factor(cols_j$time_point_hr)
    cols_j$infection <- as.factor(cols_j$infection)
    cols_j$flow_cell <- as.factor(cols_j$flow_cell)
    cols_j$infection = factor(cols_j$infection, levels=c("NI","Mtb_MOI_5"))

    ## subset on reads
    reads_j <- reads[colnames(reads) %in% rownames(cols_j)]

    ## rest of limma/relaimpo pipeline
    dge <- DGEList(counts = reads_j)
    dge <- calcNormFactors(dge)

    design = model.matrix(~ perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + infection + African_admixture, data = cols_j)
    v <- voom(dge, design, plot = FALSE)

    exp = v$E
    weights = v$weights

    gene_model = lm(exp[1,] ~ perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + infection + African_admixture, data = cols_j, weights = weights[1,])
    rel_impo = suppressWarnings(calc.relimp(gene_model))

    ## matrix to store
    importances = data.frame(matrix(nrow = nrow(reads_j), ncol = length(rel_impo$lmg)))
    rownames(importances) <- rownames(exp)
    colnames(importances) <- names(rel_impo$lmg)


    ## loop over all genes
    for(i in 1:nrow(exp))
    {
        if(i%%1000 == 0) print(i)
        gene_model = lm(exp[i,] ~ perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + infection + African_admixture, data = cols_j, weights = weights[i,])
        rel_impo = suppressWarnings(calc.relimp(gene_model, type = "lmg"))
        importances[i,] = rel_impo$lmg
    }

    write.table(importances, paste0("tp_subsets_importances/importances_",tp,".txt"))

    ## calculate col medians (median of all genes for each covariate -- collapses matrix)
    medians <- as.data.frame(colMedians(as.matrix(importances)))
    colnames(medians) <- c(paste0("medians_",tp))

    rownames(rela_impo_outs_ALL_time_points) <- rownames(medians)
    colnames(rela_impo_outs_ALL_time_points)[j] <- colnames(medians)[1]
    rela_impo_outs_ALL_time_points[,j] <- medians[,1]

    ## get the total R_squared of the model in an additional column
    importances$R_squared = apply(importances, 1, sum)
    ## normalize things dividing importances by the values of R_squared
    relative_importances = importances

    for(i in 1:ncol(relative_importances))
    {
        relative_importances[, i] = relative_importances[, i]/relative_importances$R_squared
    }

    ## remove the last column because it is all ones
    relative_importances = relative_importances[1:(ncol(relative_importances)-1)]

    medians_NORM <- as.data.frame(colMedians(as.matrix(relative_importances)))
    colnames(medians_NORM) <- c(paste0("medians_",tp))

    rownames(NORMALIZED_rela_impo_outs_ALL_time_points) <- rownames(medians_NORM)
    colnames(NORMALIZED_rela_impo_outs_ALL_time_points)[j] <- colnames(medians_NORM)[1]
    NORMALIZED_rela_impo_outs_ALL_time_points[,j] <- medians_NORM[,1]

    print(tp)
}

write.table(rela_impo_outs_ALL_time_points, "rela_impo_outs_ALL_time_points.txt", sep = ",", quote = FALSE)
write.table(NORMALIZED_rela_impo_outs_ALL_time_points, "NORMALIZED_rela_impo_outs_ALL_time_points.txt", sep = ",", quote = FALSE)

#######################################
### infection vs ancestry PVE plots ###
#######################################
rela_impo_outs_ALL_time_points <- read.table("rela_impo_outs_ALL_time_points.txt", header = TRUE, sep = ",")
time_points <- c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48)

rela_impo_subset <- rela_impo_outs_ALL_time_points[c(4,5),]
rela_impo_subset_t <- as.data.frame(t(rela_impo_subset))
rela_impo_subset_t$time_point <- time_points
rela_impo_subset_to_plot <- melt(rela_impo_subset_t, id.vars = "time_point")

pdf("infection_vs_ancestry_PVE_importances.pdf", width = 8, height = 6)
ggplot(rela_impo_subset_to_plot) +
    geom_point(aes(x = time_point, y = value, color = variable)) + 
    geom_line(aes(x = time_point, y = value, color = variable)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          legend.title = element_blank()) +
    ylab("percent variance explained") + 
    xlab("hrs post infection") +
    scale_x_continuous(breaks = time_points) 
dev.off()


## with normalized data
NORMALIZED_rela_impo_outs_ALL_time_points <- read.table("NORMALIZED_rela_impo_outs_ALL_time_points.txt", header = TRUE, sep = ",")
time_points <- c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,24,30,36,42,48)

NORM_rela_impo_subset <- NORMALIZED_rela_impo_outs_ALL_time_points[c(4,5),]
NORM_rela_impo_subset_t <- as.data.frame(t(NORM_rela_impo_subset))
NORM_rela_impo_subset_t$time_point <- time_points
NORM_rela_impo_subset_to_plot <- melt(NORM_rela_impo_subset_t, id.vars = "time_point")

pdf("NORMALIZED_infection_vs_ancestry_PVE_importances.pdf", width = 8, height = 6)
ggplot(NORM_rela_impo_subset_to_plot) +
    geom_point(aes(x = time_point, y = value, color = variable)) + 
    geom_line(aes(x = time_point, y = value, color = variable)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
          legend.title = element_blank()) +
    ylab("percent variance explained") + 
    xlab("hrs post infection") +
    scale_x_continuous(breaks = time_points) 
dev.off()









