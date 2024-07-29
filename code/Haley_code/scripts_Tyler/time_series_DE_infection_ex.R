library(limma)
library(dplyr)
library(edgeR)
library(qvalue)
library(stats)
library(reshape2)
library(gridExtra)

#setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")

## read in your data
## meta data is normal format (samples are rows, columns are meta data)
meta_data = read.table("meta_data_GOOD_SAMPLES_with_ccScores.txt", header = T, sep = ",", check.names = FALSE)

## read in uncorrected gene expression matrix (counts)
reads = read.table("GE_uncorrected_raw_counts_GOOD_SAMPLES.txt", header = T, sep = ",")

#length(which(colnames(reads)!=rownames(meta_data)))

#meta_data$internal_ID <- as.factor(meta_data$internal_ID)
#eta_data$sample_number <- as.factor(meta_data$sample_number)
#eta_data$flow_cell <- as.factor(meta_data$flow_cell)
#meta_data$lane <- as.factor(meta_data$lane)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
#meta_data$infection = factor(meta_data$infection, levels=c("NI","Mtb_MOI_5"))

#################################
## begin actual analysis here ###
#################################

## read in your raw counts ad a DGE object
dge <- DGEList(counts = reads)

## calculate normalization factors
dge <- calcNormFactors(dge)

## create your design matrix
## here, flow cell, perc_aligned_scale, perc_GC_scale, and perc_Dups_scale are my technical covariates
## time point hr is the column in my design matrix with my time points in hours -- make sure that this is a FACTOR!! or you will get strange results
## infection:time_point_hr is the term that evaluates my infection effect (i.e. your effect of interest) nested in time point (so, what is the effect of the infection at each time point)
design <- model.matrix(~ flow_cell + perc_Aligned_scale + perc_GC_scale + perc_Dups_scale + individual + time_point_hr + infection:time_point_hr, data = meta_data)

## this is the standard limma pipeline
## voom will apply the normalization factors as well as fit the mean-variance trend in your data and transform it to logCPM
v <- voom(dge, design, plot = FALSE)
## fit to your model
vfit <-lmFit(v, design)
vfit <- eBayes(vfit)

## effect sizes for coefficients of interest, your columns of interest will look like a combination of the two variables in the nested term (if you're looking at time point hr 1 and the noninfected condition, the column in the vfit matrix will be time_point_hr1:infectionNI)
betas = vfit$coefficients[, 19:56]
p_values = vfit$p.value[, 19:56]
fdrs = p_values
   
for(i in 1:ncol(fdrs))
    {
        fdrs[,i] = p.adjust(p_values[, i], method = "BH")
    }


#### below is just plotting the number of DE genes at each time point ####

## let's subset on the MTB infection fdrs because these are what we care about right now
fdrs_inf = data.frame(fdrs[ , c(20:38)])
    
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






