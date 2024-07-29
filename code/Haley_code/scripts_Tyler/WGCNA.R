library(matrixStats)
library(Hmisc)
library(splines)
library(foreach)
library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
library(survival)
library(GO.db)
library(preprocessCore)
library(impute)
library(WGCNA)
library(reshape2)
library(dplyr)
library(annotate)
library(org.Hs.eg.db)
library(anRichment)
# org.Hs.eg.db for go enrichment

options(stringsAsFactors = FALSE)
setwd("/Users/Haley/Desktop/lab/code/time_course/GOOD_SAMPLES/DATA_single_ctls")

## read in log2CPM corrected expression data as input into WGCNA
expression <- read.table("corrected_expression_log2CPM_voom.txt", header = TRUE, sep = ",")
meta_data <- read.table("meta_data_GOOD_SAMPLES.txt", header = T, sep = ",", check.names = FALSE)
length(which(colnames(expression)!=rownames(meta_data)))

## only include MTB infected samples for now 
infection <- c("Mtb_MOI_5")

meta_data <- meta_data[meta_data$infection %in% infection,]
meta_data$flow_cell <- as.factor(meta_data$flow_cell)
meta_data$time_point_hr <- as.factor(meta_data$time_point_hr)
meta_data$infection <- as.factor(meta_data$infection)
meta_data$infection = factor(meta_data$infection)

meta_data_clean <- meta_data[, c(6, 9, 11)]

## subset on raw reads
expression_Mtb <- expression[colnames(expression) %in% rownames(meta_data_clean)]
length(which(colnames(expression_Mtb)!=rownames(meta_data_clean)))

# transpose to match tutorial
expression_Mtb <- t(expression_Mtb)

## check for genes and samples with too many missing values
gsg = goodSamplesGenes(expression_Mtb, verbose = 3)
gsg$allOK
# TRUE

## next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
setwd("/Users/Haley/Desktop/lab/code/time_course/WGCNA")
sampleTree = hclust(dist(expression_Mtb), method = "average")

# plot the sample tree: open a graphic output window of size 12 by 9 inches
pdf("sample_outliers.pdf")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

# it appears there are 2 outliers (EU148_T36, EU144_T18)
# choose a height cut that will remove the sample (140) and use a branch cut at that height
# plot a line to show the cut
abline(h = 140, col = "red")
dev.off()

# determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 140, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = expression_Mtb[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# variable datExpr now contains the expression data ready for network analysis

# subset on meta data
datTraits <- meta_data_clean[rownames(meta_data_clean) %in% rownames(datExpr),]
length(which(rownames(datExpr)!=rownames(datTraits)))



##############################################################
#### step-by-step network construction & module detection ####
##############################################################
## constructing a weighted gene network entails the choice of the soft thresholding power β to which co-expression similarity is raised to calculate adjacency
# choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# scale-free topology fit index as a function of the soft-thresholding power
pdf("soft_thresholding_power.pdf", width = 10, height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
## choose the power 14, which is the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value (in this case, roughly 0.90)

## calculate co-expression similarity and adjacency 
softPower = 14
adjacency = adjacency(datExpr, power = softPower)

## topological overlap matrix (TOM)
# to minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity
# turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

## clustering using TOM
# we now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes
geneTree = hclust(as.dist(dissTOM), method = "average")
# plot the resulting clustering tree (dendrogram)
pdf("gene_clust_based_on_TOM_dissimilarity.pdf")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)
dev.off()
## each leaf, that is a short vertical line, corresponds to a gene
## branches of the dendrogram group together densely interconnected, highly co-expressed genes


## module identification amounts to the identification of individual branches ("cutting the branches off the dendrogram”)
# we like large modules, so we set the minimum module size relatively high
minModuleSize = 30
# module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
table(dynamicMods)

##  function returned 43 modules labeled 1–43 largest to smallest, label 0 is reserved for unassigned genes

## the above command lists the sizes of the modules, we now plot the module assignment under the gene dendrogram
# convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# plot the dendrogram and colors underneath
pdf("dynamic_cut_tree.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()


## merging of modules whose expression profiles are very similar
# Dynamic Tree Cut may identify modules whose expression profiles are very similar, it may be prudent to merge such modules since their genes are highly co-expressed
# to quantify co-expression similarity of entire modules, we calculate their eigengenes and cluster them on their correlation

# calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# plot the result
pdf("clustering_module_eigengenes.pdf")
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
## we choose a height cut of 0.25, corresponding to correlation of 0.75, to merge 
MEDissThres = 0.25
# plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# the merged module colors
mergedColors = merge$colors
# eigengenes of the new merged modules
mergedMEs = merge$newMEs

## number of new clusters (12)
length(table(mergedColors))

## to see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
pdf("original_vs_merged_clusters.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
# construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


################################################################################
### relating modules to external information and identifying important genes ###
################################################################################
## quantifying module–trait associations
# we would like to identify modules that are significantly associated with the measured clinical traits
# since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations

# define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table, we color code each association by the correlation value
# will display correlations and their p-values
pdf("module_trait_relationship.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()


## gene relationship to trait and important modules: Gene Significance and Module Membership ##
# we quantify associations of individual genes with our trait of interest by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait
# for each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile, this allows us to quantify the similarity of all genes on the array to every module

# define variable ancestry containing the ancestry column of datTrait
ancestry = as.data.frame(datTraits$African_admixture)
names(ancestry) = "ancestry"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, ancestry, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(ancestry), sep="")
names(GSPvalue) = paste("p.GS.", names(ancestry), sep="")

## to get signs right
# negative in GS.ancestry table
pdf("ex_negative_GS.pdf")
plot(as.vector(datExpr[,"UEVLD"]), as.numeric(unlist(ancestry)), xlab = "expression", ylab = "african admixture")
dev.off()

pdf("ex_positive_GS.pdf")
plot(as.vector(datExpr[,"EIF3F"]), as.numeric(unlist(ancestry)), xlab = "expression", ylab = "african admixture")
dev.off()


## intramodular analysis: identifying genes with high GS and MM
# using the GS and MM measures, we can identify genes that have a high significance for admixture as well as high module membership in interesting modules
# we plot a scatterplot of Gene Significance vs. Module Membership in modules of interest
module_list = c("ivory","black","brown","paleturquoise","green","saddlebrown")
color_list = c("purple","black","brown","paleturquoise","green","saddlebrown")

for(i in i:length(module_list)){
	module = module_list[i]
	color = color_list[i]
	column = match(module, modNames)
	moduleGenes = moduleColors==module

	pdf(paste0(module,"_module_membership_vs_gene_sig.pdf"))
	par(mfrow = c(1,1))
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	abs(geneTraitSignificance[moduleGenes, 1]),
	xlab = paste("Module Membership in", module, "module"),
	ylab = "Gene significance for african_admixture",
	main = paste("Module membership vs. gene significance\n"),
	cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = color)
	dev.off()

	print(module)
}

# clearly, GS and MM are highly correlated, illustrating that genes highly significantly associated with a trait are often also the most important (central) elements of modules associated with the trait


## summary output of network analysis results
# we have found modules with high association with our trait of interest, and have identified their central players by the Module Membership measure
# we now merge this statistical information with gene annotation and write out a file that summarizes the most important results 
# our expression data are annotated by hugo ids
# to look at genes in specific modules
colnames(datExpr)[moduleColors=="black"]

# gene symbol, module color, gene significance for weight, and module membership and p-values in all modules
# the modules will be ordered by their significance for weight, with the most significant ones to the left
geneInfo0 = data.frame(geneSymbol = colnames(datExpr),
						moduleColor = moduleColors,
						geneTraitSignificance,
						GSPvalue)
# order modules by their significance for ancestry
modOrder = order(-abs(cor(MEs, ancestry, use = "p")))

# add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
	oldNames = names(geneInfo0)
	geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
	MMPvalue[, modOrder[mod]]);
	names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
	paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.ancestry))
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo_african_admixture.csv", row.names = FALSE)



####################################################################################################
### interfacing network analysis with other data such as functional annotation and gene ontology ###
####################################################################################################
# our previous analysis has identified several modules (labeled brown, red, and salmon) that are highly associated with weight
# to facilitate a biological interpretation, we would like to know the gene ontologies of the genes in the modules, whether they are significantly enriched in certain functional categories etc

# export gene lists
# get all genes
allHugoIDS = colnames(datExpr)
# choose interesting modules
intModules = c("ivory", "black", "brown", "lightcyan1", "darkmagenta", "darkolivegreen", "paleturquoise", "darkgrey", "skyblue", "green", "saddlebrown", "grey")

for (module in intModules)
{
	modGenes = (moduleColors==module)
	modLLIDs = allHugoIDS[modGenes]
	write.table(as.data.frame(modLLIDs), paste0("module_genesets/module_",module,"_geneset.txt"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# as background in the enrichment analysis, we will use all probes in the analysis
write.table(as.data.frame(allHugoIDS), paste("module_genesets/background_geneset.txt"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)



### EXPERIMENTAL -- USE R package anRichment ALSO
## enrichment analysis directly within R
# WGCNA package now contains a function to perform GO enrichment analysis using a simple, single step
# org.Hs.eg.db

# convert symbols to Entrez IDs
entrez = convert2entrez(organism = "human", symbol = allHugoIDS)
# how many conversions were successful?
table(is.finite(entrez))

# the function takes a vector of module labels, and the Entrez (a.k.a. Locus Link) codes for the genes whose labels are given
GOenr = GOenrichmentAnalysis(moduleColors, entrez, organism = "human", nBestP = 10)

# this is an enrichment table containing the 10 best terms for each module present in moduleColors
tab = GOenr$bestPTerms[[4]]$enrichment

# names of the columns within the table can be accessed by
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

## abridged table version 
keepCols = c(1, 2, 5, 6, 7, 12, 13)
screenTab = tab[, keepCols]
# round the numeric columns to 2 decimal places:
numCols = c(3, 4)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)

# truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)

# shorten the column names
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name")
rownames(screenTab) = NULL




##########################################
### visualization of networks within R ###
##########################################
## one way to visualize a weighted network is to plot its heatmap
# each row and column of the heatmap correspond to a single gene
# the heatmap can depict adjacencies or topological overlaps, with light colors denoting low adjacency (overlap) and darker colors higher adjacency (overlap)
# the gene dendrograms and module colors are plotted along the top and left side of the heatmap
# this code can be executed only if the network was calculated using a single-block approach

# transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# set diagonal to NA for a nicer plot
diag(plotTOM) = NA

pdf("network_heatmap_plot.pdf")
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()
# note that the generating the heatmap plot may take a substantial amount of time



### visualizing the network of eigengenes ###
# it is often interesting to study the relationships among the found modules
# one can use the eigengenes as representative profiles and quantify module similarity by eigengene correlation
# the package contains a convenient function plotEigengeneNetworks that generates a summary plot of the eigengene network
# it is usually informative to add a clinical trait (or multiple traits) to the eigengenes to see how the traits fit into the eigengene network

# recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes

# isolate admixture 
admixture = as.data.frame(datTraits$African_admixture)
names(admixture) = "admixture"

# add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, admixture))

# plot the relationships among the eigengenes and the trait
pdf("relationship_eigengenes_admixture.pdf")
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
dev.off()

# the function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationships
# visualizing the gene network using a heatmap plot
# the heatmap depicts the Topological Overlap Matrix (TOM) among all genes in the analysis
# light color represents low overlap and progressively darker red color represents higher overlap
# blocks of darker colors along the diagonal are the modules

## the eigengene dendrogram and heatmap identify groups of correlated eigengenes termed meta-modules

##  Panel (a) shows a hierarchical clustering dendrogram of the eigengenes in which the dissimilarity of eigengenes EI , EJ is given by 1 − cor(EI , EJ ). The heatmap in panel (b) shows the eigengene adjacency AIJ = (1 + cor(EI , EJ ))/2


##################
### ANRICHMENT ###
##################
# moduleColors, entrez are the active variables 
# download multiple collection sets and run separately
## gene ontology
GOcollection = buildGOcollection(organism = "human")

## kegg, reactome
biosysCollection = BioSystemsCollection("human")

## genomic positions in 5Mb windows
genomicPosCollection = genomicPositionCollection(
	organism = "human",
	spacings = 5e6,
	overlapFactor = 2)

## MSigDB
setwd("/Users/Haley/Desktop/lab/code/time_course/WGCNA/anRichment")
MSigDBCollection = MSigDBCollection(file = paste0("msigdb_v7.0.xml"), organism = "human")

setwd("/Users/Haley/Desktop/lab/code/time_course/WGCNA")
## popDR genes (up and downreg)
# a gene set contains the Entrez identifiers of the genes that belong to it, and, for every gene, an evidence code for the evidence that the gene belongs to the gene set, as well as source of that evidence (an article reference, web site, etc)
popDR_up <- as.data.frame(read.table("/Users/Haley/Desktop/lab/code/time_course/DE/ABSOLUTE_closest_time_point/popDR/genelists/upreg_UNION_all_timepoints.txt", sep = ","))
popDR_up <- as.character(popDR_up$V1)
popDR_down <- as.data.frame(read.table("/Users/Haley/Desktop/lab/code/time_course/DE/ABSOLUTE_closest_time_point/popDR/genelists/downreg_UNION_all_timepoints.txt", sep = ","))
popDR_down <- as.character(popDR_down$V1)

entrez_popDR_up = unique(convert2entrez(organism = "human", symbol = popDR_up))
table(is.finite(entrez_popDR_up))
entrez_popDR_down = unique(convert2entrez(organism = "human", symbol = popDR_down))
table(is.finite(entrez_popDR_down))

## upregulated first
popDR_up_geneset = newGeneSet(
	geneEntrez = entrez_popDR_up,
	geneEvidence = "IEP",
	geneSource = paste0("MLS_time_course_2019"),
	ID = "dummyid",
	name = "time_course_popDR_up",
	description = "popDR genes from MLS time course",
	source = paste0("xx"),
	organism = "human",
	internalClassification = c("popDR_up", "dummy"),
	groups = "popDR_up",
	lastModified = "2020-01-08")

popDR_upGroup = newGroup(name = "popDR_up", description = "popDR up genes",
source = "time course data")

popDR_up_Collection = newCollection(dataSets = list(popDR_up_geneset), groups = list(popDR_upGroup))

## downregulated next
popDR_down_geneset = newGeneSet(
	geneEntrez = entrez_popDR_down,
	geneEvidence = "IEP",
	geneSource = paste0("MLS_time_course_2019"),
	ID = "dummyid",
	name = "time_course_popDR_down",
	description = "popDR genes from MLS time course",
	source = paste0("xx"),
	organism = "human",
	internalClassification = c("popDR_down", "dummy"),
	groups = "popDR_down",
	lastModified = "2020-01-08")

popDR_downGroup = newGroup(name = "popDR_down", description = "popDR down genes",
source = "time course data")

popDR_down_Collection = newCollection(dataSets = list(popDR_down_geneset), groups = list(popDR_downGroup))

### PERFORM ENRICHMENTS ###
## make list of all the collections of interest 
# GOcollection, 
# biosysCollection, 
# genomicPosCollection, 
# MSigDBCollection, 
# popDR_up_Collection, 
# popDR_down_Collection

# go
enrichment = enrichmentAnalysis(
	classLabels = moduleColors, identifiers = entrez,
	refCollection = GOcollection,
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")
write.csv(enrichment$enrichmentTable, file = paste0("anRichment/GOcollection_enrichmentTable.csv"), row.names = FALSE)

# biosys
enrichment = enrichmentAnalysis(
	classLabels = moduleColors, identifiers = entrez,
	refCollection = biosysCollection,
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")
write.csv(enrichment$enrichmentTable, file = paste0("anRichment/biosysCollection_enrichmentTable.csv"), row.names = FALSE)

# genomic pos
enrichment = enrichmentAnalysis(
	classLabels = moduleColors, identifiers = entrez,
	refCollection = genomicPosCollection,
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")
write.csv(enrichment$enrichmentTable, file = paste0("anRichment/genomicPosCollection_enrichmentTable.csv"), row.names = FALSE)

# msigdb
enrichment = enrichmentAnalysis(
	classLabels = moduleColors, identifiers = entrez,
	refCollection = MSigDBCollection,
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")
write.csv(enrichment$enrichmentTable, file = paste0("anRichment/MSigDBCollection_enrichmentTable.csv"), row.names = FALSE)

# popDR up
enrichment = enrichmentAnalysis(
	classLabels = moduleColors, identifiers = entrez,
	refCollection = popDR_up_Collection,
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")
write.csv(enrichment$enrichmentTable, file = paste0("anRichment/popDR_up_Collection_enrichmentTable.csv"), row.names = FALSE)

# popDR down
enrichment = enrichmentAnalysis(
	classLabels = moduleColors, identifiers = entrez,
	refCollection = popDR_down_Collection,
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")
write.csv(enrichment$enrichmentTable, file = paste0("anRichment/popDR_down_Collection_enrichmentTable.csv"), row.names = FALSE)

## this one separately
internalColl = internalCollection(organism = "human")

bloodAtlasEnrichment = enrichmentAnalysis(classLabels = moduleColors, identifiers = entrez,
	refCollection = internalColl,
	useGroups = c("BloodAtlases", "ImmunePathways"),
	useBackground = "given",
	threshold = 5e-2,
	thresholdType = "Bonferroni",
	getOverlapEntrez = TRUE,
	getOverlapSymbols = TRUE,
	geneSeparator = ",",
	ignoreLabels = "grey")

write.csv(bloodAtlasEnrichment$enrichmentTable, file = paste0("anRichment/bloodAtlas_enrichmentTable.csv"), row.names = FALSE)




#######################################
### permute african admixture label ###
#######################################
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

## permute african admixture
datTraits$permut_African_admixture <- sample(datTraits$African_admixture)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# since we have a moderately large number of modules and traits, a suitable graphical representation will help in reading the table, we color code each association by the correlation value
# will display correlations and their p-values
pdf("module_trait_rel_permuted_admix.pdf")
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()
















