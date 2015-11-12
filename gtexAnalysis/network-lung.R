# setwd("C:/Users/yzhang9/Dropbox/GTEx_batch_effect_2015/gtexAnalysis")

load("lung-rld-filted.rda")
library(DESeq2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(WGCNA)
options(stringsAsFactors = FALSE) #this setting is important, do not omit!
enableWGCNAThreads() #not working under RStudio?
allowWGCNAThreads()

#############################
## data input and cleaning ##
#############################

## set data into desired format
datExpr <- as.data.frame(t(assay(rld1))) # NOTE: rows=samples, columns=genes
colnames(datExpr) <- gtab1$Description
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

sampleTree <- hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
pdf(file = "lung-sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
cut <- 35
abline(h = cut, col = "red")
dev.off()

# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr = datExpr[keepSamples, ]
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)

########################################
## choose the soft-thresholding power ##
########################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf(file = "lung-softThreshold.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

## heatmaps of each individual soft power
pdf("lung-cor-heatmaps.pdf")
for (i in 1:6){
  softPower <- i
  adjacency <- adjacency(datExpr, power = softPower)
  dissAdj <- 1-adjacency
  heatmap.2(dissAdj, trace="none", key=FALSE, cexRow=0.25, cexCol=0.25, main=paste("Adjacency with soft-power =",i))
  ## Turn adjacency matrix into topological overlap matrix (TOM)
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1-TOM
  rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
  heatmap.2(dissTOM, trace="none", key=FALSE, cexRow=0.25, cexCol=0.25, main=paste("TOM with soft-power =",i))
}
dev.off()

##########################
## network construction ##
##########################

## set soft power
softPower <- 1

##----------##
## one-step ##
##----------##

net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned", minModuleSize = 10,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "lungTOM", 
                        verbose = 3)

## number of modules
table(net$colors); sum(table(net$colors))
moduleColors <- labels2colors(net$colors)

# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf(file = "lung-net.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##--------------##
## step-by-step ##
##--------------##

adjacency <- adjacency(datExpr, power = softPower)
dissAdj <- 1-adjacency
heatmap.2(dissAdj, trace="none", key=FALSE, cexRow=0.25, cexCol=0.25, main=paste("Adjacency matrix with soft-power =", softPower))
## Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
rownames(TOM) <- colnames(TOM) <- colnames(datExpr)
dissTOM <- 1-TOM
heatmap.2(dissTOM, trace="none", key=FALSE, cexRow=0.25, cexCol=0.25, main=paste("TOM with soft-power =",softPower))

#############################
## eigengenes and heatmaps ##
#############################

## calcuate eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
rownames(MEs) <- rownames(datExpr)
MEDiss <- 1-cor(MEs)
METree <- hclust(as.dist(MEDiss), method="average")
pdf("lung-METree.pdf", width=7, height=6)
plot(METree, main="Clustering of module eigengenes", xlab="", sub="")
dev.off()

## heatmap: eigengenes vs. samples
pdf("lung-heatmap-eigengenes.pdf", width=15, height=5)
heatmap.eigengenes <- heatmap.2(t(MEs), trace="none", key=FALSE, col=blueWhiteRed(75), xlab="samples", cexCol=0.4, cexRow=0.8, main="Eigengenes vs. samples")
dev.off()

## heatmap: genes vs. samples
datExpr.centered <- sweep(datExpr, 2, colMeans(datExpr), "-")
norm <- sqrt(colSums(datExpr.centered^2))
datExpr.normed <- sweep(datExpr.centered, 2, norm, "/")
pdf("lung-heatmap-genes.pdf", width=15, height=10)
heatmap.2(t(datExpr.normed), Colv=heatmap.eigengenes$colDendrogram, Rowv=as.dendrogram(net$dendrograms[[1]]), trace="none", key=FALSE, col=blueWhiteRed(75), xlab="samples", cexCol=0.4, cexRow=0.3, main="Actual genes vs. samples")
dev.off()

#########################
## clinical trait data ##
#########################

clinical <- colData(rld1)

##-------------------------##
## clinical trait: DTHHRDY ##
##-------------------------##

## distribution of death classification
death <- as.data.frame(as.factor(clinical$DTHHRDY))
hist(death)

# Isolate DTHHRDY from the clinical traits
ventilator <- as.data.frame(as.numeric(clinical$DTHHRDY==0))
rownames(ventilator) <- rownames(clinical)
names(ventilator) <- "Ventilator"
# Add to existing module eigengenes
MET <- orderMEs(cbind(MEs, ventilator))
# Plot the relationships among the eigengenes and the trait
pdf("lung-death-ventilator.pdf", width=6, height=12)
par(cex = 0.9)
plotEigengeneNetworks(MET, "Ventilator and eigengenes", marDendro = c(0,4,2,2), marHeatmap = c(6,6,2,2), cex.lab = 0.8, xLabelsAngle = 90, printAdjacency = TRUE)
dev.off()

## MEblue has high correlation with ventilator
MEblue <- MEs$MEblue
ventilator[,1] <- as.factor(ventilator[,1])
annotation_col <- data.frame(MEblue,ventilator)
ann_colors <- list(MEblue=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
pheatmap(t(datExpr.normed[order(MEblue),mergedColors=="blue"]), scale="none", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=FALSE, annotation_col=annotation_col, annotation_colors=ann_colors, filename="lung-pheatmap-ventilator-blue.pdf", width=6, height=6)
# pheatmap(t(MEs))

# # Plot the dendrogram
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
#                       plotHeatmaps = FALSE)
# # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
# par(cex = 1.0)
# plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
#                       plotDendrograms = FALSE, xLabelsAngle = 0)

pdf("lung-death-other.pdf", width=6, height=6)
otherDeath <- c("Violent","Fast","Intermediate","Slow")
for (i in 1:4){
  death <- as.data.frame(as.numeric(clinical$DTHHRDY==i))
  rownames(death) <- rownames(clinical)
  names(death) <- otherDeath[i]
  # Add the death to existing module eigengenes
  MET <- orderMEs(cbind(MEs, death))
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  plotEigengeneNetworks(MET, "Cause of death and eigengenes", marHeatmap = c(6,6,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90, printAdjacency = TRUE)
}
dev.off()

# ## investigate MEblue and fast death relation
# fast <- as.data.frame(as.numeric(clinical$DTHHRDY==2))
# rownames(fast) <- rownames(clinical)
# names(fast) <- "Fast"
# MEblue <- MEs$MEblue
# fast[,1] <- as.factor(fast[,1])
# annotation_col <- data.frame(MEblue,fast)
# ann_colors <- list(MEblue=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
# pheatmap(t(datExpr.normed[order(MEblue),mergedColors=="blue"]), scale="none", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=FALSE, annotation_col=annotation_col, annotation_colors=ann_colors, filename="lung-pheatmap-fast-blue.pdf", width=6, height=6)

##------------------------##
## clinical trait: GENDER ##
##------------------------##

# Isolate GENDER from the clinical traits
gender <- as.data.frame(clinical$GENDER)
gender[gender==2] <- 0
rownames(gender) <- rownames(clinical)
names(gender) <- "Gender"
# Add to existing module eigengenes
MET <- orderMEs(cbind(MEs, gender))
# Plot the relationships among the eigengenes and the trait
pdf("lung-gender.pdf", width=6, height=12)
par(cex = 0.9)
plotEigengeneNetworks(MET, "Gender and eigengenes", marDendro = c(0,4,2,2), marHeatmap = c(6,6,2,2), cex.lab = 0.8, xLabelsAngle = 90, printAdjacency = TRUE)
dev.off()

## MEbrown contains mostly genes on Y chromosome
MEbrown <- MEs$MEbrown
gender[gender==1] <- "M"
gender[gender==0] <- "F"
annotation_col<- data.frame(MEbrown,gender)
ann_colors <- list(MEbrown=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
pheatmap(t(datExpr.normed[order(MEbrown),mergedColors=="brown"]), scale="none", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=FALSE, annotation_col=annotation_col, annotation_colors=ann_colors, filename="lung-pheatmap-gender-brown.pdf", width=6, height=6)

##---------------------##
## clinical trait: AGE ##
##---------------------##

# Isolate AGE from the clinical traits
age <- as.data.frame(clinical$AGE)
rownames(age) <- rownames(clinical)
names(age) <- "Age"

## distribution of age
pdf("lung-hist-age.pdf")
hist(age)
dev.off()

## a) make age as oridal variable
age.a <- rep(NA, nrow(age))
age.a[age=="20-29 years"] <- 2
age.a[age=="30-39 years"] <- 3
age.a[age=="40-49 years"] <- 4
age.a[age=="50-59 years"] <- 5
age.a[age=="60-69 years"] <- 6
age.a[age=="70-79 years"] <- 7
# Add to existing module eigengenes
MET <- orderMEs(cbind(MEs, age.a))
# Plot the relationships among the eigengenes and the trait
pdf("lung-age.pdf", width=6, height=12)
par(cex = 0.9)
plotEigengeneNetworks(MET, "Age and eigengenes", marDendro = c(0,4,2,2), marHeatmap = c(6,6,2,2), cex.lab = 0.8, xLabelsAngle = 90, printAdjacency = TRUE)
dev.off()

## b) make age as binary variable: age < 50, age >=50
age.b <- rep(NA, nrow(age))
age.b[age=="20-29 years"] <- 0
age.b[age=="30-39 years"] <- 0
age.b[age=="40-49 years"] <- 0
age.b[age=="50-59 years"] <- 1
age.b[age=="60-69 years"] <- 1
age.b[age=="70-79 years"] <- 1
# Add to existing module eigengenes
MET <- orderMEs(cbind(MEs, age.b))
# Plot the relationships among the eigengenes and the trait
pdf("lung-age-binary.pdf", width=6, height=12)
par(cex = 0.9)
plotEigengeneNetworks(MET, "Age and eigengenes", marDendro = c(0,4,2,2), marHeatmap = c(6,6,2,2), cex.lab = 0.8, xLabelsAngle = 90, printAdjacency = TRUE)
dev.off()


##################
## network plot ##
##################

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 1);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("brown");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);



############################################################################

## draft ##

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
rownames(dissTOM) <- colnames(dissTOM) <- colnames(datExpr)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^9;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
# sizeGrWindow(9,9)
pdf("lung-network-heatmap.pdf", height=9, width=9)
heatmap(plotTOM)
# TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()


nSelect = 100
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectAdj = adjacency[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
# selectTree = hclust(as.dist(selectTOM), method = "average")
# selectColors = moduleColors[select];
# Open a graphical window
# sizeGrWindow(9,9)
pdf("lung-network-heatmap-selected.pdf")
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^10;
plotAdj = selectAdj^10
diag(plotDiss) = NA;
heatmap(plotDiss, cexRow=0.4, cexCol=0.4)
# TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

#######################
## blockwise network ##
#######################

# t <- system.time(bwnet <- blockwiseModules(datExpr, power = 5,
#                        TOMType = "unsigned", minModuleSize = 30, maxBlockSize = 10000,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "lungTOM-blockwise", 
#                        verbose = 3))

###########################################
## save
save(t,bwnet,datExpr,file="lung-bwnet.rda")
###########################################

## reload
load("lung-bwnet.rda")

## number of modules
table(bwnet$colors); sum(table(bwnet$colors))
bwModuleColors = labels2colors(bwnet$colors)

# open a graphics window
# sizeGrWindow(6,6)
pdf(file = "lung-bwnet.pdf", width = 12, height = 6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
