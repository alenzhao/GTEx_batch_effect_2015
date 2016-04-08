library(bnlearn)
library(minet)
library(WGCNA)
options(stringsAsFactors = FALSE) #this setting is important, do not omit!
enableWGCNAThreads() #not working under RStudio?
allowWGCNAThreads()
## set soft power
softPower <- 1

#########################
## Gaussian simulation ##
#########################
library(mnormt)

## settings and parameters
nGenes <- 20
nSamples <- 100

nClusters <- 4
clusterSize <- rep(5,nClusters)
clusterCor <- (6:3)/10

bgCor <- 0.05 #background correlation
bin.thred <- 0.2 #threshold for binary weights

## correlation matrix
sigma0 <- matrix(bgCor, ncol=nGenes, nrow=nGenes)
diag(sigma0) <- 1
for (i in 1:nClusters){
  temp <- matrix(clusterCor[i], ncol=clusterSize[i], nrow=clusterSize[i])
  diag(temp) <- 1
  index0 <- sum(clusterSize[1:i])-clusterSize[i]+1
  index1 <- sum(clusterSize[1:i])
  index <- index0:index1
  sigma0[index,index] <- temp
}

## true network
true.net <- sigma0*0
true.net[sigma0!=bgCor] <- 1
diag(true.net) <- 0

## multivariate normal
set.seed(999)
datExpr <- rmnorm(n=nSamples, varcov=sigma0) # NOTE: rows=samples, columns=genes
datExpr <- as.data.frame(datExpr)

###########
## WGCNA ##
###########

## one-step network construction
net.wgcna <- blockwiseModules(datExpr, power = softPower,
                              TOMType = "unsigned", minModuleSize = 1,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = FALSE,
                              verbose = 3)

## plot the dendrogram and the module colors underneath
moduleColors <- labels2colors(net.wgcna$colors)
pdf(file = "wgcna-net-sim0.pdf", width = 12, height = 6)
plotDendroAndColors(net.wgcna$dendrograms[[1]], moduleColors[net.wgcna$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = NULL, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

## estimated weights on edges
adjMat <- adjacency(datExpr, power = softPower) #i.e. correlation if softPower=1
est.wgcna <- adjMat
diag(est.wgcna) <- 0 #no edge between itself
## binary weights
est.wgcna.bin <- est.wgcna>bin.thred

## MAE = mean absolute error
##     = sum(abs(estimated edge weights - true edge weights (i.e. sigma0)))/N, where N=nGenes^2
mae.wgcna <- sum(abs(est.wgcna-sigma0))/(nGenes^2)
mae.wgcna
## binary
mae.wgcna.bin <- sum(abs(est.wgcna.bin-true.net))/(nGenes^2)
mae.wgcna.bin

############
## ARACNE ##
############

## fit MI network, which gives estimated weights on edges
net.aracne <- minet(datExpr, method = "aracne")
est.aracne <- net.aracne
## binary weights
est.aracne.bin <- est.aracne>bin.thred

## MAE
mae.aracne <- sum(abs(est.aracne-sigma0))/(nGenes^2)
mae.aracne
## binary
mae.aracne.bin <- sum(abs(est.aracne.bin-true.net))/(nGenes^2)
mae.aracne.bin

################################################
## Bayesian network - hill climbing algorithm ##
################################################

## Hill climbing greedy search
bn.hc <- hc(datExpr)

## find linear coeficients from bn.hc
bn.hc.fit <- bn.fit(bn.hc, datExpr)
coefs <- coef(bn.hc.fit)

## edge
child <- vector()
parent <- vector()
weight <- vector()
for (i in 1:length(coefs)){
  z <- coefs[[i]][-1]
  child <- c(child, rep(names(coefs)[i],length(z)))
  parent <- c(parent, names(z))
  weight <- c(weight, z)
}
child <- gsub("V","C",child)
parent <- gsub("V","P",parent)
edge <- data.frame(child, parent, weight)

## estimated weights on edges
est.bn0 <- sigma0*0
rownames(est.bn0) <- paste0("C",1:nGenes) #row=child
colnames(est.bn0) <- paste0("P",1:nGenes) #col=parent
for (i in 1:nrow(edge)){
  est.bn0[child[i],parent[i]] <- weight[i] 
}
est.bn <- est.bn0 + t(est.bn0)
## binary weights
est.bn.bin <- est.bn>bin.thred

## MAE
mae.bn <- sum(abs(est.bn-sigma0))/(nGenes^2)
mae.bn
## binary
mae.bn.bin <- sum(abs(est.bn.bin-true.net))/(nGenes^2)
mae.bn.bin

################################################################################################

## repeat for Niter times
Niter <- 1000
mae.wgcna <- mae.wgcna.bin <- mae.aracne <- mae.aracne.bin <- mae.bn <- mae.bn.bin <- vector()
est.wgcna.all <- est.aracne.all <- est.bn.all <- true.net.all <- NULL

set.seed(999)
for (iter in 1:Niter){
  ## multivariate normal
  datExpr <- rmnorm(n=nSamples, varcov=sigma0) # NOTE: rows=samples, columns=genes
  datExpr <- as.data.frame(datExpr)
  
  ## WGCNA: estimated weights on edges
  adjMat <- adjacency(datExpr, power = softPower) #i.e. correlation if softPower=1
  est.wgcna <- adjMat
  diag(est.wgcna) <- 0 #no edge between itself
  est.wgcna.all <- cbind(est.wgcna.all,est.wgcna)
  ## binary weights
  est.wgcna.bin <- est.wgcna>bin.thred
  ## MAE
  mae.wgcna[iter] <- sum(abs(est.wgcna-sigma0))/(nGenes^2)
  ## binary
  mae.wgcna.bin[iter] <- sum(abs(est.wgcna.bin-true.net))/(nGenes^2)
  
  ## ARACNE: estimated weights on edges
  net.aracne <- minet(datExpr, method = "aracne")
  est.aracne <- net.aracne
  est.aracne.all <- cbind(est.aracne.all,est.aracne)
  ## binary weights
  est.aracne.bin <- est.aracne>bin.thred
  ## MEA
  mae.aracne[iter] <- sum(abs(est.aracne-sigma0))/(nGenes^2)
  ## binary
  mae.aracne.bin[iter] <- sum(abs(est.aracne.bin-true.net))/(nGenes^2)
  
  ## BAYESIAN: Hill climbing greedy search
  bn.hc <- hc(datExpr)
  ## find linear coeficients from bn.hc
  bn.hc.fit <- bn.fit(bn.hc, datExpr)
  coefs <- coef(bn.hc.fit)
  ## edge
  child <- vector()
  parent <- vector()
  weight <- vector()
  for (i in 1:length(coefs)){
    z <- coefs[[i]][-1]
    child <- c(child, rep(names(coefs)[i],length(z)))
    parent <- c(parent, names(z))
    weight <- c(weight, z)
  }
  child <- gsub("V","C",child)
  parent <- gsub("V","P",parent)
  edge <- data.frame(child, parent, weight)
  ## estimated weights on edges
  est.bn0 <- sigma0*0
  rownames(est.bn0) <- paste0("C",1:nGenes) #row=child
  colnames(est.bn0) <- paste0("P",1:nGenes) #col=parent
  for (i in 1:nrow(edge)){
    est.bn0[child[i],parent[i]] <- weight[i] 
  }
  est.bn <- est.bn0 + t(est.bn0)
  est.bn.all <- cbind(est.bn.all,est.bn)
  ## binary weights
  est.bn.bin <- est.bn>bin.thred
  ## MAE
  mae.bn[iter] <- sum(abs(est.bn-sigma0))/(nGenes^2)
  ## binary
  mae.bn.bin[iter] <- sum(abs(est.bn.bin-true.net))/(nGenes^2)
  
  ## true.net record
  true.net.all <- cbind(true.net.all,true.net)
  
  if((iter/100)%%1==0){cat(iter,"\n")}
}

mean(mae.wgcna);mean(mae.wgcna.bin)
mean(mae.aracne);mean(mae.aracne.bin)
mean(mae.bn);mean(mae.bn.bin)

###############
## ROC curve ##
###############
library(ROCR)

## WGCNA
predictions <- as.vector(est.wgcna.all)
labels <- as.vector(true.net.all)
pred <- prediction(predictions, labels)
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,"auc")
slot(auc,"y.values")

pdf("roc_wgcna_sim0.pdf", width=6, height=6)
plot(perf, lwd=2, main="ROC for WGCNA", colorize=TRUE, print.cutoffs.at=c(0.2))
legend("bottom", paste("AUC =", round(as.numeric(slot(auc,"y.values")),4)), bty="n")
dev.off()

## ARECNE
predictions <- as.vector(est.aracne.all)
unique.p <- length(unique(predictions))
labels <- as.vector(true.net.all)
pred <- prediction(predictions, labels)
perf <- performance(pred,"tpr","fpr")
xval <- slot(perf,"x.values")[[1]][unique.p-1]
yval <- slot(perf,"y.values")[[1]][unique.p-1]
auc <- performance(pred,"auc")
slot(auc,"y.values")

pdf("roc_aracne_sim0.pdf", width=6, height=6)
plot(perf, lwd=2, main="ROC for ARACNE", colorize=TRUE, print.cutoffs.at=c(0.2), xlim=c(0,xval), ylim=c(0,yval))
legend("bottom", paste("AUC =", round(as.numeric(slot(auc,"y.values")),4)), bty="n")
dev.off()

## BAYESIAN
predictions <- as.vector(abs(est.bn.all))
unique.p <- length(unique(predictions))
labels <- as.vector(true.net.all)
pred <- prediction(predictions, labels)
perf <- performance(pred,"tpr","fpr")
xval <- slot(perf,"x.values")[[1]][unique.p-1]
yval <- slot(perf,"y.values")[[1]][unique.p-1]
auc <- performance(pred,"auc")
slot(auc,"y.values")

pdf("roc_bn_sim0.pdf", width=6, height=6)
plot(perf, lwd=2, main="ROC for BN", colorize=TRUE, print.cutoffs.at=c(0.2), xlim=c(0,xval), ylim=c(0,yval))
legend("bottom", paste("AUC =", round(as.numeric(slot(auc,"y.values")),4)), bty="n")
dev.off()
