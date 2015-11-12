## preprocess and organize GTEx lung data 
load(file="gtex-gene-counts-lung.rda")

## add subject specific data
subjtab <- read.table(file="GTEx_Data_V4_Annotations_SubjectPhenotypes_DS.txt", fill=T,
                      stringsAsFactors=FALSE, sep="\t", header=T)
subjIDs <- sapply(strsplit(colnames(lungdat),"-"),function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])
ctab <- data.frame(subjtab[map,], lungtab)

### preprocess: rlog transformation and filtering ###

## rlog transformation
library(DESeq2)
lungDES <- DESeqDataSetFromMatrix(countData = lungdat,
                                  colData = ctab,
                                  design = ~ 1)
rownames(lungDES) <- gtab$Name
rld <- rlog(lungDES)

####################################################
## save data
save(rld, gtab, file="lung-rld.rda")
####################################################

load("lung-rld.rda")
expr <- assay(rld)
dim(expr)

## filtering based on low expression and small variance
gmean <- rowMeans(expr)
summary(gmean)
gsd <- apply(expr,1,sd)
summary(gsd)

pdf("lung-rld-mean-sd.pdf")
hist(gmean, breaks=50)
hist(gsd,breaks=50)
dev.off()

## rld1: harsh filtering
i.keep1 <- which(gmean>5&gsd>1.5)
rld1 <- rld[i.keep1,]
gtab1 <- gtab[i.keep1,]

####################################################
## save data
save(rld1, gtab1, file="lung-rld-filted.rda")
####################################################

load("lung-rld.rda")
expr <- assay(rld)
dim(expr)

## filtering based on low expression and small variance
gmean <- rowMeans(expr)
summary(gmean)
gsd <- apply(expr,1,sd)
summary(gsd)

pdf("lung-rld-mean-sd.pdf")
hist(gmean, breaks=50)
hist(gsd,breaks=50)
dev.off()

## rld1: harsh filtering
i.keep2 <- which(gmean>5&gsd>1)
rld2 <- rld[i.keep2,]
gtab2 <- gtab[i.keep2,]

####################################################
## save data
save(rld2, gtab2, file="lung-rld-filted2.rda")
####################################################
