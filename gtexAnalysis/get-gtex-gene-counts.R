## read in GTEx gene read count data and perform preliminary analysis
dat <- read.table("GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz",
                  skip=2, sep="\t",header=TRUE,stringsAsFactors=FALSE)

## extract gene annotation
gtab <- dat[,1:2]

## format count data
colnames(dat) <- gsub(".","-",colnames(dat),fixed=T)
dat <- dat[,-c(1,2)]
dat <- as.matrix(dat)

## load sample annotation
stab <- read.table("GTEx_Data_V4_Annotations_SampleAttributesDS.txt",
                   header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
                   quote="")
ind <- which(stab$SAMPID%in%colnames(dat))
stab <- stab[ind,]

## check that count data match sample annotation
identical(stab$SAMPID,colnames(dat))

## save objects
save(gtab, stab, dat, file="gtex-gene-counts.rda")

## save just lung data
ind <- which(stab$SMTSD=="Lung")
lungdat <- dat[,ind]
lungtab <- stab[ind,]

save(gtab, lungdat, lungtab, file="gtex-gene-counts-lung.rda")
