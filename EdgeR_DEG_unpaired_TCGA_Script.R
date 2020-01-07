##############################################################
#######        RNASeq gene expression analysis       #########        
##############        COMPLETE DATASET          ##############
##############################################################

dir <- #WORKING_DIRECTORY
setwd(dir)

library(edgeR)

rawdata <- read.delim("TCGA_RNASeq_counts.txt", check.names=FALSE, stringsAsFactors=FALSE)

# gene symbol as rownames
rownames(rawdata)<-rawdata$SYMBOL

# Manually load the 'Tissue' vector T for tumor sample N for non-tumor sample
Tissue <- c("T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","N","T","N","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","N","T","N","T","T","T","N","T","N","T","N","T","T","N","T","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","T","N","T","T","N","T","N","T","N","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","N","T","T","T","N","T","N","T","T","T","N","T","T","N","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","N","T","T","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","N","T","N","T","N","T","N","T","T","N","T","N","T","T","T","T","N","T","T","T","N","T","T","T","T","T","T","T","N","T","N","T","N","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","N","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T")

# Obtaining the size of each gene to normalize according to the size of the genes
gene.length<-read.table("his-Size.tab",header=T)
idx<-match(rawdata$SYMBOL,gene.length$Gene)
results_counts<-gene.length[idx,]
results_counts[is.na(results_counts$Length),"Length"]<-0
nrow(results_counts)

# DGEList object construction using 'Tissue' as group
  yET <- DGEList(counts=rawdata[,2:550], genes=results_counts, group=Tissue)

# low expression genes deletion
nrow(yET$counts)
keep<-rowSums(cpm(yET)>1) >= 52 #number of non-tumoral samples
yET<-yET[keep, ]
yET$samples$lib.size <-colSums(yET$counts)
nrow(yET$counts)

# Data normalization
dgenorm<-calcNormFactors(yET)

#RPKM
RPKM_general <- rpkm(yET)

# Common dispersion estimation 
dgenorm <- estimateCommonDisp(dgenorm,verbose = T)

# Tagwise dispersion estimation
dgenorm <- estimateTagwiseDisp(dgenorm)

# DE analysis with 'exactTest' 
et <- exactTest(dgenorm,pair=c("N","T"))

# keep results in 'TopTags' table (Benjamini–Hochberg adjust method equal to FDR)
top<-topTags(et, n=nrow(et), adjust.method="BH",sort.by="PValue")


head(top$table)

# under-expressed and over-expressed genes summary
summary(decideTests(et))


##############################################################################################

