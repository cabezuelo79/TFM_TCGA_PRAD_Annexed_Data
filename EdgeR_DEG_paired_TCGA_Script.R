##############################################################
#######        RNASeq gene expression analysis       #########        
##############        PAIRED   DATASET          ##############
##############################################################

dir <- #WORKING_DIRECTORY
setwd(dir)

library(edgeR)

rawdataPAR <- read.delim("TCGA_RNASeq_PAIRED_counts.txt", check.names=FALSE, stringsAsFactors=FALSE)

# gene symbol as rownames
rownames(rawdataPAR)<-rawdataPAR$SYMBOL

# Obtaining the size of each gene to normalize according to the size of the genes
gene.length_PAR<-read.table("his-Size.tab",header=T)
idx_PAR<-match(rawdataPAR$SYMBOL,gene.length_PAR$Gene)
results_counts_PAR<-gene.length_PAR[idx_PAR,]
results_counts_PAR[is.na(results_counts_PAR$Length),"Length"]<-0
nrow(results_counts_PAR)

# DGEList object construction
yPAR <- DGEList(counts=rawdataPAR[,2:105], genes=results_counts_PAR)

# low expression genes deletion
nrow(yPAR$counts)
keep_par<-rowSums(cpm(yPAR)>1) >= 52 #number of NT samples
yPAR<-yPAR[keep_par, ]
yPAR$samples$lib.size <-colSums(yPAR$counts)
nrow(yPAR$counts)

# Data normalization
yPAR <- calcNormFactors(yPAR)

# RPKM
RPKM_PAR <- rpkm(yPAR)

# Manually load 'Tissue' & 'Patient' vectors
TissuePAR <- factor(c("N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T","N","T"))
PatientPAR <- factor(c("CH-5761","CH-5761","CH-5767","CH-5767","CH-5768","CH-5768","CH-5769","CH-5769","EJ-7115","EJ-7115","EJ-7123","EJ-7123","EJ-7125","EJ-7125","EJ-7314","EJ-7314","EJ-7315","EJ-7315","EJ-7317","EJ-7317","EJ-7321","EJ-7321","EJ-7327","EJ-7327","EJ-7328","EJ-7328","EJ-7330","EJ-7330","EJ-7331","EJ-7331","EJ-7781","EJ-7781","EJ-7782","EJ-7782","EJ-7783","EJ-7783","EJ-7784","EJ-7784","EJ-7785","EJ-7785","EJ-7786","EJ-7786","EJ-7789","EJ-7789","EJ-7792","EJ-7792","EJ-7793","EJ-7793","EJ-7794","EJ-7794","EJ-7797","EJ-7797","EJ-A8FO","EJ-A8FO","G9-6333","G9-6333","G9-6342","G9-6342","G9-6348","G9-6348","G9-6351","G9-6351","G9-6356","G9-6356","G9-6362","G9-6362","G9-6363","G9-6363","G9-6365","G9-6365","G9-6384","G9-6384","G9-6496","G9-6496","G9-6499","G9-6499","HC-7211","HC-7211","HC-7737","HC-7737","HC-7738","HC-7738","HC-7740","HC-7740","HC-7742","HC-7742","HC-7745","HC-7745","HC-7747","HC-7747","HC-7752","HC-7752","HC-7819","HC-7819","HC-8258","HC-8258","HC-8259","HC-8259","HC-8260","HC-8260","HC-8262","HC-8262","J4-A83J","J4-A83J"))

# Data frame construction including 'Patient' & 'Tissue' factors
data.frame(Sample=colnames(yPAR),PatientPAR,TissuePAR)

# Design matrix based on the two factors ('Tissue' and 'Patient')
designPAR <- model.matrix(~PatientPAR+TissuePAR)
rownames(designPAR) <- colnames(yPAR)

# Estimated dispersion calculation and visualization of the common dispersion
yPAR <- estimateDisp(yPAR, designPAR, robust=TRUE) 
yPAR$common.dispersion

# Likelihood-ratio test
fitPAR <- glmFit(yPAR, designPAR)
lrtPAR <- glmLRT(fitPAR)

# 'topTags' element with the result obtained (n=Inf to include all the results)
top_PAR <- topTags(lrtPAR,n=Inf)

head(top_PAR)

# under-expressed and over-expressed genes summary
summary(decideTests(lrtPAR))


##############################################################################################


