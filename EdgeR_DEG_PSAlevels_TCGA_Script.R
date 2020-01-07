##############################################################
#######        RNASeq gene expression analysis       #########        
##############        PSA LEVELS FACTOR         ##############
##############################################################

dir <- #WORKING_DIRECTORY
setwd(dir)

library(edgeR)

rawdata_psa <- read.delim("PSAlevels_FILTER_TCGA_RNASeq_counts.txt", check.names=FALSE, stringsAsFactors=FALSE)

# gene symbol as rownames
rownames(rawdata_psa)<-rawdata_psa$SYMBOL

# Manually load 'psa' vector (PSA0 = PSA levels under 10 ng/mL, PSA1 = PSA levels from 10 to 20 ng/mL, PSA2 = PSA levels over 20 ng/mL)
psa <- factor(c("PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA1","PSA0","PSA1","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA1","PSA0","PSA2","PSA0","PSA1","PSA1","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA1","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA1","PSA1","PSA0","PSA0","PSA0","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA2","PSA0","PSA1","PSA1","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA0","PSA2","PSA0","PSA0","PSA0"))

# Obtaining the size of each gene to normalize according to the size of the genes
gene.length_psa <-read.table("his-Size.tab",header=T)
idx_psa <-match(rawdata_psa$SYMBOL,gene.length_psa$Gene)
results_counts_psa <-gene.length_psa[idx_psa,]
results_counts_psa[is.na(results_counts_psa$Length),"Length"]<-0
nrow(results_counts_psa)

# DGEList object construction using 'psa' as group
yPSA <- DGEList(counts=rawdata_psa[,2:441], genes=results_counts_psa, group=psa)

# low expression genes deletion
nrow(yPSA$counts)
keep<-rowSums(cpm(yD)>1) >= 11  #number of 'PSA2' samples (lowest sample group)
yPSA<-yPSA[keep, ]
yPSA$samples$lib.size <-colSums(yPSA$counts)
nrow(yPSA$counts)

# Data normalization
yPSA <- calcNormFactors(yPSA)

# GLM matrix design
designPSA <- model.matrix(~0+psa, data=yPSA$samples)
colnames(designPSA) <- levels(yPSA$samples$group)

designPSA

# Dispersion estimation (Common, Trended, Tagwise)
yPSA <- estimateGLMCommonDisp(yPSA, designPSA)
yPSA <- estimateGLMTrendedDisp(yPSA, designPSA)
yPSA <- estimateGLMTagwiseDisp(yPSA, designPSA)


# Quasi-likelihood test
fitPSA <- glmQLFit(yPSA, designPSA)

# CONTRASTS and topTags tables creation
my.contrasts_PSA <- makeContrasts(PSA1vsPSA0=PSA1-PSA0, PSA2vsPSA1=PSA2-PSA1, PSA2vsPSA0=PSA2-PSA0, levels=designPSA)
qlf.PSA1vsPSA0 <- glmQLFTest(fitPSA, contrast=my.contrasts_PSA[,"PSA1vsPSA0"])
top_PSA1vsPSA0 <- topTags(qlf.PSA1vsPSA0, n=Inf)
qlf.PSA2vsPSA1 <- glmQLFTest(fitPSA, contrast=my.contrasts_PSA[,"PSA2vsPSA1"])
top_PSA2vsPSA1 <- topTags(qlf.PSA2vsPSA1, n=Inf)
qlf.PSA2vsPSA0 <- glmQLFTest(fitPSA, contrast=my.contrasts_PSA[,"PSA2vsPSA0"])
top_PSA2vsPSA0 <- topTags(qlf.PSA2vsPSA0, n=Inf)


 # under-expressed and over-expressed genes summaries
summary(decideTests(qlf.PSA1vsPSA0))
summary(decideTests(qlf.PSA2vsPSA1))
summary(decideTests(qlf.PSA2vsPSA0))

##############################################################
