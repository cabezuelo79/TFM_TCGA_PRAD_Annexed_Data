##############################################################
#######        RNASeq gene expression analysis       #########        
##############        TREATMENT FACTOR          ##############
##############################################################

dir <- #WORKING_DIRECTORY
setwd(dir)

library(edgeR)

  rawdata_drugs <- read.delim("Treatment_FILTER_TCGA_RNASeq_counts.txt", check.names=FALSE, stringsAsFactors=FALSE)
  
  # gene symbol as rownames
  rownames(rawdata_drugs)<-rawdata_drugs$SYMBOL
  
  # Manually load 'Drugs' vector (DR1 = treatment based on a single drug target, DR2 = treatment in which the drug target has been changed, DR3 = chemotherapy)
  Drugs <- factor(c("DR1","DR1","DR1","DR2","DR1","DR1","DR1","DR1","DR2","DR3","DR1","DR1","DR1","DR1","DR1","DR1","DR1","DR2","DR1","DR3","DR1","DR1","DR1","DR1","DR3","DR1","DR2","DR1","DR2","DR2","DR1","DR2","DR1","DR1","DR1","DR1","DR2","DR1","DR3","DR2","DR1","DR1","DR1","DR1","DR1","DR1","DR2","DR1","DR1","DR3","DR1","DR2","DR2","DR2","DR2","DR2","DR2","DR2","DR2","DR1","DR1","DR1","DR2","DR2","DR2","DR2","DR2","DR2","DR1","DR2","DR1","DR1","DR1"
  ))
  
  # Obtaining the size of each gene to normalize according to the size of the genes
  gene.length_D <-read.table("his-Size.tab",header=T)
  idx_D <-match(rawdata_drugs$SYMBOL,gene.length_D$Gene)
  results_counts_D <-gene.length_D[idx_D,]
  results_counts_D[is.na(results_counts_D$Length),"Length"]<-0
  nrow(results_counts_D)
  
  # DGEList object construction using 'Drugs' as group
  yD <- DGEList(counts=rawdata_drugs[,2:74], genes=results_counts_D, group=Drugs)
  
  # low expression genes deletion
  nrow(yD$counts)
  keep<-rowSums(cpm(yD)>1) >= 5  #number of 'DR3' samples (lowest sample group)
  yD<-yD[keep, ]
  yD$samples$lib.size <-colSums(yD$counts)
  nrow(yD$counts)
  
  # Data normalization
  yD <- calcNormFactors(yD)
  
  # GLM matrix design
  designD <- model.matrix(~0+Drugs, data=yD$samples)
  colnames(designD) <- levels(yD$samples$group)

  designD
  
  # Dispersion estimation (Common, Trended, Tagwise)
  yD <- estimateGLMCommonDisp(yD, designD)
  yD <- estimateGLMTrendedDisp(yD, designD)
  yD <- estimateGLMTagwiseDisp(yD, designD)
  

  # Quasi-likelihood test
  fitD <- glmQLFit(yD, designD)
  
  # CONTRASTS and topTags tables creation
  my.contrasts_D <- makeContrasts(DR2vsDR1=DR2-DR1, DR3vsDR2=DR3-DR2, DR3vsDR1=DR3-DR1, levels=designD)
  qlf.DR2vsDR1 <- glmQLFTest(fitD, contrast=my.contrasts_D[,"DR2vsDR1"])
  top_DR2vsDR1 <- topTags(qlf.DR2vsDR1, n=Inf)
  qlf.DR3vsDR2 <- glmQLFTest(fitD, contrast=my.contrasts_D[,"DR3vsDR2"])
  top_DR3vsDR2 <- topTags(qlf.DR3vsDR2, n=Inf)
  qlf.DR3vsDR1 <- glmQLFTest(fitD, contrast=my.contrasts_D[,"DR3vsDR1"])
  top_DR3vsDR1 <- topTags(qlf.DR3vsDR1, n=Inf)
 

 # under-expressed and over-expressed genes summaries
summary(decideTests(qlf.DR2vsDR1))
summary(decideTests(qlf.DR3vsDR2))
summary(decideTests(qlf.DR3vsDR1))

##############################################################

