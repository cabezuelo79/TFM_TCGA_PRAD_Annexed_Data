##############################################################
#######        RNASeq gene expression analysis       #########        
##############        GLEASON SCORE FACTOR      ##############
##############################################################

dir <- #WORKING_DIRECTORY
setwd(dir)

library(edgeR)

rawdataG <- read.delim("TCGA_RNASeq_counts.txt", check.names=FALSE, stringsAsFactors=FALSE)

# gene symbol as rownames
rownames(rawdataG)<-rawdataG$SYMBOL
  
  
# Manually load 'Gleason' vector (G0 = Non-tumoral samples, G1 = Gleason score under 7, G2 = Gleason score of 7 or more)
Gleason <- factor(c("G1","G1","G2","G1","G2","G2","G2","G2","G1","G1","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G2","G2","G2","G2","G2","G0","G2","G0","G1","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G0","G2","G0","G2","G2","G2","G0","G2","G0","G2","G0","G2","G2","G0","G1","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G0","G2","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G1","G2","G2","G0","G2","G2","G2","G2","G0","G1","G2","G1","G0","G2","G0","G2","G2","G2","G0","G2","G2","G0","G2","G0","G2","G2","G0","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G0","G2","G2","G2","G0","G2","G2","G0","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G1","G1","G2","G2","G2","G2","G1","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G0","G2","G0","G2","G0","G2","G0","G2","G2","G0","G2","G0","G2","G1","G2","G2","G0","G2","G2","G2","G0","G2","G2","G2","G1","G2","G2","G2","G0","G1","G0","G1","G0","G2","G2","G0","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G0","G2","G1","G2","G2","G2","G1","G2","G1","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G1","G2","G2","G1","G2","G2","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G1","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2"
  ))
  
 # Obtaining the size of each gene to normalize according to the size of the genes
  gene.length_G <-read.table("his-Size.tab",header=T)
  idx_G <-match(rawdataG$SYMBOL,gene.length_G$Gene)
  results_counts_G <-gene.length_G[idx_G,]
  results_counts_G[is.na(results_counts_G$Length),"Length"]<-0
  nrow(results_counts_G)
  
 # DGEList object construction using 'Gleason' as group
  yG <- DGEList(counts=rawdataG[,2:550], genes=results_counts_G, group=Gleason)
  
 # low expression genes deletion
  nrow(yG$counts)
  keep<-rowSums(cpm(yG)>1) >= 45 #number of 'G1' samples (lowest sample group)
  yG<-yG[keep, ]
  yG$samples$lib.size <-colSums(yG$counts)
  nrow(yG$counts)
  
  # Data normalization
  yG <- calcNormFactors(yG)
  
  # GLM matrix design
  designG <- model.matrix(~0+Gleason, data=yG$samples)
  colnames(designG) <- levels(yG$samples$group)

  designG
  
  # Dispersion estimation (Common, Trended, Tagwise)
  yG <- estimateGLMCommonDisp(yG, designG)
  yG <- estimateGLMTrendedDisp(yG, designG)
  yG <- estimateGLMTagwiseDisp(yG, designG)
  
  
  # Quasi-likelihood test
  fitG <- glmQLFit(yG, designG)
  
  # CONTRASTS and topTags tables creation
  my.contrasts <- makeContrasts(G1vsG0=G1-G0, G2vsG1=G2-G1, G2vsG0=G2-G0, levels=designG)
  qlf.G1vsG0 <- glmQLFTest(fitG, contrast=my.contrasts[,"G1vsG0"])
  top_G1vsG0 <- topTags(qlf.G1vsG0, n=Inf)
  qlf.G2vsG1 <- glmQLFTest(fitG, contrast=my.contrasts[,"G2vsG1"])
  top_G2vsG1 <- topTags(qlf.G2vsG1, n=Inf)
  qlf.G2vsG0 <- glmQLFTest(fitG, contrast=my.contrasts[,"G2vsG0"])
  top_G2vsG0 <- topTags(qlf.G2vsG0, n=Inf)
  
 # under-expressed and over-expressed genes summaries
  summary(decideTests(qlf.G1vsG0))
  summary(decideTests(qlf.G2vsG1))
  summary(decideTests(qlf.G2vsG0))
  
 ############################################################## 
