library(edgeR)
library(goseq)
library(GO.db)

DEG <- rownames(DEG_file)

universe <- rownames(GENE_BACKGROUND_file)

mytable <- as.integer(unique(universe)%in%DEG)

names(mytable) <- unique(universe)
table(mytable)

pwf=nullp(mytable,"hg19","geneSymbol")



############################ GO TERMS  ###########################################

GO.wall=goseq(pwf,"hg19","geneSymbol", test.cats = c("GO:BP","GO:CC","GO:MF"))

GO.BP <- GO.wall[GO.wall$ontology=="BP",]
GO.BP$FDRunder <- p.adjust(GO.BP$under_represented_pvalue, n=nrow(GO.BP))
GO.BP$FDRover <- p.adjust(GO.BP$over_represented_pvalue, n=nrow(GO.BP))

GO.CC <- GO.wall[GO.wall$ontology=="CC",]
GO.CC$FDRunder <- p.adjust(GO.CC$under_represented_pvalue, n=nrow(GO.CC))
GO.CC$FDRover <- p.adjust(GO.CC$over_represented_pvalue, n=nrow(GO.CC))

GO.MF <- GO.wall[GO.wall$ontology=="MF",]
GO.MF$FDRunder <- p.adjust(GO.MF$under_represented_pvalue, n=nrow(GO.MF))
GO.MF$FDRover <- p.adjust(GO.MF$over_represented_pvalue, n=nrow(GO.MF))

GO.BP <- GO.BP[GO.BP$FDRover<=0.05,]
GO.CC <- GO.CC[GO.CC$FDRover<=0.05,]
GO.MF <- GO.MF[GO.MF$FDRover<=0.05,]


######################### KEGG PATHWAYS ###########################################

mapPathwayToName <- function(organism) {
  KEGG_PATHWAY_LIST_BASE <- "http://rest.kegg.jp/list/pathway/"
  pathway_list_REST_url <- paste(KEGG_PATHWAY_LIST_BASE, organism, sep="")
  
  pathway_id_name <- data.frame()
  cont<-0
  for (line in readLines(pathway_list_REST_url)) {
    cont<-cont+1
    tmp <- strsplit(line, "\t")[[1]]
    pathway_id <- strsplit(tmp[1], organism)[[1]][2]
    pathway_name <- tmp[2]
    pathway_name <- strsplit(pathway_name, "\\s+-\\s+")[[1]][1]
    pathway_id_name[cont, 1] = pathway_id
    pathway_id_name[cont, 2] = pathway_name
    
  }
  names(pathway_id_name) <- c("path","pathway_name")
  pathway_id_name
}

  kegg_DE <- goseq(pwf,'hg19','geneSymbol',test.cats="KEGG")
  kegg_DE$FDRunder <- p.adjust(kegg_DE$under_represented_pvalue, n=nrow(kegg_DE))
  kegg_DE$FDRover <- p.adjust(kegg_DE$over_represented_pvalue, n=nrow(kegg_DE))
  
  kegg_path<-mapPathwayToName("hsa")
  idx<-match(kegg_DE$category,kegg_path$path)
  pathway_name<-kegg_path[idx,2]
  data<-cbind(kegg_DE,pathway_name)
  
  KEGG <- data[data$FDRover<=0.05,]







