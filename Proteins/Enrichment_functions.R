#########################################################
#  wrapper functions for enrichment analysis
#########################################################

#########################################################
####
# wrapper function to call enrichGO..
WrapEnrichGO <- function(GeneList, Universe, ont="BP"){
  enrichGO(gene = GeneList,
           OrgDb = org.Mm.eg.db,
           keyType = 'ENSEMBL',
           ont = ont,
           universe = Universe,
           qvalueCutoff = 0.1,
           pvalueCutoff = 0.1,
           pAdjustMethod = "BH",
           readable=TRUE)
}

#########################################################
####
# KEGG enrichment analysis requires some extra functions

# a function to replace a string of / separated entrez ids with gene symbols
KeggParse <- function(entrez.string){
  gene.id <- as_vector(str_split(entrez.string,"/"))
  symbol.string <- NULL
  for(i in 1:length(gene.id)){
    symbol.string <- str_c(symbol.string,
                           annot.gene[which(annot.gene$entrez_id==gene.id[i])[1],]$symbol,"/")
  }
  str_sub(symbol.string,1,-2)[1]
}
# # test
#  entrez.string <- "70358/22337/20531/64177/15203/17750/71775/105243/12309/320718"
# KeggParse(entrez.string)
# rm(entrez.string)

# function to apply kegg.parse() to enrichKEGG output @results slot
KeggReplace <- function(kegg.result){
  for(i in 1:dim(kegg.result)[1]){
    kegg.result[i,"geneID"] <- KeggParse(kegg.result[i,"geneID"])
  }
  kegg.result
}

WrapEnrichKegg <- function(GeneList, Universe){
  kegg.out <- enrichKEGG(gene = GeneList,
                         organism = "mmu",
                         universe = Universe,
                         qvalueCutoff = 0.1,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
  
  kegg.out@result <- KeggReplace(kegg.out@result)
  
  kegg.out
}
# qvals <- ProteinTest$q.Age.M7
# thresh <- 0.01

#########################################################
# filters for cleaning up the enrichments

###
# filter enrichment results with q-value threshold
enrich.filter <- function(enrichment, thresh){
  
  # keep only significant enrichments
  enrichment@result <- enrichment@result[enrichment@result$qvalue < thresh,]
  
  enrichment
}
# enrichments <- lapply(enrichments, enrich.filter, 0.05)


####
# trim long category descriptions for cleaner plotting
enrich.trim <- function(enrichment, len=42){
  
  
  # trim category names to help with plotting
  enrichment@result$Description <- str_sub(enrichment@result$Description, 1, len)
  
  enrichment
}
# enrichments <- lapply(enrichments, enrich.trim, 42)


