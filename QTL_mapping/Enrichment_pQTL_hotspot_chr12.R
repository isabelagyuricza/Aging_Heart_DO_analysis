#################### Subsetting pQTL hotspot on chromosome 12 and doing enrichment analysis #########################

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(gplots)
library(RColorBrewer)
library(corrplot)
library(qtl2)
library(pcaMethods)
library(tidyverse)
library(stringr)
library(broom)

################################################################################
############ load data

# load QTLViewer version of data
load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")

################################################################################

### for transbands of protein interactive scan for age

### Q x Age interactions on Chr 12

#####################
###  extract the lod peaks and start trimming them to just the hotspot

lod.peaks <- as_data_frame(dataset.DOheart.protein.debatch.nopoly$lod.peaks$age_int) %>%
  filter(substr(marker.id,1,3) == "12_") %>%
  separate(marker.id, into=c("qtl.chr","qtl.pos")) %>%
  mutate(qtl.chr=as.integer(qtl.chr), qtl.pos=as.double(qtl.pos)/(10^6)) %>%
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.protein)

# use this graphic to fine tune the edges of the hotspot
quartz()
ggplot(lod.peaks, aes(x=qtl.pos, y=lod)) +
  geom_point() +
  xlim(98,102) #Region with the highest lod scores ~98Mb to 102Mb

lod.peaks <- dplyr::filter(lod.peaks, qtl.pos>98 & qtl.pos<102)
dim(lod.peaks)
# 177 genes

############################
# grab normalized mRNA expr data, regress out Sex (adjusting for sex), then rankz it.

# expression data
expr.data <- dataset.DOheart.protein.debatch.nopoly$data$norm[,lod.peaks$protein.id]

# compute residuals to adjust x wrt covariate
myresids <- function(x, cov){
  residuals(lm(x~cov,na.action=na.exclude))
}

#
expr.data <- apply(expr.data,2,myresids,dataset.DOheart.protein.debatch.nopoly$annot.samples$Sex)

# Apply a normal scores transform to the data.
rankZ <- function(x) {
  x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))
}

expr.data <- apply(expr.data,2,rankZ) #Across mice

rm(myresids, rankZ)

###
# order the data by QTL peak location
lod.peaks <- arrange(lod.peaks, by=qtl.pos)
expr.data <- expr.data[,lod.peaks$protein.id]
###

############################
# look at correlations
test <- expr.data

colnames(test) <- lod.peaks$symbol

col <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7",
                          "#ef8a62", "#b2182b"))

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/correlation_all_proteins_chr12.pdf",
    height = 12,width = 12)
#quartz()
corrplot(cor(test, use = "pairwise.complete.obs"),
         method="color", outline=FALSE, order="hclust",
         tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1))
dev.off()
#
# grab average absolute correlations and filter the lod peaks
tmp.m <- apply(abs(cor(expr.data, use = "pairwise.complete.obs")),2,mean)
sum(tmp.m > 0.3)  #161
sum(tmp.m > 0.35) #152
sum(tmp.m > 0.4)  #136
#
lod.peaks <- mutate(lod.peaks, mean.corr = tmp.m)

write.csv(lod.peaks,file = "GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/cor_table_proteins_chr12.csv")
#
lod.peaks <- filter(lod.peaks, mean.corr>0.3)
expr.data <- expr.data[,lod.peaks$protein.id]


# hierarchical clustering
quartz()
hm <- heatmap.2(cor(expr.data, use = "pairwise.complete.obs"),
                Rowv = TRUE,Colv=TRUE,dendrogram ="column",
                trace="none",#cexCol = 0.5,
                #colsep = c(64,128),sepcolor = "black",
                keysize = 0.9, key.title = "Color_key",margins = c(8,10),
                col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

test <- expr.data

colnames(test) <- lod.peaks$symbol


# hierarchical clustering

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/protein_chr12/cor_heatmap_protein_age_int_chr12.pdf")

#quartz()
heatmap.2(cor(test, use = "pairwise.complete.obs"),
          Rowv = TRUE,Colv=TRUE,dendrogram ="column",
          trace="none",cexCol = 0.5,cexRow = 0.5,key = TRUE,
          keysize = 0.9, key.title = "Color_key",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

dev.off()

## There are a group of genes that are positvely correlated to each other, and there 
## are some of them that also mapped to the hotspot on chromosome 3 (myh, col, and histones)
## I'm going to check if they are the same isoforms. 
## In adittion to that, there are other extracellular matrix components, such as Fbn1, 
## Hspg2 adnd other isoforms of collagens that are positively correlated. 

# reorder genes in lod.peaks and expr.data
lod.peaks <- lod.peaks[hm$colInd,]
expr.data <- expr.data[,hm$colInd]

#####################
# GeneSet Enrichment analysis only on the proteins that are correlated (abs_cor > 3)

###
# load wrapper functions
source("Protein_enrich_analysis/Protein_GK/Scripts/Enrichment_Functions.R")
library(org.Mm.eg.db)
library(clusterProfiler)

# library loads this will mask some dplyr functions

###
# compute enrichment results and capture as a named list
# so we can lapply functions to all at once or attach()

annots <- dataset.DOheart.protein.debatch.nopoly$annot.protein

enrichments <- list(
  BP = WrapEnrichGO(lod.peaks$gene.id,annots$gene.id, ont="BP"),
  CC = WrapEnrichGO(lod.peaks$gene.id,annots$gene.id, ont="CC"),
  MF = WrapEnrichGO(lod.peaks$gene.id,annots$gene.id, ont="MF")  )

###
# filter the enrichments to retain only the most significant
enrichments <- lapply(enrichments, enrich.filter, 0.05)

save(enrichments, file = "GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/protein_chr12/enrichments.RData")

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/protein_chr12/enrich_proteinscor_age_int_chr12.pdf")

categories <- enrichments$BP[enrichments$BP$ID %in% c("GO:0006936","GO:0098781")]  

tmp <- enrichments$BP

tmp@result <- categories
  
#quartz()
cnetplot(tmp,colorEdge=TRUE) + guides(size = F) +
  ggtitle("Chr12 BP, q < 0.05")

#quartz()
cnetplot(enrichments$CC,colorEdge=TRUE,
         showCategory = 5) + guides(size = F) +
  ggtitle("Chr12 CC, q < 0.05")

dev.off()

# This enrichment analysis using only the correlated proteins also show us 
# pathways associated with cardiac muscle contraction and contractile fibers. In
# contrast with the previous enrichment analysis that showed us only mitochondria 
# respiratory chain components. 
