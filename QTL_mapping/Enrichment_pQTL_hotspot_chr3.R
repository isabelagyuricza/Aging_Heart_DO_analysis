#################### Subsetting pQTL hotspot on chromosome 3 and doing enrichment analysis #########################

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

### Q x Age interactions on Chr 3

#####################
###  extract the lod peaks and start trimming them to just the hotspot

lod.peaks <- as_data_frame(dataset.DOheart.protein.debatch.nopoly$lod.peaks$age_int) %>%
  filter(substr(marker.id,1,2) == "3_") %>%
  separate(marker.id, into=c("qtl.chr","qtl.pos")) %>%
  mutate(qtl.chr=as.integer(qtl.chr), qtl.pos=as.double(qtl.pos)/(10^6)) %>%
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.protein)

# use this graphic to fine tune the edges of the hotspot
#quartz()
ggplot(lod.peaks, aes(x=qtl.pos, y=lod)) +
  geom_point() +
  xlim(145,150) #Region with the highest lod scores ~145Mb to 149Mb
#
lod.peaks <- dplyr::filter(lod.peaks, qtl.pos>145 & qtl.pos<149)
dim(lod.peaks)
# 167 genes

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

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/correlation_all_proteins_chr3.pdf",
    height = 12,width = 12)
#quartz()
corrplot(cor(test, use = "pairwise.complete.obs"),
         method="color", outline=FALSE, order="hclust",
         tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1))
dev.off()
#
# grab average absolute correlations and filter the lod peaks
tmp.m <- apply(abs(cor(expr.data, use = "pairwise.complete.obs")),2,mean)
sum(tmp.m > 0.3)  #130
sum(tmp.m > 0.35) #110
sum(tmp.m > 0.4)  #74
#
lod.peaks <- mutate(lod.peaks, mean.corr = tmp.m)

write.csv(lod.peaks,file = "GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/cor_table_proteins_chr3.csv")

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

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/protein_chr3/cor_heatmap_protein_age_int_chr3.pdf")

#quartz()
heatmap.2(cor(test, use = "pairwise.complete.obs"),
          Rowv = TRUE,Colv=TRUE,dendrogram ="column",
          trace="none",cexCol = 0.5,cexRow = 0.5,key = TRUE,
          keysize = 0.9, key.title = "Color_key",
          col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

dev.off()

#There are a group of genes that seem to be positively correlated, such as: 
# histones and proteins responsible for fiber/sarcomere organization and function
# (myosins, colagen)

# reorder genes in lod.peaks and expr.data

lod.peaks <- lod.peaks[hm$colInd,]
expr.data <- expr.data[,hm$colInd]

###################################

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

save(enrichments,file = "GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/protein_chr3/enrichments.RData")

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Results/protein_chr3/enrich_proteinscor_age_int_chr3.pdf")

categories <- enrichments$CC[enrichments$CC$ID %in% c("GO:0000502","GO:0032982","GO:0000786")]

tmp <- enrichments$CC

tmp@result <- categories

#quartz()
cnetplot(tmp,colorEdge=TRUE) + guides(size = F) +
  ggtitle("Chr3 CC, q < 0.05")

dev.off()

# Different from the previous enrichment analysis with all the protein in the transbands,
# in which we found zero enriched pathways,
# when using just the proteins that have a negative/positive correlation with each other
# we see enriched cellular components terms for myosin filaments, nucleosome and proteosome complexes.
