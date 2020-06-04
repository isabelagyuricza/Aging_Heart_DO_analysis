############################################### DE analysis heart - DO 

####### This script is for running the DE analysis using DESeq2 for heart data - DO

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(hexbin)
library(biomaRt)

load("data_setup/RNA_GBRS/rna_heart_expr.RData")

annot_samples$Age <- annot_samples$Age/6

annot_samples$Generation <- as.factor(annot_samples$Generation)

dds <- DESeqDataSetFromMatrix(countData = t(floor(mrna.counts_filtered)),
                              colData = annot_samples,
                              design = ~ Sex + Generation + Age)

rm(mrna.counts_filtered)

dds_vst <- vst(dds,blind = TRUE)

dds_de <- DESeq(dds)

result <- results(dds_de)

####################### Diagnostic plots

pdf("DE_analysis/diagnostic_plots/plotPCA.pdf")

#quartz()
plotPCA(dds_vst,intgroup=c("Sex","Generation","Age"))+
  theme(legend.position = "bottom")
#quartz()
plotPCA(dds_vst,intgroup=c("Sex","Age"))+
  theme(legend.position = "bottom")

dev.off()

pdf("DE_analysis/diagnostic_plots/MAplots.pdf")
#quartz()
plotMA(result)
dev.off()

pdf("DE_analysis/diagnostic_plots/histogram_plots.pdf")
#quartz()
hist(result$pvalue,main = "Histogram of pvalues",xlab = "p-values") 

dev.off()

######################################## Subsetting significant genes (padj < 0.1)

result_signif <- result[which(result$padj < 0.1),] # We have 2287 genes that change with age..

write.csv(result,"GBRS/DE_analysis/Results/results_DE.csv")

save(dds_de,file = "GBRS/DE_analysis/Results/DE_object.RData")

write.csv(result_signif,"GBRS/DE_analysis/Results/results_signif_DE.csv")

############################## For continuous variables, the reported LFC means the log2 fold change is per unit of change of that variable!
############################## So, a LFC of -1 for one gene, means that for each unit of age, the gene increases in 2x.
############################### A LFC of -1 for one gene, means that for each unit of age, the gene decreases 1/2.

