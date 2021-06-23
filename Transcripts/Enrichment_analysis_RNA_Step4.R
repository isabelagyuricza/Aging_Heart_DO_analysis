############################################## Gene_enrichment analysis using ClusterProfile

######################## DO Heart - data

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(DOSE)
library(GO.db)
library(topGO)
library(org.Mm.eg.db)
library(GSEABase)
library(clusterProfiler)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(DESeq2)
library(tibble)
library(tidyverse)

######################################################## 1) Setting up the data

result_signif <- read.csv("GBRS/DE_analysis_heart/Results/results_signif_DE.csv",header = TRUE, row.names = 1) #Significant genes from DE analysis

result <- read.csv("GBRS/DE_analysis_heart/Results/results_DE.csv",header = TRUE, row.names = 1) #Whole genes from DE analysis - Going to use as background

genes_signif <- rownames(result_signif)

background <- rownames(result)


############### 1.1) Getting entrez id for input

gene_entrez <- bitr(genes_signif, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db) #There are some duplicated ENSEMBL IDs that match different entrez_ID (different transcripts or gene models), 
#######################################I'm going to keep only the unique genes.

temp_tibble <- tibble(ENSEMBL = rownames(result_signif),
                      log2FoldChange = result_signif$log2FoldChange)

gene_entrez <- gene_entrez %>% 
  left_join(temp_tibble, by = 'ENSEMBL')

rm(temp_tibble)

tb <- table(gene_entrez$ENSEMBL) #Taking the duplicated

tb[tb>1]

dup <- which(duplicated(gene_entrez$ENSEMBL))

gene_entrez <- gene_entrez[-dup,] #Removing the duplicated


############ 1.2) Creating the input to the function (vector with LFC from Deseq2 and entrez names)

gene_entrez_input <- gene_entrez$log2FoldChange


names(gene_entrez_input) <- gene_entrez$ENTREZID


############ 1.3) Doing the same for the background

background_entrez <- bitr(background, fromType = "ENSEMBL",
                   toType = c("ENTREZID", "SYMBOL"),
                   OrgDb = org.Mm.eg.db)


temp_tibble <- tibble(ENSEMBL = rownames(result),
                      log2FoldChange = result$log2FoldChange)

background_entrez <- background_entrez %>% 
  left_join(temp_tibble, by = 'ENSEMBL')

rm(temp_tibble)

tb <- table(background_entrez$ENSEMBL) #Taking the duplicated

tb[tb>1]

dup <- which(duplicated(background_entrez$ENSEMBL))

background_entrez <- background_entrez[-dup,] #Removing duplicated

########### 1.4) Creating the input to the function

background_entrez_input <- background_entrez$log2FoldChange


names(background_entrez_input) <- background_entrez$ENTREZID

#_____________________________________________________________________________________________________________


###################################################### 2) Gene ontology enrichment

go_enrich_output <- enrichGO(gene = names(gene_entrez_input),
                             OrgDb = org.Mm.eg.db,
                             ont = "BP", #Choosing BP, I tryed other ones and had the same results
                             universe = names(background_entrez_input),
                             qvalueCutoff = 0.1,
                             pvalueCutoff = 0.1,
                             pAdjustMethod = "BH",
                             readable = TRUE)

go_enrich_output_resul <- go_enrich_output@result[go_enrich_output@result$qvalue < 0.05,] #Taking only the significant pathways


write.csv(go_enrich_output_resul,"GBRS/gene_enrichment/Results/GO_enrich_clusterProfile.csv")


######### 2.1) Making some nice plots

pdf("GBRS/gene_enrichment/Results/dotplot_GO.pdf")

#quartz()
dotplot(go_enrich_output,x="Count",color="qvalue",title="Enriched GO terms")

dev.off()


pdf("GBRS/gene_enrichment_heart/Results/cnetplot_GO.pdf",width = 13,height = 10)

# Selecting from the result just the pathways I want.

go_enrich_output_resul_filt <- go_enrich_output_resul[c("GO:0014812","GO:0070838","GO:0006953","GO:0051147"),8]

# Subseting just the genes within each pathway

list <- list()
for (i in 1:length(go_enrich_output_resul_filt)) {
  list[i] <- str_split(go_enrich_output_resul_filt[i],"/",simplify = FALSE)}

# Keeping just the unique ones

symbols <- unlist(list) %>% unique()

# On my data.frame containing the LOD fold change (gene_entrez), keep just the information
# for the genes within the pathways.
gene_entrez_filt <- gene_entrez %>% 
  filter(SYMBOL %in% symbols)

# Creating fold change vector for cnetplot

gene_entrez_input <- gene_entrez_filt$log2FoldChange


names(gene_entrez_input) <- gene_entrez_filt$ENTREZID

# Scalling based on the fold change distribution for plotting

my_scaler <- function(x){
  x <- ifelse(x > 0.5,0.5,x)
  x <- ifelse(x < 0.5 & x > 0,0,x)
  x <- ifelse(x < 0,-0.5,x)
  return(x)
}

gene_entrez_input <- my_scaler(gene_entrez_input)

# Choosing the colors from http://colorbrewer2.org/#type=diverging&scheme=PiYG&n=5

colpal <- c('#2166ac','#ef8a62', '#b2182b')

# Subsetting the categories I want for plotting

categories <-go_enrich_output[go_enrich_output$ID %in% c("GO:0014812","GO:0070838","GO:0006953","GO:0051147")]  

tmp <- go_enrich_output

tmp@result <- categories

pdf("GBRS/gene_enrichment_heart/Results/cnetplot_GO_2.pdf", width = 12, height = 10)

#quartz()
cnetplot(tmp, foldChange = gene_entrez_input, colorEdge=F) +
  scale_colour_gradientn(colours = colpal, name = "Log Fold Change (LFC)", 
                         breaks = c(0.5,0,-0.5),labels = c("LFC > 0.5","0.5 < LFC < 0","LFC < 0")) +
  labs(colour = 'Log Fold-Change') +
  guides(size = F)

dev.off()
