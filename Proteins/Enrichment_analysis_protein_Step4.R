############################################## Enrichment analysis - heart DO proteins dataset ######################################

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
library(tidyverse)
library(stringr)
library(broom)
library(ensimplR)

# Load wrapper functions for enrichment analysis using ClusterProfiler (Gary Churchill)

source("Protein_enrich_analysis/Protein_GK/Scripts/Enrichment_functions.R")
# library loads this will mask some dplyr functions


# Load anova results

load("Protein_enrich_analysis/Protein_GK/Results/result_anova_protein.RData")

enrichments_age <- list(
  protein.age.BP = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Age<0.05)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="BP"),
  protein.age.CC = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Age<0.05)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="CC"),
  protein.age.MF = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Age<0.05)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="MF"))

enrichments_age_filtered <- lapply(enrichments_age, enrich.filter,0.05) 

dim(enrichments_age_filtered$protein.age.BP)
# [1] 2 9

dim(enrichments_age_filtered$protein.age.CC)
#[1] 30  9

# Using this padj cutoff we found positive regulation of cellular proliferation, 
# intracellular protein transport and some mitochondria pathways in the CC.


# Enrichment for sex:age interaction

  enrichments_int <- list(protein.int.BP = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Int<0.05)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="BP"),
  protein.int.CC = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Int<0.05)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="CC"),
  protein.int.MF = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Int<0.05)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="MF"))
  
  enrichments_int_filtered <- lapply(enrichments_int, enrich.filter,0.1) 
  
# dim(enrichments_int_filtered$protein.int.BP)
# [1] 16  9
  
  
# Enrichment for sex
  
  enrichments_sex <- list(protein.sex.BP = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Sex<0.001)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="BP"),
  protein.sex.CC = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Sex<0.001)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="CC"),
  protein.sex.MF = WrapEnrichGO(
    unique(filter(df_final_wth, padj.Sex<0.001)$Gene_ID),
    unique(df_final_wth$Gene_ID), ont="MF")   )

enrichments_sex_filtered <- lapply(enrichments_sex, enrich.filter,0.05)

# dim(enrichments_sex_filtered$protein.sex.BP)
# [1] 30  9

attach(enrichments_age_filtered)

pdf("Protein_enrich_analysis/Protein_GK/Results/emapplot_age_1.pdf",width = 8,height = 8)
#quartz()
emapplot(protein.age.BP) + ggtitle("Age BP, q < 0.05")
emapplot(protein.age.CC) + ggtitle("Age CC, q < 0.05")
emapplot(protein.age.MF) + ggtitle("Age MF, q < 0.05")

dev.off()


############################### BP cnetplot ####################################

fc <- tibble(estimate = df_final_wth$estimate.Age,
             id = df_final_wth$Gene_ID)

all_symbols <- batchGenes(unique(fc$id)) %>% 
  select(id,symbol) %>% 
  left_join(fc, by = "id")

result_filt <- protein.age.BP@result[,8]

# Subseting just the genes within each pathway

list <- list()
for (i in 1:length(result_filt)) {
  list[i] <- str_split(result_filt[i],"/",simplify = FALSE)}

# Keeping just the unique ones

cat_symbols <- unlist(list) %>% unique()

# On my data.frame containing the LOD fold change (gene_entrez), keep just the information
# for the genes within the pathways.

all_symbols_filt <- all_symbols %>% 
  filter(symbol %in% cat_symbols)

# Creating fold change vector for cnetplot

input <- all_symbols_filt$estimate

names(input) <- all_symbols_filt$id

# Scalling based on the fold change distribution for plotting

my_scaler <- function(x){
  x <- ifelse(x > 0.5,0.5,x)
  x <- ifelse(x < 0.5 & x > 0,0,x)
  x <- ifelse(x < 0,-0.5,x)
  return(x)
}

input <- my_scaler(input)

colpal <- c('#2166ac','#ef8a62', '#b2182b')

pdf("Protein_enrich_analysis/Protein_GK/Results/cnetplot_age_BP_2.pdf",width = 12,height = 10)

cnetplot(protein.age.BP, foldChange = input, colorEdge=F) +
  scale_colour_gradientn(colours = colpal, name = "Slope", 
                         breaks = c(0.5,0,-0.5),labels = c("slope > 0.5","0.5 < slope < 0","slope < 0")) +
  guides(size = F)+
  ggtitle("Protein Age BP, q < 0.05")

dev.off()

############################### CC cnetplot ####################################


result_filt <- protein.age.CC@result[c("GO:0098800","GO:0098803"),8]

# Subseting just the genes within each pathway

list <- list()
for (i in 1:length(result_filt)) {
  list[i] <- str_split(result_filt[i],"/",simplify = FALSE)}

# Keeping just the unique ones

cat_symbols <- unlist(list) %>% unique()

# On my data.frame containing the LOD fold change (gene_entrez), keep just the information
# for the genes within the pathways.

all_symbols_filt <- all_symbols %>% 
  filter(symbol %in% cat_symbols)

# Creating fold change vector for cnetplot

input <- all_symbols_filt$estimate

names(input) <- all_symbols_filt$id

# Scalling based on the fold change distribution for plotting

my_scaler <- function(x){
  x <- ifelse(x > 0.5,0.5,x)
  x <- ifelse(x < 0.5 & x > 0,0,x)
  x <- ifelse(x < 0,-0.5,x)
  return(x)
}

input <- my_scaler(input)

colpal <- c('#2166ac','#ef8a62', '#b2182b')

categories <-protein.age.CC[protein.age.CC@result$ID %in% 
                              c("GO:0098800","GO:0098803")]  

tmp <- protein.age.CC

tmp@result <- categories

pdf("Protein_enrich_analysis/Protein_GK/Results/cnetplot_age_CC_2.pdf",width = 12,height = 10)

cnetplot(tmp, foldChange = input, colorEdge=F) +
  scale_colour_gradientn(colours = colpal, name = "Slope", 
                         breaks = c(0.5,0,-0.5),labels = c("slope > 0.5","0.5 < slope < 0","slope < 0")) +
  guides(size = F)+
  ggtitle("Protein Age CC, q < 0.05")

dev.off()


########################### BOTH ON THE SAME PLOT ##############################


fc <- tibble(estimate = df_final_wth$estimate.Age,
             id = df_final_wth$Gene_ID)

all_symbols <- batchGenes(unique(fc$id)) %>% 
  select(id,symbol) %>% 
  left_join(fc, by = "id")

result_filt <- c(protein.age.BP@result[,8],protein.age.CC@result[c("GO:0098800","GO:0098803"),8])

# Subseting just the genes within each pathway

list <- list()
for (i in 1:length(result_filt)) {
  list[i] <- str_split(result_filt[i],"/",simplify = FALSE)}

# Keeping just the unique ones

cat_symbols <- unlist(list) %>% unique()

# On my data.frame containing the LOD fold change (gene_entrez), keep just the information
# for the genes within the pathways.

all_symbols_filt <- all_symbols %>% 
  filter(symbol %in% cat_symbols)

# Creating fold change vector for cnetplot

input <- all_symbols_filt$estimate

names(input) <- all_symbols_filt$id

# Scalling based on the fold change distribution for plotting

my_scaler <- function(x){
  x <- ifelse(x > 0.5,0.5,x)
  x <- ifelse(x < 0.5 & x > 0,0,x)
  x <- ifelse(x < 0,-0.5,x)
  return(x)
}

input <- my_scaler(input)

colpal <- c('#2166ac','#ef8a62', '#b2182b')


categories <-rbind(protein.age.BP@result,protein.age.CC@result[protein.age.CC@result$ID %in% 
                                                                 c("GO:0098800","GO:0098803"),])  

tmp <- protein.age.CC

tmp@result <- categories

pdf("Protein_enrich_analysis/Protein_GK/Results/cnetplot_age_all.pdf",width = 20,height = 16)

#quartz()
cnetplot(tmp, foldChange = input, colorEdge=F) +
  scale_colour_gradientn(colours = colpal, name = "Slope", 
                         breaks = c(0.5,0,-0.5),labels = c("slope > 0.5","0.5 < slope < 0","slope < 0")) +
  guides(size = F)+
  ggtitle("Protein Age BP and CC, q < 0.05")

dev.off()


#############################################################################

attach(enrichments_int)

pdf("Protein_enrich_analysis/Protein_GK/Results/emapplot_int.pdf",width = 8,height = 8)
#quartz()
emapplot(protein.int.BP) + ggtitle("Int BP, q < 0.1")
#emapplot(protein.int.CC) + ggtitle("Int CC, q < 0.1")
emapplot(protein.int.MF) + ggtitle("Int MF, q < 0.1")

dev.off()


fc <- df_final_wth$estimate.Int
names(fc) <- df_final_wth$Gene_ID

colpal <- c('#d7b5d8','#df65b0','#ce1256')

pdf("Protein_enrich_analysis/Protein_GK/Results/cnetplot_int.pdf",width = 8,height = 6)

#quartz()
cnetplot(enrich.trim(protein.int.BP,32), foldChange = fc, colorEdge=F) +
  scale_colour_gradientn(colours = colpal) +
  labs(colour = 'Slope') +
  guides(size = F)+
  ggtitle("Protein Age:Sex BP, q < 0.1")

dev.off()

detach(enrichments_int)


#############################################################################

attach(enrichments_sex)

pdf("Protein_enrich_analysis/Protein_GK/Results/emapplot_sex.pdf",width = 10,height = 8)

quartz()
emapplot(protein.sex.BP) + ggtitle("Sex BP, q < 0.05")
emapplot(protein.sex.CC) + ggtitle("Sex CC, q < 0.05")
emapplot(protein.sex.MF) + ggtitle("Sex MF, q < 0.05")

dev.off()

fc <- df_final_wth$estimate.Sex
names(fc) <- df_final_wth$Gene_ID

colpal <- c('#4dac26','#b8e186','#f7f7f7','#f1b6da','#d01c8b')

pdf("Protein_enrich_analysis/Protein_GK/Results/cnetplot_sex.pdf",width = 10,height = 8)

#quartz()
cnetplot(enrich.trim(protein.sex.BP,32), foldChange = fc, colorEdge=F) +
  scale_colour_gradientn(colours = colpal,limits=c(-1,2)) +
  labs(colour = 'Slope') +
  guides(size = F) +
  ggtitle("Protein Sex BP, q < 0.05")

dev.off()

################## Saving enrichment data

save(enrichments_age_filtered,file = "Protein_enrich_analysis/Protein_GK/Results/enrichment_age.RData")

save(enrichments_int_filtered,file = "Protein_enrich_analysis/Protein_GK/Results/enrichment_int.RData")

save(enrichments_sex_filtered,file = "Protein_enrich_analysis/Protein_GK/Results/enrichment_sex.RData")


