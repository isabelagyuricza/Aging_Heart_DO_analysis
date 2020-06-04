######################## Protein complexes analysis ############################

# This script uses the information on the table 
# to gather the genes and proteins of each complex
# and check how their correlation is affected by age

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(tidyverse)
library(corrplot)
library(beeswarm)
library(ggbeeswarm)
library(ensimplR)

################################################################################
############ load data

load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")

complex_info <- readRDS(
"Protein_enrich_analysis/Protein_GK/Protein_complexes/ori_protein_complex_long_dat.RDS")

rm(genoprobs,K,map,markers,ensembl.version)

complex_info <- complex_info %>% 
  select(complex_name,mouse_ensembl_id,symbol) %>% 
  unique()

complex_info <- complex_info %>% 
  rename(gene.id = mouse_ensembl_id) %>% 
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.protein, by = "gene.id") %>% 
  select(complex_name,gene.id,protein.id,symbol.x) %>% 
  mutate(complex_name = gsub(" ","_",complex_name)) %>% 
  rename(symbol = symbol.x)

match <- c("Proteasome","NPC","cytoplasmic_ribosomal_small","cytoplasmic_ribosomal_large",
           "mitochondrial_ribosomal_large","Nucleosomal_protein","COPII",
           "mitochondrial_ribosomal_small","ATP_synthase","Cytochrome_c",
           "mitochondrial_pyruvate","COPI","Cytochrome_bc1","mitochondrial_inner_membrane",
           "Fatty_Acid_Import","mitochondrial_outer","ERAD","Respiratory_chain_complex_I")

complex_info <-  complex_info %>% 
  filter(grepl(paste(match,collapse = "|"),complex_name))
  
complex_info <- na.exclude(complex_info)

complex_info$complex_name[which(complex_info$complex_name == "F0/F1_ATP_synthase_(complex_V)")] <- "Mitochondrial_complex_V"

complex_info$complex_name[which(complex_info$complex_name == "Cytochrome_bc1_complex_(Ubiquinol-cytochrome_c_reductase_complex,_complex_III)")] <- "Mitochondrial_complex_III"

complex_info$complex_name[which(complex_info$complex_name == "Cytochrome_c_oxidase_(complex_IV)")] <- "Mitochondrial_complex_IV"

complex_info$complex_name[which(complex_info$complex_name == "mitochondrial_inner_membrane_presequence_translocase_complex")] <- "Mitochondrial_inner_membrane_translocase"

complex_info$complex_name[which(complex_info$complex_name == "mitochondrial_outer_membrane_translocase_complex")] <- "Mitochondrial_outer_membrane_translocase"

complex_info$complex_name[which(complex_info$complex_name == "Respiratory_chain_complex_I_(early_intermediate_NDUFAF1_assembly),__mitochondrial")] <- "Mitochondrial_complex_I"


genes_info <- batchGenes(c("ENSMUSG00000021577","ENSMUSG00000009863","ENSMUSG00000058076","ENSMUSG00000000171"))


more_info <- data_frame(complex_name = "Mitochondrial_complex_II",
                        gene.id = genes_info$id)

more_info <- more_info %>% 
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.protein, by = "gene.id") %>% 
  select(complex_name,gene.id,protein.id,symbol)

complex_info <- rbind(complex_info,more_info)

protein_expression <- data.frame(dataset.DOheart.protein.debatch.nopoly$data$norm[,
  which(colnames(dataset.DOheart.protein.debatch.nopoly$data$norm)
        %in% unique(complex_info$protein.id))]) %>% 
  rownames_to_column("mouse.id") %>% 
  gather("protein.id","protein_expression",-mouse.id)

gene_expression <- data.frame(dataset.DOheart.mrna$data$norm[,
    which(colnames(dataset.DOheart.mrna$data$norm)
        %in% unique(complex_info$gene.id))]) %>% 
  rownames_to_column("mouse.id") %>% 
  gather("gene.id","gene_expression",-mouse.id)

df <- complex_info %>% 
  left_join(protein_expression, by = "protein.id") %>% 
  left_join(gene_expression, by = c("gene.id","mouse.id")) %>% 
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.sample %>% 
              select(mouse.id, Age), by = "mouse.id")


####################### Testing RNA-protein correlation ########################

cor.SE <- function(x) {
  unname(sqrt((1 - x$estimate^2)/x$parameter))
}

df_cor <- df %>% 
  group_by(complex_name,Age) %>% 
  summarise(
    cor = cor.test(protein_expression,gene_expression)$estimate,
    SE = cor.SE(cor.test(protein_expression,gene_expression))
    ) %>% 
  ungroup()

pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/RNA_protein_correlation.pdf",
    width = 12, height = 8)

#quartz()
df_cor %>%
  mutate(Age = factor(Age, levels = sort(unique(Age)))) %>% 
  ggplot(aes(x = Age, y = cor)) +
  geom_point()+
  geom_errorbar(aes(ymin = cor - SE, ymax = cor + SE))+
  facet_wrap(~ complex_name, ncol = 4,scales = "free") +
  theme(strip.text.x = element_text(size = 7))

dev.off()

rm(df_cor)
# It doesn't seem to be a difference. The error bars are huge. 

######## Testing change in correlation among proteins from the same complex and 
######## how they change with age 

col <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7",
                          "#ef8a62", "#b2182b"))

pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/corrplots_protein.pdf",
    width = 14, height = 5.5)

order_rule <- vector('list', length = length(unique(df$complex_name)))
names(order_rule) <- unique(df$complex_name)

for (i in unique(df$complex_name)){

df_test_6 <- df %>% 
  filter(complex_name == i,
         Age == 6) %>% 
  select(mouse.id,protein.id,protein_expression) %>% 
  spread(protein.id,protein_expression) %>% 
  column_to_rownames("mouse.id")

df_test_12 <- df %>% 
  filter(complex_name == i,
         Age == 12) %>% 
  select(mouse.id, protein.id,protein_expression) %>% 
  spread(protein.id,protein_expression) %>% 
  column_to_rownames("mouse.id")

df_test_18 <- df %>% 
  filter(complex_name == i,
         Age == 18) %>% 
  select(mouse.id, protein.id,protein_expression) %>% 
  spread(protein.id,protein_expression) %>% 
  column_to_rownames("mouse.id")

par(mfrow=c(1,3))

COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")

CLUST <- hclust(dist(abs(COR_06)), method = 'ward.D2')
ORDER <- CLUST$order

order_rule[[i]] <- CLUST$labels[ORDER]

corrplot(COR_06[ORDER, ORDER],
         method="color", outline=FALSE, order="original",
         tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
         title = paste0(i, " - Age 6"),mar=c(0,0,1.5,0))

corrplot(COR_12[ORDER, ORDER],
         method="color", outline=FALSE, order="original",
         tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
         title = paste0(i, " - Age 12"),mar=c(0,0,1.5,0))

corrplot(COR_18[ORDER, ORDER],
         method="color", outline=FALSE, order="original",
         tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
         title = paste0(i, " - Age 18"),mar=c(0,0,1.5,0))
}

dev.off()


pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/corrplots_genes.pdf",
    width = 14, height = 5.5)

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == 6) %>% 
    select(mouse.id,gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6 <- df_test_6[,!is.na(colSums(df_test_6))]
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 12) %>% 
    select(mouse.id, gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df_test_12[,!is.na(colSums(df_test_12))]
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 18) %>% 
    select(mouse.id, gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df_test_18[,!is.na(colSums(df_test_18))]
  
  par(mfrow=c(1,3))
  #quartz()
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  CLUST <- hclust(dist(abs(COR_06)), method = 'ward.D2')
  
  ORDER <- CLUST$order
  
  corrplot(COR_06[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
           title = paste0(i, " - Age 6"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
           title = paste0(i, " - Age 12"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
           title = paste0(i, " - Age 18"),mar=c(0,0,1.5,0))
  
}

dev.off()


######## DIFFERENT PLOT FOR PAPER !


pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/corrplots_genes_paper.pdf",
    width = 14, height = 5.5)

for (i in c("26S_Proteasome","COPII", "Mitochondrial_complex_V")){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == 6) %>% 
    select(mouse.id,gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_6 <- df_test_6[,!is.na(colSums(df_test_6))]
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 12) %>% 
    select(mouse.id, gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df_test_12[,!is.na(colSums(df_test_12))]
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 18) %>% 
    select(mouse.id, gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df_test_18[,!is.na(colSums(df_test_18))]
  
  par(mfrow=c(1,3))
  #quartz()
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  CLUST <- hclust(dist(abs(COR_06)), method = 'ward.D2')
  
  df_filt <- df %>% 
    filter(complex_name == i) %>% 
    select(complex_name,gene.id,protein.id) %>% 
    unique() %>% 
    mutate(protein.id = factor(protein.id, levels = order_rule[[i]])) %>% 
    arrange(protein.id)
  
  ORDER <- df_filt$gene.id
  
  corrplot(COR_06[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
           title = paste0(i, " - Age 6"),mar=c(0,0,1.5,0))
  
  corrplot(COR_12[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
           title = paste0(i, " - Age 12"),mar=c(0,0,1.5,0))
  
  corrplot(COR_18[ORDER, ORDER],
           method="color", outline=FALSE, order="original",
           tl.cex=0.3,tl.col = "black",col = col(200),cl.lim=c(-1,1),
           title = paste0(i, " - Age 18"),mar=c(0,0,1.5,0))
  
}

dev.off()


################ Checking expression proteins/genes complexes ##################

pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/expression_across_age.pdf",
    width = 12, height = 10)

for (i in unique(df$complex_name)){

  #quartz()
 plot <-  df %>% 
    filter(complex_name == i) %>% 
    gather("expression_type","expression_value",c(protein_expression,gene_expression)) %>% 
    ggplot(aes(x = factor(Age, levels = c(6,12,18)), y = expression_value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(aes(color = symbol), priority = 'density', alpha = 0.5, groupOnX = TRUE, cex = .4) +
    labs(title = i, x = "Age", color = "Gene symbol") +
    facet_wrap(~ expression_type, scales = "free")
 
 print(plot)
  
}
  
dev.off()  


