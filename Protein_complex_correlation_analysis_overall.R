##################### Protein complexes analysis overall #######################

# This script uses the information on Greg's table 
# to gather the genes and proteins of each complex
# and check how their correlation is affected by age.

# Gathering the correlation values into protein complexes, applying the regression
# and permutting..

# Author: Isabela Gerdes Gyuricza
# Date: 07_01_2020

################################################################################
############ clear workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(tidyverse)
library(corrplot)
library(ensimplR)
library(DGCA)

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
           "mitochondrial_ribosomal_large","COPII",
           "mitochondrial_ribosomal_small","ATP_synthase","Cytochrome_c",
           "mitochondrial_pyruvate","COPI","Cytochrome_bc1","mitochondrial_inner_membrane"
           ,"mitochondrial_outer","Respiratory_chain_complex_I")

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

################## Fitting regression for each complex #########################

# 1) Proteins

age_trend_real_protein_list <- list()

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
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
  
  df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                   cor_12 = COR_12[upper.tri(COR_12)],
                   cor_18 = COR_18[upper.tri(COR_18)],
                   comparing_this = dimnames(COR_06)[[1]][names[,1]],
                   to_this = dimnames(COR_06)[[2]][names[,2]])
  
  df_cor <- df_cor %>% 
    unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
    gather("Age","cor",-pair_comparison) %>% 
    mutate(Age = gsub("cor_","",Age)) 
  
  MODEL <- suppressMessages(lmerTest::lmer(
    cor ~ (1 + Age | pair_comparison) + Age,
    data = df_cor %>% mutate(Age = as.numeric(Age))
  ))
  
  print(lmerTest::ranova(MODEL))
  
  FIXEF <- lme4::fixef(MODEL) 
  
  coefs <- summary(MODEL)$coefficients
  
  coefs <- coefs["Age",which(colnames(coefs) %in% c("Std. Error","Pr(>|t|)"))]
  
  names(coefs) <- c("SE","pvalue_model")
  
  AGE_EFFECT <- FIXEF[2]
  
  age_trend_real_protein_list[[i]] <- c(AGE_EFFECT, coefs)
  
}

rm(i,df_test_12,df_test_18,df_test_6,COR_06,COR_12,COR_18,df_cor,
   names,match,more_info)

# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

permut_age_trend <- function(df, times, seed = 19940418) {
  
  set.seed(seed = seed)
  age_trend_permut_protein_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- vector(mode = 'numeric', length = times)
    
    for (j in 1:times) {
      
      df_test_6 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id),64)) %>% 
        select(mouse.id,protein.id,protein_expression) %>% 
        spread(protein.id,protein_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_12 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id[
                 !df$mouse.id %in% rownames(df_test_6)
                 ]),62)) %>% 
        select(mouse.id, protein.id,protein_expression) %>% 
        spread(protein.id,protein_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_18 <- df %>% 
        filter(complex_name == i,
               !mouse.id %in% rownames(df_test_6) &
                 !mouse.id %in% rownames(df_test_12)) %>% 
        select(mouse.id, protein.id,protein_expression) %>% 
        spread(protein.id,protein_expression) %>% 
        column_to_rownames("mouse.id")
      
      COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
      COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
      COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
      
      names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
      
      df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                       cor_12 = COR_12[upper.tri(COR_12)],
                       cor_18 = COR_18[upper.tri(COR_18)],
                       comparing_this = dimnames(COR_06)[[1]][names[,1]],
                       to_this = dimnames(COR_06)[[2]][names[,2]])
      
      df_cor <- df_cor %>% 
        unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
        gather("Age","cor",-pair_comparison) %>% 
        mutate(Age = gsub("cor_","",Age)) 
      
      MODEL <- suppressMessages(lme4::lmer(
        cor ~ (1 + Age | pair_comparison) + Age,
        data = df_cor %>% mutate(Age = as.numeric(Age))
      ))
      
      FIXEF <- lme4::fixef(MODEL)
      
      AGE_EFFECT <- FIXEF[2]
      
      age_trend_perm[j] <- as.numeric(AGE_EFFECT)
      
    }
    
    age_trend_permut_protein_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_protein_list)
}

age_trend_permut_protein_results <- permut_age_trend(df,1000)


save(age_trend_permut_protein_results, 
     file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_overall_cor_method/age_trend_perms_protein.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_protein_results)){
  
  complex_df_perm <- age_trend_permut_protein_results[[i]]
  
  complex_df_real <- age_trend_real_protein_list[[i]][[1]]
  
    
    pval <- bigEmpPVals(stat = abs(complex_df_real), 
                        stat0 = abs(complex_df_perm))
    
  empirical_pvals[[i]] <- pval
  
}

empirical_pvals <- bind_rows(empirical_pvals) %>% 
  gather("complex_name","pvalue_perm_protein")


final_df <- data.frame(age_trend_real_protein_list)

colnames(final_df)[which(colnames(final_df) == "X26S_Proteasome")] <- "26S_Proteasome"

colnames(final_df)[which(colnames(final_df) == "Nuclear_pore_complex_.NPC.")] <- "Nuclear_pore_complex_(NPC)"


final_df <- final_df %>% 
  rownames_to_column("Index") %>% 
  gather("complex_name","Value",-Index) %>% 
  spread("Index","Value") %>% 
  rename(CorTrendAge_Protein = Age,
         pvalue_model_Protein = pvalue_model,
         SE_model_Protein = SE) %>% 
  left_join(empirical_pvals, by = "complex_name")



################################################################################


# 1) Genes

age_trend_real_gene_list <- list()

for (i in unique(df$complex_name)){
  
  df_test_6 <- df %>% 
    filter(complex_name == i,
           Age == 6) %>% 
    select(mouse.id,gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_12 <- df %>% 
    filter(complex_name == i,
           Age == 12) %>% 
    select(mouse.id, gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  df_test_18 <- df %>% 
    filter(complex_name == i,
           Age == 18) %>% 
    select(mouse.id, gene.id,gene_expression) %>% 
    unique() %>% 
    spread(gene.id,gene_expression) %>% 
    column_to_rownames("mouse.id")
  
  COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
  COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
  COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
  
  names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
  
  df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                   cor_12 = COR_12[upper.tri(COR_12)],
                   cor_18 = COR_18[upper.tri(COR_18)],
                   comparing_this = dimnames(COR_06)[[1]][names[,1]],
                   to_this = dimnames(COR_06)[[2]][names[,2]])
  
  df_cor <- df_cor %>% 
    unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
    gather("Age","cor",-pair_comparison) %>% 
    mutate(Age = gsub("cor_","",Age)) 
  
  MODEL <- suppressMessages(lmerTest::lmer(                       
    cor ~ (1 + Age | pair_comparison) + Age,
    data = df_cor %>% mutate(Age = as.numeric(Age))
  ))
  
  print(lmerTest::ranova(MODEL))
  
  FIXEF <- lme4::fixef(MODEL)
  
  coefs <- summary(MODEL)$coefficients
  
  coefs <- coefs["Age",which(colnames(coefs) %in% c("Std. Error","Pr(>|t|)"))]
  
  names(coefs) <- c("SE","pvalue_model")
  
  
  AGE_EFFECT <- FIXEF[2]
  
  age_trend_real_gene_list[[i]] <- c(AGE_EFFECT, coefs) 
}

rm(i,df_test_12,df_test_18,df_test_6,COR_06,COR_12,COR_18,df_cor,
   names)

# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

permut_age_trend <- function(df, times, seed = 20200701) {
  
  set.seed(seed = seed)
  age_trend_permut_gene_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- vector(mode = 'numeric', length = times)
    
    for (j in 1:times) {
      
      df_test_6 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id),64)) %>% 
        select(mouse.id,gene.id,gene_expression) %>% 
        unique() %>% 
        spread(gene.id,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_12 <- df %>% 
        filter(complex_name == i,
               mouse.id %in% sample(unique(df$mouse.id[
                 !df$mouse.id %in% rownames(df_test_6)
                 ]),62)) %>% 
        select(mouse.id, gene.id,gene_expression) %>% 
        unique() %>% 
        spread(gene.id,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      df_test_18 <- df %>% 
        filter(complex_name == i,
               !mouse.id %in% rownames(df_test_6) &
                 !mouse.id %in% rownames(df_test_12)) %>% 
        select(mouse.id, gene.id,gene_expression) %>% 
        unique() %>% 
        spread(gene.id,gene_expression) %>% 
        column_to_rownames("mouse.id")
      
      COR_06 <- cor(df_test_6, use = "pairwise.complete.obs")
      COR_12 <- cor(df_test_12, use = "pairwise.complete.obs")
      COR_18 <- cor(df_test_18, use = "pairwise.complete.obs")
      
      names <- which( upper.tri(COR_06,diag=F) , arr.ind = TRUE )
      
      df_cor <- tibble(cor_6 = COR_06[upper.tri(COR_06)],
                       cor_12 = COR_12[upper.tri(COR_12)],
                       cor_18 = COR_18[upper.tri(COR_18)],
                       comparing_this = dimnames(COR_06)[[1]][names[,1]],
                       to_this = dimnames(COR_06)[[2]][names[,2]])
      
      df_cor <- df_cor %>% 
        unite(col = "pair_comparison", c("comparing_this","to_this"),sep = "_") %>% 
        gather("Age","cor",-pair_comparison) %>% 
        mutate(Age = gsub("cor_","",Age)) 
      
      MODEL <- suppressMessages(lme4::lmer(
        cor ~ (1 + Age | pair_comparison) + Age,
        data = df_cor %>% mutate(Age = as.numeric(Age))
      ))
      
      FIXEF <- lme4::fixef(MODEL)
      
      AGE_EFFECT <- FIXEF[2]
      
      age_trend_perm[j] <- as.numeric(AGE_EFFECT)
      
    }
    
    age_trend_permut_gene_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_gene_list)
}

age_trend_permut_gene_results <- permut_age_trend(df,1000)


save(age_trend_permut_gene_results, 
     file = "protein_enrich_analysis/protein_GK/protein_complexes/regression_overall_cor_method/age_trend_perms_gene.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_gene_results)){
  
  complex_df_perm <- age_trend_permut_gene_results[[i]]
  
  complex_df_real <- age_trend_real_gene_list[[i]][[1]]
  
  
  pval <- bigEmpPVals(stat = abs(complex_df_real), 
                      stat0 = abs(complex_df_perm))
  
  empirical_pvals[[i]] <- pval
  
}

empirical_pvals <- bind_rows(empirical_pvals) %>% 
  gather("complex_name","pvalue_perm_gene")


final_df_2 <- data.frame(age_trend_real_gene_list)

colnames(final_df_2)[which(colnames(final_df_2) == "X26S_Proteasome")] <- "26S_Proteasome"

colnames(final_df_2)[which(colnames(final_df_2) == "Nuclear_pore_complex_.NPC.")] <- "Nuclear_pore_complex_(NPC)"


final_df <- final_df_2 %>% 
  rownames_to_column("Index") %>% 
  gather("complex_name","Value",-Index) %>% 
  spread("Index","Value") %>% 
  rename(CorTrendAge_Gene = Age,
         pvalue_model_Gene = pvalue_model,
         SE_model_Gene = SE) %>% 
  left_join(empirical_pvals, by = "complex_name") %>% 
  left_join(final_df, by = "complex_name") %>% 
  select(complex_name, CorTrendAge_Gene,pvalue_model_Gene,pvalue_perm_gene,SE_model_Gene,
         CorTrendAge_Protein,pvalue_model_Protein,pvalue_perm_protein,SE_model_Protein)

write.csv(final_df, file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_overall_cor_method/age_trend_overall.csv")

plot_df <- final_df %>% 
  mutate(STD_CorTrengAge_Gene = CorTrendAge_Gene/SE_model_Gene,
         STD_CorTrengAge_Protein = CorTrendAge_Protein/SE_model_Protein,
         sig_perm_protein = ifelse(pvalue_perm_protein < 0.05, TRUE,FALSE),
         sig_perm_gene = ifelse(pvalue_perm_gene < 0.05, TRUE,FALSE),
         sig_model_protein = ifelse(pvalue_model_Protein < 0.05, TRUE,FALSE),
         sig_model_gene = ifelse(pvalue_model_Gene < 0.05, TRUE,FALSE)) %>% 
  select(complex_name,STD_CorTrengAge_Gene,sig_perm_gene,sig_model_gene,
         STD_CorTrengAge_Protein,sig_perm_protein, sig_model_protein)


# After talking to Greg he said that we should be more conservative and consider
# the pvalues permutations more reliable.. So, for that reason: ..

pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_overall_cor_method/STD_coef_scatterplot_complexes2.pdf",
    width = 13,height = 10)
#quartz()
plot_df %>%
  mutate(
    `Significance (p-value < 0.05)` = ifelse(
      plot_df$sig_perm_gene == TRUE & plot_df$sig_perm_protein == FALSE, 'Transcript Only',
      ifelse(
        plot_df$sig_perm_gene == FALSE & plot_df$sig_perm_protein == TRUE, 'Protein Only',
        ifelse(
          plot_df$sig_perm_gene == TRUE & plot_df$sig_perm_protein == TRUE, 'Both',
          'Neither'
        )
      )
    ),
    `Significance (p-value < 0.05)` = factor(
      `Significance (p-value < 0.05)`,
      levels = c('Transcript Only', 'Protein Only', 'Both', 'Neither')
      )
  ) %>% 
  ggplot(aes(x = STD_CorTrengAge_Gene, y = STD_CorTrengAge_Protein)) +
  geom_point(aes(colour = `Significance (p-value < 0.05)`), size = 5,alpha = 0.7) + 
  geom_text(aes(label = complex_name),hjust = 0, vjust = 0, size = 3) +
  scale_x_continuous(
    lim = c(-34,34),
    breaks = seq(-200, 200, 20),
    minor_breaks =  seq(-200, 200, 10)
  ) +
  scale_y_continuous(
    lim = c(-34, 34),
    breaks = seq(-200, 200, 20),
    minor_breaks = seq(-200, 200, 10)
  ) +
    scale_colour_manual(breaks = c('Transcript Only', 'Protein Only', 'Both', 'Neither'),
                      values = c("firebrick","royalblue","mediumorchid4","gray")) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_hline(yintercept = 0,linetype = "dashed") +
theme_minimal() 
  
dev.off()








