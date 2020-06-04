######################## Protein complexes analysis ############################

# This script uses the information on table 
# to gather the genes and proteins of each complex
# and check how their correlation is affected by age.

# Testing the regression analysis and permutation.

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

match <- c("Proteasome","ATP_synthase","Cytochrome_c",
           "COPI","Cytochrome_bc1",
           "Respiratory_chain_complex_I")

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


######## Testing change in correlation among proteins from the same complex and 
######## how they change with age 


################ Fitting regression line for each pair within each complex ################

# 1) Proteins

est_age_trend <- function(x, age){
  MODEL <- lm(x ~ age)
  MODEL <- summary(MODEL)
  AGE_TREND <- MODEL$coefficients[2,1]
  return(AGE_TREND)
}

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
  
  age_trend_real <- df_cor %>% 
    mutate(Age = as.numeric(Age)) %>% 
    group_by(pair_comparison) %>% 
    summarise(
      CorrTrend = est_age_trend(cor, Age)) %>% 
    ungroup()
  
  age_trend_real_protein_list[[i]] <- age_trend_real
  
}

#### Checking if direction of correlated is separate into groups within complex

list_trend <- list()

for(i in 1:length(age_trend_real_protein_list)) {
  
  t <- age_trend_real_protein_list[[i]] %>% 
    separate(pair_comparison, into = c("protein.id","protein_2"),sep = "_") %>% 
    left_join(df %>% select(symbol,protein.id) %>% unique(), 
              by ="protein.id") %>% 
    select(symbol,protein_2,CorrTrend) %>% 
    rename(protein_1 = symbol,
           protein.id = protein_2) %>% 
    left_join(df %>% select(protein.id,symbol) %>% unique()
              , by = "protein.id") %>% 
    select(protein_1, symbol, CorrTrend) %>% 
    rename(protein_2 = symbol)
  
  list_trend[[i]] <- t
  
}

rm(i,df_test_12,df_test_18,df_test_6,age_trend_real,COR_06,COR_12,COR_18,df_cor,
   names,match,more_info)

# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

permut_age_trend <- function(df, times) {
  
  age_trend_permut_protein_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- list()
    
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
      
      age_trend <- df_cor %>% 
        mutate(Age = as.numeric(Age)) %>% 
        group_by(pair_comparison) %>% 
        summarise(
          CorrTrend = est_age_trend(cor, Age)) %>% 
        ungroup()
      
      age_trend_perm[[j]] <- age_trend
      
    }
    
    age_trend_permut_protein_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_protein_list)
}

age_trend_permut_protein_results <- permut_age_trend(df,1000)

# For each complex and each pair, compute the 95% quantile confidence

age_trend_permut_protein_results <- sapply(age_trend_permut_protein_results,bind_rows, 
                                           USE.NAMES = TRUE, simplify = FALSE)

save(age_trend_permut_protein_results, 
     file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_cor_method/age_trend_perms_protein.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_protein_results)){
  
  complex_df_perm <- age_trend_permut_protein_results[[i]]
  
  complex_df_real <- age_trend_real_protein_list[[i]]
  
  complex_pair_comparison <- tibble(
    pair_comparison = complex_df_real$pair_comparison,
    CorrTrend = complex_df_real$CorrTrend,
    pval = 0
  )
  
  for (j in unique(complex_df_real$pair_comparison)) {
    
    complex_df_perm_subset <- complex_df_perm %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    complex_df_real_subset <- complex_df_real %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    pval <- bigEmpPVals(stat = abs(complex_df_real_subset), 
                        stat0 = abs(complex_df_perm_subset))
    
    complex_pair_comparison$pval[
      complex_pair_comparison$pair_comparison == j
      ] <- pval
    
  }
  empirical_pvals[[i]] <- complex_pair_comparison
  
}


adjust <- function (x){
  
  padj <- adjustPVals(x$pval, "BH")
  
  x$pval.adj <- padj
  
  return(x)
  
}


complex_protein_age_trend_final <- sapply(empirical_pvals, adjust, simplify = FALSE)

complex_protein_age_trend_final <- bind_rows(complex_protein_age_trend_final,.id = "complex_name")

complex_protein_age_trend_final_significant <- complex_protein_age_trend_final %>% 
  filter(pval.adj < 0.1)

save(complex_protein_age_trend_final_significant, 
     file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_cor_method/significant_protein_pair.RData")

write.csv(complex_protein_age_trend_final_significant, 
     file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_cor_method/significant_protein_pair.csv")


# We have significant changes in correlation with age for basically all the complexes: 
# proteasome, COPII and all the mitochondrial complexes. For all the proteasome pairs 
# that significantly change with age, all of them decrease with age. THe other complexes 
# have some pairs increasing correlation with age, but the majority decrease (103 out of 112)

# Plot the heatmaps for the age trend for each pair

age_trend_real_protein <- bind_rows(age_trend_real_protein_list,
                                    .id = "complex_name")

col <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7",
                          "#ef8a62", "#b2182b"))

order_rule <- vector('list', length = length(unique(age_trend_real_protein$complex_name)))
names(order_rule) <- unique(age_trend_real_protein$complex_name)

pdf("Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_cor_method/age_trend_heatmap_protein.pdf", width = 8)

for (i in unique(age_trend_real_protein$complex_name)){
  
  d <- age_trend_real_protein %>% 
    mutate(signif = ifelse(age_trend_real_protein$pair_comparison %in% 
                             complex_protein_age_trend_final_significant$pair_comparison,
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    filter(complex_name == i) %>% 
    select(pair_comparison,CorrTrend, signif) %>% 
    separate(pair_comparison, into = c("protein_1","protein_2"),sep = "_")
  
  # Define order
  IDs <- df %>% 
    filter(complex_name == i) %>% 
    group_by(complex_name, gene.id, protein.id, symbol) %>% 
    summarise() %>% 
    ungroup()
  
  TREND_MAT <- matrix(
    data = 0,
    nrow = length(unique(c(d$protein_1, d$protein_2))),
    ncol = length(unique(c(d$protein_1, d$protein_2))),
    dimnames = list(unique(c(d$protein_1, d$protein_2)), unique(c(d$protein_1, d$protein_2)))
  )
  for(r in 1:nrow(d)){
    TREND_MAT[d$protein_1[r], d$protein_2[r]] <- d$CorrTrend[r]
    TREND_MAT[d$protein_2[r], d$protein_1[r]] <- d$CorrTrend[r]
  }
  CLUSTERING <- hclust(dist(TREND_MAT), method = 'ward.D2')
  ORDER <- colnames(TREND_MAT)[CLUSTERING$order]
  order_rule[[i]] <- IDs %>% 
    mutate(
      protein.id = factor(protein.id, levels = ORDER)
    ) %>% 
    arrange(protein.id) %>% 
    mutate(
      gene.id = factor(gene.id, levels = unique(gene.id))
    ) %>% 
    arrange(gene.id)
  
  d <- d %>% 
    mutate(
      protein_1_old = factor(protein_1, levels = levels(order_rule[[i]]$protein.id)),
      protein_2_old = factor(protein_2, levels = levels(order_rule[[i]]$protein.id))
    ) %>% 
    mutate(
      protein_1 = ifelse(
        as.numeric(protein_1_old) < as.numeric(protein_2_old),
        as.character(protein_1_old), as.character(protein_2_old)
      ),
      protein_2 = ifelse(
        as.numeric(protein_1_old) < as.numeric(protein_2_old),
        as.character(protein_2_old), as.character(protein_1_old)
      )
    ) %>% 
    mutate(
      protein_1 = factor(protein_1, levels = levels(order_rule[[i]]$protein.id)),
      protein_2 = factor(protein_2, levels = levels(order_rule[[i]]$protein.id))
    )
  
  SIZE <- 32/length(unique(d$protein_1))
  
  #quartz()
  plot <- d %>% 
    ggplot(aes(protein_1, protein_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$protein_2))) - .5, color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$protein_1))) - .5, color="#eeeeee", size = .25) +
    ggtitle(paste0('Age trend estimates - ',i)) +
    theme_bw() +
    xlab('protein.id') +
    ylab('protein.id') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age estimate", limits = c(-0.065, 0.065), breaks = seq(-0.06, 0.06, 0.03)) +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$protein_1)[-length(levels(d$protein_1))]
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$protein_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 7, colour = "black"),
      axis.text.y=element_text(size = 7, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )
  
  print(plot)
  
}

dev.off()


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
  
  age_trend_real <- df_cor %>% 
    mutate(Age = as.numeric(Age)) %>% 
    group_by(pair_comparison) %>% 
    summarise(
      CorrTrend = est_age_trend(cor, Age)) %>% 
    ungroup()
  
  age_trend_real_gene_list[[i]] <- age_trend_real
  
}

rm(i,df_test_12,df_test_18,df_test_6,age_trend_real,COR_06,COR_12,COR_18,df_cor)


# Creating a function to do that permuting mice 1000 times, disrupting only 
#the mouse - age association.

permut_age_trend <- function(df, times) {
  
  age_trend_permut_gene_list <- list()
  
  for (i in unique(df$complex_name)){
    
    age_trend_perm <- list()
    
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
      
      age_trend <- df_cor %>% 
        mutate(Age = as.numeric(Age)) %>% 
        group_by(pair_comparison) %>% 
        summarise(
          CorrTrend = est_age_trend(cor, Age)) %>% 
        ungroup()
      
      age_trend_perm[[j]] <- age_trend
      
    }
    
    age_trend_permut_gene_list[[i]] <- age_trend_perm
    
    
  }
  
  return(age_trend_permut_gene_list)
}

age_trend_permut_gene_results <- permut_age_trend(df,1000)

age_trend_permut_gene_results <- sapply(age_trend_permut_gene_results,bind_rows, 
                                        USE.NAMES = TRUE, simplify = FALSE)

save(age_trend_permut_gene_results, 
     file = "protein_enrich_analysis/protein_GK/protein_complexes/regression_cor_method/age_trend_perms_gene.RData")


# Using bigEmpPVals function from DGCA to compute the pvals from the permutation

empirical_pvals <- list() 

for (i in names(age_trend_permut_gene_results)){
  
  complex_df_perm <- age_trend_permut_gene_results[[i]]
  
  complex_df_real <- age_trend_real_gene_list[[i]]
  
  complex_pair_comparison <- tibble(
    pair_comparison = complex_df_real$pair_comparison,
    CorrTrend = complex_df_real$CorrTrend,
    pval = 0
  )
  
  for (j in unique(complex_df_real$pair_comparison)) {
    
    complex_df_perm_subset <- complex_df_perm %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    complex_df_real_subset <- complex_df_real %>% 
      filter(pair_comparison == j) %>% 
      select(CorrTrend) %>% 
      unlist()
    
    pval <- bigEmpPVals(stat = abs(complex_df_real_subset), 
                        stat0 = abs(complex_df_perm_subset))
    
    complex_pair_comparison$pval[
      complex_pair_comparison$pair_comparison == j
      ] <- pval
    
  }
  empirical_pvals[[i]] <- complex_pair_comparison
  
}


complex_gene_age_trend_final <- sapply(empirical_pvals, adjust, simplify = FALSE)

complex_gene_age_trend_final <- bind_rows(complex_gene_age_trend_final,.id = "complex_name")

complex_gene_age_trend_final_significant <- complex_gene_age_trend_final %>% 
  filter(pval.adj < 0.1)

save(complex_gene_age_trend_final_significant, 
     file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_cor_method/significant_gene_pair.RData")

write.csv(complex_gene_age_trend_final_significant, 
          file = "Protein_enrich_analysis/Protein_GK/Protein_complexes/regression_cor_method/significant_gene_pair.csv")

# All the changes in RNA are positive. And they only appear for some mitochondria
# complexes. 

# Plot the heatmaps for the age trend for each pair

age_trend_real_gene <- bind_rows(age_trend_real_gene_list,
                                 .id = "complex_name")

col <- colorRampPalette(c("#2166ac", "#67a9cf", "#d1e5f0", "#f7f7f7", "#fddbc7",
                          "#ef8a62", "#b2182b"))

pdf("protein_enrich_analysis/protein_GK/protein_complexes/regression_cor_method/age_trend_heatmap_gene.pdf",width = 8)

for (i in unique(age_trend_real_gene$complex_name)){
  
  d <- age_trend_real_gene %>% 
    mutate(signif = ifelse(age_trend_real_gene$pair_comparison %in% 
                             complex_gene_age_trend_final_significant$pair_comparison,
                           "TRUE", "FALSE")) %>% 
    mutate(signif = as.logical(signif)) %>% 
    filter(complex_name == i) %>% 
    select(pair_comparison,CorrTrend,signif) %>% 
    separate(pair_comparison, into = c("gene_1","gene_2"),sep = "_")
  
  d <- d %>% 
    mutate(
      gene_1_old = factor(gene_1, levels = levels(order_rule[[i]]$gene.id)),
      gene_2_old = factor(gene_2, levels = levels(order_rule[[i]]$gene.id))
    ) %>% 
    mutate(
      gene_1 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_1_old), as.character(gene_2_old)
      ),
      gene_2 = ifelse(
        as.numeric(gene_1_old) < as.numeric(gene_2_old),
        as.character(gene_2_old), as.character(gene_1_old)
      )
    ) %>% 
    mutate(
      gene_1 = factor(gene_1, levels = levels(order_rule[[i]]$gene.id)),
      gene_2 = factor(gene_2, levels = levels(order_rule[[i]]$gene.id))
    )
  
  SIZE <- 32/length(unique(d$gene_1))
  
  #quartz()
  plot <- d %>% 
    ggplot(aes(gene_1, gene_2)) +
    geom_hline(yintercept = seq(2, length(levels(d$gene_2))) - .5, color="#eeeeee", size = .25) + 
    geom_vline(xintercept = seq(2, length(levels(d$gene_1))) - .5, color="#eeeeee", size = .25) +
    ggtitle(paste0('Age trend estimates - ',i)) +
    theme_bw() +
    xlab('gene.id') +
    ylab('gene.id') +
    geom_tile(aes(fill = CorrTrend), color='#eeeeee') +
    geom_point(data = d %>% filter(signif), color='black', shape = 8, size = SIZE, stroke = SIZE/2) +
    scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b" , 
                         name = "Age estimate", limits = c(-0.065, 0.065), breaks = seq(-0.06, 0.06, 0.03)) +
    scale_x_discrete(
      drop = FALSE, 
      limits = levels(d$gene_1)[-length(levels(d$gene_1))]
    ) +
    scale_y_discrete(
      drop = FALSE, 
      limits = levels(d$gene_2)[-1]
    ) +
    theme(
      axis.text.x=element_text(angle=90, size = 7, colour = "black"),
      axis.text.y=element_text(size = 7, colour = 'black'),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      panel.border=element_blank(),
      panel.grid.major = element_blank()
    )
  
  print(plot)
  
}

dev.off()


test_gene <- complex_gene_age_trend_final_significant %>% 
  separate(pair_comparison, into = c("gene.id","gene_2"),sep = "_") %>% 
  left_join(df %>% select(gene.id,symbol,complex_name, protein.id) %>% unique(), 
            by =c("gene.id","complex_name")) %>% 
  select(complex_name,symbol,gene_2,CorrTrend,protein.id) %>% 
  rename(gene_1 = symbol,
         gene.id = gene_2) %>% 
  left_join(df %>% select(gene.id,symbol,complex_name) %>% unique()
            , by = c("gene.id","complex_name")) %>% 
  select(complex_name, gene_1, symbol, CorrTrend, protein.id) %>% 
  rename(gene_2 = symbol)

test_protein <- complex_protein_age_trend_final_significant %>% 
  separate(pair_comparison, into = c("protein.id","protein_2"),sep = "_")

test <- test_gene %>% 
  right_join(test_protein, by = c("protein.id","complex_name")) %>%
  left_join(df %>% select(protein.id,symbol,complex_name) %>% unique(), 
            by =c("protein.id","complex_name")) %>%  
  select(complex_name,gene_1,gene_2,CorrTrend.x,symbol,protein_2,
         CorrTrend.y) %>% 
  rename(protein.id = protein_2,
         protein_1 = symbol) %>% 
  left_join(df %>% select(protein.id,symbol,complex_name) %>% unique(), 
            by =c("protein.id","complex_name")) %>%  
  select(complex_name, gene_1, gene_2,CorrTrend.x,protein_1,symbol, CorrTrend.y) %>% 
  rename(protein_2 = symbol)

# There is just one pair comparison, in which RNA cor increase with age and protein
# cor decreasew with age: 

#Mitochondrial_complex_IV	Cox5a	Cox7c	0.02462800	Cox5a	Cox7c	-0.05137750
