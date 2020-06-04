###################### DO Heart mRNA permutations ############################

# Analyzing permutations for mRNA obtained by Sumner jobs (using rslurm) 

# loading data

library(tidyverse)

load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")

# Reading the results scped from Sumner:

result.list <- list.files("_rslurm_permutation_apply_gene", pattern = "results_")

genes_list <- list()
for (i in result.list){
  genes_list[[i]] <- readRDS(paste0("_rslurm_permutation_apply_gene/",i))
}; rm(i)

# Checking the number of genes, it should be 9694 in total.

check_genes <- c()
for (i in 1:length(genes_list)){
  check_genes[i] <- length(genes_list[[i]])
}; rm(i)

sum(check_genes) #[1] 9694 Perfect, all genes are here.

# Setting up the list:

lods_list <- unlist(genes_list,recursive = FALSE)

names(lods_list) <- str_split(names(lods_list),"RDS.",simplify = TRUE)[,2]

# Taking the 95% quantile for each gene.

for (i in 1:length(lods_list)){
  if (class(lods_list[[i]]) != "numeric") {print(i)}
}

#Nothing, no problem.

# Taking the lod threshold based on the 95% quantile

gene_quantile <- lapply(lods_list, function(x) quantile(sort(x), 0.95))

lods_genes <- bind_rows(gene_quantile) %>% 
  gather("gene.id", "lod_threshold") %>% 
  left_join(dataset.DOheart.mrna$lod.peaks$age_int, by = "gene.id") %>% 
  select(gene.id, lod_threshold,lod, qtl.chr, qtl.pos, cis) %>% 
  rename(lod_scan = lod) %>% 
  mutate(signif = lod_scan >= lod_threshold)


summary(lods_genes$lod_threshold)

# summary(lods_genes$lod_threshold)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.053   7.289   7.601   7.641   7.936  11.519 

pdf("GBRS/QTL_viewer/heart_newversion/Results/permutation/Garys_method/allgenes_permutations_hist.pdf")
hist(lods_genes$lod_threshold, main = "Histogram of lod threshold for all genes")
dev.off()


# Checking the p-values for each gene.

gene_pval <- c()

for (i in 1:nrow(lods_genes)){
  
  gene <- lods_genes$gene.id[i]
  
  lods_perm <- lods_list[[gene]]
  
  lod_scan <- lods_genes$lod_scan[i]
  
  gene_pval[i] <- mean(lods_perm >= lod_scan)
  
  names(gene_pval)[i] <- gene
} 

lods_genes <- lods_genes %>% 
  mutate(pval = gene_pval)

dim(lods_genes[lods_genes$pval < 0.05,]) #[1] 1302    9

dim(lods_genes[lods_genes$pval < 0.05 & lods_genes$cis == TRUE,]) #[1] 6 9

# Focus on that significant CIS QTLs to describe some of the QTLs.
