###################### DO Heart protein permutations ############################

# Analyzing permutations for proteins obtained by Sumner jobs (using rslurm) 

# loading data

library(tidyverse)

load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")

# Reading the results scped from Sumner:

result.list <- list.files("_rslurm_permutation_apply_protein/", pattern = "results_")

proteins_list <- list()
for (i in result.list){
  proteins_list[[i]] <- readRDS(paste0("_rslurm_permutation_apply_protein/",i))
}; rm(i)

# Checking the number of proteins, it should be 2731 in total.

check_proteins <- c()
for (i in 1:length(proteins_list)){
  check_proteins[i] <- length(proteins_list[[i]])
}; rm(i)

sum(check_proteins) #[1] 2731 Perfect, all proteins are here.

# Setting up the list:

lods_list <- unlist(proteins_list,recursive = FALSE)

names(lods_list) <- str_split(names(lods_list),"RDS.",simplify = TRUE)[,2]

# Taking the 95% quantile for each gene.

for (i in 1:length(lods_list)){
  if (class(lods_list[[i]]) != "numeric") {print(i)}
}

#Nothing, no problem.

# Taking the lod threshold based on the 95% quantile

protein_quantile <- lapply(lods_list, function(x) quantile(sort(x), 0.95))

lods_proteins <- bind_rows(protein_quantile) %>% 
  gather("protein.id", "lod_threshold") %>% 
  left_join(dataset.DOheart.protein.debatch.nopoly$lod.peaks$age_int, by = "protein.id") %>% 
  select(protein.id, lod_threshold,lod, qtl.chr, qtl.pos, cis) %>% 
  rename(lod_scan = lod) %>% 
  mutate(signif = lod_scan >= lod_threshold)


summary(lods_proteins$lod_threshold)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.251   7.322   7.622   7.654   7.959  10.071 

pdf("GBRS/QTL_viewer/heart_newversion/Results/permutation/Garys_method/allproteins_permutations_hist.pdf")
hist(lods_proteins$lod_threshold, main = "Histogram of lod threshold for all proteins")
dev.off()


# Checking the p-values for each gene.

protein_pval <- c()

for (i in 1:nrow(lods_proteins)){
  
  protein <- lods_proteins$protein.id[i]
  
  lods_perm <- lods_list[[protein]]
  
  lod_scan <- lods_proteins$lod_scan[i]
  
  protein_pval[i] <- mean(lods_perm >= lod_scan)
  
  names(protein_pval)[i] <- protein
} 

lods_proteins <- lods_proteins %>% 
  mutate(pval = protein_pval)

dim(lods_proteins[lods_proteins$pval < 0.05,]) #[1] 628    8

dim(lods_proteins[lods_proteins$pval < 0.05 & lods_proteins$cis == TRUE,]) #[1] 0 8

# Checking the ones that are in the hotspots 

lods_proteins <- lods_proteins %>% 
  filter(pval < 0.05,
         qtl.chr %in% c(3,5,12))

lods_proteins_filt <- lods_proteins %>% 
  rowwise() %>% 
  mutate(hotspot = if (qtl.chr == 3 & qtl.pos > 145 & qtl.pos < 149) { TRUE }
   else if (qtl.chr == 5 & qtl.pos > 63 & qtl.pos < 67) { TRUE }
   else if (qtl.chr == 12 & qtl.pos > 98 & qtl.pos < 102) { TRUE }
   else { FALSE }) %>% 
  filter(hotspot)
