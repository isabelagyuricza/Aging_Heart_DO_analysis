#################################### RNA filter - heart DO #################################

####### This script is to filter RNA counts from GBRS. 

load("GBRS/heart/rna_gbrs_expr_heart.RData")
load("data_setup/annot_samples_mod.RData")

###### Filtering in a way to include Cdk2na gene, 
###### The gene is a marker of senescence and I's associated with aging process.

rna_gbrs_expr_heart[,"ENSMUSG00000044303"] # - Cdk2na gene! Very low expression

mrna.counts <- rna_gbrs_expr_heart

rm(rna_gbrs_expr_heart)

# compute gene-level summary statistics
# for each gene, the median and total read counts
# and the number of samples with at least one read

col.med <- apply(mrna.counts,2,median)
col.tot <- apply(mrna.counts,2,sum)
col.pass192 <- apply(mrna.counts>=1,2,sum)  

# index for subsetting genes
indx.keep <- which(col.pass192 >= 96 & col.med >=1) #Keeping only genes that have more than 1 read for at least half of the samples
                                                    #And a median of at least 1.

mrna.counts_filtered <- mrna.counts[annot_samples$Mouse.ID, indx.keep] #reordering mice to match annots and keeping only the filtered genes
dim(mrna.counts_filtered)
# 192 x 21018
# this is the size of the filtered data

##Check if the Cdk2na is in the dataset! 

mrna.counts_filtered[,"ENSMUSG00000044303"] # - great! It is there :)

###############################################################################

###### Creating the filtered matrix RData with annots

save(mrna.counts_filtered,annot_samples,file = "data_setup/RNA_GBRS/rna_heart_expr.RData")

