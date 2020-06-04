###################### RESIDUAL PERMUTATIONS ############################

# This script uses the functions below and rslurm package for running 
# 1000 permutations for each gene/protein on Sumner. 


Get_expression <- function(id, data,rankz=TRUE){
  
data <- get(data)
  # assumes a QTLViewer formatted environment
  
  # fetch the expression data
  if(rankz){expr <- data$data$rz[,id]
  
  } else {
    
    data$data$norm[,id]
  }
  
  # convert to matrix for qtl2
  # note this preserves rownames
  expr <- as.matrix(expr, ncol=1)
  
  expr
}

Scan_expression <- function(data, interactive_term, expr, rankz=TRUE){
  # assumes a QTLViewer formatted environment
  
  data <- get(data)
  
  # define covars
  acov <- data$covar.matrix
  icov <- data$covar.matrix[,interactive_term]
  
  # additive scan
  scan.add <- scan1(genoprobs, expr, addcovar = acov, cores=4)
  
  # interactive scan
  scan.int <- scan1(genoprobs, expr, addcovar = acov,
                    intcovar= icov, cores=4)
  
  # delta
  scan.diff <- scan.int - scan.add
  
  # combine them
  scans <- cbind(scan.add, scan.int, scan.diff)
  
  # add scan type to column names
  colnames(scans) <- c("add", "full", "int")
  
  #return scans
  scans
}


######
# Function to run residual permutation of interaction LOD
# of the four options for methods, "kinship" seems best
Permute_expression <- function(id, data, interactive_term, nperms=1, 
                               peak_type="kinship", rankz=TRUE){
  #######
  # get the trait
  pheno <- Get_expression(id, data)
  # scan the trait
  scan.out <- Scan_expression(data, interactive_term, pheno)

  data_1 <- get(data)
  
  if(peak_type=="kinship"){
    
    K_Q <- calc_kinship(genoprobs, type = "overall")
    
    
    geno <- matrix(rep(1,nrow(data_1$covar.matrix)),ncol=1)
    row.names(geno) <- data_1$annot.samples$mouse.id
    
  }else{
    # pull genotypes at peak of int LOD
    peak_marker <- switch(peak_type,
                          "add" = names(which.max(scan.out[,"add"])),
                          "full" = names(which.max(scan.out[,"full"])),
                          
                          "int" = names(which.max(scan.out[,"int"]))  )
    
    geno <- pull_genoprobpos(genoprobs, peak_marker)
    
    chr <- str_split(peak_marker, pattern="_")[[1]][1]
    K_Q <- K[[chr]]
  }
  
  # define covars
  acov <- data_1$covar.matrix
  icov <- data_1$covar.matrix[,interactive_term]
  
  # get fitted values
  fit <- fit1(geno, pheno, K_Q, addcovar = acov)
  # fit1.scrib$fitted
  
  # get residuals
  resids <- pheno - fit$fitted
  
  # permutation loop
  perm_lod <- rep(0, nperms)
  for(i in 1:nperms){
    # compute new pheno from fit + permuted residuals
    pheno.perm <- fit$fitted + sample(resids)
    
    # scan permuted trait
    scan.perm <- Scan_expression(data, interactive_term, pheno.perm)
    
    ## Grab maximum LOD per permutation
    perm_lod[i] <- max(scan.perm[,"int"])
    print(i)
  }
  
  perm_lod
}
# end function definitions
################################################################################

#Creating arguments dataframe for running permutations for each gene/protein
#using rslurm. Testing with RNA.

# library(rslurm)
# library(qtl2)
# library(tidyverse)
# 
# load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")
# 
# Gene <- dataset.DOheart.mrna$lod.peaks$age_int$gene.id %>% unique
# 
# # For each row (gene/protein) run the permutation's functions using the specific
# # arguments (columns).
# 
# df <- data.frame(id = Gene, data = "dataset.DOheart.mrna",interactive_term = "Age",
#                  nperms = 1000, peak_type = "kinship",rankz = TRUE, 
#                  stringsAsFactors = FALSE)
# 
# rm(Gene)
# 
# ### Testing arguments:
# 
# test <- Get_expression(df[1,1],df[1,2])
# 
# test <- Scan_expression(df[1,2],df[1,3],Get_expression(df[1,1],df[1,2]))
# 
# test <- Permute_expression(
#   data = df[1,2],
#   id = df[1,1],
#   interactive_term = df[1,3],
#   nperms =  df[1,4],
#   peak_type = df[1,5]
#   )
# 
# # All the functions working! Creating scripts for sumner (cluster)
# 
# sjob <- slurm_apply(Permute_expression, df, jobname = 'permutation_apply_gene',
#                     nodes = 9694, cpus_per_node = 4, submit = FALSE)
# 
# # I copied the files to sumner and made some modifications so the scrips 
# # run on singularity.
# 
# ##############################################################################


########################### FOR PROTEIN ########################################

# 
# load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")
# 
# Protein <- dataset.DOheart.protein.debatch.nopoly$lod.peaks$age_int$protein.id %>% unique
# 
# # For each row (gene/protein) run the permutation's functions using the specific
# # arguments (columns).
# 
# df <- data.frame(id = Protein, data = "dataset.DOheart.protein.debatch.nopoly",interactive_term = "Age",
#                  nperms = 1000, peak_type = "kinship",rankz = TRUE,
#                  stringsAsFactors = FALSE)
# 
# rm(Protein)
# 
# ### Testing arguments:
# 
# test <- Get_expression(df[1,1],df[1,2])
# 
# test <- Scan_expression(df[1,2],df[1,3],Get_expression(df[1,1],df[1,2]))
# 
# test <- Permute_expression(
#   data = df[1,2],
#   id = df[1,1],
#   interactive_term = df[1,3],
#   nperms =  df[1,4],
#   peak_type = df[1,5]
#   )
# 
# # All the functions working! Creating scripts for sumner (cluster)
# 
# sjob <- slurm_apply(Permute_expression, df, jobname = 'permutation_apply_protein',
#                     nodes = 2731, cpus_per_node = 4, submit = FALSE)
# 
# # I copied the files to sumner and made some modifications so the scrips
# # run on singularity.
# #
# # ##############################################################################
