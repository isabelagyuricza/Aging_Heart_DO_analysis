############ Gathering counts from all the samples for heart

######################################################################

library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)
library("gplots")
library("RColorBrewer")
library("genefilter")

#################################### Heart data

file_list <- list.files("GBRS/heart/expected_counts_heart", pattern="*.txt") # create list of all .gene.counts.txt files in folder

# read in each .txt file in file_list and create a data frame with the same name as the file

list_counts <- list()
for (i in file_list){
  x <- read.table(paste("GBRS/heart/expected_counts_heart/", i, sep=""),header=TRUE,sep = "\t")
  colnames(x) <- c("Target_ID","A","B","C","D","E","F","G","H","Total","Notes")
  list_counts[[i]] <- x
}; rm(x,i)

#Checking if all the genes for all the samples are the same

ls <- list()
for (i in 1:length(list_counts)){
  ls <- list_counts[[i]]$Target_ID
}
#Factor with 47643 levels = taking all the genes for all the samples there are 47643 levels (genes). There is no difference of the set of the genes among samples.

######################## Creating input data with ALL the genes (list_counts)

df_counts <- data.frame(row.names=list_counts[[1]]$Target_ID,
                        stringsAsFactors = FALSE) #I took the genes from the first sample, because they are all the same.

for (i in names(list_counts)){
  x <- str_split(i,pattern = "-",simplify = TRUE)
  x <- x[1]
  df_counts[[x]] <- list_counts[[i]]$Total
  
}

##################### Changing the colnames

for (i in 1:length(colnames(df_counts))) {
  x <- str_split(colnames(df_counts)[i],"_",simplify = TRUE)
  x <- x[2]
  x <- as.character(x)
  x <- as.numeric(x)
  x <- formatC(x,width = 4,format = "d",flag = "0")
  colnames(df_counts)[i] <- x
  
}

for (i in 1:length(colnames(df_counts))){
  y <- colnames(df_counts)[i]
  x <- paste0("DO.",y)
  colnames(df_counts)[i] <- x
  
}

rna_gbrs_expr_heart <- t(df_counts)

save(rna_gbrs_expr_heart,file = "GBRS/heart/rna_gbrs_expr_heart.RData")


rm(list = ls())

