################################################# DO heart protein differential expression analysis

library(tidyverse)
library(stringr)
library(broom)

############################################### Step 1: loading/setting up the data #########################################################

#Sourcing all data from QTLViewer file

load("GBRS/QTL_viewer/QTLv_data_str_newversion.RData")


#Taking Tag and Generation from the raw data (The QTL RData object contains only Sex and Age)
do_covariates <- read_csv("ProteinData_raw/Peptide/Data/SampleAnnotations.csv") %>%
  dplyr::rename(tag = Tag) %>%
  mutate(tag = gsub(x = tag, pattern = "a", replacement = "n", fixed = TRUE),
         tag = gsub(x = tag, pattern = "b", replacement = "c", fixed = TRUE),
         tag = factor(tag),
         Sex = factor(Sex),
         Age = factor(Age),
         sample_id = paste(Batch, tag, sep = "~")) %>%
  dplyr::select(Sex, Age, tag, Batch, Generation,sample_id,Mouse.ID)

do_covariates <- do_covariates %>% 
  mutate(Mouse.ID = gsub("-",".",Mouse.ID)) %>% 
  rename(mouse.id = Mouse.ID)

do_covariates <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>% 
  select(mouse.id) %>% 
  left_join(do_covariates, by = "mouse.id") %>% 
  rename(Mouse.ID = mouse.id)

identical(do_covariates$Mouse.ID,as.character(dataset.DOheart.protein.debatch.nopoly$annot.samples$mouse.id))

# Sample annotations - a tibble 190 x 5
# standardize Age so the change will be in units of years

annot.sample <- cbind(dataset.DOheart.protein.debatch.nopoly$annot.samples,
                      Tag=do_covariates$tag,Generation=do_covariates$Generation) %>% 
  as_tibble() %>% 
  mutate(Mouse.ID=as.character(mouse.id),Tag=as.factor(Tag),Generation=as.factor(Generation)) %>% 
  select(-mouse.id)


# protein data (Batch adjusted and non-polymorphic)

protein.data <- dataset.DOheart.protein.debatch.nopoly$data$norm

dim(protein.data)
#[1]  190 4062

# protein annotations
annot.protein <- dataset.DOheart.protein.debatch.nopoly$annot.protein
# 4062 x 9 tibble

####
# clean up
rm(dataset.DOheart.mrna, dataset.DOheart.protein.debatch.nopoly, ensembl.version,
   genoprobs, K, map, markers)

save(protein.data,annot.protein,annot.sample,file = "GBRS/GC/Data/Protein.RData")

############################################# Step 2: ANOVA computation ###############################################

# To see differences by unit of years

annot.sample <- annot.sample %>% 
  mutate(Age=(Age/12)-1)

# Checking if the datasets are in the same order. Just continue if that is TRUE!
identical(rownames(protein.data),annot.sample$Mouse.ID)
#[1] TRUE

# Using PCA to test for Tag and Generation effect

df_comput <- protein.data %>% 
  as_tibble() %>% 
  dplyr::mutate(Sex=as.factor(annot.sample$Sex), Age=as.numeric(annot.sample$Age),Tag=as.factor(annot.sample$Tag),Gen=as.factor(annot.sample$Generation))


aggregate(df_comput$Age,by=list(df_comput$Tag),table) 

#      Group.1 x.-0.5 x.0 x.0.5
# 1      126      4  10     5
# 2     127c      7   7     5
# 3     127n      6   6     7
# 4     128c      7   9     3
# 5     128n      6   9     4
# 6     129c      9   3     7
# 7     129n      4   8     7
# 8     130c     10   2     7
# 9     130n      6   4     9
# 10     131      5   4    10  It doesn't seems to be a coffound effect between age and tag like in the kidney dataset

aggregate(df_comput$Sex,by=list(df_comput$Tag),table) #It doesn't seem to be coffounded with Sex either

#     Group.1 x.F x.M
# 1      126  11   8
# 2     127c  11   8
# 3     127n   8  11
# 4     128c   8  11
# 5     128n  10   9
# 6     129c  10   9
# 7     129n  12   7
# 8     130c   8  11
# 9     130n  11   8
# 10     131   5  14

pca <- prcomp(df_comput[,grep(pattern = "ENSMUSP",colnames(df_comput))])

# 1) For tag

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=df_comput$Tag)

#quartz()
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3)

# 2) For Generation


d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=df_comput$Gen)

#quartz()
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3)

# Creating anova function

anova_comput <- function (data)
  
{
  
  #Creating a df (df_comput) containing the expression data + annotations. This is why is important the annot.sample to be in the same order than protein.data
  df_comput <- data %>% 
    as_tibble() %>% 
    dplyr::mutate(Sex=as.factor(annot.sample$Sex), Age=as.numeric(annot.sample$Age),Tag=as.factor(annot.sample$Tag),Gen=as.factor(annot.sample$Generation))
  
  df_value <- df_comput %>% dplyr::select(contains("ENSMUS"))
  
  add <- function(x){
    
    fit_additive <- summary(lm(x ~ Sex+Age, data=df_comput))$coefficients
    
    df_final_additive <- tibble(pvalue.Age = fit_additive["Age","Pr(>|t|)"],
                                estimate.Age = fit_additive["Age","Estimate"],
                                se.Age = fit_additive["Age","Std. Error"],
                                pvalue.Sex = fit_additive["SexM","Pr(>|t|)"],
                                estimate.Sex = fit_additive["SexM","Estimate"],
                                se.Sex = fit_additive["SexM","Std. Error"])
    
    return(df_final_additive)}
  
  result_add <- apply(df_value, 2, add)
  
  int <- function(x){
    
    fit_interactive <- summary(lm(x ~ Sex+Age+Sex:Age, data=df_comput))$coefficients
    
    df_final_interactive <- tibble(pvalue.Int = fit_interactive["SexM:Age","Pr(>|t|)"],
                                   estimate.Int = fit_interactive["SexM:Age","Estimate"],
                                   se.Int = fit_interactive["SexM:Age","Std. Error"])
    
    return(df_final_interactive)}
  
  result_int <- apply(df_value, 2, int)
  
  result <- list(result_add,result_int)
  
  names(result) <- c("additive","interactive")
  
  return(result)
}

comput_test_withoutrandomv <- anova_comput(protein.data)


df_additive_wth <- comput_test_withoutrandomv$additive %>% 
  bind_rows(.id = "Protein_ID")

df_interactive_wth <- comput_test_withoutrandomv$interactive %>% 
  bind_rows(.id = "Protein_ID")

df_final_wth <- df_interactive_wth %>% 
  dplyr::select(Protein_ID,pvalue.Int,estimate.Int,se.Int) %>% 
  left_join(df_additive_wth, by="Protein_ID") %>% 
  dplyr::left_join(annot.protein,by="Protein_ID") %>% 
  dplyr::select(Protein_ID,Gene_ID,pvalue.Int,estimate.Int,se.Int,pvalue.Age,estimate.Age,se.Age,pvalue.Sex,estimate.Sex,se.Sex)


pdf("Protein_enrich_analysis/Protein_GK/Results/hist_pvalues_lm.pdf",width = 12)
#quartz()
df_final_wth %>% 
  dplyr::select(starts_with("pvalue")) %>%
  gather(key=test.name, value=pvalue,pvalue.Int:pvalue.Sex) %>%
  ggplot(aes(x=pvalue)) +
  geom_histogram(bins = 100) +
  facet_wrap(~test.name)

dev.off()

############################################# Step 3: FDR adjustment ###############################################


df_final_wth <- df_final_wth %>% 
  mutate(padj.Int = p.adjust(pvalue.Int),
         padj.Age = p.adjust(pvalue.Age),
         padj.Sex = p.adjust(pvalue.Sex)) %>% 
  dplyr::select(Protein_ID,Gene_ID,pvalue.Age,padj.Age,estimate.Age,se.Age,pvalue.Sex,padj.Sex,estimate.Sex,se.Sex,pvalue.Int,padj.Int,estimate.Int,se.Int)

# count significant gene list sizes
apply(dplyr::select(df_final_wth, starts_with("padj"))<0.001, 2, sum)

#              padj.Age padj.Sex padj.Int 
#padj < 0.001   810       79        5 
#padj < 0.01    985      105       22
#padj < 0.1    1246      168       51


save(df_final_wth,file = "Protein_enrich_analysis/Protein_GK/Results/result_anova_protein.RData")

