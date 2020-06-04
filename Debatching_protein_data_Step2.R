############ Debatching protein heart data #############

########### Using linear mixed model - Greg Keele's script


#################################################################################################################


######### Isabela Gerdes Gyuricza - 11-25-2019

library(tidyverse)
library(lme4)
library(tidyverse)


do_scaled_protein_nopoly_dat <- read.csv("data_setup/processed_data_GK/scaled_protein_nopoly.csv", header = TRUE)

####### Using PCA to check for batch effects

#Checking for batch effect

df <- do_scaled_protein_nopoly_dat %>% 
  spread(key = protein_id, value = Intensity)

pca <- prcomp(df[,7:500])

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=df$Batch)

#quartz()
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3)


#Checking for tag effect

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=df$tag)

#quartz()
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3)

# It's ok!

#### Function to regress out batch (and possibly sex)
debatch_protein_lmer <- function(long_dat, 
                                 regress_out_sex = FALSE) {
  
  wide_dat <- long_dat %>%
    dplyr::select(protein_id, Intensity, Sex, Age, mouse_id, Batch, tag) %>%
    spread(key = protein_id, value = Intensity)
  
  proteins <- grep(names(wide_dat), pattern = "ENSMUS", value = TRUE)
  
  resid_dat <- matrix(NA, nrow = nrow(wide_dat), ncol = length(proteins))
  rownames(resid_dat) <- wide_dat$mouse_id
  colnames(resid_dat) <- proteins
  
  for (i in 1:length(proteins)) {
    lmer_fit <- suppressWarnings(lmer(formula(paste(proteins[i], "~ Sex + Age + (1 | Batch)")), data = wide_dat))
    
    resid_dat[,i] <- wide_dat[,proteins[i]] - ranef(lmer_fit)$Batch[wide_dat$Batch, 1]
    if (regress_out_sex) {
      resid_dat[,i] <- resid_dat[,proteins[i]] - fixef(lmer_fit)["SexM"] * model.matrix(~Sex, data = wide_dat)[, -1, drop = FALSE]
    }
  }
  
  resid_dat <- resid_dat %>% 
    as.data.frame %>%
    rownames_to_column("mouse_id") %>%
    left_join(wide_dat %>%
                dplyr::select(mouse_id, Sex, Age, Batch, tag) %>%
                mutate(mouse_id = as.character(mouse_id)))
  resid_dat
}


debatch <- debatch_protein_lmer(do_scaled_protein_nopoly_dat)

####### Using PCA to check for batch effects

pca <- prcomp(debatch[,2:500])

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=debatch$Batch)

#quartz()
ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3)

# Data is debatched!

 debatch <- debatch %>% 
  dplyr::select(mouse_id, Sex, Age, contains("ENSMUS"),tag) %>%
  gather(key = protein_id, value = Intensity, -c(mouse_id, Sex, Age, tag))

write.csv(debatch, file = "data_setup/processed_data_GK/Results/debatch_nopoly_protein_heart.csv", quote = FALSE, row.names = FALSE)
