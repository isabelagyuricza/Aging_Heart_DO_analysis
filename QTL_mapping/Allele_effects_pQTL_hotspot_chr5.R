#################### Investigating allele effects on pQTL hotspot on chromosome 5 #########################

################################################################################
############ clear workspace

rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

################################################################################
############ libraries

library(gplots)
library(RColorBrewer)
library(corrplot)
library(qtl2)
library(pcaMethods)
library(tidyverse)
library(stringr)
library(broom)
library(grid)
library(gridExtra)

################################################################################
############ load data

load("GBRS/QTL_viewer/heart_newversion/QTLs_data_newversion_corrected.RData")


################################################################################

### for transbands of protein interactive scan for age

### Q x Age interactions on Chr 5

###  extract the lod peaks and start trimming them to just the hotspot

lod.peaks <- as_data_frame(dataset.DOheart.protein.debatch.nopoly$lod.peaks$age_int) %>%
  filter(substr(marker.id,1,2) == "5_") %>%
  separate(marker.id, into=c("qtl.chr","qtl.pos")) %>%
  mutate(qtl.chr=as.integer(qtl.chr), qtl.pos=as.double(qtl.pos)/(10^6)) %>%
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.protein)

# use this graphic to fine tune the edges of the hotspot
#quartz()
ggplot(lod.peaks, aes(x=qtl.pos, y=lod)) +
  geom_point() +
  xlim(63,67) #Region with the highest lod scores ~63Mb to 67Mb

# There are two "peaks" on chromosome 5, but I'm going to stay with the one with 
# higher number of hits for now.
#
lod.peaks <- dplyr::filter(lod.peaks, qtl.pos>63 & qtl.pos<67)
dim(lod.peaks)
# 130 genes

############################
# grab normalized mRNA expr data, regress out Sex (adjusting for sex), then rankz it.

# expression data
expr.data <- dataset.DOheart.protein.debatch.nopoly$data$norm[,lod.peaks$protein.id]

# compute residuals to adjust x wrt covariate
myresids <- function(x, cov){
  residuals(lm(x~cov,na.action=na.exclude))
}

#
expr.data <- apply(expr.data,2,myresids,dataset.DOheart.protein.debatch.nopoly$annot.samples$Sex)

# Apply a normal scores transform to the data.
rankZ <- function(x) {
  x <- rank(x, na.last = "keep") / (length(x) + 1)
  return(qnorm(x))
} # rankZ()
#
expr.data <- apply(expr.data,2,rankZ) #Across mice

rm(myresids, rankZ)

###
# order the data by QTL peak location
lod.peaks <- arrange(lod.peaks, by=qtl.pos)
expr.data <- expr.data[,lod.peaks$protein.id]
###
#
# grab average absolute correlations and filter the lod peaks
tmp.m <- apply(abs(cor(expr.data, use = "pairwise.complete.obs")),2,mean)
sum(tmp.m > 0.3)  #110
sum(tmp.m > 0.35) #104
sum(tmp.m > 0.4)  #94
#
lod.peaks <- mutate(lod.peaks, mean.corr = tmp.m)
#
lod.peaks <- filter(lod.peaks, mean.corr>0.3)
expr.data <- expr.data[,lod.peaks$protein.id]


# hierarchical clustering
quartz()
hm <- heatmap.2(cor(expr.data, use = "pairwise.complete.obs"),
                Rowv = TRUE,Colv=TRUE,dendrogram ="column",
                trace="none",#cexCol = 0.5,
                #colsep = c(64,128),sepcolor = "black",
                keysize = 0.9, key.title = "Color_key",margins = c(8,10),
                col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))


#All the genes mapping on that transband seem to be positively correlated,
# which exception of two genes (Actg2 and Ttn), that are positively correlated to each other,
# but negatively correlated to the rest.

# reorder genes in lod.peaks and expr.data
lod.peaks <- lod.peaks[hm$colInd,]
expr.data <- expr.data[,hm$colInd]

####################################

# compute a PCA summary of the hotspot genes

# Each mouse will have one value for the hotspot (protein expression)
pca.chr5 <- pca(expr.data[,lod.peaks$protein.id], scale="uv", center=TRUE, nPcs=8)
summary(pca.chr5)
# first PC explains 55%

# save PC1 loadings
tmp <- loadings(pca.chr5)
tmp2 <- data_frame(protein.id=row.names(tmp),PC1 = tmp[,"PC1"])
lod.peaks <- left_join(lod.peaks, tmp2)
rm(tmp, tmp2)

# create a datframe for phenotypes and add PCs
pheno <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>%
  dplyr::select(mouse.id, Sex, Age) %>%
  mutate(PC1 = scores(pca.chr5)[,1])

###
# map PC1, to confirm the signal on the hotspot!

scan.additive <- scan1(genoprobs, pheno$PC1, kinship = K,
                       addcovar = dataset.DOheart.protein.debatch.nopoly$covar.matrix)
scan.interaction <- scan1(genoprobs, pheno$PC1, kinship = K,
                          addcovar = dataset.DOheart.protein.debatch.nopoly$covar.matrix,
                          intcovar = dataset.DOheart.protein.debatch.nopoly$covar.matrix[,"Age"])
scan.out <- cbind(scan.additive, scan.interaction, scan.interaction-scan.additive)

which.max(scan.out[,3])

# 5_66518572 
# 18577 

# Using fit_1 to get the allele effects

# 1) Using age 6 as baseline and getting the combined coefficients dataframe

annot.samples <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>% 
  mutate(fAge = as.factor(Age))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit <- fit1(genoprobs$`5`[,,"5_66518572"], pheno$PC1, kinship = K[5],
            addcovar = ac, intcovar = ic,zerosum = FALSE)

# 1) Using age 6 as baseline

fit_6 <- data.frame(estimate = fit$coef[LETTERS[1:8]], 
                    se = fit$SE[LETTERS[1:8]], Age = 6) %>% 
  rownames_to_column("Founder")


# 2) Using age 12 as baseline 

annot.samples$fAge <- factor(annot.samples$fAge, levels = c(12,6,18))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit_2 <- fit1(genoprobs$`5`[,,"5_66518572"], pheno$PC1, kinship = K[5],
              addcovar = ac, intcovar = ic, zerosum = FALSE)

fit_12 <- data.frame(estimate = fit_2$coef[LETTERS[1:8]], 
                     se = fit_2$SE[LETTERS[1:8]], Age = 12) %>% 
  rownames_to_column("Founder")


# 3) Using age 18 as baseline 

annot.samples$fAge <- factor(annot.samples$fAge, levels = c(18,6,12))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit_3 <- fit1(genoprobs$`5`[,,"5_66518572"], pheno$PC1, kinship = K[5],
              addcovar = ac, intcovar = ic, zerosum = FALSE)

fit_18 <- data.frame(estimate = fit_3$coef[LETTERS[1:8]], 
                     se = fit_3$SE[LETTERS[1:8]], Age = 18) %>% 
  rownames_to_column("Founder")


# Gathering them all in the same dataframe.

final_df <- bind_rows(fit_6,fit_12,fit_18)

final_df$Age <- factor(final_df$Age, levels = c(6, 12, 18))

# Comparing with summing coefficients method.


###
# plot the regression coefs

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Changing_baseline_method/Results/Hotspot_chr5/PC1_effects.pdf")

data(CCcolors)
final_df %>% 
  mutate(Founder = factor(Founder, levels = c(LETTERS[1:8]), 
                          labels = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/LtJ", 
                                     "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"))) %>% 
  ggplot(aes(x=Age, y=estimate)) +
  geom_point() +
  geom_line(aes(group=Founder),linetype = "dashed") +
  facet_wrap(~Founder, ncol=4) +
  geom_errorbar(aes(x=Age, ymin=estimate-se, ymax=estimate+se), width=0.2) +
  theme_bw() + ggtitle("PC1") -> p

g <- ggplot_gtable(ggplot_build(p))

stript <- which(grepl('strip-t', g$layout$name))
stript <- c(46, 47, 48, 49, 42, 43, 44, 45)

fills <- CCcolors

k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#quartz()
grid.draw(g)

dev.off()

rm(final_df,sum_df,fit,fit_12,fit_18,fit_2,fit_3,fit_6,ac,ic)

rm(annot.samples,j,e)

################################################################################

############ Testing models for three ER that are positively correlated.

source("GBRS/GC/Scripts/GetGene.R")

dataset.DOheart.mrna$annot.samples <- dataset.DOheart.mrna$annot.samples %>% 
  mutate(fAge = as.factor(Age)) %>% 
  rename(Mouse.ID = mouse.id)

dataset.DOheart.protein.debatch.nopoly$annot.samples <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>% 
  mutate(fAge = as.factor(Age)) %>% 
  rename(Mouse.ID = mouse.id)

expression <- left_join(GetGene(c("Hspa5","Pdia4","Ttn"),dataset.DOheart.mrna$annot.samples,
                                dataset.DOheart.mrna$annot.mrna,
                                dataset.DOheart.mrna$data$norm),
                        GetGene(c("Hspa5","Pdia4","Ttn"),dataset.DOheart.protein.debatch.nopoly$annot.samples,
                                dataset.DOheart.protein.debatch.nopoly$annot.protein,
                                dataset.DOheart.protein.debatch.nopoly$data$norm),
                        by=c("Mouse.ID","Sex","Age"), suffix=c(".mrna",".prot"))


expression <- expression %>% 
  select(-Sex, -Age, -fAge.mrna,-fAge.prot, -Ttn.prot, -Ttn.1, -Ttn.3) %>% 
  column_to_rownames("Mouse.ID")

colnames(expression)[which(colnames(expression)=="Ttn.2")] <- "Ttn.prot"

fit_all <- list()

for (i in colnames(expression)){
  
  data <- expression[,i]
  names(data) <- rownames(expression)
  
  
  # 1) Using age 6 as baseline
  
  annot.samples <- dataset.DOheart.protein.debatch.nopoly$annot.samples
  
  ac <- model.matrix(~fAge+Sex, 
                     data=annot.samples)[,-1]
  
  row.names(ac) <- annot.samples$Mouse.ID
  
  ic <- model.matrix(~fAge, data=annot.samples)[,-1]
  
  row.names(ic) <- annot.samples$Mouse.ID
  
  
  fit <- fit1(genoprobs$`5`[,,"5_66518572"], data, kinship = K[5],
              addcovar = ac, intcovar = ic,zerosum = FALSE)
  
  fit_6 <- data.frame(estimate = fit$coef[LETTERS[1:8]], 
                      se = fit$SE[LETTERS[1:8]], Age = 6) %>% 
    rownames_to_column("Founder")
  
  
  # 2) Using age 12 as baseline 
  
  annot.samples$fAge <- factor(annot.samples$fAge, levels = c(12,6,18))
  
  ac <- model.matrix(~fAge+Sex, 
                     data=annot.samples)[,-1]
  
  row.names(ac) <- annot.samples$Mouse.ID
  
  ic <- model.matrix(~fAge, data=annot.samples)[,-1]
  
  row.names(ic) <- annot.samples$Mouse.ID
  
  fit_2 <- fit1(genoprobs$`5`[,,"5_66518572"], data, kinship = K[5],
                addcovar = ac, intcovar = ic,zerosum = FALSE)
  
  fit_12 <- data.frame(estimate = fit_2$coef[LETTERS[1:8]], 
                       se = fit_2$SE[LETTERS[1:8]], Age = 12) %>% 
    rownames_to_column("Founder")
  
  
  # 3) Using age 18 as baseline 
  
  annot.samples$fAge <- factor(annot.samples$fAge, levels = c(18,6,12))
  
  ac <- model.matrix(~fAge+Sex, 
                     data=annot.samples)[,-1]
  
  row.names(ac) <- annot.samples$Mouse.ID
  
  ic <- model.matrix(~fAge, data=annot.samples)[,-1]
  
  row.names(ic) <- annot.samples$Mouse.ID
  
  fit_3 <- fit1(genoprobs$`5`[,,"5_66518572"],data, kinship = K[5],
                addcovar = ac, intcovar = ic, zerosum = FALSE)
  
  fit_18 <- data.frame(estimate = fit_3$coef[LETTERS[1:8]], 
                       se = fit_3$SE[LETTERS[1:8]], Age = 18) %>% 
    rownames_to_column("Founder")
  
  
  # Gathering them all in the same dataframe.
  
  final_df <- bind_rows(fit_6,fit_12,fit_18)
  
  
  final_df$Age <- factor(final_df$Age, levels = c(6, 12, 18))
  
  
  fit_all[[i]] <- final_df
  
}

Coef <- bind_rows(fit_all, .id = "biotype")

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Changing_baseline_method/Results/Hotspot_chr5/ER_TTN_effects.pdf")

#quartz()
Coef %>% 
  filter(grepl(".prot",biotype)) %>% 
  mutate(Founder = factor(Founder, levels = c(LETTERS[1:8]), 
                          labels = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/LtJ", 
                                     "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")),
         biotype = gsub(".prot","", biotype)) %>% 
  ggplot(aes(x=Age, y=estimate,color=biotype)) +
  geom_point() +
  geom_line(aes(group=biotype)) +
  facet_wrap(~Founder, ncol=4) +
  geom_errorbar(aes(x=Age, ymin=estimate-se, ymax=estimate+se), width=0.2) +
  theme_bw() + ggtitle("ER-TTN proteins expression")

dev.off()


pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Changing_baseline_method/Results/Hotspot_chr5/ER_TTN_effects_2.pdf",
    width = 8,height = 6)

#quartz()
Coef %>% 
  filter(grepl(".prot",biotype)) %>% 
  mutate(Founder = factor(Founder, levels = c(LETTERS[1:8]), 
                          labels = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/LtJ", 
                                     "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ")),
         biotype = gsub(".prot","", biotype)) %>% 
  ggplot(aes(x=Age, y=estimate)) +
  geom_point() +
  geom_line(aes(group=biotype)) +
  facet_grid(biotype ~ Founder, scales = "free")+
  geom_errorbar(aes(x=Age, ymin=estimate-se, ymax=estimate+se), width=0.2) +
  theme_bw() + ggtitle("ER-TTN proteins expression") -> p

g <- ggplot_gtable(ggplot_build(p))

stript <- which(grepl('strip-t', g$layout$name))

fills <- CCcolors

k <- 1
for (i in stript) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#quartz()
grid.draw(g)

dev.off()
