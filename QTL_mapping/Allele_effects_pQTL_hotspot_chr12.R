#################### Investigating allele effects on pQTL hotspot on chromosome 12 #########################

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

### Q x Age interactions on Chr 12

#####################
###  extract the lod peaks and start trimming them to just the hotspot

lod.peaks <- as_data_frame(dataset.DOheart.protein.debatch.nopoly$lod.peaks$age_int) %>%
  filter(substr(marker.id,1,3) == "12_") %>%
  separate(marker.id, into=c("qtl.chr","qtl.pos")) %>%
  mutate(qtl.chr=as.integer(qtl.chr), qtl.pos=as.double(qtl.pos)/(10^6)) %>%
  left_join(dataset.DOheart.protein.debatch.nopoly$annot.protein)

# use this graphic to fine tune the edges of the hotspot
quartz()
ggplot(lod.peaks, aes(x=qtl.pos, y=lod)) +
  geom_point() +
  xlim(98,102) #Region with the highest lod scores ~98Mb to 102Mb

lod.peaks <- dplyr::filter(lod.peaks, qtl.pos>98 & qtl.pos<102)
dim(lod.peaks)
# 177 genes

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
  x <- rank(x, na.last = "keep") / (length(x) - sum(is.na(x)) + 1)
  return(qnorm(x))
}

expr.data <- apply(expr.data,2,rankZ) #Across mice

rm(myresids, rankZ)

###
# order the data by QTL peak location
lod.peaks <- arrange(lod.peaks, by=qtl.pos)
expr.data <- expr.data[,lod.peaks$protein.id]
###

############################
# look at correlations

# grab average absolute correlations and filter the lod peaks
tmp.m <- apply(abs(cor(expr.data, use = "pairwise.complete.obs")),2,mean)
sum(tmp.m > 0.3)  #161
sum(tmp.m > 0.35) #152
sum(tmp.m > 0.4)  #136
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


## There are a group of genes that are positvely correlated to each other, and there 
## are some of them that also mapped to the hotspot on chromosome 3 (myh, col, and histones)
## I'm going to check if they are the same isoforms. 
## In adittion to that, there are other extracellular matrix components, such as Fbn1, 
## Hspg2 adnd other isoforms of collagens that are positively correlated. 

# reorder genes in lod.peaks and expr.data
lod.peaks <- lod.peaks[hm$colInd,]
expr.data <- expr.data[,hm$colInd]



####################################

# compute a PCA summary of the hotspot genes

# Each mouse will have one value for the hotspot (protein expression)
pca.chr12 <- pca(expr.data[,lod.peaks$protein.id], scale="uv", center=TRUE, nPcs=8)
summary(pca.chr12)
# first PC explains 51%

# save PC1 loadings
tmp <- loadings(pca.chr12)
tmp2 <- data_frame(protein.id=row.names(tmp),PC1 = tmp[,"PC1"])
lod.peaks <- left_join(lod.peaks, tmp2)
rm(tmp, tmp2)

# create a datframe for phenotypes and add PCs
pheno <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>%
  dplyr::select(mouse.id, Sex, Age) %>%
  mutate(PC1 = scores(pca.chr12)[,1])

###
# map PC1, to confirm the signal on the hotspot!

scan.additive <- scan1(genoprobs, pheno$PC1, kinship = K,
                       addcovar = dataset.DOheart.protein.debatch.nopoly$covar.matrix)
scan.interaction <- scan1(genoprobs, pheno$PC1, kinship = K,
                          addcovar = dataset.DOheart.protein.debatch.nopoly$covar.matrix,
                          intcovar = dataset.DOheart.protein.debatch.nopoly$covar.matrix[,"Age"])
scan.out <- cbind(scan.additive, scan.interaction, scan.interaction-scan.additive)
colnames(scan.out) <- c("add","full","diff")
rm(scan.additive, scan.interaction)

which.max(scan.out[,3])

# 12_99780940 
# 44951

# Using fit_1 to get the allele effects

# 1) Using age 6 as baseline and getting the combined coefficients dataframe

annot.samples <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>% 
  mutate(fAge = as.factor(Age))

ac <- model.matrix(~fAge+Sex, 
                   data=annot.samples)[,-1]

row.names(ac) <- annot.samples$mouse.id

ic <- model.matrix(~fAge, data=annot.samples)[,-1]

row.names(ic) <- annot.samples$mouse.id

fit <- fit1(genoprobs$`12`[,,"12_99780940"], pheno$PC1, kinship = K[12],
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

fit_2 <- fit1(genoprobs$`12`[,,"12_99780940"], pheno$PC1, kinship = K[12],
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

fit_3 <- fit1(genoprobs$`12`[,,"12_99780940"], pheno$PC1, kinship = K[12],
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

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Changing_baseline_method/Results/Hotspot_chr12/PC1_effects.pdf")

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

############ Testing models for three histones proteins that are correlated.

source("GBRS/GC/Scripts/GetGene.R")

dataset.DOheart.mrna$annot.samples <- dataset.DOheart.mrna$annot.samples %>% 
  mutate(fAge = as.factor(Age)) %>% 
  rename(Mouse.ID = mouse.id)

dataset.DOheart.protein.debatch.nopoly$annot.samples <- dataset.DOheart.protein.debatch.nopoly$annot.samples %>% 
  mutate(fAge = as.factor(Age)) %>% 
  rename(Mouse.ID = mouse.id)

expression <- left_join(GetGene(c("H2afy","H2afy2","Npm1"),dataset.DOheart.mrna$annot.samples,
                                dataset.DOheart.mrna$annot.mrna,
                                dataset.DOheart.mrna$data$norm),
                        GetGene(c("H2afy","H2afy2","Npm1"),dataset.DOheart.protein.debatch.nopoly$annot.samples,
                                dataset.DOheart.protein.debatch.nopoly$annot.protein,
                                dataset.DOheart.protein.debatch.nopoly$data$norm),
                        by=c("Mouse.ID","Sex","Age"), suffix=c(".mrna",".prot"))


expression <- expression %>% 
  select(-Sex, -Age, -fAge.mrna,-fAge.prot) %>% 
  column_to_rownames("Mouse.ID")

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
  
  
  fit <- fit1(genoprobs$`12`[,,"12_99780940"], data, kinship = K[12],
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
  
  fit_2 <- fit1(genoprobs$`12`[,,"12_99780940"], data, kinship = K[12],
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
  
  fit_3 <- fit1(genoprobs$`12`[,,"12_99780940"],data, kinship = K[12],
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

pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Changing_baseline_method/Results/Hotspot_chr12/Hist_effects.pdf")

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
  theme_bw() + ggtitle("Nucleosome proteins expression")

dev.off()


pdf("GBRS/QTL_viewer/heart_newversion/Results/Transbands/Drivers_allele_effects/Changing_baseline_method/Results/Hotspot_chr12/Hist_effects_2.pdf",
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
  theme_bw() + ggtitle("Nucleosome proteins expression") -> p

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
