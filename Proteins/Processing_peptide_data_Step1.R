####################################################################################################################
#
#   This script reads in the summarized peptides quantifications for the DO
#     heart samples from cross-sectional aging study design

options(stringsAsFactors = FALSE)
library(tidyverse)

######################################################### 
## 
##                    Data inputs
##
########################################################
## Covariate data
do_covariates <- read_csv("ProteinData_raw/Peptide/Data/SampleAnnotations.csv") %>%
  #separate(Sample.Name, c("Batch", "RQ", "TMT")) %>%
  rownames_to_column("mouse_id") %>% 
  dplyr::rename(tag = Tag) %>%
  mutate(tag = gsub(x = tag, pattern = "a", replacement = "n", fixed = TRUE),
         tag = gsub(x = tag, pattern = "b", replacement = "c", fixed = TRUE),
         tag = factor(tag),
         Sex = factor(Sex),
         Age = factor(Age),
         sample_id = paste(Batch, tag, sep = "~")) %>%
  dplyr::select(mouse_id, Sex, Age, tag, Batch, sample_id)

## Raw peptide data
raw_peptide_dat <- read_tsv("ProteinData_raw/Peptide/Data/heart_samples_peptides.tsv")

## Polymorphic peptide data from the CC
polymorphic_peptides <- read.csv("ProteinData_raw/Peptide/Data/peptides_w_polymorphisms.csv", header = TRUE) %>%
  dplyr::rename(peptide_id = X,
                A = aj.e,
                C = s129.e,
                D = nod.e,
                E = nzo.e,
                F = cast.e,
                G = pwk.e,
                H = wsb.e) %>%
  mutate(B = 1,
         peptide_id = gsub(x = peptide_id, pattern = "*", replacement = "", fixed = TRUE)) %>%
  select(peptide_id, LETTERS[1:8]) 
polymorphic_peptides$balance <- rowSums(polymorphic_peptides[,LETTERS[1:8]])

# Note: Reference BL6 is 1 and all the strains that contain peptides different than 1 (0) are considered polymorphic.

######################################################### 
## 
##                    Data processing
##
########################################################
peptide_dat <- raw_peptide_dat %>%
  dplyr::select(ProteinId, PeptideSequence, starts_with("rq"), Class, -Description) %>%
  # Filter out reverse proteins
  filter(!grepl(x = ProteinId, pattern = "##", fixed = TRUE),
         # Filter out contaminants
         !grepl(x = ProteinId, pattern = "contaminant", fixed = TRUE),
         # Filter to only things in Ensembl
         grepl(x = ProteinId, pattern = "ENSMUSP", fixed = TRUE)) %>%
  # Convert TMT-channels from columns to rows
  gather(starts_with("rq"), key = "sample_id", value = "Intensity") %>%
  # Process peptide sequence
  ## Remove first period and pre-string
  mutate(PeptideSequence = gsub(x = PeptideSequence, pattern = "^[-|A-Z]+\\.", replacement = "", perl = TRUE),
         ## Remove second period and post-string
         PeptideSequence = gsub(x = PeptideSequence, pattern = "\\.[A-Z|-]+$", replacement = "", perl = TRUE),
         ## Remove * representing differently charged methyl groups
         PeptideSequence = gsub(x = PeptideSequence, pattern = "*", replacement = "", fixed = TRUE),
         ## Remove -, always at beginning or end of AA sequence
         PeptideSequence = gsub(x = PeptideSequence, pattern = "-", replacement = "", fixed = TRUE),
         ## Remove "." and number from ProteinId
         protein_id = gsub(ProteinId, pattern = "\\.[0-9]*$", replacement = "", perl = TRUE),
         ## Pulling TMT tag
         tag = gsub(x = sample_id, pattern = "rq_", replacement = "", perl = TRUE),
         tag = gsub(x = tag, pattern = "_sn", replacement = "", fixed = TRUE),
         # Switch a to n, and b to c in sample_id to match mouse annotation file
         tag = gsub(x = tag, pattern = "a", replacement = "n", fixed = TRUE),
         tag = gsub(x = tag, pattern = "b", replacement = "c", fixed = TRUE),
         sample_id = paste0(Class, "~", tag)) %>%
  ## Rename PeptideSequence as peptide_id and Class as Batch
  dplyr::rename(peptide_id = PeptideSequence,
                Batch = Class) %>%
  dplyr::mutate(Batch = gsub(Batch, pattern = "HS", replacement = "Set")) %>%
  # Sum up multiple observations of a peptide within an experiment (per Tian)
  group_by(protein_id, peptide_id, Batch, tag, sample_id) %>%
  summarize(Intensity = sum(Intensity, na.rm = TRUE)) %>%
  ungroup

## devtools::install_github("churchill-lab/ensimplR")
library(ensimplR)
gene_annot <- batchGenes(ids = unique(peptide_dat$protein_id))

## Merge into peptide data
peptide_dat <- peptide_dat %>%
  left_join(gene_annot %>%
              dplyr::rename(protein_id = id) %>%
              dplyr::select(protein_id, gene_id, symbol))

## Looking the polymorphic peptides in the data
polymorphic_peptides <- polymorphic_peptides %>%
  filter(peptide_id %in% unique(peptide_dat$peptide_id))

## 7,076 polymorphic peptides out of 347,742 total peptides are being analyzed for mapping

## Merging in polymorphic data
## adding some summaries of polymorphisms
peptide_dat <- peptide_dat %>%
  left_join(polymorphic_peptides %>%
              dplyr::select(peptide_id, balance) %>%
              mutate(polymorphic = TRUE)) %>%
  mutate(balance = ifelse(is.na(balance), 8, balance),
         polymorphic = ifelse(is.na(polymorphic), FALSE, polymorphic))

# Note: Changing balance for 8 means that there is no polymorphisms in any of the strains for the specific peptide.

# Summaries of polymorphic peptides per protein
peptides_per_protein_dat <- peptide_dat %>%
  dplyr::select(protein_id, symbol, peptide_id, balance, polymorphic) %>%
  distinct %>%
  dplyr::count(protein_id)

# Note: There are duplicated rows due to the samples, we just want to know about the proteins so we removed the "samples"  IDs.

polymorphic_peptides_per_protein_dat <- peptide_dat %>%
  dplyr::select(protein_id, symbol, peptide_id, balance, polymorphic) %>%
  distinct %>%
  filter(balance < 8) %>%
  dplyr::count(protein_id)
peptide_to_protein_summaries <- inner_join(peptides_per_protein_dat %>% 
                                             dplyr::rename(num_peptides = n),
                                           polymorphic_peptides_per_protein_dat %>%
                                             dplyr::rename(num_poly_peptides = n)) %>%
  mutate(poly_ratio = num_poly_peptides/num_peptides)



####################################################
##
##            Filter out polymorphic peptides
##
####################################################
## Filtering out polymorphic peptides
peptide_dat_nopoly <- peptide_dat %>%
  filter(!polymorphic)

# dim(peptide_dat_nopoly)
# [1] 4616300      10

## Derive normalization factors
## Ratio of max cumulative intensities for a tag to each cumulative tag intensity 

# Note: the normalization consists in divide the maximum of the sum of intensities (expression) for each batch for the sum of intensities for each mouse.

tmt_max_nopoly <- peptide_dat_nopoly %>%
  group_by(Batch, tag) %>%
  summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
  ungroup %>%
  group_by(Batch) %>%
  summarise(max_sum_inten = max(sum_inten))


tmt_norm_nopoly <- peptide_dat_nopoly %>%
  group_by(sample_id) %>%
  summarise(sum_inten = sum(Intensity, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(Batch = gsub(x = sample_id, pattern = "~.*$", replacement = "", perl = TRUE)) %>%
  mutate(Batch = gsub(x=Batch,pattern = "HS",replacement = "Set",perl = TRUE)) %>% 
  left_join(tmt_max_nopoly) %>%
  mutate(norm_factor = max_sum_inten/sum_inten) %>%
  select(-c(sum_inten, max_sum_inten))

## Adding in zero for unobserved proteins
scaled_protein_dat_nopoly <- peptide_dat_nopoly %>%
  
  ## Gathering information for PROTEIN not peptide!
  
  ## Summing peptides for a protein
  group_by(protein_id, sample_id) %>%
  summarize(Intensity = sum(Intensity)) %>%
  ungroup %>%
  # Remove to add NAs
  complete(protein_id, sample_id) %>%
  # Add back Batch and tag
  mutate(Batch = gsub(x = sample_id, pattern = "~.*$", replacement = "", perl = TRUE),
         tag = gsub(x = sample_id, pattern = "^[A-Z0-9]*~", replacement = "", perl = TRUE),
         Batch = gsub(x = Batch, pattern = "HS", replacement = "Set", perl = TRUE)) %>%
  ## Merge in tag-specific normalizations
  left_join(tmt_norm_nopoly) %>%
  ## Adding in 0 for unobserved proteins
  mutate(Intensity = ifelse(is.na(Intensity), 0, Intensity)) %>%
  # Standardize across samples in a Batch
  mutate(Intensity = Intensity * norm_factor) %>%
  # Remove norm_factor
  select(-norm_factor)

# Take log2
scaled_protein_dat_nopoly <- scaled_protein_dat_nopoly %>%
  mutate(Intensity = log2(Intensity + 1))

## Save the data

scaled_protein_dat_nopoly <- scaled_protein_dat_nopoly %>% 
mutate(sample_id = gsub(x = sample_id, pattern = "HS", replacement = "Set", perl = TRUE))  

# Scaled data

scaled_protein_dat_nopoly <-  scaled_protein_dat_nopoly %>%
  # Filter to the sample used in Chick & Munger
  inner_join(do_covariates %>%
               select(sample_id, mouse_id, Sex, Age))

  write.csv(scaled_protein_dat_nopoly, file = "data_setup/processed_data_GK/scaled_protein_nopoly.csv", row.names = FALSE, quote = FALSE)

