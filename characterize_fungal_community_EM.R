# -----------------------------------------------------------------------------#
# Gather and Format Data for WFDP for EM community 
# Original Author: L. McKinley Nevins 
# February 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.2.1
#                     fundiversity v 1.1.1
#                     vegan v 2.6.10
#                     cluster v 2.1.8
#                     FD v 1.0.12.3
#                     ade4 v 1.7.23
#                     phyloseq v 1.48.0
#                     ape v 5.8.1
#                     agricolae v 1.3.7
#                     rstatix v 0.7.2
#                     ggpubr v 0.6.0
#                     compositions v 2.0.8
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(fundiversity); packageVersion("fundiversity")
library(vegan); packageVersion("vegan")
library(cluster); packageVersion("cluster")
library(FD); packageVersion("FD")
library(ade4); packageVersion("ade4")
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
library(agricolae); packageVersion("agricolae")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(compositions); packageVersion("compositions")
#################################################################################
#                               Main workflow                                   #
#  Subset fungal community data to just EM in WFDP. Merge with exploration type #
#  functional assignments. Perform clr transformation to generate comparable    #
#  relative abundance values.                                                   #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/"
setwd(wd)


# load in original FunGuild object after EM assignments
# this has all of the host taxa for all of the sites 

ps_EM <- readRDS("./Phylogeny_Outputs/EM_funguild_final.RDS")

# Select just WFDP as the site 
WFDP <- subset_samples(ps_EM, Site=="WFDP")

# get otu table 
WFDP_otu <- otu_table(WFDP)

WFDP_otu <- as.data.frame(WFDP_otu)

# save this as a csv file 
write.csv(WFDP_otu, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_otus_EM.csv")

# This was pared down to just the trees that we have both community and trait measurements for 
# which results in 61 trees 

# load back in 
WFDP_otu_sub <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_otus_EM_trait_trees.csv")

#need the samples to be row names
WFDP_otu_sub <- data.frame(WFDP_otu_sub[,-1], row.names=WFDP_otu_sub[,1])

#convert to phyloseq compatible object 
WFDP_otu_sub <- phyloseq::otu_table(as.matrix(WFDP_otu_sub), taxa_are_rows = F)


#####
#preparation of the tax_table from the full WFDP phyloseq object
WFDP_tax <- tax_table(WFDP)

WFDP_tax <- as.data.frame(WFDP_tax)

#convert to phyloseq compatible object 
WFDP_tax <- phyloseq::tax_table(as.matrix(WFDP_tax))


###create subset phyloseq object 

# Load in the TEMPORARY sample data file for WFDP 
env_WFDP <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_data/WFDP_enviro_data_all.csv", row.names = 1)


## create object 
ps_WFDP_final <- phyloseq(otu_table(WFDP_otu_sub), tax_table(WFDP_tax), sample_data(env_WFDP))


#### Trim to reflect just the community in WFDP
# number of taxa - 2,542
ntaxa(ps_WFDP_final)

# number of samples - 60 host trees 
nsamples(ps_WFDP_final)

asv <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

sample_count <- rowSums(asv)

sample_count <- as.data.frame(sample_count) #all trees still have some ASV's

# some of these look pretty tiny, can be revisited later 

asv_count <- colSums(asv)

asv_count <- as.data.frame(asv_count) #a lot of ASVs no longer present in the samples 

#remove taxa that aren't present in this new subset 
ps_WFDP_final <- subset_taxa(ps_WFDP_final, taxa_sums(ps_WFDP_final) > 0)

# number of taxa - now 404
ntaxa(ps_WFDP_final)


# Get taxon table for EM community in WFDP
WFDP_tax <- tax_table(ps_WFDP_final)

WFDP_tax <- as.data.frame(WFDP_tax) %>% rownames_to_column(var = "X") # 404 taxa

# pull out the genus and family levels, as well as the original ASVs. This will be necessary 
# to match them up with the site x environment matrix 

WFDP_tax <- select(WFDP_tax, X, Family, Genus)

WFDP_tax$Genus <- as.factor(WFDP_tax$Genus)

summary(WFDP_tax)

# need to fix name format for the Genus column to get rid of the 'g_'
WFDP_tax$Genus <- gsub("g__", "", WFDP_tax$Genus)


# Load in table of Genera with functional classifications 
exp <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Functional_analyses/em_functions_all.csv")

exp$Genus <- gsub("g__", "", exp$Genus)

exp$Genus <- as.factor(exp$Genus)

summary(exp)

# this has a categorical hydrophobicity column, and then two columns of hydrophylic and hydrophobic that are 
# binary coded for what the trait is 

# The Exploration Types have all been one-hot encoded as separate binary columns - there were 14 categories

# Merge datasets together 
matrix <- merge(WFDP_tax, exp, by = "Genus", all.y = TRUE)

# remove rows that have NA's - these are genera that are not present in WFDP
matrix <- na.omit(matrix)

#get column name set for name merging down below
matrix$ASV <- matrix$X

# the number is correct - 358 ASVs with genus and functional assignments 


#start names file with the original ASVs 
ASV <- phyloseq::taxa_names(ps_WFDP_final)
ASV_long <- as.data.frame(ASV)

#change ASV names to something nicer to work with
taxa_names(ps_WFDP_final)
n_seqs <- seq(ntaxa(ps_WFDP_final))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_WFDP_final) <- paste("ASV", formatC(n_seqs, 
                                                  width = len_n_seqs, 
                                                  flag = "0"), sep = "_")
taxa_names(ps_WFDP_final)

# get shortened names 
ASV2 <- taxa_names(ps_WFDP_final) 
ASV_short <- as.data.frame(ASV2)

# join two dataframes
all_names <- cbind(ASV_long, ASV_short)


#merge taxa_names with the traits file to get the updated ASV names 
matrix <- merge(matrix, all_names, by = "ASV")

# reformat a tiny bit just to get ASVs as the species in rownames, and my traits only 
traits_WFDP_EM <- matrix %>% select(ASV2, hydrophilic, hydrophobic, ET_contact, ET_contact_short, ET_short, ET_contact_medium,
                                    ET_contact_medium_fringe, ET_contact_medium_smooth, ET_medium_smooth, 
                                    ET_medium_fringe, ET_medium_mat, ET_medium_long, ET_medium_long_smooth,
                                    ET_medium_long_fringe, ET_contact_long_smooth, ET_long) %>% column_to_rownames(var = "ASV2") 

#ASVs are now the row names and there are columns for the binary coding of each trait level 

# Save traits file for later analyses 

write.csv(traits_WFDP_EM, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_EM_trait_mat.csv")



# We want these to be the ones that are contained within the final phyloseq object 

taxa_to_keep <- matrix$ASV2

# Prune taxa from the phyloseq object that are NOT in the taxa_to_keep list 
ps_WFDP_final <- prune_taxa(taxa_to_keep, ps_WFDP_final)


#SAVE THE FINAL RAW DATA PHYLOSEQ OBJECT 
# For the taxa that have functional assignments in the community 
saveRDS(ps_WFDP_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_final.RDS")


# Get to know the phyloseq data ####

# number of taxa - 358
ntaxa(ps_WFDP_final)

# number of samples - 60
nsamples(ps_WFDP_final)

# sample names
sample_names(ps_WFDP_final)
rank_names(ps_WFDP_final)
# taxa names
taxa_names(ps_WFDP_final)

# ASV table
otu_table(ps_WFDP_final) %>% View()

# how many sequences observed in each sample?
seq_counts <- otu_table(ps_WFDP_final) %>% rowSums() %>% as.data.frame()
# No trees with no sequences 


# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_WFDP_final) %>% colSums()


# how many different samples was each taxon found in?
asv <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

#inspect tax_table 
WFDP_tax_table <- as.data.frame(tax_table(ps_WFDP_final))

#count number of classifications in each column to determine coverage to taxonomic levels 

colSums(!is.na(WFDP_tax_table))

#  Kingdom  Phylum    Class    Order   Family   Genus    Species 
#  358       358       358     358     358      358      187 
#  100%     100%      100%     100%    100%     100%     52.2% 



#################################################################################

########################### -- 
# (2) DATA TRANSFORMATION
########################### -- 

# Pursuing a transformation because sequence data are inherently compositional, but standard rarefaction 
# is not repeatable and leads to the exclusion of a huge fraction of your data 

# trying a centered log-ratio transformation per
# https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full

# Also: Van den Boogaart, K. G., and Tolosana-Delgado, R. (2013). Analyzing Compositional Data with R, 
# London, UK: Springer.

## Using robust clr transformation, which only transforms the non-zero values in the dataset, and is thus 
# less affected by high zero counts than standard clr is 

asv <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame()

clr_WFDP_EM <- decostand(asv, method = "rclr")

clr_WFDP_EM <- as.data.frame(clr_WFDP_EM)


# load the clr data back into a phyloseq object 
ps_WFDP_final <- phyloseq(tax_table(tax_table(ps_WFDP_final)),
                      otu_table(otu_table(clr_WFDP_EM, taxa_are_rows=FALSE)),
                      sample_data(sample_data(ps_WFDP_final)))


#SAVE THE FINAL PHYLOSEQ OBJECT 
saveRDS(ps_WFDP_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_clr_final.RDS")

## This now only contains the fungal taxa that we have functional information for 

# -- #### END ##### -- 
