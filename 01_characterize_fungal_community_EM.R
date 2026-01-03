# -----------------------------------------------------------------------------#
# Gather and Format Data for WFDP for EM community 
# Original Author: L. McKinley Nevins 
# January 1, 2026
# Software versions:  R v 4.5.2
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.3.0
#                     vegan v 2.7.2
#                     cluster v 2.1.8.1
#                     ade4 v 1.7.23
#                     phyloseq v 1.54.0
#                     ape v 5.8.1
#                     agricolae v 1.3.7
#                     rstatix v 0.7.3
#                     ggpubr v 0.6.2
#                     compositions v 2.0.9
#                     stringr v 1.6.0
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(vegan); packageVersion("vegan")
library(cluster); packageVersion("cluster")
library(ade4); packageVersion("ade4")
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
library(agricolae); packageVersion("agricolae")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(compositions); packageVersion("compositions")
library(stringr); packageVersion("stringr")

#################################################################################
#                               Main workflow                                   #
#  Subset fungal community data to just EM in WFDP. Merge with exploration type #
#  functional assignments. Perform clr transformation to generate comparable    #
#  relative abundance values.                                                   #
#                                                                               #
#################################################################################

############### --- 
# (1) DATA PREP
############### --- 

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/"
setwd(wd)


# load in the phyloseq that has had guild assignments made using both FungalTraits and FunGuild
ps_EM <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_phyloseq_final_2025.RDS")


# Select just WFDP as the site 
WFDP <- subset_samples(ps_EM, Site=="WFDP")


# get otu table 
WFDP_otu <- otu_table(WFDP)

WFDP_otu <- as.data.frame(WFDP_otu)
# Contains 85 trees because these are the total number of trees in WFDP that we have fungal community data for 

# save this as a csv file 
write.csv(WFDP_otu, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_otus_EM.csv")

# This was pared down to just the EM or DUAL trees that we have both community and trait measurements
# which results in 40 trees 

# load back in 
WFDP_otu_sub <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_otus_EM_trait_trees.csv")

#need the samples to be row names
WFDP_otu_sub <- data.frame(WFDP_otu_sub[,-1], row.names=WFDP_otu_sub[,1])

#convert to phyloseq compatible object 
WFDP_otu_sub <- phyloseq::otu_table(as.matrix(WFDP_otu_sub), taxa_are_rows = F)


# Load in table of Genera with functional classifications 
## This should align with the genera present in WFDP without any additional work 
exp <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/em_functions_2025.csv")

exp$Genus <- as.factor(exp$Genus)

# this has a categorical hydrophobicity column, and then two columns of hydrophylic and hydrophobic that are 
# binary coded for what the trait is 

# The Exploration Types have all been one-hot encoded as separate binary columns - there were 14 categories


# Load in the TEMPORARY sample data file for WFDP 
# This contains all studied trees, not just the EM and DUAL hosts 
env_WFDP <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_data/WFDP_enviro_data_all.csv", row.names = 1)

#####
#preparation of the tax_table from the full WFDP phyloseq object
WFDP_tax <- tax_table(WFDP)

WFDP_tax <- as.data.frame(WFDP_tax)

#convert to phyloseq compatible object 
WFDP_tax <- phyloseq::tax_table(as.matrix(WFDP_tax))
# 1,870 ASVs


#################################################################################

####################################### --- 
# (2) CREATE AND SAVE PHYLOSEQ OBJECT
####################################### --- 

## create object 
ps_WFDP_final <- phyloseq(otu_table(WFDP_otu_sub), tax_table(WFDP_tax), sample_data(env_WFDP))


#### Trim to reflect just the community in WFDP
# number of taxa - 1,870 ASVs - not all actually present, just haven't been trimmed yet 
ntaxa(ps_WFDP_final)

# number of samples - 40 host trees that are either EM or DUAL
nsamples(ps_WFDP_final)

asv <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

sample_count <- rowSums(asv)

sample_count <- as.data.frame(sample_count) #all trees still have some reads

asv_count <- colSums(asv)

asv_count <- as.data.frame(asv_count) #a lot of ASVs no longer present in the samples 

#remove taxa that aren't present in this new subset 
ps_WFDP_final <- subset_taxa(ps_WFDP_final, taxa_sums(ps_WFDP_final) > 0)

# number of taxa - now 289 ASVs
ntaxa(ps_WFDP_final)


# Get taxon table for EM community in WFDP
WFDP_tax <- tax_table(ps_WFDP_final)

WFDP_tax <- as.data.frame(WFDP_tax) %>% rownames_to_column(var = "X") # 289 taxa


# Save as the file to assess the total ASV's and figure out proportions that have trait assignments, etc. 

write.csv(WFDP_tax, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/all_EM_tax.csv")


# pull out the genus and family levels, as well as the original ASVs. This will be necessary 
# to match them up with the site x environment matrix 

WFDP_tax <- dplyr::select(WFDP_tax, X, Family, Genus)

WFDP_tax$Genus <- as.factor(WFDP_tax$Genus)

summary(WFDP_tax)

# Merge datasets together 
matrix <- merge(WFDP_tax, exp, by = "Genus", all.y = TRUE)

# remove rows that have NA's - these are genera that are not present in WFDP
matrix <- na.omit(matrix)

#get column name set for name merging down below
matrix$ASV <- matrix$X


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

# save matrix file as key to the long and short ASV assignments 
write.csv(matrix, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/all_ASV_names_WFDP_EM_2025.csv")


# reformat a tiny bit just to get ASVs as the species in rownames, and my traits only 
traits_WFDP_EM <- matrix %>% dplyr::select(ASV2, hydrophilic, hydrophobic, ET_contact, ET_contact_short, ET_short, ET_contact_medium,
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
# Now contains the 260 ASVs that have functional assignments in the community of EM and DUAL tree hosts 
saveRDS(ps_WFDP_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_EM_final.RDS")

#################################################################################

# Get to know the phyloseq data ####

# number of taxa - 260
ntaxa(ps_WFDP_final)

# number of samples - 40
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
# No trees with no reads


# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_WFDP_final) %>% colSums()


# how many different samples was each taxon found in?
asv <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

#inspect tax_table 
WFDP_tax_table <- as.data.frame(tax_table(ps_WFDP_final))

#count number of classifications in each column to determine coverage to taxonomic levels 

colSums(!is.na(WFDP_tax_table))

#  Kingdom  Phylum    Class    Order   Family   Genus    Species 
#  260      260       260      260     260      260      144 
#  100%     100%      100%     100%    100%     100%     55.4% 


#################################################################################

########################### -- 
# (3) DATA TRANSFORMATION
########################### -- 

# Pursuing a transformation because sequence data are inherently compositional, but standard rarefaction 
# is not repeatable and leads to the exclusion of a huge fraction of your data 

# Doing a centered log-ratio transformation per
# https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full

# Also: Van den Boogaart, K. G., and Tolosana-Delgado, R. (2013). Analyzing Compositional Data with R, 
# London, UK: Springer.

## Using clr transformation with an added pseudocount, to 

asv <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame()

clr_WFDP_EM <- decostand(asv, method = "clr", pseudocount = 1e-06)

clr_WFDP_EM <- as.data.frame(clr_WFDP_EM)


# load the clr data back into a phyloseq object 
ps_WFDP_final <- phyloseq(tax_table(tax_table(ps_WFDP_final)),
                      otu_table(otu_table(clr_WFDP_EM, taxa_are_rows=FALSE)),
                      sample_data(sample_data(ps_WFDP_final)))


#SAVE THE FINAL PHYLOSEQ OBJECT 
saveRDS(ps_WFDP_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_clr_EM_final.RDS")

## This now only contains the fungal taxa that we have functional information for 

#################################################################################

##################################### --- 
# (4) VISUALIZE COMMUNITY VARIATION 
##################################### --- 

## Look at ASV taxonomic richness across the samples 
# Doesn't factor in relative abundance, so using the raw phyloseq object. This is just showing the presence or 
# absence of an ASV


# Load in the raw reads for ASVs
ps_WFDP_raw <- readRDS("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_final.RDS")

# Look at the spread of ASVs present in the samples 
tax <- tax_table(ps_WFDP_raw) %>% as.data.frame()

#move the rownames into a column 
tax <- rownames_to_column(tax, var = "ASV")

# Make genus a factor 
tax$Genus <- as.factor(tax$Genus)
summary(tax$Genus)


# Grab OTU table 
otu <- otu_table(ps_WFDP_raw) %>% as.data.frame()

otu_long <- otu %>%
  rownames_to_column("Tree_ID") %>%
  pivot_longer(
    cols = -Tree_ID,
    names_to = "ASV",
    values_to = "count"
  ) %>%
  filter(count > 0)


# Merge otu and tax tables together 
otu_tax <- merge(otu_long, tax, by = "ASV")


# Calculate genus richness for each tree 
genus_richness <- otu_tax %>%
  group_by(Tree_ID, Genus) %>%
  summarise(
    n_asvs = n_distinct(ASV),
    .groups = "drop"
  )


# Calculate relative ASV richness to account for the fact that the trees have different numbers of ASVS
genus_rel <- genus_richness %>%
  group_by(Tree_ID) %>%
  mutate(rel_asvs = n_asvs / sum(n_asvs)) %>%
  ungroup()


# add column of Host_ID
genus_rel <- genus_rel %>%
  mutate(
    Host_ID = str_extract(Tree_ID, "(?<=-)[A-Z]{4}(?=-)")
  )

# add column to shorten the tree_id for plotting 
genus_rel <- genus_rel %>%
  mutate(
    tree_short = str_extract(Tree_ID, "[A-Z]{4}-\\d{2}")
  )





# Visualize 
richness_plot_EM <- ggplot(genus_rel,
                           aes(x = tree_short, y = rel_asvs, fill = Genus)) +
  geom_col(width = 0.9) +
  facet_wrap(~ Host_ID, scales = "free_x") +
  labs(
    x = "",
    y = "Relative EM fungal ASV richness",
    fill = "Genus"
  ) +
  theme_bw() +
theme(legend.position = "bottom")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) + 
  theme(
    axis.text.x = element_text(size = 0, colour="black", angle = 90, hjust = 1),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 11, colour="black"), 
    panel.spacing = unit(1, "lines"))


richness_plot_EM


## Beta-diversity can be calculated using Aitchison Distance, which is essentially the euclidean
# distance calculated between pairs of samples that have been transformed by CLR

# calculate Aitchison distance using dist() from base R 

aitchison_WFDP_EM <- dist(clr_WFDP_EM, method = "euclidean")

peek <- as.matrix(aitchison_WFDP_EM) %>% as.data.frame()

# save Aitchison Distance matrix 
save(aitchison_WFDP_EM, file="~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/aitchison_dist_EM.Rdata")


## visualize ##

# PCA is appropriate for euclidean distances 

# get a few columns to load back into the clr dataframe 

sample_data <-  data.frame(sample_data(ps_WFDP_final))

sites <- dplyr::select(sample_data, EM_Sample_Name, WFDP_Code, Host_ID, Cell)

clr_sites_EM <- merge(sites, clr_WFDP_EM, by = 'row.names')

#PCA of differences in composition for HOSTS
# Scaling because these CLR values are not scaled 
pca_clr_EM = prcomp(clr_sites_EM[6:265], center = T, scale = T)

sd.pca_clr_EM = pca_clr_EM$sdev
loadings.pca_clr_EM = pca_clr_EM$rotation
names.pca_clr_EM = colnames(clr_sites_EM[6:265])
scores.pca_clr_EM = as.data.frame(pca_clr_EM$x)
scores.pca_clr_EM$Cell = clr_sites_EM$Cell
scores.pca_clr_EM$Host_ID = clr_sites_EM$Host_ID
summary(pca_clr_EM)


# PCA scores are 'scores.pca_clr_EM' with column for Host_ID

loadings.pca_clr_EM <- as.data.frame(loadings.pca_clr_EM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_clr_EM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


# set colors for hosts 
                # ABAM      ABGR      ALRU       TABR        TSHE        
EM_hosts <- c("#9b5fe0", "#16a4d8", "#60dbe8", "#efdf48", "#d64e12")


# Plot the Results by host alone
PCA_plot_host <- ggplot(scores.pca_clr_EM, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_colour_manual(values=EM_hosts, 
                      name="Focal Tree Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "TABR", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU","TABR", "TSHE")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Focal Tree Species") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)  
  )

PCA_plot_host

# Looks pretty wacky 


# Outcome: Just a bit of exploratory work for now. Next step is to represent the traits of the fungal 
# community in a numerical way, to be able to add it to the models. 


# -- #### END ##### -- 
