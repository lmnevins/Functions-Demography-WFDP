# -----------------------------------------------------------------------------#
# Gather and Format Data for WFDP for AM community 
# Original Author: L. McKinley Nevins 
# January 2, 2026
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
#  Subset fungal community data to just AM in WFDP. Merge with the table of     #
#  trait assignments. Perform clr transformation to generate comparable         #
#  relative abundance values.                                                   #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/"
setwd(wd)


# Load phyloseq object produced from filtered FUNGuild guilds to target just AM fungi ####
ps_AM <- readRDS("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/Phylogeny_Outputs/AM_funguild_final_2025.RDS")

# Select just WFDP as the site 
WFDP_AM <- subset_samples(ps_AM, Site=="WFDP")

# get otu table 
WFDP_AM_asv <- otu_table(WFDP_AM)

WFDP_AM_asv <- as.data.frame(WFDP_AM_asv)

# save this as a csv file 
write.csv(WFDP_AM_asv, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_asvs_AM.csv")

# This was pared down to just the trees that we have both community and host tree trait measurements
# which results in 28 trees that are either AM or dual 

# load back in 
WFDP_AM_asv_sub <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_asvs_AM_trait_trees.csv")

#need the samples to be row names
WFDP_AM_asv_sub <- data.frame(WFDP_AM_asv_sub[,-1], row.names=WFDP_AM_asv_sub[,1])

#convert to phyloseq compatible object 
WFDP_AM_asv_sub <- phyloseq::otu_table(as.matrix(WFDP_AM_asv_sub), taxa_are_rows = F)


## Would load in table of families with trait assignments here, once it exists 
traits <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_functions_2026.csv")

traits$Family <- as.factor(traits$Family)


# Load in the TEMPORARY sample data file for WFDP 
env_WFDP <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_data/WFDP_enviro_data_AM_subset.csv", row.names = 1)


#####
#preparation of the tax_table from the full WFDP phyloseq object
WFDP_AM_tax <- tax_table(WFDP_AM)

WFDP_AM_tax <- as.data.frame(WFDP_AM_tax)

#convert to phyloseq compatible object 
WFDP_AM_tax <- phyloseq::tax_table(as.matrix(WFDP_AM_tax))


#################################################################################

####################################### --- 
# (2) CREATE AND SAVE PHYLOSEQ OBJECT
####################################### --- 


## create object 
ps_WFDP_AM_final <- phyloseq(otu_table(WFDP_AM_asv_sub), tax_table(WFDP_AM_tax), sample_data(env_WFDP))


#### Trim to reflect just the community in WFDP
# number of taxa - 438
ntaxa(ps_WFDP_AM_final)

# number of samples - 28 host trees - that are either AM or dual 
nsamples(ps_WFDP_AM_final)

asv <- otu_table(ps_WFDP_AM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

sample_count <- rowSums(asv)

sample_count <- as.data.frame(sample_count) #all trees still have some ASV's

asv_count <- colSums(asv)

asv_count <- as.data.frame(asv_count) #a lot of ASVs no longer present in the samples 

#remove taxa that aren't present in this new subset 
ps_WFDP_AM_final <- subset_taxa(ps_WFDP_AM_final, taxa_sums(ps_WFDP_AM_final) > 0)

# number of taxa - now 125
ntaxa(ps_WFDP_AM_final)


# Get taxon table for AM community in WFDP
WFDP_AM_tax <- tax_table(ps_WFDP_AM_final)

WFDP_AM_tax <- as.data.frame(WFDP_AM_tax) %>% rownames_to_column(var = "X") # 125 taxa


# Save as the file to assess the total ASV's and figure out proportions that have trait assignments, etc. 

write.csv(WFDP_AM_tax, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/all_AM_tax.csv")


# pull out the genus and family levels, as well as the original ASVs. This will be necessary 
# to match them up with the site x environment matrix 

WFDP_AM_tax <- dplyr::select(WFDP_AM_tax, X, Family, Genus, Species)

WFDP_AM_tax$Genus <- as.factor(WFDP_AM_tax$Genus)
WFDP_AM_tax$Family <- as.factor(WFDP_AM_tax$Family)

summary(WFDP_AM_tax)



# Merge datasets together 
matrix <- merge(WFDP_AM_tax, traits, by = "Family", all.y = TRUE)

#get column name set for name merging down below
matrix$ASV <- matrix$X

# the number is correct - 115 ASVs with Family and trait assignments 


#start names file with the original ASVs 
ASV <- phyloseq::taxa_names(ps_WFDP_AM_final)
ASV_long <- as.data.frame(ASV)

#change ASV names to something nicer to work with
taxa_names(ps_WFDP_AM_final)
n_seqs <- seq(ntaxa(ps_WFDP_AM_final))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_WFDP_AM_final) <- paste("ASV", formatC(n_seqs, 
                                                  width = len_n_seqs, 
                                                  flag = "0"), sep = "_")
taxa_names(ps_WFDP_AM_final)

# get shortened names 
ASV2 <- taxa_names(ps_WFDP_AM_final) 
ASV_short <- as.data.frame(ASV2)

# join two dataframes
all_names <- cbind(ASV_long, ASV_short)


#merge taxa_names with the traits file to get the updated ASV names 
matrix <- merge(matrix, all_names, by = "ASV")

# save matrix file as key to the long and short ASV assignments 
write.csv(matrix, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/all_ASV_names_WFDP_AM_2025.csv")


## !! update for whatever AM trait columns are 

# reformat a tiny bit just to get ASVs as the species in rownames, and my traits only 
traits_WFDP_AM <- matrix %>% dplyr::select(ASV2, Weber_Guild, intraradical_hyphae, extraradical_hyphae, Ancestral,
                                    Edaphophilic, Rhizophilic, spore_vol, C_S_R, C_S_R__C, C_S_R__S, 
                                    C_S_R__R, C_S_R__S_R, C_S_R__R_C, Y_A_S, Y_A_S__Y, Y_A_S__A, 
                                    Y_A_S__S, Y_A_S__Y_S) %>% column_to_rownames(var = "ASV2") 

#ASVs are now the row names and there are columns for the binary coding of each trait level 

# Save traits file for later analyses 

write.csv(traits_WFDP_AM, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_AM_trait_mat.csv")



# We want these to be the ones that are contained within the final phyloseq object 
taxa_to_keep <- matrix$ASV2

# Prune taxa from the phyloseq object that are NOT in the taxa_to_keep list 
ps_WFDP_AM_final <- prune_taxa(taxa_to_keep, ps_WFDP_AM_final)


#SAVE THE FINAL RAW DATA PHYLOSEQ OBJECT 
# For the taxa that have functional assignments in the community 
saveRDS(ps_WFDP_AM_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_AM_final.RDS")
#################################################################################

# Get to know the phyloseq data ####

# number of taxa - 115
ntaxa(ps_WFDP_AM_final)

# number of samples - 28
nsamples(ps_WFDP_AM_final)

# sample names
sample_names(ps_WFDP_AM_final)
rank_names(ps_WFDP_AM_final)
# taxa names
taxa_names(ps_WFDP_AM_final)

# ASV table
otu_table(ps_WFDP_AM_final) %>% View()

# how many sequences observed in each sample?
seq_counts <- otu_table(ps_WFDP_AM_final) %>% rowSums() %>% as.data.frame()
# No trees with no sequences 


# how many times was each taxon observed across the samples?
otu_table <- otu_table(ps_WFDP_AM_final) %>% colSums()


# how many different samples was each taxon found in?
asv <- otu_table(ps_WFDP_AM_final) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame

#inspect tax_table 
WFDP_tax_table <- as.data.frame(tax_table(ps_WFDP_AM_final))

#count number of classifications in each column to determine coverage to taxonomic levels 

colSums(!is.na(WFDP_tax_table))

#  Kingdom  Phylum    Class    Order   Family   Genus    Species 
#  115       115       115      115     115      93       60 
#  100%     100%      100%     100%    100%     80.8%     52.2% 


#################################################################################

########################### -- 
# (3) DATA TRANSFORMATION
########################### -- 

# Pursuing a transformation because sequence data are inherently compositional, but standard rarefaction 
# is not repeatable and leads to the exclusion of a huge fraction of your data 

# trying a centered log-ratio transformation per
# https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full

# Also: Van den Boogaart, K. G., and Tolosana-Delgado, R. (2013). Analyzing Compositional Data with R, 
# London, UK: Springer.

## Using clr transformation with an added pseudocount, to 

asv <- otu_table(ps_WFDP_AM_final) %>% as("matrix") %>% as.data.frame()

clr_WFDP_AM <- decostand(asv, method = "clr", pseudocount = 1e-06)

clr_WFDP_AM <- as.data.frame(clr_WFDP_AM)


# load the clr data back into a phyloseq object 
ps_WFDP_AM_final <- phyloseq(tax_table(tax_table(ps_WFDP_AM_final)),
                          otu_table(otu_table(clr_WFDP_AM, taxa_are_rows=FALSE)),
                          sample_data(sample_data(ps_WFDP_AM_final)))


#SAVE THE FINAL PHYLOSEQ OBJECT 
saveRDS(ps_WFDP_AM_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_clr_AM_final.RDS")

## This now only contains the fungal taxa that we have functional information for 


#################################################################################

##################################### --- 
# (4) VISUALIZE COMMUNITY VARIATION 
##################################### --- 


## Update all of this to be for the AM data 



## Look at ASV taxonomic richness across the samples 
# Doesn't factor in relative abundance, so using the raw phyloseq object. This is just showing the presence or 
# absence of an ASV


# Load in the raw reads for ASVs
ps_WFDP_raw <- readRDS("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_AM_final.RDS")

# Look at the spread of ASVs present in the samples 
tax <- tax_table(ps_WFDP_raw) %>% as.data.frame()

#move the rownames into a column 
tax <- rownames_to_column(tax, var = "ASV")

# Make Family a factor 
tax$Family <- as.factor(tax$Family)
summary(tax$Family)


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
fam_richness <- otu_tax %>%
  group_by(Tree_ID, Family) %>%
  summarise(
    n_asvs = n_distinct(ASV),
    .groups = "drop"
  )


# Calculate relative ASV richness to account for the fact that the trees have different numbers of ASVS
fam_rel <- fam_richness %>%
  group_by(Tree_ID) %>%
  mutate(rel_asvs = n_asvs / sum(n_asvs)) %>%
  ungroup()


# add column of Host_ID
fam_rel <- fam_rel %>%
  mutate(
    Host_ID = str_extract(Tree_ID, "(?<=-)[A-Z]{4}(?=-)")
  )

# add column to shorten the tree_id for plotting 
fam_rel <- fam_rel %>%
  mutate(
    tree_short = str_extract(Tree_ID, "[A-Z]{4}-\\d{2}")
  )





# Visualize 
richness_plot_EM <- ggplot(fam_rel,
                           aes(x = tree_short, y = rel_asvs, fill = Family)) +
  geom_col(width = 0.9) +
  facet_wrap(~ Host_ID, scales = "free_x") +
  labs(
    x = "",
    y = "Relative AM Fungal ASV Richness",
    fill = "Family"
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

aitchison_WFDP_AM <- dist(clr_WFDP_AM, method = "euclidean")

peek <- as.matrix(aitchison_WFDP_AM) %>% as.data.frame()

# save Aitchison Distance matrix 
save(aitchison_WFDP_AM, file="~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/aitchison_dist_AM.Rdata")


## visualize ##

# PCA is appropriate for euclidean distances 

# get a few columns to load back into the clr dataframe 

sample_data <-  data.frame(sample_data(ps_WFDP_AM_final))

sites <- dplyr::select(sample_data, AM.Sample.Name, WFDP_Code, Host_ID, Cell)

clr_sites_AM <- merge(sites, clr_WFDP_AM, by = 'row.names')

#PCA of differences in composition for HOSTS
# Scaling because these CLR values are not scaled 
pca_clr_AM = prcomp(clr_sites_AM[6:115], center = T, scale = T)

sd.pca_clr_AM = pca_clr_AM$sdev
loadings.pca_clr_AM = pca_clr_AM$rotation
names.pca_clr_AM = colnames(clr_sites_AM[6:115])
scores.pca_clr_AM = as.data.frame(pca_clr_AM$x)
scores.pca_clr_AM$Cell = clr_sites_AM$Cell
scores.pca_clr_AM$Host_ID = clr_sites_AM$Host_ID
summary(pca_clr_AM)


# PCA scores are 'scores.pca_clr_AM' with column for Host_ID

loadings.pca_clr_AM <- as.data.frame(loadings.pca_clr_AM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_clr_AM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


# set colors for hosts 
               # ALRU      CONU     TABR        THPL        
AM_hosts <- c("#60dbe8", "#8bd346","#efdf48", "#f9a52F")


# Plot the Results by host alone
PCA_plot_host_AM <- ggplot(scores.pca_clr_AM, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal(base_size = 11) +
  scale_colour_manual(values=AM_hosts, 
                      name="Focal Tree Species",
                      breaks=c("ALRU", "CONU", "TABR", "THPL"),
                      labels=c("ALRU", "CONU", "TABR", "THPL")) +
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

PCA_plot_host_AM

# Looks pretty wacky 


# Outcome: Just a bit of exploratory work for now. Next step is to represent the traits of the fungal 
# community in a numerical way, to be able to add it to the models. 


# -- #### END ##### -- 
