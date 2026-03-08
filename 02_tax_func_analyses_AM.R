# -----------------------------------------------------------------------------#
# Exploring fungal taxonomic and functional variation in WFDP for AM community 
# Original Author: L. McKinley Nevins 
# January 22, 2026
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.2.1
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
#                     cowplot v 1.2.0
#                     ggfortify v 0.4.19
#                     gginnards v 0.2.0.2
#                     ggrepel v 0.9.6
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
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
library(cowplot); packageVersion("cowplot")
library(ggfortify); packageVersion("ggfortify")
library(gginnards); packageVersion("gginnards")
library(ggrepel); packageVersion("ggrepel")

#################################################################################
#                               Main workflow                                   #
#  Generate dataframe to characterize functions of the AM community of each     #
#  host in the plot, and prepare the functional data in an appropriate format   #
#  for the Bayesian models.                                                     #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/"
setwd(wd)

# Read in WFDP EM phyloseq object with the raw ASV abundances
ps_WFDP_raw <- readRDS("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_AM_final.RDS")

# grab final dataframe of raw ASVs 
trees_AM_raw <- otu_table(ps_WFDP_raw) %>% as("matrix") %>% as.data.frame() # 27 AM or dual trees 


# Load in the hot-coded trait values for each ASV
traits_WFDP_AM <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_AM_trait_mat.csv")

# Set ASVs as rownames 
traits_AM <- column_to_rownames(traits_WFDP_AM, var = "X")

# Create separate dataframes for each of the frameworks, and right now exclude spore size 
traits_REA <- dplyr::select(traits_AM, Ancestral, Edaphophilic, Rhizophilic)

traits_CSR <- dplyr::select(traits_AM, C_S_R__C, C_S_R__S, C_S_R__R, C_S_R__S_R, C_S_R__R_C)

traits_YAS <- dplyr::select(traits_AM, Y_A_S__Y, Y_A_S__A, Y_A_S__S, Y_A_S__Y_S)

#################################################################################

#################################### -- 
# (2) ANALYZE FUNCTIONAL COMPOSITION
#################################### -- 

# I want to look at the fungal functions for each tree community and compare the 
# relative portions of traits that are present using the clr values. 

## Want to use the raw ASV abundances first, then do clr transformation of the traits after they have been 
# aggregated across the ASVs


## Lots of steps here that need to be repeated for the three different frameworks 


# --- ## REA ## --- ## 

# check structure 
str(trees_AM_raw)
str(traits_REA)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(trees_AM_raw)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_AM <- as.data.frame(trees_AM_raw)
tree_matrix_AM[] <- lapply(tree_matrix_AM, as.numeric)  
tree_matrix_AM <- as.matrix(tree_matrix_AM)
rownames(tree_matrix_AM) <- tree_ids 

# Keep rownames
asv_ids <- rownames(traits_REA)

traits_matrix_REA <- as.data.frame(traits_REA)
traits_matrix_REA[] <- lapply(traits_matrix_REA, as.numeric)
traits_matrix_REA <- as.matrix(traits_matrix_REA)
rownames(traits_matrix_REA) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_AM)) # None        
sum(is.na(traits_matrix_REA))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_AM), rownames(traits_matrix_REA))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_AM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_REA <- traits_matrix_REA[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_AM) == rownames(traits_matrix_REA)) #TRUE

# All match, good to go here 

############# 

# Trait composition: trees × traits using raw ASV counts 

trait_abund_per_tree_REA <- tree_matrix_AM %*% traits_matrix_REA
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion
trait_prop_per_tree_REA <- trait_abund_per_tree_REA / rowSums(trait_abund_per_tree_REA)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
trait_clr_per_tree_REA <- decostand(trait_prop_per_tree_REA, method = "clr", pseudocount = 1e-06)


## Save file of CLR abundance for each trait 
write.csv(trait_clr_per_tree_REA, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_trait_clr_per_tree_REA.csv")

# Get dataframe for plotting 
trait_sums_df_REA <- as.data.frame(trait_clr_per_tree_REA)
trait_sums_df_REA$Sample_code <- rownames(trait_sums_df_REA)


# Reshape to long format
trait_sums_long_REA <- trait_sums_df_REA %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Abund")


# Join table of tree environmental data 
sample_data <- sample_data(ps_WFDP_raw) %>% as("matrix") %>% as.data.frame()

sample_data <- sample_data %>% rownames_to_column(var = "Sample_code")

trait_sums_long_REA <- trait_sums_long_REA %>%
  left_join(sample_data, by = "Sample_code")

################# -- 

# Can now explore the clr-weighted trait profiles of each host tree! 

# The 'compatible name' column aligns with the mycorrhizal community sample code and can be used to label 
# individual trees 

# Set guild order 
REA_order <- c("Rhizophilic", "Edaphophilic", "Ancestral")

# apply order to the trait column 
trait_sums_long_REA <- trait_sums_long_REA %>%
  mutate(Trait = factor(Trait, levels = REA_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(3, option = "D", direction = -1)
print(viridis_colors)

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos_REA <- trait_sums_long_REA %>%
  filter(CLR_Abund > 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

traits_neg_REA <- trait_sums_long_REA %>%
  filter(CLR_Abund < 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

# Merge for plotting 

# Combine
traits_diverging_REA <- bind_rows(traits_pos_REA, traits_neg_REA)

# Make a couple variables factors 
traits_diverging_REA$Host_ID <- as.factor(traits_diverging_REA$Host_ID)
traits_diverging_REA$Compatible_Name <- as.factor(traits_diverging_REA$Compatible_Name)


# Clean Compatible_Name for less visual clutter 
traits_diverging_REA <- traits_diverging_REA %>%
  mutate(Tree_ID = Compatible_Name %>%
           str_remove("^W."))


# Diverging bar plot for trait representation in individual trees 

REA_diverging <- ggplot(traits_diverging_REA, aes(x = Tree_ID, y = CLR_Abund_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 2) +
  scale_fill_manual(
    values = setNames(viridis(3, option = "D", direction = -1), REA_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Representation of REA Functional Guilds",
    fill = "REA Functional Guild") +
  geom_hline(yintercept = 0, color = "red3", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    legend.text = element_text(size = 14, colour="black"),
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.position = "bottom") 

REA_diverging

ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_traits_REA_plot.png", 
       plot = REA_diverging, width = 7, height = 8, units = "in", dpi = 300)

## Done with visualization of Rhizophilic-Edaphophilic-Ancestral Guilds 


# --- ## CSR ## --- ## 

# check structure 
str(trees_AM_raw)
str(traits_CSR)

# trees matrix was already processed above, so just need to do for CSR

# make numeric matrix

# Keep rownames
asv_ids <- rownames(traits_CSR)

traits_matrix_CSR <- as.data.frame(traits_CSR)
traits_matrix_CSR[] <- lapply(traits_matrix_CSR, as.numeric)
traits_matrix_CSR <- as.matrix(traits_matrix_CSR)
rownames(traits_matrix_CSR) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_AM)) # None        
sum(is.na(traits_matrix_CSR))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_AM), rownames(traits_matrix_CSR))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_AM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_CSR <- traits_matrix_CSR[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_AM) == rownames(traits_matrix_CSR)) #TRUE

# All match, good to go here 

############# 

# Trait composition: trees × traits using raw ASV counts 

trait_abund_per_tree_CSR <- tree_matrix_AM %*% traits_matrix_CSR
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion
trait_prop_per_tree_CSR <- trait_abund_per_tree_CSR / rowSums(trait_abund_per_tree_CSR)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
trait_clr_per_tree_CSR <- decostand(trait_prop_per_tree_CSR, method = "clr", pseudocount = 1e-06)


## Save file of CSR abundance for each trait 
write.csv(trait_clr_per_tree_CSR, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_trait_clr_per_tree_CSR.csv")

# Get dataframe for plotting 
trait_sums_df_CSR <- as.data.frame(trait_clr_per_tree_CSR)
trait_sums_df_CSR$Sample_code <- rownames(trait_sums_df_CSR)

# Clean up the trait names 
trait_sums_df_CSR <- trait_sums_df_CSR %>%
  mutate(Competitive = C_S_R__C, 
         `Stress Tolerant`= C_S_R__S, 
         Ruderal = C_S_R__R, 
         `Stress Tolerant-Ruderal` = C_S_R__S_R, 
         `Competitive-Ruderal` = C_S_R__R_C)

# Subset 
trait_sums_df_CSR <- dplyr::select(trait_sums_df_CSR, Competitive, `Stress Tolerant`, Ruderal, `Stress Tolerant-Ruderal`, 
                                   `Competitive-Ruderal`, Sample_code)

# Reshape to long format
trait_sums_long_CSR <- trait_sums_df_CSR %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Abund")


# Join table of tree environmental data 
trait_sums_long_CSR <- trait_sums_long_CSR %>%
  left_join(sample_data, by = "Sample_code")

################# -- 

# Can now explore the clr-weighted trait profiles of each host tree! 

# The 'compatible name' column aligns with the mycorrhizal community sample code and can be used to label 
# individual trees 


# Set guild order 
CSR_order <- c("Competitive", "Stress Tolerant", "Ruderal", "Stress Tolerant-Ruderal", "Competitive-Ruderal")

# apply order to the trait column 
trait_sums_long_CSR <- trait_sums_long_CSR %>%
  mutate(Trait = factor(Trait, levels = CSR_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(5, option = "C", direction = -1)
print(viridis_colors)

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos_CSR <- trait_sums_long_CSR %>%
  filter(CLR_Abund > 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

traits_neg_CSR <- trait_sums_long_CSR %>%
  filter(CLR_Abund < 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

# Merge for plotting 

# Combine
traits_diverging_CSR <- bind_rows(traits_pos_CSR, traits_neg_CSR)

# Make a couple variables factors 
traits_diverging_CSR$Host_ID <- as.factor(traits_diverging_CSR$Host_ID)
traits_diverging_CSR$Compatible_Name <- as.factor(traits_diverging_CSR$Compatible_Name)


# Clean Compatible_Name for less visual clutter 
traits_diverging_CSR <- traits_diverging_CSR %>%
  mutate(Tree_ID = Compatible_Name %>%
           str_remove("^W."))


# Diverging bar plot for trait representation in individual trees 

CSR_diverging <- ggplot(traits_diverging_CSR, aes(x = Tree_ID, y = CLR_Abund_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 2) +
  scale_fill_manual(
    values = setNames(viridis(5, option = "C", direction = -1), CSR_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Representation of CSR Strategies",
    fill = "CSR Strategies") +
  geom_hline(yintercept = 0, color = "red3", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    legend.text = element_text(size = 14, colour="black"),
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

CSR_diverging

ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_traits_CSR_plot.png", 
       plot = CSR_diverging, width = 7, height = 8, units = "in", dpi = 300)

## Done with visualization of Competitive-Stress Tolerant-Ruderal Strategies 


# --- ## YAS ## --- ## 

# check structure 
str(trees_AM_raw)
str(traits_YAS)

# trees matrix was already processed above, so just need to do for YAS

# make numeric matrix

# Keep rownames
asv_ids <- rownames(traits_YAS)

traits_matrix_YAS <- as.data.frame(traits_YAS)
traits_matrix_YAS[] <- lapply(traits_matrix_YAS, as.numeric)
traits_matrix_YAS <- as.matrix(traits_matrix_YAS)
rownames(traits_matrix_YAS) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_AM)) # None        
sum(is.na(traits_matrix_YAS))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_AM), rownames(traits_matrix_YAS))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_AM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_YAS <- traits_matrix_YAS[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_AM) == rownames(traits_matrix_YAS)) #TRUE

# All match, good to go here 

############# --

# Trait composition: trees × traits using raw ASV counts 

trait_abund_per_tree_YAS <- tree_matrix_AM %*% traits_matrix_YAS
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion
trait_prop_per_tree_YAS <- trait_abund_per_tree_YAS / rowSums(trait_abund_per_tree_YAS)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
trait_clr_per_tree_YAS <- decostand(trait_prop_per_tree_YAS, method = "clr", pseudocount = 1e-06)


## Save file of YAS abundance for each trait 
write.csv(trait_clr_per_tree_YAS, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_trait_clr_per_tree_YAS.csv")

# Get dataframe for plotting 
trait_sums_df_YAS <- as.data.frame(trait_clr_per_tree_YAS)
trait_sums_df_YAS$Sample_code <- rownames(trait_sums_df_YAS)

# Clean up the trait names 
trait_sums_df_YAS <- trait_sums_df_YAS %>%
  mutate(`Growth Yield` = Y_A_S__Y, 
         `Resource Acquisition`= Y_A_S__A, 
         `Stress Tolerance` = Y_A_S__S, 
         `Growth Yield-Stress Tolerance` = Y_A_S__Y_S)

# Subset 
trait_sums_df_YAS <- dplyr::select(trait_sums_df_YAS, `Growth Yield`, `Resource Acquisition`, `Stress Tolerance`, 
                                   `Growth Yield-Stress Tolerance`, Sample_code)

# Reshape to long format
trait_sums_long_YAS <- trait_sums_df_YAS %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Abund")


# Join table of tree environmental data 
trait_sums_long_YAS <- trait_sums_long_YAS %>%
  left_join(sample_data, by = "Sample_code")

################# -- 

# Can now explore the clr-weighted trait profiles of each host tree! 

# The 'compatible name' column aligns with the mycorrhizal community sample code and can be used to label 
# individual trees 


# Set guild order 
YAS_order <- c("Growth Yield", "Resource Acquisition", "Stress Tolerance", "Growth Yield-Stress Tolerance")

# apply order to the trait column 
trait_sums_long_YAS <- trait_sums_long_YAS %>%
  mutate(Trait = factor(Trait, levels = YAS_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(4, option = "G", direction = -1)
print(viridis_colors)

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos_YAS <- trait_sums_long_YAS %>%
  filter(CLR_Abund > 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

traits_neg_YAS <- trait_sums_long_YAS %>%
  filter(CLR_Abund < 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

# Merge for plotting 

# Combine
traits_diverging_YAS <- bind_rows(traits_pos_YAS, traits_neg_YAS)

# Make a couple variables factors 
traits_diverging_YAS$Host_ID <- as.factor(traits_diverging_YAS$Host_ID)
traits_diverging_YAS$Compatible_Name <- as.factor(traits_diverging_YAS$Compatible_Name)


# Clean Compatible_Name for less visual clutter 
traits_diverging_YAS <- traits_diverging_YAS %>%
  mutate(Tree_ID = Compatible_Name %>%
           str_remove("^W."))


# Diverging bar plot for trait representation in individual trees 

YAS_diverging <- ggplot(traits_diverging_YAS, aes(x = Tree_ID, y = CLR_Abund_scaled, fill = Trait)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 2) +
  scale_fill_manual(
    values = setNames(viridis(4, option = "G", direction = -1), YAS_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Representation of YAS Strategies",
    fill = "YAS Strategies") +
  geom_hline(yintercept = 0, color = "red3", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    legend.text = element_text(size = 14, colour="black"),
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

YAS_diverging

ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_traits_YAS_plot.png", 
       plot = YAS_diverging, width = 7, height = 8, units = "in", dpi = 300)

## Done with visualization of Growth Yield-Resource Acquisition-Stress Tolerance Strategies 


#################################################################################

############################### -- 
# (3) FORMAT VALUES FOR MODEL  
############################### -- 

# !! Can skip this if only going to the PCA

# Scale and check the three AM trait dataframes to get them ready for the model 

# --- ## REA ## --- ## 

trait_REA_df <- as.data.frame(trait_clr_per_tree_REA)

# Scale across all trees to make things comparable 
trait_REA_df_scaled <- trait_REA_df %>%
  mutate(across(
    everything(),
    ~ as.numeric(scale(.x))
  ))

# Drop Ancestral (could be any of them) to stop guilds from being perfectly 
# colinear (this is what is done automatically in linear regression and other models)

trait_REA_df_scaled <- trait_REA_df_scaled %>%
  select(-Ancestral)

# Check colinearity 
pairs(trait_REA_df_scaled %>% select(everything()))

# Looks good, nothing is perfectly related 

# Interpretation: Negative value = this tree emphasizes this strategy less than average
# Positive value = more than average
# Magnitude = strength of deviation across trees

# Save dataframe of scaled REA functional guild relative abundance values to use in the model 

write.csv(trait_REA_df_scaled, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_traits_REA_scaled.csv")


# --- ## CSR ## --- ## 

trait_CSR_df <- as.data.frame(trait_clr_per_tree_CSR)

# Scale across all trees to make things comparable 
trait_CSR_df_scaled <- trait_CSR_df %>%
  mutate(across(
    everything(),
    ~ as.numeric(scale(.x))
  ))


# Fix trait names again 
trait_CSR_df_scaled <- trait_CSR_df_scaled %>%
  mutate(Competitive = C_S_R__C, 
         `Stress Tolerant`= C_S_R__S, 
         Ruderal = C_S_R__R, 
         `Stress Tolerant-Ruderal` = C_S_R__S_R, 
         `Competitive-Ruderal` = C_S_R__R_C)

# Drop Competitive-Ruderal (could be any of them) to stop guilds from being perfectly 
# colinear (this is what is done automatically in linear regression and other models)
# Also drop the columns with the old name format 

trait_CSR_df_scaled <- trait_CSR_df_scaled %>%
  select(Competitive, `Stress Tolerant`, `Ruderal`, `Stress Tolerant-Ruderal`)

# Check colinearity 
pairs(trait_CSR_df_scaled %>% select(everything()))

# Looks good, nothing is perfectly related 

# Save dataframe of scaled CSR functional guild relative abundance values to use in the model 

write.csv(trait_CSR_df_scaled, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_traits_CSR_scaled.csv")


# --- ## YAS ## --- ## 

trait_YAS_df <- as.data.frame(trait_clr_per_tree_YAS)

# Scale across all trees to make things comparable 
trait_YAS_df_scaled <- trait_YAS_df %>%
  mutate(across(
    everything(),
    ~ as.numeric(scale(.x))
  ))

# Fix trait names again 
trait_YAS_df_scaled <- trait_YAS_df_scaled %>%
  mutate(`Growth Yield` = Y_A_S__Y, 
         `Resource Acquisition`= Y_A_S__A, 
         `Stress Tolerance` = Y_A_S__S, 
         `Growth Yield-Stress Tolerance` = Y_A_S__Y_S)

# Drop Growth Yield-Stress Tolerance (could be any of them) to stop guilds from being perfectly 
# colinear (this is what is done automatically in linear regression and other models)
# Also drop the columns with the old name format 

trait_YAS_df_scaled <- trait_YAS_df_scaled %>%
  select(`Growth Yield`, `Resource Acquisition`, `Stress Tolerance`)

# Check colinearity 
pairs(trait_YAS_df_scaled %>% select(everything()))

# Looks good, nothing is perfectly related 

# Save dataframe of scaled YAS functional guild relative abundance values to use in the model 

write.csv(trait_YAS_df_scaled, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/AM_traits_YAS_scaled.csv")

# Traits plotted and scaled dataframes saved for modeling 

#################################################################################

############################### -- 
# (4) FUNGAL TRAIT PCA  
############################### -- 

# Perform PCA on AM traits to reduce dimensionality and find 1-2 main axes of variation that can 
# be loaded into the model 

# Going to merge together the three different trait frameworks. They are compatible so it could 
# be interesting to see what patterns emerge across the frameworks on the two axes 

# Using the raw trait dataframes because the PCA will center and scale the traits during the analysis 

trait_REA_df <- as.data.frame(trait_clr_per_tree_REA) %>% rownames_to_column(var = "Sample_code")

trait_CSR_df <- as.data.frame(trait_clr_per_tree_CSR) %>% rownames_to_column(var = "Sample_code")

trait_YAS_df <- as.data.frame(trait_clr_per_tree_YAS) %>% rownames_to_column(var = "Sample_code")


# Take sample_data and merge in some columns for later 
env <- dplyr::select(sample_data, Sample_code, Host_ID, WFDP_Code)


# Merge all three frameworks together, plus the categorical variables 
trait_df <- merge(trait_REA_df, trait_CSR_df, by = "Sample_code")

trait_df <- merge(trait_df, trait_YAS_df, by = "Sample_code")

trait_df <- merge(trait_df, env, by = "Sample_code")

# All ready to go 


# Perform PCA on AM functional frameworks
am.pca = prcomp(trait_df[2:13], center = T, scale = T)

sd.am = am.pca$sdev
loadings.am = am.pca$rotation
trait.names.am = colnames(trait_df[2:13])
scores.am = as.data.frame(am.pca$x)
scores.am$WFDP_Code = trait_df$WFDP_Code
scores.am$Host_ID = trait_df$Host_ID
summary(am.pca)

# Save loadings for AM functional groups 
write.csv(loadings.am, "./PCA_loadings_AM_traits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.am, "./PCA_scores_AM_traits.csv")


# PCA scores are 'scores.am' with column for Host_ID

loadings.am <- as.data.frame(loadings.am)

# get proportion of variance explained to add to each axis label 
pca_var <- am.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

# Change loadings names to something cleaner
new_loadings <- c("Ancestral (REA)", "Edaphophilic (REA)", "Rhizophilic (REA)", "Competitor (CSR)", 
                  "Stress Tolerator (CSR)", "Ruderal (CSR)", "Stress Tolerator-Ruderal (CSR)", 
                  "Competitor-Ruderal (CSR)", "Growth Yield (YAS)", "Resource Acquisition (YAS)", 
                  "Stress Tolerance (YAS)", "Growth Yield-Stress Tolerance (YAS)")

rownames(loadings.am) <- new_loadings

# set colors for hosts 
                # ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")

# Visualize
PCA_plot_am <- ggplot(scores.am, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  geom_segment(data = loadings.am, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.am, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.am)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values=all_hosts, name="Focal Species",
                     breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                     labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Focal Species") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black"))

PCA_plot_am


# Get top 3 traits for PC1
top_PC1_am <- loadings.am[order(abs(loadings.am$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_am <- loadings.am[order(abs(loadings.am$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in Rhizophilic (REA), Growth Yield (YAS), and Ruderal (CSR) strategies. This is 
# good because these are one from each of the functional frameworks, but they can be relatable in 
# thinking about fast growing Glomeraceae taxa that produce a lot of hyphal biomass in the roots 

# PC2 is being driven by Growth Yield-Stress Tolerance (YAS), Stress Tolerator-Ruderal (CSR), and Competitor (CSR)

## Together the axes explain 78.4% of the variation 

# Not that many species points showing because a lot of them had the same scores, especially for the 
# REA framework. 


##Broken-Stick test for the significance of the loadings
print(am.pca)

plot(am.pca, type = "l")

ev = am.pca$sdev^2

evplot = function(ev) {
  # Broken stick model (MacArthur 1957)
  n = length(ev)
  bsm = data.frame(j=seq(1:n), p=0)
  bsm$p[1] = 1/n
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p = 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

evplot(ev)


summary(am.pca)

# 3 axes retained in the broken stick, so can definitely confidently use the first two


# -- #### END ##### -- 

