# -----------------------------------------------------------------------------#
# Exploring fungal taxonomic and functional variation in WFDP for EM community 
# Original Author: L. McKinley Nevins 
# February 9, 2025
# Software versions:  R v 4.5.2
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
library(ggfortify); packageVersion("ggfortify")
library(gginnards); packageVersion("gginnards")
library(ggrepel); packageVersion("ggrepel")

#################################################################################
#                               Main workflow                                   #
#  Generate dataframe to characterize functions of the EM community of each     #
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
ps_WFDP_raw <- readRDS("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_EM_final.RDS")

# grab final dataframe of raw ASVs 
trees_EM_raw <- otu_table(ps_WFDP_raw) %>% as("matrix") %>% as.data.frame()

# Load in the hot-coded trait values for each ASV
traits_WFDP_EM <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_EM_trait_mat.csv")

# Set ASVs as rownames 
traits_EM <- column_to_rownames(traits_WFDP_EM, var = "X")

# Pull out hydrophobicity, as I only want to consider the exploration types 
traits_EM_ET <- dplyr::select(traits_EM, -hydrophilic, -hydrophobic)

#################################################################################

#################################### -- 
# (2) ANALYZE FUNCTIONAL COMPOSITION
#################################### -- 

# I want to look at the fungal functions for each tree community and compare the 
# relative portions of traits that are present using the clr values. 

## Want to use the raw ASV abundances first, then do clr transformation of the traits after they have been 
# aggregated across the ASVs

# check structure 
str(trees_EM_raw)
str(traits_EM_ET)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(trees_EM_raw)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_EM <- as.data.frame(trees_EM_raw)
tree_matrix_EM[] <- lapply(tree_matrix_EM, as.numeric)  
tree_matrix_EM <- as.matrix(tree_matrix_EM)
rownames(tree_matrix_EM) <- tree_ids 


# Keep rownames
asv_ids <- rownames(traits_EM_ET)

traits_matrix_EM_ET <- as.data.frame(traits_EM_ET)
traits_matrix_EM_ET[] <- lapply(traits_matrix_EM_ET, as.numeric)
traits_matrix_EM_ET <- as.matrix(traits_matrix_EM_ET)
rownames(traits_matrix_EM_ET) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_EM)) # None        
sum(is.na(traits_matrix_EM_ET))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM_ET))
# all are shared

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM_ET <- traits_matrix_EM_ET[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM_ET)) #TRUE

# All match, good to go here 

############# 

# Trait composition: trees × traits using raw ASV counts 

trait_abund_per_tree_ET <- tree_matrix_EM %*% traits_matrix_EM_ET
# Output is a matrix of totals of how much each trait is represented in a particular
# tree’s fungal community, using raw abundance 

# convert the abundance of each trait for each tree into a proportion
trait_prop_per_tree_ET <- trait_abund_per_tree_ET / rowSums(trait_abund_per_tree_ET)


# clr transform this dataset to get the CLR transformed abundance for each trait and tree 
trait_clr_per_tree_ET <- decostand(trait_prop_per_tree_ET, method = "clr", pseudocount = 1e-06)


## Save file of CLR abundance for each trait 
write.csv(trait_clr_per_tree_ET, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/EM_trait_clr_per_tree.csv")

# Get dataframe for plotting 
trait_sums_df_ET <- as.data.frame(trait_clr_per_tree_ET)
trait_sums_df_ET$Sample_code <- rownames(trait_sums_df_ET)


# Reshape to long format
trait_sums_long_ET <- trait_sums_df_ET %>%
  pivot_longer(-Sample_code, names_to = "Trait", values_to = "CLR_Abund")


# Join table of tree environmental data 
sample_data <- sample_data(ps_WFDP_raw) %>% as("matrix") %>% as.data.frame()

sample_data <- sample_data %>% rownames_to_column(var = "Sample_code")

trait_sums_long_ET <- trait_sums_long_ET %>%
  left_join(sample_data, by = "Sample_code")

################# -- 

# Can now explore the clr-weighted trait profiles of each host tree! 

# Clean up the trait names 
trait_sums_long_ET <- trait_sums_long_ET %>%
  mutate(Trait_clean = Trait %>%
           str_remove("^ET_") %>%                # remove prefix
           str_replace_all("_", "-") %>%         # underscores → dashes
           str_to_title()                        # capitalize each word
  )


# The 'compatible name' column aligns with the mycorrhizal community sample code and can be used to label 
# individual trees 


# Set exploration type order to reflect range of short to long distance investment 
exploration_order <- c("Contact", "Contact-Short", "Short", "Contact-Medium", "Contact-Medium-Smooth",
                       "Contact-Medium-Fringe", "Contact-Long-Smooth", "Medium-Smooth", "Medium-Fringe", 
                       "Medium-Mat", "Medium-Long", "Medium-Long-Smooth", "Medium-Long-Fringe", "Long")

# apply order to the trait column 
trait_sums_long_ET <- trait_sums_long_ET %>%
  mutate(Trait_clean = factor(Trait_clean, levels = exploration_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(14, option = "E", direction = -1)
print(viridis_colors)


## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(trait_sums_long_ET$Compatible_Name)
all_traits <- unique(trait_sums_long_ET$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
trait_sums_long_ET_full <- trait_sums_long_ET %>%
  dplyr::select(Compatible_Name, Trait_clean, CLR_Abund, Host_ID) %>%
  complete(Compatible_Name = all_trees, Trait_clean = all_traits, fill = list(CLR_Abund = 0))

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos_ET <- trait_sums_long_ET_full %>%
  filter(CLR_Abund > 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

traits_neg_ET <- trait_sums_long_ET_full %>%
  filter(CLR_Abund < 0) %>%
  group_by(Compatible_Name, Host_ID) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup()

# Merge for plotting 

# Combine
traits_diverging_ET <- bind_rows(traits_pos_ET, traits_neg_ET)

# Make a couple variables factors 
traits_diverging_ET$Host_ID <- as.factor(traits_diverging_ET$Host_ID)
traits_diverging_ET$Compatible_Name <- as.factor(traits_diverging_ET$Compatible_Name)


# Clean Compatible_Name for less visual clutter 
traits_diverging_ET <- traits_diverging_ET %>%
  mutate(Tree_ID = Compatible_Name %>%
           str_remove("^W."))


# Diverging bar plot for trait representation in individual trees 

ET_diverging <- ggplot(traits_diverging_ET, aes(x = Tree_ID, y = CLR_Abund_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Host_ID, scales = "free_x", nrow = 2) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "E", direction = -1), exploration_order)) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  labs(
    x = "",
    y = "Relative Representation of Exploration Types",
    fill = "Exploration\n Type") +
  geom_hline(yintercept = 0, color = "red3", linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 0, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    legend.text = element_text(size = 14, colour="black"),
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=5,byrow=TRUE))

ET_diverging

ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/EM_traits_ET_plot.png", 
       plot = ET_diverging, width = 10, height = 11, units = "in", dpi = 300)


## Done with visualization of exploration types 


#################################################################################

############################### -- 
# (3) FORMAT VALUES FOR MODEL  
############################### -- 

# trait_clr_per_tree_ET is the dataframe of the CLR weighted trait values for each tree 

trait_df <- as.data.frame(trait_clr_per_tree_ET)

# Scale across all trees to make things comparable 
trait_df_scaled <- trait_df %>%
  mutate(across(
    starts_with("ET_"),
    ~ as.numeric(scale(.x))
  ))


# Drop long exploration type (could be any of them) to stop ET's from being perfectly 
# colinear (this is what is done automatically in linear regression and other models)

trait_df_scaled <- trait_df_scaled %>%
  select(-ET_long)


# Check colinearity 
pairs(trait_df_scaled %>% select(starts_with("ET_")))

# Looks good, nothing is perfectly related 


# Interpretation: Negative value = this tree emphasizes this strategy less than average
# Positive value = more than average
# Magnitude = strength of deviation across trees


# Save dataframe of scaled exploration type relative abundance values to use in the model - if you want 
# to load these directly in 

write.csv(trait_df_scaled, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/EM_traits_scaled.csv")


#################################################################################

############################### -- 
# (4) FUNGAL TRAIT PCA  
############################### -- 

# Perform PCA on EMF traits to reduce dimensionality and find 1-2 main axes of variation that can 
# be loaded into the model 

# Using the raw trait_df object because the PCA will center and scale the traits during the analysis 
trait_df <- as.data.frame(trait_clr_per_tree_ET)

# Take sample_data and merge in some columns for later 

env <- dplyr::select(sample_data, Sample_code, Host_ID, WFDP_Code)

trait_df <- rownames_to_column(trait_df, var = "Sample_code")

# Merge 
trait_df <- merge(trait_df, env, by = "Sample_code")


# Perform PCA on exploration types 
emf.pca = prcomp(trait_df[2:15], center = T, scale = T)

sd.emf = emf.pca$sdev
loadings.emf = emf.pca$rotation
trait.names.emf = colnames(trait_df[2:15])
scores.emf = as.data.frame(emf.pca$x)
scores.emf$WFDP_Code = trait_df$WFDP_Code
scores.emf$Host_ID = trait_df$Host_ID
summary(emf.pca)

# Save loadings for EMF exploration type traits
write.csv(loadings.emf, "./PCA_loadings_EMF_traits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.emf, "./PCA_scores_EMF_traits.csv")


# PCA scores are 'scores.emf' with column for Host_ID

loadings.emf <- as.data.frame(loadings.emf)

# get proportion of variance explained to add to each axis label 
pca_var <- emf.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

# Change loadings names to something cleaner 
new_loadings <- c("Contact", "Contact-Short", "Short", "Contact-Medium", "Contact-Medium Fringe", 
                  "Contact-Medium Smooth", "Medium Smooth", "Medium Fringe", "Medium Mat", 
                  "Medium-Long", "Medium-Long Smooth", "Medium-Long Fringe", "Contact-Long Smooth", "Long")

rownames(loadings.emf) <- new_loadings

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")

# Visualize
PCA_plot_emf <- ggplot(scores.emf, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  geom_segment(data = loadings.emf, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.emf, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.emf)),
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

PCA_plot_emf


# Get top 3 traits for PC1
top_PC1_emf <- loadings.emf[order(abs(loadings.emf$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_emf <- loadings.emf[order(abs(loadings.emf$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in Medium-Smooth, Long, and Contact-Medium Fringe types 

# PC2 is being driven by Contact-Short, Contact-Medium Smooth, and Medium-Long types

## Together the axes explain 40% of the variation 


# Perform broken stick analyses to check how many axes are capturing the majority of 
# the variation 

##Broken-Stick test for the significance of the loadings
print(emf.pca)

plot(emf.pca, type = "l")

ev = emf.pca$sdev^2

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


summary(emf.pca)


# 6 axes retained in the broken stick, so can definitely confidently use the first two

# -- #### END ##### -- 

