# -----------------------------------------------------------------------------#
# Exploring fungal taxonomic and functional variation in WFDP for EM community 
# Original Author: L. McKinley Nevins 
# February 9, 2025
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

#################################################################################
#                               Main workflow                                   #
#  Generate dataframes to characterize functions and perform analyses of the    #
#  taxonomic and functional composition of the fungal community of each host    #
#  tree across the plot.                                                        #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/"
setwd(wd)

# Read in WFDP EM phyloseq object 
ps_WFDP_final <- readRDS("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_clr_final.RDS")

# grab final dataframe of clr transformed ASVs 
clr_WFDP_EM <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame()

# Load in the hot-coded trait values for each ASV
traits_WFDP_EM <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_EM_trait_mat.csv")

# Set ASVs as rownames 

traits_WFDP_EM <- column_to_rownames(traits_WFDP_EM, var = "X")

#################################################################################

#################################### -- 
# (2) ANALYZE TAXONOMIC COMPOSITION
#################################### -- 

## Beta-diversity can be calculated using Aitchison Distance, which is essentially the 
# euclidean distance calculated between pairs of samples that have been transformed by CLR

# calculate Aitchison distance using dist() from base R 

aitchison_WFDP_EM <- dist(clr_WFDP_EM, method = "euclidean")

aitchison_WFDP_EM <- as.matrix(aitchison_WFDP_EM)

# save Aitchison Distance matrix 
save(aitchison_WFDP_EM, file="~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/aitchison_dist_WFDP_EM.Rdata")


######## Explore variation in the fungal communities between host tree species ## -- 

# PCA is appropriate for euclidean distances 

# get a few metdata columns to load back into the clr dataframe 

sample_data <-  data.frame(sample_data(ps_WFDP_final)) %>% rownames_to_column("Sample_ID")

sites <- select(sample_data, Sample_ID, Host_ID) %>% column_to_rownames("Sample_ID")

clr_sites_EM <- merge(sites, clr_WFDP_EM, by = 'row.names')

#PCA of differences in composition for host taxa  
pca_clr_EM = prcomp(clr_sites_EM[3:360], center = T, scale = F)

sd.pca_clr_EM = pca_clr_EM$sdev
loadings.pca_clr_EM = pca_clr_EM$rotation
names.pca_clr_EM = colnames(clr_sites_EM[3:360])
scores.pca_clr_EM = as.data.frame(pca_clr_EM$x)
scores.pca_clr_EM$Host_ID = clr_sites_EM$Host_ID
summary(pca_clr_EM)


# PCA scores are 'scores.pca_clr_EM' with column for Host 

loadings.pca_clr_EM <- as.data.frame(loadings.pca_clr_EM)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_clr_EM$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


#set colors for hosts 
                  # ABAM        ABGR        ALRU          CONU         TABR        THPL        TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66")


# Plot the Results by host
PCA_plot_host <- ggplot(scores.pca_clr_EM, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Host_ID), type = "norm", linewidth = 1, size = 1) +
  theme_minimal() +
  scale_colour_manual(values=all_hosts, 
                      name="Host_ID",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(title = "PCA Biplot: Community Variation Across Host Tree Taxa",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host_ID") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12))

PCA_plot_host

# PCA axes don't explain a lot of the variation 

#####ASSESS MULTIVARIATE HOMOGENEITY OF HOSTS########### -- 
#betadisper function in vegan 

#first object needs to be a dist object of Aitchison distances 

# rerun to get back in a distance format 
aitchison_WFDP_EM <- dist(clr_WFDP_EM, method = "euclidean")

#second object is the groups of interest, as a vector
betadisper.host <- vegan::betadisper(aitchison_WFDP_EM, sample_data$Host_ID, type = "median", sqrt.dist = FALSE)

betadisper.host

permutest(betadisper.host)
#groups dispersions are not different p = 0.103

#if the dispersion is different between groups, then examine
scores(betadisper.host, display = c("sites", "centroids"),
       choices = c(1,2))


#visualize 
plot(betadisper.host, axes = c(1,2), ellipse = FALSE, segments = FALSE, lty = "solid", label = TRUE, 
     label.cex = 0.5, col = c("#0D0887FF", "#5402A3FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66"))



boxplot(betadisper.host)
mod.HSD <- TukeyHSD(betadisper.host)
mod.HSD
plot(mod.HSD)

# no significant differences. The centroids are all very similar, there are just big visual 
# differences for the distances 

# Nicer boxplot 
distances_EM <- data.frame(
  Host_ID = betadisper.host$group,
  DistanceToCentroid = betadisper.host$distances
)


centroid_plot_EM <- ggplot(distances_EM, aes(x = Host_ID, y = DistanceToCentroid, fill = Host_ID)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  theme_minimal() +
  labs(title = "Beta Dispersion by Tree Host",
       x = "Host",
       y = "Distance to Centroid") +
  scale_fill_manual(values=all_hosts, 
                    name="Host_ID",
                    breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                    labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  theme(legend.position = "right")

centroid_plot_EM


# Test for significant differences between centroids 
centroid_aov <- aov(DistanceToCentroid ~ Host_ID,
                    data = distances_EM
)

summary(centroid_aov)

# Not significant here either p = 0.911


permanova.host <- vegan::adonis2(aitchison_WFDP_EM ~ Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.host

# vegan::adonis2(formula = aitchison_WFDP_EM ~ Host_ID, data = sample_data, permutations = 999, method = "euclidean")
#         Df SumOfSqs      R2      F Pr(>F)
# Model     6    356.5 0.09379 0.9142  0.911
# Residual 53   3444.4 0.90621              
# Total    59   3800.8 1.00000  

# Not significant 


##################################### -- 
x # Come back to this when I have the kriged envrionmental data for each host tree 

# Take distance to centroid values and regress with key environmental variables to assess if there are 
# relationships with environmental variation across the sites

# need to make dataset that has distance to centroid and the environmental variables 

# Have distance to centroid and environmental data for each individual tree 
distances_EM <- tibble::rownames_to_column(distances_EM, "Sample_ID")

# use environmental data that was pulled in earlier 
view(sample_data)

# pick some specific environmental variables of interest to compare 
enviro <- select(sample_data, Sample_ID, slope, aspect, elevation_m)

# combine with the distance to centroid summary table
centroid_enviro_EM <- merge(distances_EM, enviro, by = "Sample_ID")

### Right now just looking at this for individual tree info of slope, aspect, and elevation 
# Once I have the kriged data I can look at this across many more environmental variables 


### Testing these using the individual tree data

## Slope
slope <- ggplot(centroid_enviro_EM, aes(x = slope, y = DistanceToCentroid)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope

# test relationships 
lm_slope <- lm(DistanceToCentroid ~ slope, data = centroid_enviro_EM)
summary(lm_slope) # NOT SIGNIFICANT


## Aspect
aspect <- ggplot(centroid_enviro_EM, aes(x = aspect, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect

# test relationships 
lm_aspect <- lm(DistanceToCentroid ~ aspect, data = centroid_enviro_EM)
summary(lm_aspect) # NOT SIGNIFICANT


## Elevation
elev <- ggplot(centroid_enviro_EM, aes(x = elevation_m, y = DistanceToCentroid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
elev

# test relationships 
lm_elev <- lm(DistanceToCentroid ~ elevation_m, data = centroid_enviro_EM)
summary(lm_elev) # NOT SIGNIFICANT


### RESULT: No significant relationship between distance to centroid and 
# slope, aspect, or elevation for individual trees 


#################################################################################

#################################### -- 
# (3) ANALYZE FUNCTIONAL COMPOSITION
#################################### -- 

# I want to look at the fungal functions for each tree community and compare the 
# relative portions of traits that are present using the clr values. 

# check structure 
str(clr_WFDP_EM)
str(traits_WFDP_EM)

# make both numeric matrices
# Keep rownames
tree_ids <- rownames(clr_WFDP_EM)

# Convert to numeric matrix and add back in the rownames 
tree_matrix_EM <- as.data.frame(clr_WFDP_EM)
tree_matrix_EM[] <- lapply(tree_matrix_EM, as.numeric)  
tree_matrix_EM <- as.matrix(tree_matrix_EM)
rownames(tree_matrix_EM) <- tree_ids 

# Keep rownames
asv_ids <- rownames(traits_WFDP_EM)

traits_matrix_EM <- as.data.frame(traits_WFDP_EM)
traits_matrix_EM[] <- lapply(traits_matrix_EM, as.numeric)
traits_matrix_EM <- as.matrix(traits_matrix_EM)
rownames(traits_matrix_EM) <- asv_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_EM)) # None        
sum(is.na(traits_matrix_EM))  # None   

# Make sure ASVs match between traits and trees
# Find shared ASVs between the two datasets 
shared_asvs <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM))

# all are shared 

# set the ASVs to be in the same order 
asv_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM <- traits_matrix_EM[asv_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM)) #TRUE

# All match, good to go here 

############# -- 

# CLR-weighted trait composition: trees × traits
trait_sums_per_tree <- tree_matrix_EM %*% traits_matrix_EM
# Output is a matrix of CLR-weighted totals of how much each trait is represented in a particular
# tree’s fungal community.


# Get dataframe for plotting 
trait_sums_df <- as.data.frame(trait_sums_per_tree)
trait_sums_df$Sample_ID <- rownames(trait_sums_df)


# Reshape to long format
trait_sums_long <- trait_sums_df %>%
  pivot_longer(-Sample_ID, names_to = "Trait", values_to = "CLR_Sum")


# Join table of tree environmental data 
x # Pull in the imputed environmental data when I have it 

trait_sums_long <- trait_sums_long %>%
  left_join(sample_data, by = "Sample_ID")


################# -- 

# Can now explore the clr-weighted trait profiles of each host tree! 

# Interested in exploration type and hydrophobicity separately, so pull these out of 
# the bigger dataset 

exploration_traits <- trait_sums_long %>%
  filter(str_starts(Trait, "ET_"))

hydro_traits <- trait_sums_long %>%
  filter(str_starts(Trait, "hydro"))

# Clean up the trait names 
exploration_traits <- exploration_traits %>%
  mutate(Trait_clean = Trait %>%
           str_remove("^ET_") %>%                # remove prefix
           str_replace_all("_", "-") %>%         # underscores → dashes
           str_to_title()                        # capitalize each word
  )

hydro_traits <- hydro_traits %>%
  mutate(Trait_clean = Trait %>%
           str_to_title()                       
  )

# Clean up label names 
exploration_traits <- exploration_traits %>%
  mutate(Tree = Compatible_Name %>%
           str_remove("^W.")                
  )


hydro_traits <- hydro_traits %>%
  mutate(Tree = Compatible_Name %>%
           str_remove("^W.")                
  )

# This creates a 'Tree' label which aligns with the numbering system for the fungal communities 


# Set exploration type order to reflect range of short to long distance investment 
exploration_order <- c("Contact", "Contact-Short", "Short", "Contact-Medium", "Contact-Medium-Smooth",
                       "Contact-Medium-Fringe", "Contact-Long-Smooth", "Medium-Smooth", "Medium-Fringe", 
                       "Medium-Mat", "Medium-Long", "Medium-Long-Smooth", "Medium-Long-Fringe", "Long"
)

# apply order to the trait column 
exploration_traits <- exploration_traits %>%
  mutate(Trait_clean = factor(Trait_clean, levels = exploration_order))

# get viridis colors for the traits 
library(viridis)

viridis_colors <- viridis(14, option = "D", direction = -1)
print(viridis_colors)


# plot per tree 
ET_tree_plot <- ggplot(exploration_traits, aes(x = Tree, y = CLR_Sum, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative Exploration Type Abundance",
       fill = "Exploration Type"
  ) +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

ET_tree_plot



## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(exploration_traits$Tree)
all_traits <- unique(exploration_traits$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
exploration_traits_full <- exploration_traits %>%
  select(Tree, Trait_clean, CLR_Sum) %>%
  complete(Tree = all_trees, Trait_clean = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 

# Make a 'Tree' column in sample_data
sample_data <- sample_data %>%
  mutate(Tree = Compatible_Name %>%
           str_remove("^W.")
  )

# Grab Host_ID from sample_data
hosts <- select(sample_data, Tree, Host_ID)

# merge hosts to data 
exploration_traits_full <- merge(exploration_traits_full, hosts, by = "Tree")


## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
traits_pos <- exploration_traits_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(CLR_Sum)) %>%
  ungroup()

traits_neg <- exploration_traits_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup() %>%
  mutate(CLR_Sum_scaled = -abs(CLR_Sum_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
traits_diverging <- bind_rows(traits_pos, traits_neg)


# Diverging bar plot for trait representation in individual trees 

ET_diverging_WFDP <- ggplot(traits_diverging, aes(x = Tree, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = setNames(viridis(14, option = "D", direction = -1), exploration_order)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Tree",
    y = "Exploration Type Abundance",
    fill = "Exploration Type") +
  geom_hline(yintercept = 0, color = "red", linewidth = 2) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

ET_diverging_WFDP


# set hydro palette
# hydrophilic   hydrophobic
hydro_colors <- c("#4682B4", "#FF6347")


# plot per tree 
hydro_tree_plot <- ggplot(hydro_traits, aes(x = Tree, y = CLR_Sum, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "CLR-Weighted Hydrophilic Trait Composition per Tree",
    x = "",
    y = "Relative Trait Abundance",
    fill = "Trait") +
  scale_fill_manual(values=hydro_colors, 
                    name="Trait",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

hydro_tree_plot


## Fill in some trees that have no bars 
# Get list of all trees and traits
all_trees <- unique(hydro_traits$Sample_ID)
all_traits <- unique(hydro_traits$Trait_clean)

# Expand to all Tree × Trait combinations and fill missing with 0
# For positive trait abundance 
hydro_traits_full <- hydro_traits %>%
  select(Tree, Trait_clean, CLR_Sum) %>%
  complete(Tree = all_trees, Trait_clean = all_traits, fill = list(CLR_Sum = 0))

# Add a couple other data columns for plotting later 
# Grab host data 
hosts <- select(sample_data, Tree, Host_ID)

# merge hosts to data 
hydro_traits_full <- merge(hydro_traits_full, hosts, by = "Tree")

## Here these can diverge to consider the traits that are relatively more and less abundant in the tree 
# communities separately 

# Split positive and negative and scale separately 
hydro_traits_pos <- hydro_traits_full %>%
  filter(CLR_Sum > 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(CLR_Sum)) %>%
  ungroup()

hydro_traits_neg <- hydro_traits_full %>%
  filter(CLR_Sum < 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Sum_scaled = CLR_Sum / sum(abs(CLR_Sum))) %>%
  ungroup() %>%
  mutate(CLR_Sum_scaled = -abs(CLR_Sum_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
hydro_traits_diverging <- bind_rows(hydro_traits_pos, hydro_traits_neg)


# Diverging bar plot for trait representation in individual trees 

hydro_diverging <- ggplot(hydro_traits_diverging, aes(x = Tree, y = CLR_Sum_scaled, fill = Trait_clean)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=hydro_colors, 
                    name="Hyphal Hydrophobicity",
                    breaks=c("Hydrophilic", "Hydrophobic"),
                    labels=c("Hydrophilic", "Hydrophobic")) +
  theme_minimal(base_size = 12) +
  labs(x = "Tree",
       y = "Hydrophobicity Abundance",
       fill = "Hydrophobicity") +
  geom_hline(yintercept = 0, color = "red", linewidth = 2) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

hydro_diverging

#################################################################################

#################################### -- 
# (3) ANALYZE FUNCTIONAL COMPOSITION
#################################### -- 

# Get weighted taxon abundance in tree communities 

# Get dataframe for plotting 
tree_abund_df <- as.data.frame(tree_matrix_EM) %>% rownames_to_column(var = "Sample_ID")


# Reshape to long format
tree_abund_long <- tree_abund_df %>%
  pivot_longer(-Sample_ID, names_to = "ASV", values_to = "CLR_Abund")


# Join table of tree environmental data 
x # Pull in the imputed environmental data when I have it 

tree_abund_full <- tree_abund_long %>%
  merge(sample_data, by = "Sample_ID")


# Pull out taxa table from WFDP phyloseq object 
WFDP_tax_table <- phyloseq::tax_table(ps_WFDP_final)

## Merge with taxa data 
WFDP_tax <- as.data.frame(WFDP_tax_table) %>% rownames_to_column(var = "ASV2")

# Merge with all_names file to get shortened ASV names 
WFDP_tax_trim <- merge(WFDP_tax, all_names, by = "ASV2")

WFDP_tax_trim <- select(WFDP_tax_trim, Genus, Species, ASV = ASV2)

tree_tax_full <- merge(WFDP_tax_trim, tree_abund_full, by = "ASV")


# Clean up taxa names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Genus = Genus %>%
           str_replace_all("g__", "")) %>%
  mutate(Species = Species %>%
           str_replace_all("s__", "")
  )

tree_tax_full$Genus <- as.factor(tree_tax_full$Genus)


## Test plotting taxon relative abundance per tree community

# plot per tree 
Tax_tree_plot <- ggplot(tree_tax_full, aes(x = Tree, y = CLR_Abund, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative EM Genus Abundance",
       fill = "Genus"
  ) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.text = element_text(colour="black", size = 9)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

Tax_tree_plot


# Get a species summary just to look at 
species_summary <- tree_tax_full %>%
  group_by(Genus, Species) %>%
  summarise(
    mean_CLR = mean(CLR_Abund),
    max_CLR = max(CLR_Abund),
    .groups = "drop"
  )


## Calculate positive and negative relative abundance for plotting 

# Split positive and negative and scale separately 
all_genera_pos <- tree_tax_full %>%
  filter(CLR_Abund > 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(CLR_Abund)) %>%
  ungroup()

all_genera_neg <- tree_tax_full %>%
  filter(CLR_Abund < 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup() %>%
  mutate(CLR_Abund_scaled = -abs(CLR_Abund_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
all_genera_diverging <- bind_rows(all_genera_pos, all_genera_neg)

# Diverging bar plot for trait representation in individual trees 

all_genera_diverging_plot <- ggplot(all_genera_diverging, aes(x = Tree, y = CLR_Abund_scaled, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative EM Genus Abundance",
       fill = "Genus"
  ) +
  geom_hline(yintercept = 0, color = "red", linewidth = 2) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

all_genera_diverging_plot


#### 41 genera is a lot to look at, and I'm interested in the broadly abundant ones 


## Try grouping rare genera into an "other" category 
tree_tax_full %>%
  filter(CLR_Abund > 0) %>%
  group_by(Genus) %>%
  summarise(tree_count = n_distinct(Sample_ID)) %>%
  arrange(desc(tree_count))

# Find genera to keep that are present for at least 3 host trees 
genera_to_keep <- tree_tax_full %>%
  filter(CLR_Abund > 0) %>%               
  group_by(Genus) %>%
  summarise(tree_count = n_distinct(Sample_ID)) %>% 
  filter(tree_count >= 3) %>%   # Keep genera found in ≥ 3 trees
  pull(Genus)

tree_tax_subset <- tree_tax_full %>%
  mutate(Genus = if_else(Genus %in% genera_to_keep, Genus, "Other"))


# Repeat steps and plot subset 

# Split positive and negative and scale separately 
common_genera_pos <- tree_tax_subset %>%
  filter(CLR_Abund > 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(CLR_Abund)) %>%
  ungroup()

common_genera_neg <- tree_tax_subset %>%
  filter(CLR_Abund < 0) %>%
  group_by(Tree) %>%
  mutate(CLR_Abund_scaled = CLR_Abund / sum(abs(CLR_Abund))) %>%
  ungroup() %>%
  mutate(CLR_Abund_scaled = -abs(CLR_Abund_scaled))  # ensure values are negative

# Merge for plotting 

# Combine
common_genera_diverging <- bind_rows(common_genera_pos, common_genera_neg)

# Diverging bar plot for trait representation in individual trees 

common_genera_diverging_plot <- ggplot(common_genera_diverging, aes(x = Tree, y = CLR_Abund_scaled, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative EM Genus Abundance",
       title = "Commen Genera Relative Abundance",
       fill = "Genus"
  ) +
  geom_hline(yintercept = 0, color = "red", linewidth = 2) +
  theme(legend.title = element_text(colour="black", size=12, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 12))

common_genera_diverging_plot


# -- #### END ##### -- 

