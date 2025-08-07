# -----------------------------------------------------------------------------#
# Exploring fungal taxonomic and functional variation in WFDP for EM community 
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
#  Load in EM fungal taxa table, and the table of functions for genera that     #
#  had successful assignments. Match them up and format into site x species     #
#  matrix for functional diversity analyses with fundiversity package. Also     #
#  perform analyses of the taxonomic composition of the fungal community of     #
#  each host tree across the plot.                                              #
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

asv_count <- as.data.frame(asv_count) #a lot of ASV's no longer present in the samples 

#remove taxa that aren't present in this new subset 
ps_WFDP_final <- subset_taxa(ps_WFDP_final, taxa_sums(ps_WFDP_final) > 0)

# number of taxa - now 404
ntaxa(ps_WFDP_final)


# Get taxon table for EM community in WFDP
WFDP_tax <- tax_table(ps_WFDP_final)

WFDP_tax <- as.data.frame(WFDP_tax) %>% rownames_to_column(var = "X") # 404 taxa

# pull out the genus and family levels, as well as the original OTUs. This will be necessary 
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
matrix$OTU <- matrix$X

# the number is correct - 358 OTUs with genus and functional assignments 


#start names file with the original OTUs 
OTU <- phyloseq::taxa_names(ps_WFDP_final)
OTU_long <- as.data.frame(OTU)

#change OTU names to something nicer to work with
taxa_names(ps_WFDP_final)
n_seqs <- seq(ntaxa(ps_WFDP_final))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_WFDP_final) <- paste("OTU", formatC(n_seqs, 
                                                  width = len_n_seqs, 
                                                  flag = "0"), sep = "_")
taxa_names(ps_WFDP_final)

# get shortened names 
OTU2 <- taxa_names(ps_WFDP_final) 
OTU_short <- as.data.frame(OTU2)

# join two dataframes
all_names <- cbind(OTU_long, OTU_short)


#merge taxa_names with the traits file to get the updated OTU names 
matrix <- merge(matrix, all_names, by = "OTU")

# reformat a tiny bit just to get OTU's as the species in rownames, and my traits only 
traits_WFDP_EM <- matrix %>% select(OTU2, hydrophilic, hydrophobic, ET_contact, ET_contact_short, ET_short, ET_contact_medium,
                                    ET_contact_medium_fringe, ET_contact_medium_smooth, ET_medium_smooth, 
                                    ET_medium_fringe, ET_medium_mat, ET_medium_long, ET_medium_long_smooth,
                                    ET_medium_long_fringe, ET_contact_long_smooth, ET_long) %>% column_to_rownames(var = "OTU2") 

#OTUs are now the row names and there are columns for the binary coding of each trait level 


# We want these to be the ones that are contained within the final phyloseq object 

taxa_to_keep <- matrix$OTU2

# Prune taxa from the phyloseq object that are NOT in the taxa_to_keep list 
ps_WFDP_final <- prune_taxa(taxa_to_keep, ps_WFDP_final)


#SAVE THE FINAL RAW DATA PHYLOSEQ OBJECT 
# For the taxa that have functional assignments in the community 
saveRDS(ps_WFDP_final, file = "~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/WFDP_phyloseq_funcs_final.RDS")



# Getting to know your phyloseq data ####

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



#### 
# Pursuing a transformation because sequence data are inherently compositional, but standard rarefaction 
# is not repeatable and leads to the exclusion of a huge fraction of your data 

# trying a centered log-ratio transformation per https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.02224/full

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

# grab final dataframe of clr transformed OTUs that have trait values 
clr_WFDP_EM <- otu_table(ps_WFDP_final) %>% as("matrix") %>% as.data.frame()

####################################################

####################################
# (1) ANALYZE TAXONOMIC COMPOSITION
####################################

## Beta-diversity can be calculated using Aitchison Distance, which is essentially the euclidean
# distance calculated between pairs of samples that have been transformed by CLR

# calculate Aitchison distance using dist() from base R 

aitchison_WFDP_EM <- dist(clr_WFDP_EM, method = "euclidean")

aitchison_WFDP_EM <- as.matrix(aitchison_WFDP_EM)

# save Aitchison Distance matrix 
save(aitchison_WFDP_EM, file="~/Dropbox/WSU/WFDP_Chapter_3_Project/Fungal_Communities/aitchison_dist_WFDP_EM.Rdata")


######## visualize ##

# PCA is appropriate for euclidean distances 

# get a few site-level columns to load back into the clr dataframe 

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
                  # ABAM        ABGR        ALRU         CONU        TABR          THPL        TSHE        
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

#####assess multivariate homogeneity of sites###########
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

# Not significant here either 


permanova.host <- vegan::adonis2(aitchison_WFDP_EM ~ Host_ID, data = sample_data, method = "euclidean", permutations = 999)
permanova.host

# vegan::adonis2(formula = aitchison_WFDP_EM ~ Host_ID, data = sample_data, permutations = 999, method = "euclidean")
#         Df SumOfSqs      R2      F Pr(>F)
# Model     6    356.5 0.09379 0.9142  0.911
# Residual 53   3444.4 0.90621              
# Total    59   3800.8 1.00000  

# Not significant 


#####################################
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


########################################################################

# I want to look at the traits for each tree community and compare the relative portions of traits that 
# are present using the clr values. 


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
otu_ids <- rownames(traits_WFDP_EM)

traits_matrix_EM <- as.data.frame(traits_WFDP_EM)
traits_matrix_EM[] <- lapply(traits_matrix_EM, as.numeric)
traits_matrix_EM <- as.matrix(traits_matrix_EM)
rownames(traits_matrix_EM) <- otu_ids

# check for any NAs in the data 
sum(is.na(tree_matrix_EM)) # None        
sum(is.na(traits_matrix_EM))  # None   

# Make sure OTUs match between traits and trees
# Find shared OTUs between the two datasets 
shared_otus <- intersect(colnames(tree_matrix_EM), rownames(traits_matrix_EM))

# all are shared 

# set the OTU's to be in the same order 
otu_order <- colnames(tree_matrix_EM)

# reorder the trait matrix rows to match the tree matrix columns
traits_matrix_EM <- traits_matrix_EM[otu_order, ]

# double check alignment 
all(colnames(tree_matrix_EM) == rownames(traits_matrix_EM)) #TRUE

# All match, good to go here 

############# 

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


#################
 
# Can now explore the clr-weighted trait profiles of each host tree! 

# I'm interested in exploration type and hydrophobicity separately, so I want to pull these out of 
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
  geom_hline(yintercept = 0, color = "black", linewidth = 2) +
  theme(legend.title = element_text(colour="black", size=14, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 14))

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
  geom_hline(yintercept = 0, color = "black", linewidth = 2) +
  theme(legend.title = element_text(colour="black", size=14, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11)) +
  theme(axis.text.y = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 14))

hydro_diverging


## Here if I want it



########################################################

# Get weighted taxon abundance in tree communities 

# Get dataframe for plotting 
tree_abund_df <- as.data.frame(tree_matrix_EM) %>% rownames_to_column(var = "Sample_ID")


# Reshape to long format
tree_abund_long <- tree_abund_df %>%
  pivot_longer(-Sample_ID, names_to = "OTU", values_to = "CLR_Abund")


# Join table of tree environmental data 
x # Pull in the imputed environmental data when I have it 

tree_abund_full <- tree_abund_long %>%
  merge(sample_data, by = "Sample_ID")


## Merge with taxa data 
WFDP_tax <- as.data.frame(WFDP_tax_table) %>% rownames_to_column(var = "OTU2")

# Merge with all_names file to get shortened OTU names 
WFDP_tax_trim <- merge(WFDP_tax, all_names, by = "OTU2")

WFDP_tax_trim <- select(WFDP_tax_trim, Genus, Species, OTU = OTU2)

tree_tax_full <- merge(WFDP_tax_trim, tree_abund_full, by = "OTU")


# Clean up taxa names a bit 
tree_tax_full <- tree_tax_full %>%
  mutate(Genus = Genus %>%
           str_replace_all("g__", "")) %>%
  mutate(Species = Species %>%
                    str_replace_all("s__", "")
  )

tree_tax_full$Genus <- as.factor(tree_tax_full$Genus)


## Try grouping rare genera into an "other" category, because 41 genera are too many for clear, 
# descriptive plotting 
tree_tax_full %>%
  filter(CLR_Abund > 0) %>%
  group_by(Genus) %>%
  summarise(tree_count = n_distinct(Sample_ID)) %>%
  arrange(desc(tree_count))



# Find genera to keep that are present for at least 3 host trees 
genera_to_keep <- tree_tax_full %>%
  filter(CLR_Abund > 0) %>%               # Consider genus present if CLR > 0
  group_by(Genus) %>%
  summarise(tree_count = n_distinct(Sample_ID)) %>%  # Count unique trees
  filter(tree_count >= 3) %>%                  # Keep genera found in ≥ 3 trees
  pull(Genus)

tree_tax_subset <- tree_tax_full %>%
  mutate(Genus = if_else(Genus %in% genera_to_keep, Genus, "Other"))

## Test plotting taxon relative abundance per tree community 

# plot per tree 
Tax_tree_plot <- ggplot(tree_tax_subset, aes(x = Tree, y = CLR_Abund, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "",
       y = "Relative EM Genus Abundance",
       fill = "Genus"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
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


##########
# Can come back to this and take the same diverging approach for taxon abundance later if
# it is desired 
##########




##################

########## PAUSE














##### SANDBOX #######


## visualizations of some trait variation 

# take asv_count table and merge the shortened names file 
# move asv out of rownames 

asv_count <- asv_count %>% rownames_to_column(var = "ASV") 

# merge with 'all_names' file 
common_asvs <- merge(asv_count, all_names, by = "ASV")

# Calculate the threshold for the top 25% most common ASVs
threshold <- quantile(common_asvs$asv_count, probs = 0.75)

# Retain only ASVs in the top 25%
common_asvs <- common_asvs %>% filter(asv_count >= threshold)

#this gives 102 as the top 25% most common 


# subset matrix to just most common - can merge 'common_asvs' to 'matrix' using the ASV2 column
matrix_common <- merge(common_asvs, matrix, by = 'ASV2')

# results in 95 ASVs with functional assignments 

# need to fix name format for the Family column to get rid of the 'f_'
matrix_common$Family <- gsub("f__", "", matrix_common$Family)

plotting_traits <- select(matrix_common, Genus, Family, exploration_type, hydro, hydro_binary, ASV2)

# exp type - genus
stack_plot <- ggplot(plotting_traits, aes(x = exploration_type, fill = Genus)) +
  geom_bar(position = "stack") +  # "fill" makes it proportional (100% stacked)
  labs(x = "Exploration Type", y = "Proportion", fill = "Genus") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_fill_viridis_d(name = "Genus", option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 11))

stack_plot

# very very interesting 

# hydro - genus
stack_plot2 <- ggplot(plotting_traits, aes(x = hydro, fill = Genus)) +
  geom_bar(position = "stack") +  # "fill" makes it proportional (100% stacked)
  labs(x = "Hydrophobicity", y = "Proportion", fill = "Genus") +
  scale_fill_viridis_d(name = "Genus", option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 11))

stack_plot2

# exp type - family
stack_plot3 <- ggplot(plotting_traits, aes(x = exploration_type, fill = Family)) +
  geom_bar(position = "stack") +  # "fill" makes it proportional (100% stacked)
  labs(x = "Exploration Type", y = "Proportion", fill = "Family") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_fill_viridis_d(name = "Family", option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), 
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 11))

stack_plot3

# hydro - family 
stack_plot4 <- ggplot(plotting_traits, aes(x = hydro, fill = Family)) +
  geom_bar(position = "stack") +  # "fill" makes it proportional (100% stacked)
  labs(x = "Hydrophobicity", y = "Proportion", fill = "Family") +
  scale_fill_viridis_d(name = "Family", option = "plasma") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11), 
        legend.text = element_text(size = 11))

stack_plot4





