# -----------------------------------------------------------------------------#
# Exploring WFDP tree root and leaf trait data 
# Original Author: L. McKinley Nevins 
# January 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 3.5.1
#                     rstatix v 0.7.2
#                     ggpubr v 0.6.0
#                     fundiversity v 1.1.1
#                     vegan v 2.6.10
#                     FD v 1.0.12.3
#                     rcompanion v 2.5
#                     ggfortify v 0.4.17
#                     gginnards v 0.2.0.1
#                     ggrepel v 0.9.6
#                     corrplot v 0.95
#                     car v 3.1.3
#                     FD 1.0.12.3
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(fundiversity); packageVersion("fundiversity")
library(vegan); packageVersion("vegan")
library(FD); packageVersion("FD")
library(rcompanion); packageVersion("rcompanion")
library(ggfortify); packageVersion("ggfortify")
library(gginnards); packageVersion("gginnards")
library(ggrepel); packageVersion("ggrepel")
library(corrplot); packageVersion("corrplot")
library(car); packageVersion("car")
library(FD); packageVersion("FD")

#################################################################################
#                               Main workflow                                   #
#  Explore the leaf and root trait data. Determine if there are differences     #
#  between species, and maybe try to map some of the variation across the plot? #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
############### 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/"
setwd(wd)

# Load in leaf and root trait datasets 

leaf <- read.csv("WFDP_leaf_traits.csv")

root <- read.csv("WFDP_root_traits.csv")

# combine into one dataset for the trees 
traits <- merge(leaf, root, by = 'code')

# subset just variables of interest, excluding petiole data since it is absent for
# needle-leaf species 

traits <- select(traits, code, WFDP_Code = WFDP_Code.x, sub_plot = sub_plot.x, Host_ID = Host_ID.x, SLA_leaf, LDMC_leaf, LMA_leaf, 
                 leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, leaf_13C, specific_root_length, specific_root_area,
                 root_dry_matter_cont, root_CN, root_15N, root_13C, avg_root_dia, root_pct_N, root_pct_C)

# make species a factor 
traits$Host_ID <- as.factor(traits$Host_ID)

# filter out PSME as a host since there were so few sampled 
traits <- traits %>% filter(Host_ID != 'PSME')


# filter out T-TABR-03 as a host since it had no fungal community data 
traits <- traits %>% filter(code != 'T-TABR-03')

##############################################
## create species x site matrix 

# load in all WFDP environmental data 

x # This can be updated when krieged data is in hand 

env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# grab just sub_plot (referring to the subplot in WFDP - the site level) and Host_ID

env <- select(env, Host_ID, sub_plot)


# create a presence-absence matrix from these 
sites_tree <- env %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Host_ID, values_from = presence, values_fill = list(presence = 0))

print(sites_tree)

sites_tree$sub_plot <- as.factor(sites_tree$sub_plot)

#move subplot into the rownames 
sites_tree <- sites_tree %>%
  column_to_rownames(var = "sub_plot")

##############################################
## create species x trait matrix 

# need species as the rownames and columns to only be the traits 

# remove code column 
traits_tree <- subset(traits, select = -c(code, WFDP_Code, sub_plot))

#######################
# pause and look at variation in raw traits between species 

# Convert to long format for ggplot
long_traits_tree <- pivot_longer(traits_tree, cols = -species, names_to = "trait", values_to = "value")

# Boxplot
trait_plot_OG <- ggplot(long_traits_tree, aes(x = species, y = value, fill = species)) +
  geom_boxplot() +
  facet_wrap(~trait, scales = "free_y") +  # Separate plots for each trait
  theme_minimal() +
  labs(title = "Trait Variation Across Species", y = "Trait Value") +
  theme(legend.position = "none") + # Removes legend since species are on x-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

trait_plot_OG

# correlations between raw trait values 
# can use just the 'traits' data 

# leaf C and N
leaf_CN<- ggplot(traits, aes(x = leaf_pct_N, y = leaf_pct_C, color = species)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Leaf N vs Leaf C", x = "Leaf N (%)", y = "Leaf C (%)")

leaf_CN

# shows clustering of ALRU and CONU away from other species, makes sense since they 
# are the only broadleaf species 


# root C and N
root_CN<- ggplot(traits, aes(x = root_pct_N, y = root_pct_C, color = species)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Root N vs Root C", x = "Root N (%)", y = "Root C (%)")

root_CN

# no visible trends 

######################

#check the spread of the raw trait data 
summary(traits_tree)

#need to get average traits for each species 
traits_tree <- traits_tree %>% group_by(Host_ID) %>%
  summarise(across(everything(), mean))


# traits do not have the same sample sizes, and should be standardized so they 
# are comparable 

#move species into the rownames 
traits_tree <- traits_tree %>%
  column_to_rownames(var = "Host_ID")


#scaling
traits_tree_sc <- scale(traits_tree)
summary(traits_tree_sc)

traits_tree_sc <- as.data.frame(traits_tree_sc)


# #could also scale all traits between zero and 1 
# 
# min_values <- as.numeric(lapply(as.data.frame(traits_tree), min))
# max_values <- as.numeric(lapply(as.data.frame(traits_tree), max))
# 
# traits_tree_minmax <- apply(traits_tree, 1, function(x) {
#   (x - min_values)/(max_values - min_values)
# })
# traits_tree_minmax <- t(traits_tree_minmax)
# summary(traits_tree_minmax)

################################

## Explore the scaled trait data

# make a scaled dataset that doesn't have rownames 

scale_4_plots <- traits_tree_sc %>%
  rownames_to_column(var = "Host_ID")


# backtrack formatting a bit 
# move rownames back into a column 
long_traits_tree_sc <- traits_tree_sc %>%
  rownames_to_column(var = "Host_ID")

# Convert to long format for ggplot
long_traits_tree_sc <- pivot_longer(long_traits_tree_sc, cols = -Host_ID, names_to = "trait", values_to = "value")

# Boxplot
trait_plot <- ggplot(long_traits_tree_sc, aes(x = Host_ID, y = value, fill = Host_ID)) +
  geom_boxplot() +
  facet_wrap(~trait, scales = "free_y") +  # Separate plots for each trait
  theme_minimal() +
  labs(title = "Trait Variation Across Species", y = "Trait Value") +
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

trait_plot

# means, so not showing any spread

#################################################################################

###############
# (2) PRINCIPLE COMPONENT ANALYSIS
###############

# perform a PCA to assess the axes of trait variation for the trees 
    # using unscaled data because the PCA does the scaling internally, and not using means 

# first need to check for autocorrelation in my traits, since many of them are very related

# remove code and species columns to get everything numeric
traits_corr <- subset(traits, select = -c(code, WFDP_Code, sub_plot, Host_ID))

cor_matrix <- cor(traits_corr, use = "everything")
print(cor_matrix)

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)


# dropping a few highly correlated traits would help, and would make the PCA results cleaner 
# making choices to retain certain traits that are most related to the fungal collaboration 
# gradient 

  # leaf_CN - highly correlated with leaf_pct_N
  # root_CN - highly correlated with root_pct_N
  # specific_root_area - highly correlated with specific_root_length, but SRL is a key collab trait 
  # LDMC_leaf - highly correlated with SLA and LMA, also redundant 
  # leaf_pct_C - Highly correlated with SLA and redundant in showing structural investment
  # SLA_leaf - inverse to LMA, so redundant 

# trying with just removing these 3 and we'll see how it goes 

# traits 
traits_nocorr <- subset(traits, select = -c(code, WFDP_Code, sub_plot, Host_ID, leaf_CN, root_CN, specific_root_area, 
                                            LDMC_leaf, leaf_pct_C, SLA_leaf))

#check correlations again 
cor_matrix2 <- cor(traits_nocorr, use = "everything")
print(cor_matrix2)

corrplot(cor_matrix2, method = "color", type = "upper", tl.cex = 0.7)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
# much better

# load Host_ID column back into traits_nocorr
traits_nocorr$Host_ID <- traits$Host_ID
 
###################################
# perform PCA on 'traits_nocorr' dataset
# "traits_nocorr" data has a column for Host_ID in the very last position, so need to exclude that 


#All Traits
alltraits.pca = prcomp(traits_nocorr[1:11], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(traits_nocorr[1:11])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$Host_ID = traits_nocorr$Host_ID
summary(alltraits.pca)

# Save loadings for all traits
write.csv(loadings.alltraits, "./PCA/PCA_loadings_nocorrelatedtraits_tree.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltraits, "./PCA/PCA_scores_nocorrelatedtraits_tree.csv")


# PCA scores are 'scores.alltraits' with column for Host_ID

loadings.alltraits <- as.data.frame(loadings.alltraits)

# get proportion of variance explained to add to each axis label 
pca_var <- alltraits.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage




PCA_plot <- ggplot(scores.alltraits, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  geom_segment(data = loadings.alltraits, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.alltraits, aes(x = PC1 * 10, y = PC2 * 10, label = rownames(loadings.alltraits)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_color_viridis_d(name = "Host_ID", option = "turbo") +  # Viridis colorblind-friendly
  labs(title = "PCA Biplot: Tree Leaf and Root Traits",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Species") +
  theme(legend.position = "right")

PCA_plot

## beautiful 



##Broken-Stick test for the significance of the loadings
print(alltraits.pca)

plot(alltraits.pca, type = "l")

ev = alltraits.pca$sdev^2

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


summary(alltraits.pca)

#can safely retain first two - PC1, PC2, and maybe the third one too 

# figure out which traits were most strongly retained on the first two PC's

# Get top 4 traits for PC1
top_PC1 <- loadings.alltraits[order(abs(loadings.alltraits$PC1), decreasing = TRUE), ][1:4, ]

# Get top 4 traits for PC2
top_PC2 <- loadings.alltraits[order(abs(loadings.alltraits$PC2), decreasing = TRUE), ][1:4, ]


# PC1 is showing a perfect spread of "Do It yourself" and "Outsourcing" with SRL, leaf N, 
# root N, and average root diameter being the strongest associated. 

# PC2 is showing a spread of LMA, root C, and leaf and root 15N. This could suggest that 
# there are differing amounts of investment into leaf and root tissues, and maybe an impact 
# of nitrogen nutrition too 

##### 

# Want to do the PCA to get the position of each individual host tree 

# retain tree code and species columns 
traits2 <- subset(traits, select = -c(leaf_CN, root_CN, specific_root_area, 
                                            LDMC_leaf, leaf_pct_C, SLA_leaf))


# Set up PCA for all trees 
alltrees.pca = prcomp(traits2[3:13], center = T, scale = T)

sd.alltrees = alltrees.pca$sdev
loadings.alltrees = alltrees.pca$rotation
trait.names.alltrees = colnames(traits2[3:13])
scores.alltrees = as.data.frame(alltrees.pca$x)
scores.alltrees$code = traits2$code
scores.alltrees$species = traits2$species
summary(alltrees.pca)

# Save loadings for PCA of all individual tree positions 
write.csv(loadings.alltrees, "./PCA/PCA_loadings_tree_positions.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltrees, "./PCA/PCA_scores_tree_positions.csv")


# PCA scores are 'scores.alltrees' with column for species 

loadings.alltrees <- as.data.frame(loadings.alltrees)

# get proportion of variance explained to add to each axis label 
pca_var <- alltrees.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

# Visualize
PCA_plot_trees <- ggplot(scores.alltrees, aes(x = PC1, y = PC2, color = species)) +
  geom_point(size = 3) +
  geom_segment(data = loadings.alltrees, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.alltrees, aes(x = PC1 * 10, y = PC2 * 10, label = rownames(loadings.alltrees)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_color_viridis_d(name = "Species", option = "turbo") +  # Viridis colorblind-friendly
  labs(title = "PCA Biplot: Tree Leaf and Root Traits",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Species") +
  theme(legend.position = "right")

PCA_plot_trees


## Broken stick will be the same, so don't need to do anything with that 

# Could get the eigenvectors for each host tree's position along these axes, which could 
# summarize its position in this functional space 

# Then could regress these with environmental properties in the plot. 
# Could also look at fungal community composition dispersion to see if this is 
# related to the traits of the host tree 

# Grab PC1 and PC2 for each tree 
alltrees_PCs <- select(scores.alltrees, PC1, PC2, code, species)








#################################################################################

###############
# (3) FUNCTIONAL ANALYSES
###############

## Functional analyses using the fundiversity package 
# following this vignette: https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
###############

# not using any species x site matrix, as I'm not making comparisons between sites (all data are 
# from WFDP). If I want to in the future I can make comparisons across some category of 
# environmental variable and perhaps use that as the 'site' 

# Remember that trait values are now scaled, so they can be compared in relation to eachother 
# but the values themselves do not have the same biological meaning 

# Test FD package 
disper <- dbFD(traits_tree_sc, calc.FRic = TRUE, calc.CWM = FALSE)

# using just the traits this gives a value to summarize the whole forest 






