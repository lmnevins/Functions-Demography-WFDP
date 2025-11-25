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
#                     rcompanion v 2.5.0
#                     ggfortify v 0.4.17
#                     gginnards v 0.2.0.1
#                     ggrepel v 0.9.6
#                     corrplot v 0.95
#                     car v 3.1.3
#                     FD 1.0.12.3
#                     multcomp v 1.4.28
#                     multcompView v 0.1.10
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
library(multcomp); packageVersion("multcomp")
library(multcompView); packageVersion("multcompView")

#################################################################################
#                               Main workflow                                   #
#  Explore the leaf and root trait data. Determine if there are differences     #
#  between species, and maybe try to map some of the variation across the plot? #
#                                                                               #
#################################################################################

############### -- 
# (1) DATA PREP
############### -- 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/"
setwd(wd)

# Load in trait datasets and clean a bit 

##Leaves 
leaf <- read.csv("WFDP_leaf_traits.csv")

# make species a factor 
leaf$Host_ID <- as.factor(leaf$Host_ID)

# filter out PSME as a host since there were so few sampled 
leaf <- leaf %>% filter(Host_ID != 'PSME')

# filter out T-TABR-03 as a host since it had no fungal community data 
leaf <- leaf %>% filter(code != 'T-TABR-03')

## Roots
root <- read.csv("WFDP_root_traits.csv")

root$Host_ID <- as.factor(root$Host_ID)

root <- root %>% filter(Host_ID != 'PSME')
root <- root %>% filter(code != 'T-TABR-03')


# combine into one dataset for the trees 
traits <- merge(leaf, root, by = 'code')

# subset just variables of interest, excluding petiole data since it is absent for
# needle-leaf species 

traits <- select(traits, code, WFDP_Code = WFDP_Code.x, sub_plot = sub_plot.x, Host_ID = Host_ID.x, SLA_leaf, LDMC_leaf, LMA_leaf, 
                 leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, leaf_13C, specific_root_length, specific_root_area,
                 root_dry_matter_cont, root_CN, root_15N, root_13C, avg_root_dia, root_pct_N, root_pct_C)

############################################## -- 
## create species x site matrix 

# load in all WFDP environmental data 

x # This can be updated when krieged data is in hand 

env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# grab just sub_plot (referring to the subplot in WFDP - the site level) and Host_ID

env$sub_plot <- env$Cell

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

############################################## -- 
## create species x trait matrix 

# need species as the rownames and columns to only be the traits 

# remove code column 
traits_tree <- subset(traits, select = -c(code, WFDP_Code, sub_plot))



#################################################################################

########################################## -- 
# (2) ASSESS VARIATION IN RAW TRAIT DATA
########################################## -- 

# pause and look at variation in raw traits between species 

# Convert to long format for ggplot
long_traits_tree <- pivot_longer(traits_tree, cols = -Host_ID, names_to = "trait", values_to = "value")

# Boxplot
trait_plot_OG <- ggplot(long_traits_tree, aes(x = Host_ID, y = value, fill = Host_ID)) +
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
leaf_CN<- ggplot(traits, aes(x = leaf_pct_N, y = leaf_pct_C, color = Host_ID)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Leaf N vs Leaf C", x = "Leaf N (%)", y = "Leaf C (%)")

leaf_CN

# shows clustering of ALRU and CONU away from other species, makes sense since they 
# are the only broadleaf species 


# root C and N
root_CN<- ggplot(traits, aes(x = root_pct_N, y = root_pct_C, color = Host_ID)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Root N vs Root C", x = "Root N (%)", y = "Root C (%)")

root_CN

# no visible trends 

####

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


# make a scaled dataset that doesn't have rownames 

scale_4_plots <- traits_tree_sc %>%
  rownames_to_column(var = "Host_ID")


# backtrack formatting a bit 
# move rownames back into a column 
long_traits_tree_sc <- traits_tree_sc %>%
  rownames_to_column(var = "Host_ID")

# Convert to long format for ggplot
long_traits_tree_sc <- pivot_longer(long_traits_tree_sc, cols = -Host_ID, names_to = "trait", values_to = "value")


#################################################################################

##################################### -- 
# (3) PRINCIPLE COMPONENT ANALYSIS
##################################### -- 

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

# load Host_ID and WFDP_code columns back into traits_nocorr
traits_nocorr$Host_ID <- traits$Host_ID

traits_nocorr$WFDP_Code <- traits$WFDP_Code

 
################################### -- 
# perform PCA on 'traits_nocorr' dataset
# "traits_nocorr" data has a column for Host_ID and WFDP_code in the last two positions, 
# so need to exclude those 

#All Traits and individual tree positions 
alltraits.pca = prcomp(traits_nocorr[1:11], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(traits_nocorr[1:11])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$WFDP_Code = traits_nocorr$WFDP_Code
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


# set colors for hosts 
                  # ABAM        ABGR        ALRU         CONU        TABR          THPL        TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66")


PCA_plot <- ggplot(scores.alltraits, aes(x = PC1, y = PC2, color = Host_ID)) +
  geom_point(size = 3) +
  geom_segment(data = loadings.alltraits, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.alltraits, aes(x = PC1 * 10, y = PC2 * 10, label = rownames(loadings.alltraits)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_colour_manual(values=all_hosts, 
                      name="Host Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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

###

# Could get the eigenvectors for each host tree's position along these axes, which could 
# summarize its position in this functional space 

# Then could regress these with environmental properties in the plot. 
# Could also look at fungal community composition dispersion to see if this is 
# related to the traits of the host tree 

# Grab PC1 and PC2 for each tree 
alltrees_PCs <- select(scores.alltraits, PC1, PC2, WFDP_Code, Host_ID)

# Save this file for later 
write.csv(alltrees_PCs, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/tree_PC_scores.csv")


################################### -- 

# Perform separate PCAs for leaf and root traits 

# Grab environmental data again for some association labels 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

env <- dplyr::select(env, WFDP_Code, Association)

## Use the original leaf and root trait datasets 

## Leaves first: 

# Subset to traits excluding petiole 
leaf_sub <- select(leaf, code, WFDP_Code, sub_plot, Host_ID, SLA_leaf, LDMC_leaf, LMA_leaf, 
                   leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, leaf_13C)

# remove code and species columns to get everything numeric
leaf_traits_corr <- subset(leaf_sub, select = -c(code, WFDP_Code, sub_plot, Host_ID))

leaf_cor_matrix <- cor(leaf_traits_corr, use = "everything")
print(leaf_cor_matrix)

corrplot(leaf_cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# Try again removing strongly related traits 
leaf_traits_corr2 <- subset(leaf_traits_corr, select = -c(leaf_CN, LMA_leaf, LDMC_leaf))

leaf_cor_matrix2 <- cor(leaf_traits_corr2, use = "everything")
print(leaf_cor_matrix2)
corrplot(leaf_cor_matrix2, method = "color", type = "upper", tl.cex = 0.7)

# This looks better, a couple are still related but we'll keep them for now 

# load Host_ID and WFDP_code columns back into leaf_traits_corr2
leaf_traits_corr2$Host_ID <- leaf$Host_ID
leaf_traits_corr2$WFDP_Code <- leaf$WFDP_Code

# Then roots: 
root_sub <- select(root, code, WFDP_Code, sub_plot, Host_ID, specific_root_length, 
                   specific_root_area, root_dry_matter_cont, root_CN, root_15N, 
                   root_13C, avg_root_dia, root_pct_N, root_pct_C)


# remove code and species columns to get everything numeric
root_traits_corr <- subset(root_sub, select = -c(code, WFDP_Code, sub_plot, Host_ID))

root_cor_matrix <- cor(root_traits_corr, use = "everything")
print(root_cor_matrix)

corrplot(root_cor_matrix, method = "color", type = "upper", tl.cex = 0.7)

# Try again removing strongly related traits 
root_traits_corr2 <- subset(root_traits_corr, select = -c(root_CN, specific_root_area))

root_cor_matrix2 <- cor(root_traits_corr2, use = "everything")
print(root_cor_matrix2)
corrplot(root_cor_matrix2, method = "color", type = "upper", tl.cex = 0.7)


# This looks better, a couple are still related but we'll keep them for now 

# load Host_ID and WFDP_code columns back into root_traits_corr2
root_traits_corr2$Host_ID <- root$Host_ID
root_traits_corr2$WFDP_Code <- root$WFDP_Code

################################### -- 

# PCAs

## Leaves
leaf.pca = prcomp(leaf_traits_corr2[1:5], center = T, scale = T)

sd.leaf = leaf.pca$sdev
loadings.leaf = leaf.pca$rotation
trait.names.leaf = colnames(leaf_traits_corr2[1:5])
scores.leaf = as.data.frame(leaf.pca$x)
scores.leaf$WFDP_Code = leaf_traits_corr2$WFDP_Code
scores.leaf$Host_ID = leaf_traits_corr2$Host_ID
summary(leaf.pca)

# Save loadings for leaf traits
write.csv(loadings.leaf, "./PCA/PCA_loadings_leaf_traits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.leaf, "./PCA/PCA_scores_leaf_traits.csv")


# PCA scores are 'scores.leaf' with column for Host_ID

loadings.leaf <- as.data.frame(loadings.leaf)

# get proportion of variance explained to add to each axis label 
pca_var <- leaf.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

#Merge in mycorrhizal association for plotting 
scores.leaf <- merge(scores.leaf, env, by = "WFDP_Code")

# Change loadings names to something cleaner 
new_loadings <- c("SLA", "PctN", "PctC", "d15N", "d13C")
rownames(loadings.leaf) <- new_loadings

# Visualize
PCA_plot_leaf <- ggplot(scores.leaf, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_segment(data = loadings.leaf, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.leaf, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.leaf)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("EM" = "#FFC20A", "Both" = "#0C7BDC"), name = "Mycorrhizal Association") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Host Species") +
  labs(title = "Tree Leaf Traits",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Species") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PCA_plot_leaf

  
# Get top 3 traits for PC1
top_PC1_leaf <- loadings.leaf[order(abs(loadings.leaf$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_leaf <- loadings.leaf[order(abs(loadings.leaf$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in SLA and leaf C and N content, which matches what 
# would be expected for broad leaf vs needle-leaf species 

# PC2 is being driven by the d13C, d15N, and leaf C content. 

## Together the axes explain 81.8% of the variation 

### Calculate average PC1 score for each species and plot, compare statistically 

PC1_leaf <- dplyr::select(scores.leaf, PC1, WFDP_Code, Host_ID)

PC1_leaf_summary <- PC1_leaf %>%
  group_by(Host_ID) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1   = sd(PC1, na.rm = TRUE),
    n        = n(),
    se_PC1   = sd_PC1 / sqrt(n)
  )

# Add in mycorrhizal associations 

PC1_leaf_summary$Association <- c("EM", "EM", "Both", "Both", "Both", "Both", "EM")


# Test for significant differences between species 
aov_PC1_leaf <- aov(PC1 ~ Host_ID, data = PC1_leaf)
summary(aov_PC1_leaf)

tuk_PC1_leaf <- TukeyHSD(aov_PC1_leaf)
tuk_PC1_leaf

# Significant differences between species, essentially between 
# ALRU/CONU and all others 

# extract pairwise p-values
tuk_cld <- multcompLetters4(aov_PC1_leaf, tuk_PC1_leaf)

# combine letters with summary table
PC1_leaf_summary$letters <- tuk_cld$Host_ID %>% 
  as.data.frame() %>% 
  pull(Letters)

# Visualize 
PC1_leaf_plot <- ggplot(PC1_leaf_summary, aes(x = Host_ID, y = mean_PC1, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("EM" = "#FFC20A", "Both" = "#0C7BDC"), , name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 1 - LES") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC1_leaf_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 


### Calculate average PC2 score for each species and plot, compare statistically 

PC2_leaf <- dplyr::select(scores.leaf, PC2, WFDP_Code, Host_ID)

PC2_leaf_summary <- PC2_leaf %>%
  group_by(Host_ID) %>%
  summarise(
    mean_PC2 = mean(PC2, na.rm = TRUE),
    sd_PC2   = sd(PC2, na.rm = TRUE),
    n        = n(),
    se_PC2   = sd_PC2 / sqrt(n)
  )

# Add in mycorrhizal associations 

PC2_leaf_summary$Association <- c("EM", "EM", "Both", "Both", "Both", "Both", "EM")


# Test for significant differences between species 
aov_PC2_leaf <- aov(PC2 ~ Host_ID, data = PC2_leaf)
summary(aov_PC2_leaf)

tuk_PC2_leaf <- TukeyHSD(aov_PC2_leaf)
tuk_PC2_leaf

# Significant differences between species

# TABR-ABAM < 0.001
# THPL-ABAM < 0.001
# TABR-ABGR < 0.001
# THPL-ABGR < 0.001
# TABR-ALRU < 0.001
# TABR-CONU < 0.001
# THPL-TABR 0.004
# TSHE-TABR < 0.001
# TSHE-THPL < 0.001

# extract pairwise p-values
tuk_cld_PC2_leaf <- multcompLetters4(aov_PC2_leaf, tuk_PC2_leaf)

# combine letters with summary table
PC2_leaf_summary$letters <- tuk_cld_PC2_leaf$Host_ID %>% 
  as.data.frame() %>% 
  pull(Letters)

# Visualize 
PC2_leaf_plot <- ggplot(PC2_leaf_summary, aes(x = Host_ID, y = mean_PC2, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("EM" = "#FFC20A", "Both" = "#0C7BDC"), , name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 2 - LES") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC2_leaf_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 


############### -- 

## Roots
root.pca = prcomp(root_traits_corr2[1:7], center = T, scale = T)

sd.root = root.pca$sdev
loadings.root = root.pca$rotation
trait.names.root = colnames(root_traits_corr2[1:7])
scores.root = as.data.frame(root.pca$x)
scores.root$WFDP_Code = root_traits_corr2$WFDP_Code
scores.root$Host_ID = root_traits_corr2$Host_ID
summary(root.pca)

# Save loadings for root traits
write.csv(loadings.root, "./PCA/PCA_loadings_root_traits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.root, "./PCA/PCA_scores_root_traits.csv")


# PCA scores are 'scores.root' with column for Host_ID

loadings.root <- as.data.frame(loadings.root)

# Change loadings names to something cleaner 
new_loadings_root <- c("SRL", "RDMC", "d15N", "d13C", "RD", "PctN", "PctC")
rownames(loadings.root) <- new_loadings_root

#Merge in mycorrhizal association for plotting 
scores.root <- merge(scores.root, env, by = "WFDP_Code")

# get proportion of variance explained to add to each axis label 
pca_var <- root.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

# Visualize
PCA_plot_root <- ggplot(scores.root, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_segment(data = loadings.root, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.root, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.root)),
                  color = "black", size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("EM" = "#FFC20A", "Both" = "#0C7BDC"), name = "Mycorrhizal Association") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Host Species") +
  labs(title = "Tree Root Traits",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Species") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PCA_plot_root


# Get top 3 traits for PC1
top_PC1_root <- loadings.root[order(abs(loadings.root$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_root <- loadings.root[order(abs(loadings.root$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in Specific root length, root N content, and avg root 
# diameter. 

# PC2 is being driven by the d13C, root dry matter content, and root C content. 

## Together the axes explain 59.9% of the variation 


### Calculate average PC1 score for each species and plot, compare statistically 

PC1_root <- dplyr::select(scores.root, PC1, WFDP_Code, Host_ID)

PC1_root_summary <- PC1_root %>%
  group_by(Host_ID) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1   = sd(PC1, na.rm = TRUE),
    n        = n(),
    se_PC1   = sd_PC1 / sqrt(n)
  )

# Add in mycorrhizal associations 
PC1_root_summary$Association <- c("EM", "EM", "Both", "Both", "Both", "Both", "EM")


# Test for significant differences between species 
aov_PC1_root <- aov(PC1 ~ Host_ID, data = PC1_root)
summary(aov_PC1_root)

tuk_PC1_root <- TukeyHSD(aov_PC1_root)
tuk_PC1_root

# Significant differences between species, essentially between 
# ALRU-ABAM 0.02
# ALRU-ABGR 0.006
# CONU-ALRU 0.011
# TABR-ALRU 0.0006
# TSHE-ALRU 0.0005
# THPL-TABR 0.03
# TSHE-THPL 0.028

# extract pairwise p-values
tuk_cld_root <- multcompLetters4(aov_PC1_root, tuk_PC1_root)

# combine letters with summary table
PC1_root_summary$letters <- tuk_cld_root$Host_ID %>% 
  as.data.frame() %>% 
  pull(Letters)

# Visualize 
PC1_root_plot <- ggplot(PC1_root_summary, aes(x = Host_ID, y = mean_PC1, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("EM" = "#FFC20A", "Both" = "#0C7BDC"), , name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 1 - RES") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC1_root_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 


### Calculate average PC2 score for each species and plot, compare statistically 

PC2_root <- dplyr::select(scores.root, PC2, WFDP_Code, Host_ID)

PC2_root_summary <- PC2_root %>%
  group_by(Host_ID) %>%
  summarise(
    mean_PC2 = mean(PC2, na.rm = TRUE),
    sd_PC2   = sd(PC2, na.rm = TRUE),
    n        = n(),
    se_PC2   = sd_PC2 / sqrt(n)
  )

# Add in mycorrhizal associations 
PC2_root_summary$Association <- c("EM", "EM", "Both", "Both", "Both", "Both", "EM")


# Test for significant differences between species 
aov_PC2_root <- aov(PC2 ~ Host_ID, data = PC2_root)
summary(aov_PC2_root)

tuk_PC2_root <- TukeyHSD(aov_PC2_root)
tuk_PC2_root

# NO significant differences between species

# Visualize 
PC2_root_plot <- ggplot(PC2_root_summary, aes(x = Host_ID, y = mean_PC2, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c("EM" = "#FFC20A", "Both" = "#0C7BDC"), , name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 2 - RES") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11)) +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11))

PC2_root_plot

## Need to add significance values to this but they are a bit complicated, 
# so can come back and do this 


#################################################################################

############################################ -- 
# (4) ENVIRONMENTAL RELATIONSHIPS WITH PCAs
############################################ -- 

# Read back in environmental data for the trees 

x # Right now this is just slope, aspect, and elevation, but will be able to pull in 
 # krieged data when it exists 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# subset to data of interest 

env <- select(env, Cell, WFDP_Code, slope, aspect, elevation_m)

tree_PCs_env <- merge(env, alltrees_PCs, by = "WFDP_Code")

#### Assess relationships between PC1 and PC2 values and slope, aspect, and elevation 

### PC1

## Slope
slope <- ggplot(tree_PCs_env, aes(x = slope, y = PC1)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope

# test relationships 
lm_slope <- lm(PC1 ~ slope, data = tree_PCs_env)
summary(lm_slope) # NOT SIGNIFICANT


## Aspect
aspect <- ggplot(tree_PCs_env, aes(x = aspect, y = PC1)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect

# test relationships 
lm_aspect <- lm(PC1 ~ aspect, data = tree_PCs_env)
summary(lm_aspect) # NOT SIGNIFICANT


## Elevation
elev <- ggplot(tree_PCs_env, aes(x = elevation_m, y = PC1)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
elev

# test relationships 
lm_elev <- lm(PC1 ~ elevation_m, data = tree_PCs_env)
summary(lm_elev) # NOT SIGNIFICANT


### PC2

## Slope
slope2 <- ggplot(tree_PCs_env, aes(x = slope, y = PC2)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope2

# test relationships 
lm_slope2 <- lm(PC2 ~ slope, data = tree_PCs_env)
summary(lm_slope2) # SIGNIFICANT

# Adjusted R-squared:  0.05096, p-value: 0.04574


## Aspect
aspect2 <- ggplot(tree_PCs_env, aes(x = aspect, y = PC2)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect2

# test relationships 
lm_aspect2 <- lm(PC1 ~ aspect, data = tree_PCs_env)
summary(lm_aspect2) # NOT SIGNIFICANT


## Elevation
elev2 <- ggplot(tree_PCs_env, aes(x = elevation_m, y = PC2)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
elev2

# test relationships 
lm_elev2 <- lm(PC2 ~ elevation_m, data = tree_PCs_env)
summary(lm_elev2) # NOT SIGNIFICANT


## RESULT: Only significant relationship was between slope and PC2 values 

#################################################################################

############### -- 
# (4) FUNCTIONAL ANALYSES
############### -- 

## Functional analyses using the fundiversity package 
# following this vignette: https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
############### -- 

# not using any species x site matrix, as I'm not making comparisons between sites (all data are 
# from WFDP). If I want to in the future I can make comparisons across some category of 
# environmental variable and perhaps use that as the 'site' 

# Remember that trait values are now scaled, so they can be compared in relation to eachother 
# but the values themselves do not have the same biological meaning 

# Test FD package 
disper <- dbFD(traits_tree_sc, calc.FRic = TRUE, calc.CWM = FALSE)

# using just the traits this gives a value to summarize the whole forest 






