# -----------------------------------------------------------------------------#
# Exploring WFDP tree root and leaf trait data 
# Original Author: L. McKinley Nevins 
# January 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 4.0.1
#                     rstatix v 0.7.3
#                     vegan v 2.7.2
#                     rcompanion v 2.5.1
#                     ggfortify v 0.4.19
#                     gginnards v 0.2.0.2
#                     ggrepel v 0.9.6
#                     corrplot v 0.95
#                     car v 3.1.3
#                     multcomp v 1.4.29
#                     multcompView v 0.1.10
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(vegan); packageVersion("vegan")
library(rcompanion); packageVersion("rcompanion")
library(ggfortify); packageVersion("ggfortify")
library(gginnards); packageVersion("gginnards")
library(ggrepel); packageVersion("ggrepel")
library(corrplot); packageVersion("corrplot")
library(car); packageVersion("car")
library(multcomp); packageVersion("multcomp")
library(multcompView); packageVersion("multcompView")

#################################################################################
#                               Main workflow                                   #
#  Explore the leaf and root trait data. Determine if there are differences     #
#  between species, and major PCA axes for the separate and combined traits.    #
#  Explore relationships between traits typically considered to be equivalent   #
#  above and belowground.                                                       #
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

traits <- dplyr::select(traits, code, WFDP_Code = WFDP_Code.x, sub_plot = sub_plot.x, Host_ID = Host_ID.x, SLA_leaf, LDMC_leaf, LMA_leaf, 
                 leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, leaf_13C, specific_root_length, specific_root_area,
                 root_dry_matter_cont, root_CN, root_15N, root_13C, avg_root_dia, root_pct_N, root_pct_C)

############################################## -- 
## create species x site matrix 

# load in all WFDP environmental data 

x # This can be updated when krieged data is in hand 

env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# grab just sub_plot (referring to the subplot in WFDP - the site level) and Host_ID

env$sub_plot <- env$Cell

env <- dplyr::select(env, Host_ID, sub_plot)


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

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# pause and look at variation in raw traits between species 


# Grab environmental data again for some association labels 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

env <- dplyr::select(env, WFDP_Code, Association, Host_ID)


#subset a few columns from traits 

traits_sub <- dplyr::select(traits, -c("code", "WFDP_Code", "sub_plot"))


# Convert to long format for ggplot
long_traits <- pivot_longer(traits_sub, cols = -(Host_ID), names_to = "trait", values_to = "value")

# Boxplot
trait_plot_OG <- ggplot(long_traits, aes(x = Host_ID, y = value, fill = Host_ID)) +
  geom_boxplot() +
  facet_wrap(~trait, scales = "free_y") +  # Separate plots for each trait
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Host Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(title = "Trait Variation Across Species", y = "Trait Value") +
 theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 11, colour="black"))

trait_plot_OG


# Take a closer look at just the isotopic results 
isotopes <- dplyr::select(traits_sub, Host_ID, leaf_13C, leaf_15N, root_13C, root_15N)

# Merge with env to get associations 
isotopes$Association <- env$Association

# Convert to long format for ggplot
isotopes <- pivot_longer(isotopes, cols = -c(Host_ID, Association), names_to = "trait", values_to = "value")

# Boxplot
isotopes_plot <- ggplot(isotopes, aes(x = Host_ID, y = value, fill = Association)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\n Association") + 
  facet_wrap(~trait, scales = "free_y") +  # Separate plots for each trait
  theme_bw() +
  labs(title = "", y = expression("Isotope per mil"), x = "") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black"))

isotopes_plot

## Can do some stats here if desired to test for differences between the species for 
# these values 


################################### -- 
# (2) EXPLORE TRAIT RELATIONSHIPS 
################################### -- 

# Just doing visualizations right now but could test for significant trends if desired 


# can use just the 'traits' data 
traits <- merge(traits, env, by = "WFDP_Code")

# lazy fix for naming 
traits$Host_ID <- traits$Host_ID.x


# leaf C and N
leaf_CN <- ggplot(traits, aes(x = leaf_pct_N, y = leaf_pct_C, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

leaf_CN

# shows clustering of ALRU and CONU away from other species, makes sense since they 
# are the only broadleaf species 


# root C and N
root_CN <- ggplot(traits, aes(x = root_pct_N, y = root_pct_C, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "", x = "Root N (%)", y = "Root C (%)") + 
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

root_CN

# no visible trends 


# whole plant C 

plant_C <- ggplot(traits, aes(x = leaf_pct_C, y = root_pct_C, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "", x = "Leaf C (%)", y = "Root C (%)") + 
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

plant_C


# whole plant N 

plant_N <- ggplot(traits, aes(x = leaf_pct_N, y = root_pct_N, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "", x = "Leaf N (%)", y = "Root N (%)") + 
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

plant_N


# SLA x SRL 

SLA_SRL <- ggplot(traits, aes(x = SLA_leaf, y = specific_root_length, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "", x = "SLA (mm2/mg)", y = "SRL (cm/mg)") + 
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

SLA_SRL


# LDMC x RDMC

LDMC_RDMC <- ggplot(traits, aes(x = LDMC_leaf, y = root_dry_matter_cont, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "", x = "LDMC (mg/g)", y = "RDMC (mg/g)") + 
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

LDMC_RDMC




## Can do more here, and could switch these to be elipses instead to show more species variation, 
# rather than doing it by mycorrhizal association. 



###### OPTIONAL #########

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

# "traits" data has info columns in the first four positions, so need to exclude those 

#All Traits and individual tree positions 
alltraits.pca = prcomp(traits[5:21], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(traits[5:21])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$WFDP_Code = traits$WFDP_Code
summary(alltraits.pca)

# Save loadings for all traits
write.csv(loadings.alltraits, "./PCA/PCA_loadings_alltraits_tree.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltraits, "./PCA/PCA_scores_alltraits_tree.csv")


# PCA scores are 'scores.alltraits' with column for Host_ID

loadings.alltraits <- as.data.frame(loadings.alltraits)

# get proportion of variance explained to add to each axis label 
pca_var <- alltraits.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

#Merge in mycorrhizal association for plotting 
scores.alltraits <- merge(scores.alltraits, env, by = "WFDP_Code")

# Change loadings names to something cleaner 
new_loadings <- c("SLA", "LDMC", "LMA", "Leaf_PctN", "Leaf_PctC", "Leaf_CN", "Leaf_d15N", 
                  "Leaf_d13C", "SRL", "SRA", "RDMC", "Root_CN", 
                  "Root_d15N", "Root_d13C", "RD", "Root_PctN", "Root_PctC")
rownames(loadings.alltraits) <- new_loadings

PCA_plot <- ggplot(scores.alltraits, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_segment(data = loadings.alltraits, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.alltraits, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.alltraits)),
                  color = "black", size = 5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Focal Species") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 14)) + 
  theme(
    axis.text.x = element_text(size = 14, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 14, colour="black"))

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


# PC1 is showing a spread of Leaf_PctN, Leaf_CN, SLA, and LDMC. This reflects conservative vs 
# acquisitive strategies, with more conservative needle-leaf species with more negative values, 
# and more acquisitive broadleaf species with more positive values. 

# PC2 is showing a spread of SRA, SRL, Root_CN, and Root_PctN. This is pretty much just the equivalent 
# of PC1 reflecting leaf strategies. This also shows a "Do-it-yourself" vs "Outsourcing" strategies trade-off. 

###

# Could get the eigenvectors for each host tree's position along these axes, which could 
# summarize its position in this functional space 

# Then could regress these with environmental properties in the plot. 
# Could also look at fungal community composition dispersion to see if this is 
# related to the traits of the host tree 

# Grab PC1 and PC2 for each tree 
alltrees_PCs <- dplyr::select(scores.alltraits, PC1, PC2, WFDP_Code, Host_ID)

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
leaf_sub <- dplyr::select(leaf, code, WFDP_Code, sub_plot, Host_ID, SLA_leaf, LDMC_leaf, LMA_leaf, 
                   leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, leaf_13C)

# Then roots: 
root_sub <- dplyr::select(root, code, WFDP_Code, sub_plot, Host_ID, specific_root_length, 
                   specific_root_area, root_dry_matter_cont, root_CN, root_15N, 
                   root_13C, avg_root_dia, root_pct_N, root_pct_C)

################################### -- 

# PCAs

## Leaves
leaf.pca = prcomp(leaf_sub[5:12], center = T, scale = T)

sd.leaf = leaf.pca$sdev
loadings.leaf = leaf.pca$rotation
trait.names.leaf = colnames(leaf_sub[5:12])
scores.leaf = as.data.frame(leaf.pca$x)
scores.leaf$WFDP_Code = leaf_sub$WFDP_Code
scores.leaf$Host_ID = leaf_sub$Host_ID
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
new_loadings <- c("SLA", "LDMC", "LMA", "PctN", "PctC", "C:N", "d15N", "d13C")
rownames(loadings.leaf) <- new_loadings


## Removing legend from all plots to save for formatting 


# Visualize
PCA_plot_leaf <- ggplot(scores.leaf, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3, aes(shape = Host_ID)) +
  geom_segment(data = loadings.leaf, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.leaf, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.leaf)),
                  color = "black", size = 5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Focal Species") +
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Species") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
  axis.text.x = element_text(size = 16, colour="black"),
  axis.text.y = element_text(size = 16, colour="black"),
  axis.title.y = element_text(size = 16, colour="black"),
  axis.title.x = element_text(size = 16, colour="black"), 
  strip.text = element_text(size = 16, colour="black"))

PCA_plot_leaf

  
# Get top 3 traits for PC1
top_PC1_leaf <- loadings.leaf[order(abs(loadings.leaf$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_leaf <- loadings.leaf[order(abs(loadings.leaf$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in SLA, C:N ratio, and LMA, which matches what 
# would be expected for broad leaf vs needle-leaf species 

# PC2 is being driven by the d13C, d15N, and leaf C content. 

## Together the axes explain 83% of the variation 

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

PC1_leaf_summary$Association <- c("EM", "EM", "DUAL", "AM", "DUAL", "AM", "EM")


# Test for significant differences between species 
aov_PC1_leaf <- aov(PC1 ~ Host_ID, data = PC1_leaf)
summary(aov_PC1_leaf)

tuk_PC1_leaf <- TukeyHSD(aov_PC1_leaf)
tuk_PC1_leaf

# Significant differences between species

# 
# Df Sum Sq Mean Sq F value Pr(>F)    
# Host_ID      6 278.23   46.37   73.22 <2e-16 ***
#   Residuals   53  33.57    0.63                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# diff        lwr        upr     p adj
# ABGR-ABAM  0.16166915 -1.0958287  1.4191670 0.9996811
# ALRU-ABAM  5.68781200  4.3540341  7.0215898 0.0000000
# CONU-ABAM  5.58953827  4.5005132  6.6785633 0.0000000
# TABR-ABAM  2.11737052  0.9985036  3.2362374 0.0000076
# THPL-ABAM  0.70377152 -0.3852535  1.7927966 0.4389146
# TSHE-ABAM  1.31822662  0.2292016  2.4072517 0.0084988
# ALRU-ABGR  5.52614285  4.0515960  7.0006897 0.0000000
# CONU-ABGR  5.42786911  4.1703713  6.6853669 0.0000000
# TABR-ABGR  1.95570137  0.6722730  3.2391297 0.0004008
# THPL-ABGR  0.54210237 -0.7153955  1.7996002 0.8390672
# TSHE-ABGR  1.15655747 -0.1009403  2.4140553 0.0905886
# CONU-ALRU -0.09827373 -1.4320516  1.2355041 0.9999878
# TABR-ALRU -3.57044148 -4.9286944 -2.2121886 0.0000000
# THPL-ALRU -4.98404048 -6.3178183 -3.6502626 0.0000000
# TSHE-ALRU -4.36958538 -5.7033632 -3.0358075 0.0000000
# TABR-CONU -3.47216775 -4.5910346 -2.3533009 0.0000000
# THPL-CONU -4.88576675 -5.9747918 -3.7967417 0.0000000
# TSHE-CONU -4.27131164 -5.3603367 -3.1822866 0.0000000
# THPL-TABR -1.41359900 -2.5324659 -0.2947321 0.0052253
# TSHE-TABR -0.79914390 -1.9180108  0.3197230 0.3188954
# TSHE-THPL  0.61445510 -0.4745700  1.7034802 0.6002402



# Visualize 
PC1_leaf_plot <- ggplot(PC1_leaf_summary, aes(x = Host_ID, y = mean_PC1, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  labs(title = "", x = "", y = "PCA Axis 1 - LES") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 16, colour="black"),
    axis.text.y = element_text(size = 16, colour="black"),
    axis.title.y = element_text(size = 16, colour="black"),
    axis.title.x = element_text(size = 16, colour="black"), 
    strip.text = element_text(size = 16, colour="black"))

PC1_leaf_plot


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

PC2_leaf_summary$Association <- c("EM", "EM", "DUAL", "AM", "DUAL", "AM", "EM")


# Test for significant differences between species 
aov_PC2_leaf <- aov(PC2 ~ Host_ID, data = PC2_leaf)
summary(aov_PC2_leaf)

tuk_PC2_leaf <- TukeyHSD(aov_PC2_leaf)
tuk_PC2_leaf

# Significant differences between species
# 

# Df Sum Sq Mean Sq F value   Pr(>F)    
# Host_ID      6  60.25  10.042   26.97 1.78e-14 ***
#   Residuals   53  19.74   0.372                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# diff         lwr        upr     p adj
# ABGR-ABAM  0.10842076 -0.8575181  1.0743596 0.9998550
# ALRU-ABAM -0.83219454 -1.8567274  0.1923383 0.1841661
# CONU-ABAM  0.09165061 -0.7448770  0.9281782 0.9998741
# TABR-ABAM -2.64568851 -3.5051389 -1.7862381 0.0000000
# THPL-ABAM -1.52220562 -2.3587332 -0.6856781 0.0000170
# TSHE-ABAM -0.05163324 -0.8881608  0.7848943 0.9999957
# ALRU-ABGR -0.94061530 -2.0732790  0.1920484 0.1645681
# CONU-ABGR -0.01677016 -0.9827090  0.9491687 1.0000000
# TABR-ABGR -2.75410927 -3.7399665 -1.7682521 0.0000000
# THPL-ABGR -1.63062639 -2.5965652 -0.6646876 0.0000709
# TSHE-ABGR -0.16005401 -1.1259928  0.8058848 0.9986397
# CONU-ALRU  0.92384514 -0.1006877  1.9483780 0.1027620
# TABR-ALRU -1.81349397 -2.8568271 -0.7701608 0.0000413
# THPL-ALRU -0.69001109 -1.7145439  0.3345218 0.3886374
# TSHE-ALRU  0.78056129 -0.2439715  1.8050941 0.2473032
# TABR-CONU -2.73733912 -3.5967895 -1.8778887 0.0000000
# THPL-CONU -1.61385623 -2.4503838 -0.7773287 0.0000051
# TSHE-CONU -0.14328385 -0.9798114  0.6932437 0.9983606
# THPL-TABR  1.12348289  0.2640325  1.9829333 0.0034590
# TSHE-TABR  2.59405527  1.7346049  3.4535056 0.0000000
# TSHE-THPL  1.47057238  0.6340448  2.3070999 0.0000333

# Visualize 
PC2_leaf_plot <- ggplot(PC2_leaf_summary, aes(x = Host_ID, y = mean_PC2, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  labs(title = "", x = "", y = "PCA Axis 2") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 16, colour="black"),
    axis.text.y = element_text(size = 16, colour="black"),
    axis.title.y = element_text(size = 16, colour="black"),
    axis.title.x = element_text(size = 16, colour="black"), 
    strip.text = element_text(size = 16, colour="black"))

PC2_leaf_plot

############### -- 

## Roots
root.pca = prcomp(root_sub[5:13], center = T, scale = T)

sd.root = root.pca$sdev
loadings.root = root.pca$rotation
trait.names.root = colnames(root_sub[5:13])
scores.root = as.data.frame(root.pca$x)
scores.root$WFDP_Code = root_sub$WFDP_Code
scores.root$Host_ID = root_sub$Host_ID
summary(root.pca)

# Save loadings for root traits
write.csv(loadings.root, "./PCA/PCA_loadings_root_traits.csv", row.names = TRUE)

#Save species scores
write.csv(scores.root, "./PCA/PCA_scores_root_traits.csv")


# PCA scores are 'scores.root' with column for Host_ID

# The PC1 axis is the equivalent of the leaf conservative-acquisitive axis, but right now the direction is 
# flipped so it's less intuitive to interpret them together. Going to multiply the PC1 values by -1 to 
# flip the orientation, and this doesn't do anything to the actual PCA, or the PC2 axis. 

# Flip PC1 scores
root.pca$x[, "PC1"] <- -1 * root.pca$x[, "PC1"]

# Flip PC1 loadings
root.pca$rotation[, "PC1"] <- -1 * root.pca$rotation[, "PC1"]

# Load back into items for plotting 
scores.root   <- as.data.frame(root.pca$x)
scores.root$WFDP_Code = root_sub$WFDP_Code
scores.root$Host_ID = root_sub$Host_ID

loadings.root <- as.data.frame(root.pca$rotation)

# Change loadings names to something cleaner 
new_loadings_root <- c("SRL", "SRA", "RDMC", "C:N", "d15N", "d13C", "RD", "PctN", "PctC")
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
                  color = "black", size = 5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal Association") + 
  scale_shape_manual(
    values = c("ABAM" = 21, "ABGR" = 22, "ALRU" = 23, "CONU" = 24, "TABR" = 25, "THPL" = 7, "TSHE" = 8), name = "Host Species") +
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Host Species") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 16, colour="black"),
    axis.text.y = element_text(size = 16, colour="black"),
    axis.title.y = element_text(size = 16, colour="black"),
    axis.title.x = element_text(size = 16, colour="black"), 
    strip.text = element_text(size = 16, colour="black"))

PCA_plot_root


# Get top 3 traits for PC1
top_PC1_root <- loadings.root[order(abs(loadings.root$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_root <- loadings.root[order(abs(loadings.root$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in Specific root area, specific root length, and root N content. 

# PC2 is being driven by the d13C, root C content, and d15N.  

## Together the axes explain 63.7% of the variation 


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
PC1_root_summary$Association <- c("EM", "EM", "DUAL", "AM", "DUAL", "AM", "EM")


# Test for significant differences between species 
aov_PC1_root <- aov(PC1 ~ Host_ID, data = PC1_root)
summary(aov_PC1_root)

tuk_PC1_root <- TukeyHSD(aov_PC1_root)
tuk_PC1_root


# Df Sum Sq Mean Sq F value  Pr(>F)   
# Host_ID      6  83.55  13.926    4.21 0.00155 **
#   Residuals   53 175.32   3.308                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# diff         lwr        upr     p adj
# ABGR-ABAM -0.669580162 -3.54760394  2.2084436 0.9912283
# ALRU-ABAM  3.019678077 -0.03292712  6.0722833 0.0543156
# CONU-ABAM  0.176600441 -2.31584126  2.6690421 0.9999904
# TABR-ABAM -0.966527325 -3.52726775  1.5942131 0.9069180
# THPL-ABAM  1.351030424 -1.14141128  3.8434721 0.6441673
# TSHE-ABAM -0.967531370 -3.45997308  1.5249103 0.8950243
# ALRU-ABGR  3.689258238  0.31447621  7.0640403 0.0236224
# CONU-ABGR  0.846180603 -2.03184318  3.7242044 0.9708396
# TABR-ABGR -0.296947164 -3.23431788  2.6404236 0.9999214
# THPL-ABGR  2.020610586 -0.85741319  4.8986344 0.3389428
# TSHE-ABGR -0.297951208 -3.17597499  2.5800726 0.9999096
# CONU-ALRU -2.843077635 -5.89568283  0.2095276 0.0833810
# TABR-ALRU -3.986205402 -7.09482638 -0.8775844 0.0043784
# THPL-ALRU -1.668647652 -4.72125285  1.3839575 0.6351514
# TSHE-ALRU -3.987209447 -7.03981464 -0.9346042 0.0034935
# TABR-CONU -1.143127767 -3.70386819  1.4176127 0.8159971
# THPL-CONU  1.174429983 -1.31801172  3.6668717 0.7755524
# TSHE-CONU -1.144131811 -3.63657352  1.3483099 0.7958352
# THPL-TABR  2.317557750 -0.24318268  4.8782982 0.1004278
# TSHE-TABR -0.001004045 -2.56174447  2.5597364 1.0000000
# TSHE-THPL -2.318561794 -4.81100350  0.1738799 0.0840558


# Visualize 
PC1_root_plot <- ggplot(PC1_root_summary, aes(x = Host_ID, y = mean_PC1, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 1 - RES") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 16, colour="black"),
    axis.text.y = element_text(size = 16, colour="black"),
    axis.title.y = element_text(size = 16, colour="black"),
    axis.title.x = element_text(size = 16, colour="black"), 
    strip.text = element_text(size = 16, colour="black"))

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
PC2_root_summary$Association <- c("EM", "EM", "DUAL", "AM", "DUAL", "AM", "EM")


# Test for significant differences between species 
aov_PC2_root <- aov(PC2 ~ Host_ID, data = PC2_root)
summary(aov_PC2_root)

tuk_PC2_root <- TukeyHSD(aov_PC2_root)
tuk_PC2_root

# NO significant differences between species

# 
# Df Sum Sq Mean Sq F value Pr(>F)
# Host_ID      6  13.86   2.310   1.875  0.102
# Residuals   53  65.29   1.232    

# Visualize 
PC2_root_plot <- ggplot(PC2_root_summary, aes(x = Host_ID, y = mean_PC2, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 2") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 16, colour="black"),
    axis.text.y = element_text(size = 16, colour="black"),
    axis.title.y = element_text(size = 16, colour="black"),
    axis.title.x = element_text(size = 16, colour="black"), 
    strip.text = element_text(size = 16, colour="black"))

PC2_root_plot


#################################################################################

############################################ -- 
# (4) ENVIRONMENTAL RELATIONSHIPS WITH PCAs
############################################ -- 

# Read back in environmental data for the trees 

x # Right now this is just slope, aspect, and elevation, but will be able to pull in 
    # krieged data when it exists 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# subset to data of interest 

env <- dplyr::select(env, Cell, WFDP_Code, slope, aspect, elevation_m)

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

# Adjusted R-squared:  0.087, p-value: 0.013


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

# This is just for the overall PCA, could also do it with the separate leaf and root 
# PCAs if desired 

#################################################################################

## Score csvs are saved for leaf and root traits already. 

# Save final plots

# All traits plot 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/all_traits_biplot.png", 
       plot = PCA_plot, width = 10, height = 8, units = "in", dpi = 300)

# leaf PCA
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/leaf_traits_biplot.png", 
       plot = PCA_plot_leaf, width = 8, height = 8, units = "in", dpi = 300)

# PC1 leaf bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC1_leaf_barplot.png", 
       plot = PC1_leaf_plot, width = 7, height = 8, units = "in", dpi = 300)

# PC2 leaf bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC2_leaf_barplot.png", 
       plot = PC2_leaf_plot, width = 7, height = 8, units = "in", dpi = 300)

# root PCA
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/root_traits_biplot.png", 
       plot = PCA_plot_root, width = 8, height = 8, units = "in", dpi = 300)

# PC1 root bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC1_root_barplot.png", 
       plot = PC1_root_plot, width = 7, height = 8, units = "in", dpi = 300)

# PC2 root bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC2__root_barplot.png", 
       plot = PC2_root_plot, width = 7, height = 8, units = "in", dpi = 300)


## -- END -- ## 

