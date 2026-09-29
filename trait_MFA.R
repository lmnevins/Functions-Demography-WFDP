# -----------------------------------------------------------------------------#
# Multiple factor analysis of leaf and root trait axes to assess relationships
# Original Author: L. McKinley Nevins 
# August 9, 2026
# Software versions:  R v 4.5.2
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 4.0.3
#                     rstatix v 0.7.3
#                     vegan v 2.7.3
#                     rcompanion v 2.5.1
#                     ggfortify v 0.4.19
#                     gginnards v 0.2.0.2
#                     ggrepel v 0.9.7
#                     corrplot v 0.95
#                     car v 3.1.3
#                     multcomp v 1.4.29
#                     multcompView v 0.1.10
#                     FactoMineR v 2.16
#                     factoextra v 2.2.0
#                     gtools v 3.9.5
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
library(FactoMineR); packageVersion("FactoMineR")
library(factoextra); packageVersion("factoextra")
library(gtools); packageVersion("gtools")


#################################################################################
#                               Main workflow                                   #
#  Explore the leaf and root trait data. Perform PCA analyses and then use      #
#  multiple factor analyses (MFA) to assess relationships between the axes.     #
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
                        leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, specific_root_length, specific_root_area,
                        root_dry_matter_cont, root_CN, root_15N, avg_root_dia, root_pct_N, root_pct_C)

## UPDATE: removed C isotope traits 


# make dataframe of full species names 

sci_name <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
              "T. plicata", "T. heterophylla")

Host_ID <- c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")


taxa <- data.frame(sci_name, Host_ID)

# Merge to traits data by Species 
traits <- merge(traits, taxa, by = "Host_ID")

leaf <- merge(leaf, taxa, by = "Host_ID")

root <- merge(root, taxa, by = "Host_ID")

# Load in site info 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

#################################################################################

##################################### -- 
# (2) PRINCIPLE COMPONENT ANALYSIS
##################################### -- 

# Perform separate PCAs for leaf and root traits 

# Grab environmental data again for some association labels 
env <- dplyr::select(env, WFDP_Code, Association)

## Use the original leaf and root trait datasets 

## Leaves first: 

# Subset to traits excluding petiole 
leaf_sub <- dplyr::select(leaf, code, WFDP_Code, sub_plot, Host_ID, SLA_leaf, LDMC_leaf, LMA_leaf, 
                          leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N)

# Then roots: 
root_sub <- dplyr::select(root, code, WFDP_Code, sub_plot, Host_ID, specific_root_length, 
                          specific_root_area, root_dry_matter_cont, root_CN, root_15N, 
                          avg_root_dia, root_pct_N, root_pct_C)



# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# Define shapes for species 

# ABAM, ABGR, ALRU, CONU, TABR, THPL, TSHE  
species_shapes <- c(15, 16, 17, 18, 7, 8, 9)


################################### -- 

# PCAs

## Leaves
leaf.pca = prcomp(leaf_sub[5:11], center = T, scale = T)

sd.leaf = leaf.pca$sdev
loadings.leaf = leaf.pca$rotation
trait.names.leaf = colnames(leaf_sub[5:11])
scores.leaf = as.data.frame(leaf.pca$x)
scores.leaf$WFDP_Code = leaf_sub$WFDP_Code
scores.leaf$Host_ID = leaf_sub$Host_ID
summary(leaf.pca)

# Save loadings for leaf traits
# write.csv(loadings.leaf, "./PCA/PCA_loadings_leaf_traits_no13C.csv", row.names = TRUE)

#Save species scores
# write.csv(scores.leaf, "./PCA/PCA_scores_leaf_traits_no13C.csv")



##Broken-Stick test for the significance of the loadings
print(leaf.pca)

plot(leaf.pca, type = "l")

ev = leaf.pca$sdev^2

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


## !! Broken stick only retains the first principal component 


# PCA scores are 'scores.leaf' with column for Host_ID

loadings.leaf <- as.data.frame(loadings.leaf)

# get proportion of variance explained to add to each axis label 
pca_var <- leaf.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

#Merge in mycorrhizal association for plotting 
scores.leaf <- merge(scores.leaf, env, by = "WFDP_Code")

# Change loadings names to something cleaner 
new_loadings <- c("SLA", "LDMC", "LMA", "PctN", "PctC", "C:N", "d15N")
rownames(loadings.leaf) <- new_loadings


# Merge in scientific name for plotting 
scores.leaf <- merge(scores.leaf, taxa, by = "Host_ID")


## Removing legend from all plots to save for formatting 

# Visualize
PCA_plot_leaf <- ggplot(scores.leaf, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3.5, aes(shape = sci_name)) +
  geom_segment(data = loadings.leaf, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.leaf, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.leaf)),
                  color = "black", size = 6.5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  guides(shape = guide_legend(nrow = 4, ncol = 2)) + 
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)")) +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16, face = "italic")) + 
  theme(
    axis.text.x = element_text(size = 18, colour="black"),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"))

PCA_plot_leaf


# Get top 3 traits for PC1
top_PC1_leaf <- loadings.leaf[order(abs(loadings.leaf$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_leaf <- loadings.leaf[order(abs(loadings.leaf$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in SLA, C:N ratio, and LMA, which matches what 
# would be expected for broad leaf vs needle-leaf species 

# PC2 is being driven by the leaf C content, d15N, and N content.  

## Together the axes explain 87.4% of the variation 


############### -- 

## Roots
root.pca = prcomp(root_sub[5:12], center = T, scale = T)

sd.root = root.pca$sdev
loadings.root = root.pca$rotation
trait.names.root = colnames(root_sub[5:12])
scores.root = as.data.frame(root.pca$x)
scores.root$WFDP_Code = root_sub$WFDP_Code
scores.root$Host_ID = root_sub$Host_ID
summary(root.pca)

# Save loadings for root traits
# write.csv(loadings.root, "./PCA/PCA_loadings_root_traits_no13C.csv", row.names = TRUE)

#Save species scores
# write.csv(scores.root, "./PCA/PCA_scores_root_traits_no13C.csv")



##Broken-Stick test for the significance of the loadings
print(root.pca)

plot(root.pca, type = "l")

ev = root.pca$sdev^2

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


## Broken stick retains PC1 and PC2 but only barely 



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
new_loadings_root <- c("SRL", "SRA", "RDMC", "C:N", "d15N", "RD", "PctN", "PctC")
rownames(loadings.root) <- new_loadings_root

#Merge in mycorrhizal association for plotting 
scores.root <- merge(scores.root, env, by = "WFDP_Code")


# Merge in scientific name for plotting 
scores.root <- merge(scores.root, taxa, by = "Host_ID")


# get proportion of variance explained to add to each axis label 
pca_var <- root.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


# Visualize
PCA_plot_root <- ggplot(scores.root, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3.5, aes(shape = sci_name)) +
  geom_segment(data = loadings.root, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.root, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.root)),
                  color = "black", size = 6.5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 18, colour="black"),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"), 
    strip.text = element_text(size = 18, colour="black"))

PCA_plot_root


# Get top 3 traits for PC1
top_PC1_root <- loadings.root[order(abs(loadings.root$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_root <- loadings.root[order(abs(loadings.root$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in Specific root area, specific root length, and root N content. 

# PC2 is being driven by root C content, d15N, and root diameter. 

## Together the axes explain 70.9% of the variation 



#################################################################################

##################################### -- 
# (3) MULTIPLE FACTOR ANALYSIS
##################################### -- 

# If I am first interested at the species level, I can take the average of the trait values for each species first 

species_traits <- traits %>%
  group_by(sci_name) %>%
  summarise(
    SLA = mean(SLA_leaf, na.rm = TRUE),
    LDMC = mean(LDMC_leaf, na.rm = TRUE),
    LMA = mean(LMA_leaf, na.rm = TRUE),
    leaf_C = mean(leaf_pct_C, na.rm = TRUE),
    leaf_N = mean(leaf_pct_N, na.rm = TRUE),
    leaf_CN = mean(leaf_CN, na.rm = TRUE),
    leaf_15N = mean(leaf_15N, na.rm = TRUE),
    SRL = mean(specific_root_length, na.rm = TRUE),
    SRA = mean(specific_root_area, na.rm = TRUE),
    RDMC = mean(root_dry_matter_cont, na.rm = TRUE),
    root_C = mean(root_pct_C, na.rm = TRUE),
    root_N = mean(root_pct_N, na.rm = TRUE),
    root_CN = mean(root_CN, na.rm = TRUE),
    root_15N = mean(root_15N, na.rm = TRUE),
    root_dia = mean(avg_root_dia, na.rm = TRUE))



# Create groups according to leaf and root traits 

leaf_traits <- c("SLA", "LDMC", "LMA", "leaf_C", "leaf_N", "leaf_CN", "leaf_15N")

root_traits <- c("SRL", "SRA", "RDMC", "root_C", "root_N", "root_CN", "root_15N", "root_dia")


# select traits 
mfa_data <- species_traits %>%
  dplyr::select(sci_name, all_of(leaf_traits), all_of(root_traits)) %>%
  drop_na()

# standardize traits 
mfa_traits <- mfa_data %>%
  dplyr::select(all_of(c(leaf_traits, root_traits))) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) %>%
  as.data.frame()

# add scientific names back in 
rownames(mfa_traits) <- mfa_data$sci_name



# Perform MFA using FactoMineR with the two groups being the overall leaf and root traits 

#specifying that the data is quantitative variables that have been scaled 's' 

mfa_res <- MFA(mfa_traits,
  group = c(length(leaf_traits), length(root_traits)),
  type = c("s", "s"), name.group = c("Leaf", "Root"), graph = FALSE) 


print(mfa_res)

mfa_res$eig # how much of the combined trait variation is represented by each MFA dimension

# Looks like the two components together quantitatively capture 83.47% of the variation 



mfa_res$group # look at the relationship between the two groups 

   # $coord = how strongly each trait group contributes to the MFA dimensions

#        Dim.1     Dim.2      Dim.3      Dim.4      Dim.5
# Leaf 0.7873941 0.2257535 0.06521141 0.07902611 0.05212294
# Root 0.7345892 0.3953419 0.11782137 0.03864593 0.03622368

# looks like both leaf and root traits contribute most strongly to the first dimension 

# $contrib
#        Dim.1    Dim.2    Dim.3    Dim.4    Dim.5
# Leaf 51.73474 36.34764 35.62827 67.15793 58.99823
# Root 48.26526 63.65236 64.37173 32.84207 41.00177



# Visualize the two trait groups 
factoextra::fviz_mfa_var(mfa_res, "group", repel = TRUE)


# Look at individual species 
factoextra::fviz_mfa_ind(mfa_res,repel = TRUE)


# Look at which traits are driving the MFA axes 
factoextra::fviz_mfa_var(mfa_res, "quanti.var", repel = TRUE)

# This one looks like the PCA plot when the traits were all analyzed together - the eigenvectors are 
# grouped the same way 


### Perform RV analyses 

# Create separate leaf and root matrices 

leaf_matrix <- as.matrix(
  mfa_traits[, leaf_traits]
)

root_matrix <- as.matrix(
  mfa_traits[, root_traits]
)


# Calculate the RV coefficient - using the standard Escoufier RV coefficient here 

# Create the function once 
rv_coeff <- function(X, Y) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  crossprod_XY <- crossprod(X, Y)
  
  numerator <- sum(crossprod_XY^2)
  
  denominator <- sqrt(
    sum(crossprod(X, X)^2) *
      sum(crossprod(Y, Y)^2)
  )
  
  numerator / denominator
}


rv_obs <- rv_coeff(
  leaf_matrix,
  root_matrix
)

rv_obs

# 0.3295231 not super strong, but this needs to be compared to the null via permutation 



## Calculate permutations 

set.seed(123)

n_perm <- 999 # 999 permutations, can adjust 

rv_null <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  
  permuted_root <- root_matrix[ # randomly shuffle the root traits by species 
    sample(seq_len(nrow(root_matrix))),
    ,
    drop = FALSE
  ]
  
  rv_null[i] <- rv_coeff(
    leaf_matrix, # keep the leaf matrix consistent 
    permuted_root
  )
}

summary(rv_null)

#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03296 0.14722 0.21241 0.26176 0.33749 0.86356 


quantile(rv_null, probs = c(0.025, 0.5, 0.975))

#     2.5%        50%      97.5% 
# 0.07126408 0.21241155 0.73397257 



p_value <- (sum(rv_null >= rv_obs) + 1) / (n_perm + 1)

p_value

# 0.266


# H₀: Leaf and root trait spaces are unrelated across species.

# H₁: Leaf and root trait spaces are more strongly associated than expected under independence.



# Visualize 

rv_null_df <- data.frame(rv = rv_null)

rv_plot <- ggplot(rv_null_df, aes(x = rv)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  geom_vline(
    xintercept = rv_obs,
    linewidth = 1) +
  labs(
    x = "RV coefficient",
    y = "Count",
    title = "Observed leaf–root trait association",
    subtitle = paste0(
      "Observed RV = ", round(rv_obs, 3),
      "; permutation p = ", signif(p_value, 3))) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))


rv_plot 




# OUTCOME: Observed RV and p-value suggest that there is no significant correlation between the root and leaf 
# trait spaces across the species 


#####################

# Now we want to look at the results from the separate leaf and root PC axes to see if they are correlated 
# right now the scores.leaf and scores.root dataframes have rows for each separate tree, so I want to 
# summarize the mean score for each species first 


species_leaf_scores <- scores.leaf %>%
  group_by(sci_name) %>%
  summarise(
    PC1_leaf = mean(PC1, na.rm = TRUE),
    PC2_leaf = mean(PC2, na.rm = TRUE))


species_root_scores <- scores.root %>%
  group_by(sci_name) %>%
  summarise(
    PC1_root = mean(PC1, na.rm = TRUE),
    PC2_root = mean(PC2, na.rm = TRUE))


# merge into one dataframe 
species_scores <- merge(species_leaf_scores, species_root_scores, by = "sci_name")


#Calculate the correlations 
pc_cor <- cor(
  species_scores %>%
    dplyr::select(PC1_leaf, PC2_leaf, PC1_root, PC2_root), use = "pairwise.complete.obs")

pc_cor


#            PC1_leaf    PC2_leaf   PC1_root    PC2_root
# PC1_leaf 1.00000000  0.04647534  0.5273439  0.29246870
# PC2_leaf 0.04647534  1.00000000  0.1661329 -0.01518681
# PC1_root 0.52734392  0.16613291  1.0000000 -0.43285705
# PC2_root 0.29246870 -0.01518681 -0.4328571  1.00000000


# Visualize this 

# format the data 


plot_data <- species_scores %>%
  dplyr::select(sci_name, PC1_leaf, PC2_leaf, PC1_root, PC2_root) %>%
  pivot_longer(
    cols = c(PC1_leaf, PC2_leaf),
    names_to = "leaf_axis",
    values_to = "leaf_score") %>%
  pivot_longer(
    cols = c(PC1_root, PC2_root),
    names_to = "root_axis",
    values_to = "root_score")


spp_plot <- ggplot(plot_data, aes(x = leaf_score, y = root_score)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(aes(label = sci_name), vjust = -0.7) +
  facet_grid(root_axis ~ leaf_axis) +
  labs(x = "Leaf PCA score", y = "Root PCA score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"))

spp_plot 



# Perform the same style of permutation test to see if the correlation we saw between the leaf PC1 and root PC1 scores 
# is higher or lower than what is expected at random 


set.seed(123)

# with only 7 species we can't do 9,999 permutations, so this allows us to calculate the exact number of permuations 
# we can do with the data, and then do all of them 

# using gtools 

root_perms <- permutations(
  n = 7,
  r = 7,
  v = species_scores$PC1_root
)

n_perm <- 5040


# now perform the rest of the calculation 
observed_r <- cor(
  species_scores$PC1_leaf,
  species_scores$PC1_root
)

null_r <- numeric(n_perm)

for (i in seq_len(n_perm)) {

  null_r[i] <- cor(
    species_scores$PC1_leaf, # keep the leaf scores consistent again 
    root_perms[i, ]
  )
}


# Calculate the p-value for the test 
p_perm <- (sum(abs(null_r) >= abs(observed_r)) + 1) /
  (n_perm + 1)

p_perm


# 0.1974


# Plot results 

corr_plot <- ggplot(data.frame(r = null_r), aes(x = r)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  geom_vline(xintercept = observed_r, linewidth = 1) +
  labs(x = "Pearson correlation (r)", 
       y = "Count",
    title = "Leaf PC1 × Root PC1",
    subtitle = paste0(
      "Observed r = ", round(observed_r, 3),
      "; exact permutation P = ", round(p_perm, 3))) +
  theme_classic()


corr_plot





###############

#### Make a loop to do this analysis for all of the species at one time 

# merge all tree PCA scores for leaves and roots 
scores_combined <- scores.leaf %>%
  dplyr::select(Host_ID, WFDP_Code, PC1_leaf = PC1, PC2_leaf = PC2) %>%
  inner_join(
    scores.root %>%
      dplyr::select(Host_ID, WFDP_Code, PC1_root = PC1, PC2_root = PC2),
  by = c("Host_ID", "WFDP_Code"))


# write a function to look the full correlation and permutation testing across each species 
test_species_correlation <- function(data,
                                     x = "PC1_leaf",
                                     y = "PC1_root",
                                     n_perm = 9999,
                                     seed = 123) {
  
  data <- data %>%
    dplyr::select(all_of(c(x, y))) %>%
    tidyr::drop_na()
  
  n <- nrow(data)
  
  # Not enough individuals for a correlation
  if (n < 3) {
    
    return(list(
      summary = tibble(
        n = n,
        observed_r = NA_real_,
        p_perm = NA_real_
      ),
      null = tibble(
        permutation = integer(),
        null_r = numeric()
      )
    ))
  }
  
  # Observed correlation
  observed_r <- cor(
    data[[x]],
    data[[y]]
  )
  
  # --------------------------------------------------
  # Generate null distribution
  # --------------------------------------------------
  
  if (n <= 8) {
    
    # Exact permutation test
    all_perms <- gtools::permutations(
      n = n,
      r = n,
      v = data[[y]]
    )
    
    null_r <- apply(
      all_perms,
      1,
      function(permuted_y) {
        cor(
          data[[x]],
          permuted_y
        )
      }
    )
    
  } else {
    
    # Random permutation test
    set.seed(seed)
    
    null_r <- numeric(n_perm)
    
    for (i in seq_len(n_perm)) {
      
      permuted_y <- sample(data[[y]])
      
      null_r[i] <- cor(
        data[[x]],
        permuted_y
      )
    }
  }
  
  # --------------------------------------------------
  # Calculate permutation p-value
  # --------------------------------------------------
  
  p_perm <- (sum(abs(null_r) >= abs(observed_r)) + 1) /
    (length(null_r) + 1)
  
  # --------------------------------------------------
  # Return both summary and null distribution
  # --------------------------------------------------
  
  summary_results <- tibble(
    n = n,
    observed_r = observed_r,
    p_perm = p_perm
  )
  
  null_results <- tibble(
    permutation = seq_along(null_r),
    null_r = null_r
  )
  
  return(list(
    summary = summary_results,
    null = null_results
  ))
}



species_list <- unique(scores_combined$Host_ID)

species_results_list <- list()
null_results_list <- list()


# Do an analysis for each species separately and organize the results 
for (sp in species_list) {
  
  sp_data <- scores_combined %>%
    filter(Host_ID == sp)
  
  test <- test_species_correlation(
    sp_data,
    x = "PC1_leaf",
    y = "PC1_root",
    n_perm = 9999
  )
  
  species_results_list[[sp]] <- test$summary %>%
    mutate(
      Host_ID = sp,
      .before = 1
    )
  
  null_results_list[[sp]] <- test$null %>%
    mutate(
      Host_ID = sp,
      .before = 1
    )
}


# Create two separate results dataframes that have the results from the observed r and the permutation p value
species_results <- bind_rows(species_results_list)

# and the null r results from the permutations
null_results <- bind_rows(null_results_list)


# perform analysis for all species together 
overall_test <- test_species_correlation(
  scores_combined,
  x = "PC1_leaf",
  y = "PC1_root",
  n_perm = 9999
)



# Bind the full species results to the dataframes including the separate species results 
species_results <- bind_rows(
  species_results,
  overall_test$summary %>%
    mutate(
      Host_ID = "All species",
      .before = 1
    )
)

null_results <- bind_rows(
  null_results,
  overall_test$null %>%
    mutate(
      Host_ID = "All species",
      .before = 1
    )
)



# Set the order to have the all species results appear first when faceting 
null_results <- null_results %>%
  mutate(
    Host_ID = factor(
      Host_ID,
      levels = c("All species", "ABAM", "ABGR", "ALRU", "CONU", 
                 "TABR", "THPL", "TSHE")))

species_results <- species_results %>%
  mutate(
    Host_ID = factor(
      Host_ID,
      levels = c("All species", "ABAM", "ABGR", "ALRU", "CONU", 
                 "TABR", "THPL", "TSHE")))


# Add some results labels to each panel 
species_labels <- species_results %>%
  mutate(
    label = paste0(
      "r = ", round(observed_r, 2),
      "\nP = ", round(p_perm, 3),
      "\nn = ", n))






# Plot full permuted results for all species and the separate species 
perm_results_all_spp <- ggplot(null_results, aes(x = null_r)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  geom_vline(data = species_results, aes(xintercept = observed_r), linewidth = 1, color = "red3") +
  facet_wrap(~ Host_ID, ncol = 4, scales = "free_y") +
  geom_text(data = species_labels, aes(x = -0.95, y = Inf, label = label),
    hjust = 0, vjust = 1.2, size = 4) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Pearson correlation (r)", y = "Count") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"),
    strip.text = element_text(size = 12, colour="black"))


perm_results_all_spp


# Save faceted plot 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/Figures/PC1_PC1_all_spp.png", 
       plot = perm_results_all_spp, width = 14, height = 8, units = "in", dpi = 300)












############################


# LATER

# Plot the correlation results for each species 
species_plot <- ggplot(
  scores_combined,
  aes(x = PC1_leaf, y = PC1_root)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Host_ID) +
  geom_text(data = species_results,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1) +
  labs(x = "Leaf PC1", y = "Root PC1") +
  theme_classic()

species_plot




# Generate list of axis combinations for analyses 

axis_combinations <- tibble(
  x = c(
    "PC1_leaf",
    "PC1_leaf",
    "PC2_leaf",
    "PC2_leaf"
  ),
  y = c(
    "PC1_root",
    "PC2_root",
    "PC1_root",
    "PC2_root"
  ),
  comparison = c(
    "Leaf PC1 × Root PC1",
    "Leaf PC1 × Root PC2",
    "Leaf PC2 × Root PC1",
    "Leaf PC2 × Root PC2"
  )
)
