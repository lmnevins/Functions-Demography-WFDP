# -----------------------------------------------------------------------------#
# Tree trait - growth relationships in WFDP
# Original Author: L. McKinley Nevins 
# December 22, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     sf v 1.0.19
#                     raster v 3.6.32
#                     cowplot v 1.2.0
#                     emmeans v 2.0.1
#                     lme4 v 1.1.38
#                     gstat v 2.1.4
#                     car v 3.1.3
#                     FNN v 1.1.4.1
#                     spdep v 1.4.1
#                     ade4 v 1.7.23
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(sf); packageVersion("sf")
library(raster); packageVersion("raster")
library(cowplot); packageVersion("cowplot")
library(emmeans); packageVersion("emmeans")
library(lme4); packageVersion("lme4")
library(gstat); packageVersion("gstat")
library(car); packageVersion("car")
library(FNN); packageVersion("FNN")
library(spdep); packageVersion("spdep")
library(ade4); packageVersion("ade4")

#################################################################################
#                               Main workflow                                   #
#  Explore the growth data for all of the trait trees sampled at WFDP. Relate   #
#  to trait information for the trees themselves.Use the values of the leaf and # 
#  root PC1 and PC2 values for each focal tree to explore relationships of      # 
#  trait strategies with focal trees growth.                                    #
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)


## PCA Results
# Read in the scores values for the leaf and root PCAs

scores.leaf <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PCA_scores_leaf_traits.csv")
  
scores.root <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PCA_scores_root_traits.csv")
  
  
# Put WFDP code into compatible column name for merging with growth data 
scores.leaf$focal_stem_tag <- scores.leaf$WFDP_Code
  
scores.root$focal_stem_tag <- scores.root$WFDP_Code
  

# Read in scores for PCA of all traits together 

scores.all <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/tree_PC_scores.csv")


## Tree Demography
# Read in tree demographic data 
growth <- read.csv("stems_WFDP_20250206_trimmed.csv")
  
# each tree is identified by its WFDP stem_tag. There are three census years - 2011, 
# 2016, and 2021.


# get diameters for each of the three timepoints for all 60 trees
# Some trees have 2010 or 2011 for their first census date 
rgr_intervals <- growth %>%
  arrange(Stem_Tag, DBH_DATE) %>%
  group_by(Stem_Tag, Species) %>%
  mutate(
    DBH_prev  = lag(DBH),
    year_prev = lag(DBH_DATE),
    RGR_interval = (log(DBH) - log(DBH_prev)) / (DBH_DATE - year_prev)
  ) %>%
  filter(!is.na(RGR_interval))

# calc mean RGR across both intervals 
diams <- rgr_intervals %>%
  summarise(
    mean_RGR = mean(RGR_interval, na.rm = TRUE),
    n_intervals = n(),
    .groups = "drop"
  )


# Add column for WFDP_code to match the PCA data
diams$WFDP_Code <- diams$Stem_Tag


## Merge diams file to the leaf and root scores datafiles according to the WFDP_Code
leaf_scores_growth <- merge(diams, scores.leaf, by = "WFDP_Code")

root_scores_growth <- merge(diams, scores.root, by = "WFDP_Code")

# Merge to the full traits PCA scores file 
all_traits_scores_growth <- merge(diams, scores.all, by = "WFDP_Code")


## Environmental/ Topographic Data 

# This can be updated once the krieged environmental data is in hand 

env_data <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# Trim to just relevant columns 
env_data <- dplyr::select(env_data, Cell, EM_Sample_Name, slope, aspect, elevation_m, Association, WFDP_Code)

# Make association a factor 
env_data$Association <- as.factor(env_data$Association)


# Merge enviro data to the existing files according to the WFDP_Code
leaf_growth_env <- merge(leaf_scores_growth, env_data, by = "WFDP_Code")

root_growth_env <- merge(root_scores_growth, env_data, by = "WFDP_Code")

# Merge to the full traits PCA scores file 
all_traits_growth_env <- merge(all_traits_scores_growth, env_data, by = "WFDP_Code")

# Make Species a factor 
leaf_growth_env$Species <- as.factor(leaf_growth_env$Species)

root_growth_env$Species <- as.factor(root_growth_env$Species)

all_traits_growth_env$Species <- as.factor(all_traits_growth_env$Species)


## These datafiles now have everything compiled to explore relationships of tree growth with 
# topographic variables, and the host tree traits 

######################### 

# Also load in the individual raw trait datasets and clean a bit 

##Leaves 
leaf <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/WFDP_leaf_traits.csv")

# make species a factor 
leaf$Host_ID <- as.factor(leaf$Host_ID)

# filter out PSME as a host since there were so few sampled 
leaf <- leaf %>% filter(Host_ID != 'PSME')

# filter out T-TABR-03 as a host since it had no fungal community data 
leaf <- leaf %>% filter(code != 'T-TABR-03')

## Roots
root <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/WFDP_root_traits.csv")

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


# Merge this traits data to the all_traits_growth_env df to use to explore individual trait 
# relationships to growth later on 


trait_growth_df <- merge(traits, all_traits_growth_env, by = "WFDP_Code")

# Reload Host_ID that got renamed when merging 
trait_growth_df$Host_ID <- trait_growth_df$Host_ID.x


# make dataframe of full species names 

full <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
          "T. plicata", "T. heterophylla")

Species <- c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")


taxa <- data.frame(full, Species)

# Merge to demo data by Species 
diams <- merge(diams, taxa, by = "Species")


###############################

# set colors for hosts 
                # ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# Define shapes for species 

# ABAM, ABGR, ALRU, CONU, TABR, THPL, TSHE  
species_shapes <- c(15, 16, 17, 18, 7, 8, 9)



spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
diams$full <- factor(diams$full, levels = spp_order)


#################################################################################

################################################ -- 
# (2) GROWTH RELATIONSHIPS TO TRAIT PCA VALUES
################################################ -- 


### RGR ####

## LEAF TRAITS

# Merge in full species names to generate a legend  
leaf_growth_env <- merge(leaf_growth_env, taxa, by = "Species")

# PC1 - All species together 
leaf_PC1_rgr_all <- ggplot(leaf_growth_env, aes(x = PC1, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                    name="Focal Species",
                    breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                             "T. plicata", "T. heterophylla"),
                    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                             "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Leaf Trait PC1 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black", face = "italic"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC1_rgr_all


# test relationships 
lm_leaf_PC1_rgr_all <- lm(mean_RGR ~ PC1, data = leaf_growth_env)

summary(lm_leaf_PC1_rgr_all)

# 
# Multiple R-squared:  0.006539,	Adjusted R-squared:  -0.01059 
# F-statistic: 0.3818 on 1 and 58 DF,  p-value: 0.5391

# No significant relationship for the species all together 


# PC1 - For separate species
leaf_PC1_rgr_sep <- ggplot(leaf_growth_env, aes(x = PC1, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC1 Value", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC1_rgr_sep


# test relationships 
lm_leaf_PC1_rgr_sep <- lm(mean_RGR ~ PC1 * Species, data = leaf_growth_env)

emtrends(lm_leaf_PC1_rgr_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend      SE df t.ratio p.value
# ABAM     -0.00669 0.00485 46  -1.377  0.3064
# ABGR     -0.00295 0.00598 46  -0.493  0.6970
# ALRU      0.02867 0.01070 46   2.667  0.0737
# CONU      0.00196 0.00499 46   0.392  0.6970
# TABR     -0.00152 0.00273 46  -0.558  0.6970
# THPL      0.00415 0.00299 46   1.388  0.3064
# TSHE      0.00685 0.00365 46   1.878  0.2333

# No significant relationships for any species 

anova(lm_leaf_PC1_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0000496 0.00004963  0.6515   0.42372    
# Species      6 0.0029733 0.00049555  6.5057 5.083e-05 ***
#   PC1:Species  6 0.0010622 0.00017703  2.3241   0.04822 *  
#   Residuals   46 0.0035040 0.00007617                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species,and the leaf PC1 values could be driving some of the variation 
# between the species 


## PC2 - All species together 
leaf_PC2_rgr_all <- ggplot(leaf_growth_env, aes(x = PC2, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Leaf Trait PC2 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC2_rgr_all


# test relationships 
lm_leaf_PC2_rgr_all <- lm(mean_RGR ~ PC2, data = leaf_growth_env)

summary(lm_leaf_PC2_rgr_all)


# Multiple R-squared:  0.01531,	Adjusted R-squared:  -0.001663 
# F-statistic: 0.902 on 1 and 58 DF,  p-value: 0.3462

# No significant relationship for the species all together 


# PC2 - For separate species
leaf_PC2_rgr_sep <- ggplot(leaf_growth_env, aes(x = PC2, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC2 Value", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC2_rgr_sep


# test relationships 
lm_leaf_PC2_rgr_sep <- lm(mean_RGR ~ PC2 * Species, data = leaf_growth_env)

emtrends(lm_leaf_PC2_rgr_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend      SE df t.ratio p.value
# ABAM      0.00472 0.00548 46   0.861  0.9606
# ABGR     -0.00063 0.01270 46  -0.050  0.9606
# ALRU     -0.00117 0.01020 46  -0.114  0.9606
# CONU     -0.00407 0.00629 46  -0.648  0.9606
# TABR     -0.00165 0.00482 46  -0.343  0.9606
# THPL     -0.00172 0.00458 46  -0.376  0.9606
# TSHE     -0.00236 0.00522 46  -0.452  0.9606

# No significant relationships for any species 

anova(lm_leaf_PC2_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC2          1 0.0001162 0.00011622  1.1882 0.2813680    
# Species      6 0.0028327 0.00047212  4.8269 0.0006768 ***
#   PC2:Species  6 0.0001410 0.00002349  0.2402 0.9608359    
# Residuals   46 0.0044992 0.00009781                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, but their leaf PC2 values are not what is driving this. 



## ROOT TRAITS

# Merge in full species names to generate a legend  
root_growth_env <- merge(root_growth_env, taxa, by = "Species")

# PC1 - All species together 
root_PC1_rgr_all <- ggplot(root_growth_env, aes(x = PC1, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Root Trait PC1 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC1_rgr_all


# test relationships 
lm_root_PC1_rgr_all <- lm(mean_RGR ~ PC1, data = root_growth_env)

summary(lm_root_PC1_rgr_all)


# Multiple R-squared:  0.05804,	Adjusted R-squared:  0.0418 
# F-statistic: 3.574 on 1 and 58 DF,  p-value: 0.06369

# No significant relationship for the species all together 


# PC1 - For separate species
root_PC1_rgr_sep <- ggplot(root_growth_env, aes(x = PC1, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC1 Value", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC1_rgr_sep


# test relationships 
lm_root_PC1_rgr_sep <- lm(mean_RGR ~ PC1 * Species, data = root_growth_env)

emtrends(lm_root_PC1_rgr_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend      SE df t.ratio p.value
# ABAM      0.00179 0.00189 46   0.948  0.4841
# ABGR      0.00117 0.00166 46   0.705  0.4841
# ALRU      0.00317 0.00141 46   2.251  0.2046
# CONU     -0.00296 0.00203 46  -1.456  0.4706
# TABR      0.00184 0.00217 46   0.850  0.4841
# THPL     -0.00125 0.00167 46  -0.748  0.4841
# TSHE     -0.00286 0.00221 46  -1.295  0.4706

# No significant relationships for any species 

anova(lm_root_PC1_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0004405 0.00044049  5.4395 0.0241183 *  
#   Species      6 0.0025396 0.00042326  5.2267 0.0003575 ***
#   PC1:Species  6 0.0008839 0.00014732  1.8192 0.1162116    
# Residuals   46 0.0037251 0.00008098                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, and PC1 values are also different. 


## PC2 - All species together 
root_PC2_rgr_all <- ggplot(root_growth_env, aes(x = PC2, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Root Trait PC2 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC2_rgr_all


# test relationships 
lm_root_PC2_rgr_all <- lm(mean_RGR ~ PC2, data = root_growth_env)

summary(lm_root_PC2_rgr_all)


# Multiple R-squared:  0.08471,	Adjusted R-squared:  0.06893 
# F-statistic: 5.368 on 1 and 58 DF,  p-value: 0.02406

# SIGNIFICANT RELATIONSHIP for the species all together 


# PC2 - For separate species
root_PC2_rgr_sep <- ggplot(root_growth_env, aes(x = PC2, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC2 Value", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC2_rgr_sep


# test relationships 
lm_root_PC2_rgr_sep <- lm(mean_RGR ~ PC2 * Species, data = root_growth_env)

emtrends(lm_root_PC2_rgr_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend      SE df t.ratio p.value
# ABAM     0.007089 0.00270 46   2.629  0.0811
# ABGR     0.003485 0.00444 46   0.784  0.6117
# ALRU    -0.000353 0.00216 46  -0.164  0.8708
# CONU     0.002275 0.00257 46   0.887  0.6117
# TABR    -0.001566 0.00535 46  -0.293  0.8708
# THPL     0.005287 0.00257 46   2.053  0.1601
# TSHE     0.004078 0.00306 46   1.331  0.4431

# No significant relationships for any species 

anova(lm_root_PC2_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value   Pr(>F)    
# PC2          1 0.0006429 0.00064289  8.3353 0.005906 ** 
#   Species      6 0.0029019 0.00048364  6.2706  7.2e-05 ***
#   PC2:Species  6 0.0004965 0.00008274  1.0728 0.392579    
# Residuals   46 0.0035479 0.00007713                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, but their root PC2 values are not what is driving this. 


## Outcome: Significant result of a relationship between RGR and the Root PC2 results, 
# but it is no longer significant when the species are considered independently. Leaf PC1 is 
# significant when species are considered independently. 


#################################################################################

#################################################### -- 
# (3) GROWTH RELATIONSHIPS TO ALL TRAIT PCA VALUES
#################################################### -- 

# Same analyses as above to explore the relationships of RGR with the tree traits, but using the PCA 
# scores for all traits considered together. 


# PC1 - All species together 
all_PC1_rgr_all <- ggplot(all_traits_growth_env, aes(x = PC1, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC1 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

all_PC1_rgr_all


# test relationships 
lm_all_PC1_rgr_all <- lm(mean_RGR ~ PC1, data = all_traits_growth_env)

summary(lm_all_PC1_rgr_all)


# Multiple R-squared:  0.02399,	Adjusted R-squared:  0.007163 
# F-statistic: 1.426 on 1 and 58 DF,  p-value: 0.2373

# No significant relationship for the species all together 


# PC1 - For separate species
all_PC1_rgr_sep <- ggplot(all_traits_growth_env, aes(x = PC1, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC1 Value", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

all_PC1_rgr_sep


# test relationships 
lm_all_PC1_rgr_sep <- lm(mean_RGR ~ PC1 * Species, data = all_traits_growth_env)

emtrends(lm_all_PC1_rgr_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend      SE df t.ratio p.value
# ABAM     -0.00452 0.00265 46  -1.706  0.1896
# ABGR     -0.00148 0.00209 46  -0.707  0.4831
# ALRU     -0.00504 0.00274 46  -1.842  0.1896
# CONU      0.00638 0.00379 46   1.683  0.1896
# TABR     -0.00223 0.00241 46  -0.927  0.4188
# THPL      0.00362 0.00256 46   1.414  0.2296
# TSHE      0.00357 0.00218 46   1.637  0.1896

# No significant relationships for any species 

anova(lm_all_PC1_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0001821 0.00018207  2.3904 0.1289369    
# Species      6 0.0027532 0.00045887  6.0246 0.0001042 ***
#   PC1:Species  6 0.0011501 0.00019169  2.5167 0.0343761 *  
#   Residuals   46 0.0035037 0.00007617                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# SIGNIFICANT: RGR is different between species, and their PC1 values may be part of what is driving this. 


## PC2 - All species together 
all_PC2_rgr_all <- ggplot(all_traits_growth_env, aes(x = PC2, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC2 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

all_PC2_rgr_all


# test relationships 
lm_all_PC2_rgr_all <- lm(mean_RGR ~ PC2, data = all_traits_growth_env)

summary(lm_all_PC2_rgr_all)


# Multiple R-squared:  0.03748,	Adjusted R-squared:  0.02088 
# F-statistic: 2.258 on 1 and 58 DF,  p-value: 0.1383

# No significant relationship for species all together 


# PC2 - For separate species
all_PC2_rgr_sep <- ggplot(all_traits_growth_env, aes(x = PC2, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC2 Value", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "bottom")

all_PC2_rgr_sep


# test relationships 
lm_all_PC2_rgr_sep <- lm(mean_RGR ~ PC2 * Species, data = all_traits_growth_env)

emtrends(lm_all_PC2_rgr_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend      SE df t.ratio p.value
# ABAM     0.000923 0.00235 46   0.392  0.6965
# ABGR     0.001723 0.00252 46   0.683  0.6965
# ALRU     0.003639 0.00159 46   2.292  0.1858
# CONU    -0.002769 0.00213 46  -1.300  0.6965
# TABR     0.001532 0.00261 46   0.588  0.6965
# THPL    -0.000909 0.00190 46  -0.478  0.6965
# TSHE    -0.003346 0.00367 46  -0.911  0.6965

# No significant relationships for any species 

anova(lm_all_PC2_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC2          1 0.0002844 0.00028443  3.3553 0.0734684 .  
# Species      6 0.0027029 0.00045049  5.3142 0.0003114 ***
#   PC2:Species  6 0.0007023 0.00011704  1.3807 0.2425869    
# Residuals   46 0.0038995 0.00008477                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# SIGNIFICANT: RGR is different between species, but the interaction of these is not what is driving variation in RGR


#################################################################################

############################################## -- 
# (4) INDIVIDUAL TRAIT RELATIONSHIPS TO GROWTH  
############################################## -- 

# Exploring some relationships of individual traits to tree RGR because these are often considered independently, 
# rather than by lumping traits together into broader strategies. 


# Studies find very mixed relationships of individual traits to growth, and often none at all (the root of much 
# of this study's rationale). Anticipating questions about this, so exploring these now. 

# Focusing on the leaf and root traits that were considered in the PCAs

# Merge in full species names to generate a legend  
trait_growth_df <- merge(trait_growth_df, taxa, by = "Species")


## LEAF TRAITS


# SLA
leaf_sla_rgr <- ggplot(trait_growth_df, aes(x = SLA_leaf, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = expression("Specific Leaf Area (mm"^2* "mg"^-1*")"), y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_sla_rgr


# test relationships 
lm_leaf_sla_rgr <- lm(mean_RGR ~ SLA_leaf, data = trait_growth_df)

summary(lm_leaf_sla_rgr)

# Multiple R-squared:  0.01383,	Adjusted R-squared:  -0.003172 
# F-statistic: 0.8135 on 1 and 58 DF,  p-value: 0.3708



# LMA
leaf_lma_rgr <- ggplot(trait_growth_df, aes(x = LMA_leaf, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = expression("Leaf Mass per Area (mg mm"^-2*")"), y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_lma_rgr


# test relationships 
lm_leaf_lma_rgr <- lm(mean_RGR ~ LMA_leaf, data = trait_growth_df)

summary(lm_leaf_lma_rgr)

# Multiple R-squared:  0.0003761,	Adjusted R-squared:  -0.01686 
# F-statistic: 0.02182 on 1 and 58 DF,  p-value: 0.8831



# LDMC
leaf_ldmc_rgr <- ggplot(trait_growth_df, aes(x = LDMC_leaf, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = expression("Leaf Dry Matter Content (mg g"^-1*")"), y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_ldmc_rgr


# test relationships 
lm_leaf_ldmc_rgr <- lm(mean_RGR ~ LDMC_leaf, data = trait_growth_df)

summary(lm_leaf_ldmc_rgr)


# Multiple R-squared:  0.008526,	Adjusted R-squared:  -0.008569 
# F-statistic: 0.4987 on 1 and 58 DF,  p-value: 0.4829



# Leaf Pct_N
leaf_pctN_rgr <- ggplot(trait_growth_df, aes(x = leaf_pct_N, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Percent N (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_pctN_rgr


# test relationships 
lm_leaf_pctN_rgr <- lm(mean_RGR ~ leaf_pct_N, data = trait_growth_df)

summary(lm_leaf_pctN_rgr)


# Multiple R-squared:  0.03179,	Adjusted R-squared:  0.01509 
# F-statistic: 1.904 on 1 and 58 DF,  p-value: 0.1729



# Leaf Pct_C
leaf_pctC_rgr <- ggplot(trait_growth_df, aes(x = leaf_pct_C, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Percent C (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_pctC_rgr


# test relationships 
lm_leaf_pctC_rgr <- lm(mean_RGR ~ leaf_pct_C, data = trait_growth_df)

summary(lm_leaf_pctC_rgr)


# Multiple R-squared:  0.0132,	Adjusted R-squared:  -0.00381 
# F-statistic: 0.776 on 1 and 58 DF,  p-value: 0.382



# Leaf CN
leaf_CN_rgr <- ggplot(trait_growth_df, aes(x = leaf_CN, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Carbon:Nitrogen Ratio", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_CN_rgr


# test relationships 
lm_leaf_CN_rgr <- lm(mean_RGR ~ leaf_CN, data = trait_growth_df)

summary(lm_leaf_CN_rgr)

# Multiple R-squared:  0.002476,	Adjusted R-squared:  -0.01472 
# F-statistic: 0.1439 on 1 and 58 DF,  p-value: 0.7058


# Leaf 15N
leaf_15N_rgr <- ggplot(trait_growth_df, aes(x = leaf_15N, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "15N Isotope per mil", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_15N_rgr


# test relationships 
lm_leaf_15N_rgr <- lm(mean_RGR ~ leaf_15N, data = trait_growth_df)

summary(lm_leaf_15N_rgr)


# Multiple R-squared:  0.01645,	Adjusted R-squared:  -0.0005047 
# F-statistic: 0.9702 on 1 and 58 DF,  p-value: 0.3287


# Leaf 13C
leaf_13C_rgr <- ggplot(trait_growth_df, aes(x = leaf_13C, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "13C Isotope per mil", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_13C_rgr


# test relationships 
lm_leaf_13C_rgr <- lm(mean_RGR ~ leaf_13C, data = trait_growth_df)

summary(lm_leaf_13C_rgr)

# Multiple R-squared:  0.004138,	Adjusted R-squared:  -0.01303 
# F-statistic: 0.241 on 1 and 58 DF,  p-value: 0.6253


## ROOT TRAITS


# SRA
root_sra_rgr <- ggplot(trait_growth_df, aes(x = specific_root_area, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = expression("Specific Root Area (cm"^2* "mg"^-1*")"), y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_sra_rgr


# test relationships 
lm_root_sra_rgr <- lm(mean_RGR ~ specific_root_area, data = trait_growth_df)

summary(lm_root_sra_rgr)


# Multiple R-squared:  0.02201,	Adjusted R-squared:  0.005143 
# F-statistic: 1.305 on 1 and 58 DF,  p-value: 0.258



# SRL
root_srl_rgr <- ggplot(trait_growth_df, aes(x = specific_root_length, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = expression("Specific Root Length (cm mg"^-1*")"), y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_srl_rgr


# test relationships 
lm_root_srl_rgr <- lm(mean_RGR ~ specific_root_length, data = trait_growth_df)

summary(lm_root_srl_rgr)


# Multiple R-squared:  0.04867,	Adjusted R-squared:  0.03227 
# F-statistic: 2.967 on 1 and 58 DF,  p-value: 0.0903



# Root diameter
root_dia_rgr <- ggplot(trait_growth_df, aes(x = avg_root_dia, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Average Root Diameter (mm)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_dia_rgr


# test relationships 
lm_root_dia_rgr <- lm(mean_RGR ~ avg_root_dia, data = trait_growth_df)

summary(lm_root_dia_rgr)

# Multiple R-squared:  0.05501,	Adjusted R-squared:  0.03872 
# F-statistic: 3.377 on 1 and 58 DF,  p-value: 0.07125


# RDMC
root_rdmc_rgr <- ggplot(trait_growth_df, aes(x = root_dry_matter_cont, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = expression("Root Dry Matter Content (mg g"^-1*")"), y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_rdmc_rgr


# test relationships 
lm_root_rdmc_rgr <- lm(mean_RGR ~ root_dry_matter_cont, data = trait_growth_df)

summary(lm_root_rdmc_rgr)


# Multiple R-squared:  0.0672,	Adjusted R-squared:  0.05111 
# F-statistic: 4.178 on 1 and 58 DF,  p-value: 0.0455

# Significant - higher structural investment related to lower RGR



# Root Pct_N
root_pctN_rgr <- ggplot(trait_growth_df, aes(x = root_pct_N, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Percent N (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_pctN_rgr


# test relationships 
lm_root_pctN_rgr <- lm(mean_RGR ~ root_pct_N, data = trait_growth_df)

summary(lm_root_pctN_rgr)


# Multiple R-squared:  0.01046,	Adjusted R-squared:  -0.006597 
# F-statistic: 0.6133 on 1 and 58 DF,  p-value: 0.4367



# Root Pct_C
root_pctC_rgr <- ggplot(trait_growth_df, aes(x = root_pct_C, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Percent C (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_pctC_rgr


# test relationships 
lm_root_pctC_rgr <- lm(mean_RGR ~ root_pct_C, data = trait_growth_df)

summary(lm_root_pctC_rgr)


# Multiple R-squared:  0.1582,	Adjusted R-squared:  0.1437 
# F-statistic:  10.9 on 1 and 58 DF,  p-value: 0.001649

# Significant - higher structural investment related to lower RGR



# Root CN
root_CN_rgr <- ggplot(trait_growth_df, aes(x = root_CN, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Carbon:Nitrogen Ratio", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_CN_rgr


# test relationships 
lm_root_CN_rgr <- lm(mean_RGR ~ root_CN, data = trait_growth_df)

summary(lm_root_CN_rgr)


# Multiple R-squared:  0.0192,	Adjusted R-squared:  0.00229 
# F-statistic: 1.135 on 1 and 58 DF,  p-value: 0.291




# Root 15N
root_15N_rgr <- ggplot(trait_growth_df, aes(x = root_15N, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "15N Isotope per mil", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_15N_rgr


# test relationships 
lm_root_15N_rgr <- lm(mean_RGR ~ root_15N, data = trait_growth_df)

summary(lm_root_15N_rgr)


# Multiple R-squared:  0.04595,	Adjusted R-squared:  0.02951 
# F-statistic: 2.794 on 1 and 58 DF,  p-value: 0.1


# Root 13C
root_13C_rgr <- ggplot(trait_growth_df, aes(x = root_13C, y = mean_RGR, colour = full)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = full)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  guides(colour = guide_legend(nrow = 1, ncol = 7), shape = guide_legend(nrow = 1, ncol = 7),) + 
  labs(x = "13C Isotope per mil", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black", face = "italic"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_13C_rgr


# test relationships 
lm_root_13C_rgr <- lm(mean_RGR ~ root_13C, data = trait_growth_df)

summary(lm_root_13C_rgr)

# Multiple R-squared:  0.01855,	Adjusted R-squared:  0.001628 
# F-statistic: 1.096 on 1 and 58 DF,  p-value: 0.2994


#################################################################################

############################################## -- 
# (5) GATHER SIGNIFICANT RESULTS OF INTEREST  
############################################## -- 

# Gather up the most interesting relationships and organize the plots to compare the relationships of RGR
# with the leaf and root PC1 and PC2 axes for each of the host trees 


# The only significant results were when the species were all considered together, so I will do those 
# plots instead of the facetted ones. 


# Gather RGR plots with all species together for the separate leaf and root traits:
# leaf_PC1_rgr_all, leaf_PC2_rgr_all, root_PC1_rgr_all, root_PC2_rgr_all

# And for the combined PCA: 
# all_PC1_rgr_all, all_PC2_rgr_all


RGR_plots <- plot_grid(all_PC1_rgr_all, all_PC2_rgr_all,
                       leaf_PC1_rgr_all, root_PC1_rgr_all, 
                       leaf_PC2_rgr_all, root_PC2_rgr_all,
                                  ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'))

RGR_plots


# Plot with just the leaf and root traits 

RGR_plots_2 <- plot_grid(leaf_PC1_rgr_all, root_PC1_rgr_all, 
                         leaf_PC2_rgr_all, root_PC2_rgr_all,
                         ncol = 2, nrow = 2, labels = c('(a)', '(b)', '(c)', '(d)'))

RGR_plots_2 


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/trait_RGR_plots.png", 
       plot = RGR_plots_2, width = 7, height = 6, units = "in", dpi = 300)



# Separate trait~RGR plots

leaf_trait_RGR_plots <- plot_grid(leaf_sla_rgr, leaf_lma_rgr, leaf_ldmc_rgr, leaf_pctC_rgr, leaf_pctN_rgr, 
                             leaf_CN_rgr, leaf_13C_rgr, leaf_15N_rgr,
                             ncol = 3, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'))

leaf_trait_RGR_plots


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/leaf_trait_RGR_plots.png", 
       plot = leaf_trait_RGR_plots, width = 12, height = 11, units = "in", dpi = 300)



root_trait_RGR_plots <- plot_grid(root_sra_rgr, root_srl_rgr, root_rdmc_rgr, root_pctC_rgr, root_pctN_rgr, 
                                  root_CN_rgr, root_13C_rgr, root_15N_rgr, root_dia_rgr, 
                                  ncol = 3, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'))

root_trait_RGR_plots


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/root_trait_RGR_plots.png", 
       plot = root_trait_RGR_plots, width = 12, height = 11, units = "in", dpi = 300)



## -- END -- ## 

