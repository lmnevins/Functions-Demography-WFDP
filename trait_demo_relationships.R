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

## Using updates scores for the PCAs not including 13C 

scores.leaf <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PCA_scores_leaf_traits_no13C.csv")
  
scores.root <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PCA_scores_root_traits_no13C.csv")
  
  
# Put WFDP code into compatible column name for merging with growth data 
scores.leaf$focal_stem_tag <- scores.leaf$WFDP_Code
  
scores.root$focal_stem_tag <- scores.root$WFDP_Code
  

# Read in scores for PCA of all traits together 

scores.all <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/tree_PC_no13C__scores.csv")


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

## And removing 13C here too 

traits <- dplyr::select(traits, code, WFDP_Code = WFDP_Code.x, sub_plot = sub_plot.x, Host_ID = Host_ID.x, SLA_leaf, LDMC_leaf, LMA_leaf, 
                        leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, specific_root_length, specific_root_area,
                        root_dry_matter_cont, root_CN, root_15N, avg_root_dia, root_pct_N, root_pct_C)


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

# lm(formula = mean_RGR ~ PC1, data = leaf_growth_env)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.020874 -0.008898 -0.003029  0.008704  0.030225 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.0104053  0.0014720   7.069 2.25e-09 ***
#   PC1         0.0003933  0.0006459   0.609    0.545    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0114 on 58 degrees of freedom
# Multiple R-squared:  0.006352,	Adjusted R-squared:  -0.01078 
# F-statistic: 0.3708 on 1 and 58 DF,  p-value: 0.545


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
# ABAM     -0.00650 0.00484 46  -1.343  0.3251
# ABGR     -0.00291 0.00599 46  -0.486  0.7314
# ALRU      0.02718 0.01040 46   2.608  0.0856
# CONU      0.00177 0.00512 46   0.345  0.7314
# TABR     -0.00147 0.00274 46  -0.538  0.7314
# THPL      0.00421 0.00301 46   1.399  0.3251
# TSHE      0.00676 0.00366 46   1.849  0.2479

# No significant relationships for any species 

anova(lm_leaf_PC1_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value   Pr(>F)    
# PC1          1 0.0000482 0.00004821  0.6272  0.43244    
# Species      6 0.0029758 0.00049597  6.4534 5.49e-05 ***
#   PC1:Species  6 0.0010298 0.00017164  2.2334  0.05655 .  
# Residuals   46 0.0035353 0.00007685                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, but doesn't seem like the leaf PC1 values could be driving some of the variation 
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

# lm(formula = mean_RGR ~ PC2, data = leaf_growth_env)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.017079 -0.008167 -0.002532  0.007552  0.033591 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.010405   0.001460   7.125 1.81e-09 ***
#   PC2         0.001843   0.001613   1.143    0.258    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01131 on 58 degrees of freedom
# Multiple R-squared:  0.02201,	Adjusted R-squared:  0.005152 
# F-statistic: 1.306 on 1 and 58 DF,  p-value: 0.2579

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
# ABAM      0.01197 0.00575 46   2.083  0.1783
# ABGR      0.00555 0.01000 46   0.553  0.7356
# ALRU      0.01385 0.01180 46   1.175  0.5737
# CONU     -0.01722 0.00859 46  -2.004  0.1783
# TABR      0.00279 0.00820 46   0.340  0.7356
# THPL     -0.00249 0.00549 46  -0.454  0.7356
# TSHE     -0.00400 0.00581 46  -0.688  0.7356

# No significant relationships for any species 

anova(lm_leaf_PC2_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC2          1 0.0001671 0.00016706  2.0377 0.1601906    
# Species      6 0.0027714 0.00046191  5.6341 0.0001892 ***
#   PC2:Species  6 0.0008793 0.00014656  1.7876 0.1226862    
# Residuals   46 0.0037713 0.00008198                      
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


# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.022055 -0.009030 -0.002947  0.007515  0.030726 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0104053  0.0014304   7.274 1.01e-09 ***
#   PC1         -0.0013488  0.0006904  -1.954   0.0556 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01108 on 58 degrees of freedom
# Multiple R-squared:  0.06174,	Adjusted R-squared:  0.04556 
# F-statistic: 3.817 on 1 and 58 DF,  p-value: 0.05558

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
# ABAM      0.00158 0.00189 46   0.833  0.4806
# ABGR      0.00120 0.00168 46   0.711  0.4806
# ALRU      0.00304 0.00139 46   2.189  0.2358
# CONU     -0.00295 0.00200 46  -1.476  0.3979
# TABR      0.00184 0.00216 46   0.850  0.4806
# THPL     -0.00132 0.00165 46  -0.800  0.4806
# TSHE     -0.00309 0.00222 46  -1.392  0.3979

# No significant relationships for any species 

anova(lm_root_PC1_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0004686 0.00046855  5.7792 0.0203022 *  
#   Species      6 0.0024965 0.00041608  5.1319 0.0004153 ***
#   PC1:Species  6 0.0008946 0.00014910  1.8391 0.1123067    
# Residuals   46 0.0037295 0.00008108                      
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

# lm(formula = mean_RGR ~ PC2, data = root_growth_env)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.018416 -0.008590 -0.003148  0.007428  0.033189 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.010405   0.001423   7.311  8.8e-10 ***
#   PC2         -0.002654   0.001259  -2.108   0.0394 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01102 on 58 degrees of freedom
# Multiple R-squared:  0.07117,	Adjusted R-squared:  0.05516 
# F-statistic: 4.444 on 1 and 58 DF,  p-value: 0.03935

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
# ABAM    -0.004960 0.00342 46  -1.449  0.2695
# ABGR    -0.001881 0.00240 46  -0.784  0.6038
# ALRU     0.008311 0.00422 46   1.972  0.1747
# CONU    -0.002431 0.00373 46  -0.652  0.6038
# TABR    -0.000387 0.00383 46  -0.101  0.9201
# THPL    -0.004421 0.00243 46  -1.823  0.1747
# TSHE     0.015825 0.00540 46   2.932  0.0366

# Significant relationships for only TSHE

anova(lm_root_PC2_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC2          1 0.0005401 0.00054012  7.5296 0.0086252 ** 
#   Species      6 0.0024413 0.00040689  5.6722 0.0001784 ***
#   PC2:Species  6 0.0013080 0.00021799  3.0389 0.0137615 *  
#   Residuals   46 0.0032997 0.00007173                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, and their root PC2 values could be what is driving this. 


## Outcome: Significant result of a relationship between RGR and the Root PC2 results.


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


# Multiple R-squared:  0.02707,	Adjusted R-squared:  0.0103 
# F-statistic: 1.614 on 1 and 58 DF,  p-value: 0.209

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
# ABAM     -0.00422 0.00277 46  -1.523  0.1982
# ABGR     -0.00154 0.00213 46  -0.722  0.4742
# ALRU     -0.00455 0.00266 46  -1.714  0.1982
# CONU      0.00680 0.00379 46   1.794  0.1982
# TABR     -0.00222 0.00240 46  -0.926  0.4192
# THPL      0.00362 0.00242 46   1.496  0.1982
# TSHE      0.00389 0.00218 46   1.782  0.1982

# No significant relationships for any species 

anova(lm_all_PC1_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0002054 0.00020544  2.7012 0.1070912    
# Species      6 0.0027240 0.00045400  5.9693 0.0001133 ***
#   PC1:Species  6 0.0011611 0.00019352  2.5445 0.0327407 *  
#   Residuals   46 0.0034986 0.00007606                      
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


# Multiple R-squared:  0.03697,	Adjusted R-squared:  0.02037 
# F-statistic: 2.227 on 1 and 58 DF,  p-value: 0.1411

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
# 
# Species PC2.trend      SE df t.ratio p.value
# ABAM     0.001008 0.00235 46   0.429  0.6702
# ABGR     0.001710 0.00251 46   0.680  0.6702
# ALRU     0.003696 0.00160 46   2.316  0.1753
# CONU    -0.002820 0.00215 46  -1.310  0.6702
# TABR     0.001522 0.00258 46   0.590  0.6702
# THPL    -0.000874 0.00192 46  -0.454  0.6702
# TSHE    -0.003187 0.00369 46  -0.865  0.6702

# No significant relationships for any species 

anova(lm_all_PC2_rgr_sep)

# 
# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value   Pr(>F)    
# PC2          1 0.0002806 0.00028060  3.3138 0.075211 .  
# Species      6 0.0027142 0.00045236  5.3422 0.000298 ***
#   PC2:Species  6 0.0006992 0.00011653  1.3762 0.244364    
# Residuals   46 0.0038951 0.00008468                                        
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


# RGR_plots <- plot_grid(all_PC1_rgr_all, all_PC2_rgr_all,
#                        leaf_PC1_rgr_all, root_PC1_rgr_all, 
#                        leaf_PC2_rgr_all, root_PC2_rgr_all,
#                                   ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'))
# 
# RGR_plots


# Plot with just the leaf and root traits: 

RGR_plots_2 <- plot_grid(leaf_PC1_rgr_all, root_PC1_rgr_all, 
                         leaf_PC2_rgr_all, root_PC2_rgr_all,
                         ncol = 2, nrow = 2, labels = c('(a)', '(b)', '(c)', '(d)'))

RGR_plots_2 


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/trait_RGR_plots.png", 
       plot = RGR_plots_2, width = 7, height = 6, units = "in", dpi = 300)



# Separate trait~RGR plots - EDIT: removing 13C from these 

leaf_trait_RGR_plots <- plot_grid(leaf_sla_rgr, leaf_lma_rgr, leaf_ldmc_rgr, leaf_pctC_rgr, leaf_pctN_rgr, 
                             leaf_CN_rgr, leaf_15N_rgr,
                             ncol = 3, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)'))

leaf_trait_RGR_plots


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/leaf_trait_RGR_plots.png", 
       plot = leaf_trait_RGR_plots, width = 12, height = 11, units = "in", dpi = 300)



root_trait_RGR_plots <- plot_grid(root_sra_rgr, root_srl_rgr, root_rdmc_rgr, root_pctC_rgr, root_pctN_rgr, 
                                  root_CN_rgr, root_15N_rgr, root_dia_rgr, 
                                  ncol = 3, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'))

root_trait_RGR_plots


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/root_trait_RGR_plots.png", 
       plot = root_trait_RGR_plots, width = 12, height = 11, units = "in", dpi = 300)



## -- END -- ## 

