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

#################################################################################
#                               Main workflow                                   #
#  Explore the growth data for all of the trait trees sampled at WFDP. Assess   #
#  for differences between species, or any relationships to environmental       #
#  variables so far. Then relate to trait information for the trees themselves. #
#  Use the values of the leaf and root PC1 and PC2 values for each focal tree   # 
#  to explore relationships of trait strategies with focal trees growth.        #
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
# 2016, and 2021. I want to collect the diameters from the first and last census year 
# to calculate how much they have grown. 


# get starting and ending diameters for all 60 trees
diams <- growth %>%
  group_by(Stem_Tag, Species) %>%
  summarise(
    dia_first = DBH[which.min(DBH_DATE)],
    dia_last = DBH[which.max(DBH_DATE)],
    year_first = min(DBH_DATE),
    year_last = max(DBH_DATE),
    .groups = "drop"
  )

# get diameter difference between the two time points 
diams <- diams %>%
  mutate(diam_diff = dia_last - dia_first)

# Calculate relative growth rate for each tree 
diams <- diams %>%
  mutate(RGR = (log(dia_last) - log(dia_first)) / (year_last - year_first))


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

#################################################################################

###################################### -- 
# (2) RGR RELATIONSHIP TO TOPOGRAPHY
###################################### -- 

# Considering relationships of both the diameter difference between the first and last time points, 
# and the RGR for each focal tree to slope, aspect, and elevation of each tree. 


# Just using the leaf_growth_env df here because this is just looking at the environmental data, and 
# nothing related to the PCA results yet 


# set colors for hosts 
                # ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#9b5fe0", "#16a4d8", "#60dbe8", "#8bd346","#efdf48", "#f9a52F", "#d64e12")


# Exploring separate slopes per species here, because we would expect that growth will vary between species, 
# and may also vary individually with topography. Also testing the linear regressions with an interaction between 
# the topographic variable and species, to get the individual slopes for each species for the relationship 
# between growth and the topographic variable. 


# Explore diameter difference relationships

## Slope
slope <- ggplot(leaf_growth_env, aes(x = slope, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Slope", y = expression("Diameter Difference (mm)")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

slope

# test relationships 
lm_slope <- lm(diam_diff ~ slope * Species, data = leaf_growth_env)

emtrends(lm_slope, ~ Species, var = "slope") %>% test(adjust = "fdr")
# This gets separate p-values for the relationship for each species 
# P-values were adjusted for multiple testing using the Benjamini–Hochberg 
# false discovery rate procedure.

# 
# Species slope.trend    SE df t.ratio p.value
# ABAM        0.03774 0.308 46   0.122  0.9893
# ABGR       -0.21533 1.000 46  -0.215  0.9893
# ALRU       -0.95355 0.641 46  -1.488  0.5027
# CONU       -0.00788 0.587 46  -0.013  0.9893
# TABR        0.16403 0.287 46   0.571  0.9893
# THPL        0.39178 0.219 46   1.790  0.5027
# TSHE       -0.16793 0.410 46  -0.409  0.9893

# No species has a significant relationship between slope and diam diff


anova(lm_slope)

# Response: diam_diff
# Df Sum Sq Mean Sq F value    Pr(>F)    
# slope          1  3.955  3.9549  3.5823 0.0646994 .  
# Species        6 32.372  5.3953  4.8870 0.0006143 ***
# slope:Species  6  5.622  0.9370  0.8487 0.5393464    
# Residuals     46 50.784  1.1040                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Diam diff is different between species, but slope is not what is driving this. 



## Aspect
aspect <- ggplot(leaf_growth_env, aes(x = aspect, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Aspect", y = expression("Diameter Difference (mm)")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

aspect

# test relationships 
lm_aspect <- lm(diam_diff ~ aspect * Species, data = leaf_growth_env)

emtrends(lm_aspect, ~ Species, var = "aspect") %>% test(adjust = "fdr")

# Species aspect.trend      SE df t.ratio p.value
# ABAM       -2.74e-03 0.01100 46  -0.249  0.9945
# ABGR       -1.24e-02 0.01100 46  -1.120  0.9945
# ALRU       -3.47e-03 0.02010 46  -0.173  0.9945
# CONU        4.51e-05 0.00650 46   0.007  0.9945
# TABR        9.66e-03 0.01650 46   0.584  0.9945
# THPL        2.05e-03 0.00362 46   0.567  0.9945
# TSHE       -2.57e-03 0.01390 46  -0.185  0.9945

# No species has a significant relationship between aspect and diam diff


anova(lm_aspect)


# Response: diam_diff
# Df Sum Sq Mean Sq F value    Pr(>F)    
# aspect          1  0.078  0.0775  0.0649 0.8000355    
# Species         6 35.297  5.8828  4.9257 0.0005772 ***
#   aspect:Species  6  2.421  0.4035  0.3378 0.9133168    
# Residuals      46 54.938  1.1943                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Diam diff is different between species, but aspect is not what is driving this. 


## Elevation
elev <- ggplot(leaf_growth_env, aes(x = elevation_m, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Elevation (m)", y = expression("Diameter Difference (mm)")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

elev

# test relationships 
lm_elev <- lm(diam_diff ~ elevation_m * Species, data = leaf_growth_env)

emtrends(lm_elev, ~ Species, var = "elevation_m") %>% test(adjust = "fdr")

# Species elevation_m.trend     SE df t.ratio p.value
# ABAM            -8.49e-02 0.0388 46  -2.186  0.2373
# ABGR             2.48e-02 0.0493 46   0.503  1.0000
# ALRU            -2.14e-06 0.0580 46   0.000  1.0000
# CONU             1.30e-03 0.0599 46   0.022  1.0000
# TABR             1.45e-02 0.0465 46   0.312  1.0000
# THPL             6.55e-02 0.0443 46   1.479  0.5110
# TSHE             1.42e-03 0.0363 46   0.039  1.0000

# No species has a significant relationship between elevation and diam diff

anova(lm_elev)

# Response: diam_diff
# Df Sum Sq Mean Sq F value    Pr(>F)    
# elevation_m          1  1.235  1.2345  1.1472 0.2897248    
# Species              6 34.126  5.6876  5.2853 0.0003259 ***
#   elevation_m:Species  6  7.871  1.3118  1.2190 0.3139835    
# Residuals           46 49.502  1.0761                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Diam diff is different between species, but elevation is not what is driving this. 


### Explore RGR relationships 

## Slope
slope2 <- ggplot(leaf_growth_env, aes(x = slope, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Slope", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

slope2

# test relationships 
lm_slope2 <- lm(RGR ~ slope * Species, data = leaf_growth_env)

emtrends(lm_slope2, ~ Species, var = "slope") %>% test(adjust = "fdr")


# Species slope.trend      SE df t.ratio p.value
# ABAM       0.000472 0.00275 46   0.172  0.9601
# ABGR       0.000888 0.00895 46   0.099  0.9601
# ALRU      -0.013852 0.00572 46  -2.423  0.1355
# CONU      -0.000263 0.00524 46  -0.050  0.9601
# TABR       0.000970 0.00256 46   0.379  0.9601
# THPL       0.002017 0.00195 46   1.033  0.9601
# TSHE       0.001428 0.00366 46   0.390  0.9601


anova(lm_slope2)

# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# slope          1 0.0001545 0.00015448  1.7594 0.1912414    
# Species        6 0.0027638 0.00046063  5.2464 0.0003465 ***
#   slope:Species  6 0.0006166 0.00010276  1.1704 0.3385485    
# Residuals     46 0.0040388 0.00008780                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



## Aspect
aspect2 <- ggplot(leaf_growth_env, aes(x = aspect, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Aspect", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

aspect2


# test relationships 
lm_aspect2 <- lm(RGR ~ aspect * Species, data = leaf_growth_env)

emtrends(lm_aspect2, ~ Species, var = "aspect") %>% test(adjust = "fdr")


# Species aspect.trend       SE df t.ratio p.value
# ABAM       -7.42e-05 9.97e-05 46  -0.744  0.9214
# ABGR       -4.77e-05 1.00e-04 46  -0.475  0.9214
# ALRU        1.93e-05 1.83e-04 46   0.105  0.9214
# CONU        1.46e-05 5.91e-05 46   0.247  0.9214
# TABR        6.85e-05 1.50e-04 46   0.456  0.9214
# THPL        1.82e-05 3.29e-05 46   0.551  0.9214
# TSHE        1.25e-05 1.26e-04 46   0.099  0.9214


anova(lm_aspect2)

# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# aspect          1 0.0000017 0.00000171  0.0174 0.8957652    
# Species         6 0.0029051 0.00048418  4.9033 0.0005984 ***
#   aspect:Species  6 0.0001246 0.00002076  0.2102 0.9717906    
# Residuals      46 0.0045423 0.00009875                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Elevation
elev2 <- ggplot(leaf_growth_env, aes(x = elevation_m, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Elevation (m)", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

elev2


# test relationships 
lm_elev2 <- lm(RGR ~ elevation_m * Species, data = leaf_growth_env)

emtrends(lm_elev2, ~ Species, var = "elevation_m") %>% test(adjust = "fdr")

# Species elevation_m.trend       SE df t.ratio p.value
# ABAM            -6.85e-04 0.000357 46  -1.920  0.4277
# ABGR             1.55e-04 0.000452 46   0.343  0.9876
# ALRU            -3.65e-04 0.000533 46  -0.685  0.9876
# CONU             8.63e-06 0.000550 46   0.016  0.9876
# TABR             8.41e-05 0.000427 46   0.197  0.9876
# THPL             4.49e-04 0.000406 46   1.105  0.9616
# TSHE            -5.96e-05 0.000333 46  -0.179  0.9876

anova(lm_elev2)

# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# elevation_m          1 0.0001646 0.00016458  1.8142 0.1846000    
# Species              6 0.0027549 0.00045915  5.0614 0.0004646 ***
#   elevation_m:Species  6 0.0004812 0.00008020  0.8840 0.5142559    
# Residuals           46 0.0041730 0.00009072                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Outcome: Same general trends between diam_diff and RGR, species differ in their growth rates, 
# but these are not being driven by any of the topographic variables. 


#################################################################################

################################################ -- 
# (3) GROWTH RELATIONSHIPS TO TRAIT PCA VALUES
################################################ -- 


# Explore diameter difference relationships with the PC1 and PC2 values for the leaf and 
# root trait PCAs


## LEAF TRAITS

# PC1 - All species together 
leaf_PC1_diam_all <- ggplot(leaf_growth_env, aes(x = PC1, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC1 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC1_diam_all


# test relationships 
lm_leaf_PC1_diam_all <- lm(diam_diff ~ PC1, data = leaf_growth_env)

summary(lm_leaf_PC1_diam_all)


# Multiple R-squared:  0.01721,	Adjusted R-squared:  0.0002624 
# F-statistic: 1.015 on 1 and 58 DF,  p-value: 0.3178

# No significant relationship for the species all together 


# PC1 - For separate species
leaf_PC1_diam_sep <- ggplot(leaf_growth_env, aes(x = PC1, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC1 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

leaf_PC1_diam_sep


# test relationships 
lm_leaf_PC1_diam_sep <- lm(diam_diff ~ PC1 * Species, data = leaf_growth_env)

emtrends(lm_leaf_PC1_diam_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend    SE df t.ratio p.value
# ABAM      -0.6960 0.562 46  -1.238  0.4818
# ABGR      -0.8160 0.692 46  -1.179  0.4818
# ALRU       2.8649 1.240 46   2.302  0.1813
# CONU       0.1687 0.578 46   0.292  0.8089
# TABR      -0.2465 0.316 46  -0.780  0.6150
# THPL       0.0841 0.346 46   0.243  0.8089
# TSHE       0.4662 0.422 46   1.104  0.4818

# No significant relationships for any species 

anova(lm_leaf_PC1_diam_sep)

# 
# Response: diam_diff
# Response: diam_diff
# Df Sum Sq Mean Sq F value    Pr(>F)    
# PC1          1  1.596  1.5957  1.5627 0.2175881    
# Species      6 33.773  5.6289  5.5127 0.0002283 ***
#   PC1:Species  6 10.394  1.7324  1.6967 0.1433046    
# Residuals   46 46.969  1.0211                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Diam diff is different between species, but their leaf PC1 values are not what is driving this. 


## PC2 - All species together 
leaf_PC2_diam_all <- ggplot(leaf_growth_env, aes(x = PC2, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC2 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

leaf_PC2_diam_all


# test relationships 
lm_leaf_PC2_diam_all <- lm(diam_diff ~ PC2, data = leaf_growth_env)

summary(lm_leaf_PC2_diam_all)

# Multiple R-squared:  0.03057,	Adjusted R-squared:  0.01385 
# F-statistic: 1.829 on 1 and 58 DF,  p-value: 0.1815

# No significant relationship for the species all together 


# PC2 - For separate species
leaf_PC2_diam_sep <- ggplot(leaf_growth_env, aes(x = PC2, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC2 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

leaf_PC2_diam_sep


# test relationships 
lm_leaf_PC2_diam_sep <- lm(diam_diff ~ PC2 * Species, data = leaf_growth_env)

emtrends(lm_leaf_PC2_diam_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend    SE df t.ratio p.value
# ABAM       0.4626 0.606 46   0.764  0.8833
# ABGR       0.2071 1.400 46   0.148  0.8833
# ALRU      -1.0243 1.130 46  -0.906  0.8833
# CONU      -0.3234 0.695 46  -0.465  0.8833
# TABR      -0.2015 0.533 46  -0.378  0.8833
# THPL       0.0838 0.507 46   0.165  0.8833
# TSHE      -0.2193 0.577 46  -0.380  0.8833

# No significant relationships for any species 

anova(lm_leaf_PC2_diam_sep)

# 
# Response: diam_diff
# Df Sum Sq Mean Sq F value   Pr(>F)   
# PC2          1  2.834  2.8344  2.3690 0.130615   
# Species      6 32.614  5.4357  4.5432 0.001073 **
#   PC2:Species  6  2.247  0.3745  0.3130 0.926966   
# Residuals   46 55.037  1.1964                    
# ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Diam diff is different between species, but their leaf PC2 values are not what is driving this. 



## ROOT TRAITS

# PC1 - All species together 
root_PC1_diam_all <- ggplot(root_growth_env, aes(x = PC1, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC1 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

root_PC1_diam_all


# test relationships 
lm_root_PC1_diam_all <- lm(diam_diff ~ PC1, data = root_growth_env)

summary(lm_root_PC1_diam_all)


# Multiple R-squared:  0.02163,	Adjusted R-squared:  0.004765 
# F-statistic: 1.282 on 1 and 58 DF,  p-value: 0.2621

# No significant relationship for the species all together 


# PC1 - For separate species
root_PC1_diam_sep <- ggplot(root_growth_env, aes(x = PC1, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC1 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

root_PC1_diam_sep


# test relationships 
lm_root_PC1_diam_sep <- lm(diam_diff ~ PC1 * Species, data = root_growth_env)

emtrends(lm_root_PC1_diam_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend    SE df t.ratio p.value
# ABAM       0.1512 0.216 46   0.699  0.5695
# ABGR       0.3020 0.190 46   1.586  0.5368
# ALRU       0.2278 0.162 46   1.408  0.5368
# CONU      -0.2121 0.233 46  -0.910  0.5368
# TABR       0.2694 0.249 46   1.082  0.5368
# THPL      -0.1687 0.192 46  -0.880  0.5368
# TSHE      -0.0709 0.253 46  -0.280  0.7808

# No significant relationships for any species 

anova(lm_root_PC1_diam_sep)

# 
# Response: diam_diff
# Df Sum Sq Mean Sq F value    Pr(>F)    
# PC1          1  2.006  2.0061  1.8823 0.1767234    
# Species      6 34.806  5.8011  5.4431 0.0002545 ***
#   PC1:Species  6  6.894  1.1491  1.0782 0.3894524    
# Residuals   46 49.026  1.0658                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Diam diff is different between species, but their root PC1 values are not what is driving this. 


## PC2 - All species together 
root_PC2_diam_all <- ggplot(root_growth_env, aes(x = PC2, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC2 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

root_PC2_diam_all


# test relationships 
lm_root_PC2_diam_all <- lm(diam_diff ~ PC2, data = root_growth_env)

summary(lm_root_PC2_diam_all)


# Multiple R-squared:  0.09472,	Adjusted R-squared:  0.07911 
# F-statistic: 6.069 on 1 and 58 DF,  p-value: 0.01675

# SIGNFIICANT RELATIONSHIP for the species all together 


# PC2 - For separate species
root_PC2_diam_sep <- ggplot(root_growth_env, aes(x = PC2, y = diam_diff, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC2 Value", y = "Diameter Difference (mm)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

root_PC2_diam_sep


# test relationships 
lm_root_PC2_diam_sep <- lm(diam_diff ~ PC2 * Species, data = root_growth_env)

emtrends(lm_root_PC2_diam_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend    SE df t.ratio p.value
# ABAM        0.835 0.294 46   2.840  0.0468
# ABGR        0.871 0.485 46   1.797  0.1842
# ALRU       -0.216 0.236 46  -0.917  0.5325
# CONU        0.137 0.280 46   0.491  0.7146
# TABR       -0.215 0.584 46  -0.368  0.7146
# THPL        0.508 0.281 46   1.809  0.1842
# TSHE        0.296 0.334 46   0.886  0.5325

# SIGNIFICANT RELATIONSHIP FOR ABAM for any species 

anova(lm_root_PC2_diam_sep)

# 
# Response: diam_diff
# Df Sum Sq Mean Sq F value    Pr(>F)    
# PC2          1  8.784  8.7837  9.5781 0.0033464 ** 
# Species      6 31.566  5.2610  5.7368 0.0001615 ***
# PC2:Species  6 10.198  1.6997  1.8534 0.1095739    
# Residuals   46 42.185  0.9171                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Diam diff is different between species, and PC2 values are also different between species. 



### RGR ####

## LEAF TRAITS

# PC1 - All species together 
leaf_PC1_rgr_all <- ggplot(leaf_growth_env, aes(x = PC1, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC1 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC1_rgr_all


# test relationships 
lm_leaf_PC1_rgr_all <- lm(RGR ~ PC1, data = leaf_growth_env)

summary(lm_leaf_PC1_rgr_all)


# Multiple R-squared:  0.006755,	Adjusted R-squared:  -0.01037 
# F-statistic: 0.3945 on 1 and 58 DF,  p-value: 0.5324

# No significant relationship for the species all together 


# PC1 - For separate species
leaf_PC1_rgr_sep <- ggplot(leaf_growth_env, aes(x = PC1, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC1 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC1_rgr_sep


# test relationships 
lm_leaf_PC1_rgr_sep <- lm(RGR ~ PC1 * Species, data = leaf_growth_env)

emtrends(lm_leaf_PC1_rgr_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend      SE df t.ratio p.value
# ABAM     -0.00694 0.00484 46  -1.434  0.2839
# ABGR     -0.00297 0.00597 46  -0.497  0.7068
# ALRU      0.02892 0.01070 46   2.697  0.0682
# CONU      0.00189 0.00498 46   0.378  0.7068
# TABR     -0.00156 0.00272 46  -0.573  0.7068
# THPL      0.00423 0.00298 46   1.420  0.2839
# TSHE      0.00692 0.00364 46   1.902  0.2219

# No significant relationships for any species 

anova(lm_leaf_PC1_rgr_sep)

# 
# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0000512 0.00005116  0.6746   0.41568    
# Species      6 0.0029377 0.00048961  6.4561 5.468e-05 ***
#   PC1:Species  6 0.0010963 0.00018272  2.4094   0.04151 *  
#   Residuals   46 0.0034885 0.00007584                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species,and the leaf PC1 values could be driving some of the variation 
# between the species 


## PC2 - All species together 
leaf_PC2_rgr_all <- ggplot(leaf_growth_env, aes(x = PC2, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC2_rgr_all


# test relationships 
lm_leaf_PC2_rgr_all <- lm(RGR ~ PC2, data = leaf_growth_env)

summary(lm_leaf_PC2_rgr_all)


# Multiple R-squared:  0.01543,	Adjusted R-squared:  -0.00154 
# F-statistic: 0.9093 on 1 and 58 DF,  p-value: 0.3443

# No significant relationship for the species all together 


# PC2 - For separate species
leaf_PC2_rgr_sep <- ggplot(leaf_growth_env, aes(x = PC2, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

leaf_PC2_rgr_sep


# test relationships 
lm_leaf_PC2_rgr_sep <- lm(RGR ~ PC2 * Species, data = leaf_growth_env)

emtrends(lm_leaf_PC2_rgr_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend      SE df t.ratio p.value
# ABAM     0.004801 0.00548 46   0.875  0.9610
# ABGR    -0.000625 0.01270 46  -0.049  0.9610
# ALRU    -0.001299 0.01020 46  -0.127  0.9610
# CONU    -0.004106 0.00630 46  -0.652  0.9610
# TABR    -0.001492 0.00483 46  -0.309  0.9610
# THPL    -0.001833 0.00459 46  -0.399  0.9610
# TSHE    -0.002370 0.00523 46  -0.453  0.9610

# No significant relationships for any species 

anova(lm_leaf_PC2_rgr_sep)

# 
# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC2          1 0.0001169 0.00011690  1.1913 0.2807571    
# Species      6 0.0027980 0.00046633  4.7522 0.0007637 ***
#   PC2:Species  6 0.0001448 0.00002414  0.2460 0.9585161    
# Residuals   46 0.0045140 0.00009813                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, but their leaf PC2 values are not what is driving this. 



## ROOT TRAITS

# PC1 - All species together 
root_PC1_rgr_all <- ggplot(root_growth_env, aes(x = PC1, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC1 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC1_rgr_all


# test relationships 
lm_root_PC1_rgr_all <- lm(RGR ~ PC1, data = root_growth_env)

summary(lm_root_PC1_rgr_all)


# Multiple R-squared:  0.05721,	Adjusted R-squared:  0.04095 
# F-statistic: 3.519 on 1 and 58 DF,  p-value: 0.06569

# No significant relationship for the species all together 


# PC1 - For separate species
root_PC1_rgr_sep <- ggplot(root_growth_env, aes(x = PC1, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC1 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC1_rgr_sep


# test relationships 
lm_root_PC1_rgr_sep <- lm(RGR ~ PC1 * Species, data = root_growth_env)

emtrends(lm_root_PC1_rgr_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend      SE df t.ratio p.value
# ABAM      0.00190 0.00189 46   1.009  0.4806
# ABGR      0.00118 0.00166 46   0.711  0.4806
# ALRU      0.00316 0.00141 46   2.240  0.2096
# CONU     -0.00300 0.00203 46  -1.476  0.4572
# TABR      0.00182 0.00217 46   0.841  0.4806
# THPL     -0.00126 0.00167 46  -0.756  0.4806
# TSHE     -0.00290 0.00221 46  -1.312  0.4572

# No significant relationships for any species 

anova(lm_root_PC1_rgr_sep)

# 
# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0004333 0.00043328  5.3463 0.0252946 *  
#   Species      6 0.0025134 0.00041890  5.1689 0.0003917 ***
#   PC1:Species  6 0.0008990 0.00014983  1.8488 0.1104378    
# Residuals   46 0.0037280 0.00008104                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, and PC1 values are also different. 


## PC2 - All species together 
root_PC2_rgr_all <- ggplot(root_growth_env, aes(x = PC2, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC2_rgr_all


# test relationships 
lm_root_PC2_rgr_all <- lm(RGR ~ PC2, data = root_growth_env)

summary(lm_root_PC2_rgr_all)


# Multiple R-squared:  0.08763,	Adjusted R-squared:  0.0719 
# F-statistic: 5.571 on 1 and 58 DF,  p-value: 0.02165

# SIGNIFICANT RELATIONSHIP for the species all together 


# PC2 - For separate species
root_PC2_rgr_sep <- ggplot(root_growth_env, aes(x = PC2, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

root_PC2_rgr_sep


# test relationships 
lm_root_PC2_rgr_sep <- lm(RGR ~ PC2 * Species, data = root_growth_env)

emtrends(lm_root_PC2_rgr_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend      SE df t.ratio p.value
# ABAM     0.007260 0.00269 46   2.701  0.0674
# ABGR     0.003506 0.00443 46   0.791  0.6058
# ALRU    -0.000314 0.00215 46  -0.146  0.8846
# CONU     0.002275 0.00256 46   0.890  0.6058
# TABR    -0.001554 0.00533 46  -0.291  0.8846
# THPL     0.005387 0.00257 46   2.099  0.1446
# TSHE     0.004110 0.00306 46   1.345  0.4319

# No significant relationships for any species 

anova(lm_root_PC2_rgr_sep)

# 
# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC2          1 0.0006637 0.00066369  8.6600  0.005082 ** 
#   Species      6 0.0028723 0.00047872  6.2464 7.465e-05 ***
#   PC2:Species  6 0.0005123 0.00008539  1.1142  0.368930    
# Residuals   46 0.0035254 0.00007664                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# RGR is different between species, but their root PC2 values are not what is driving this. 


## Outcome: Significant result of a relationship between RGR and the Root PC2 results, 
# but it is no longer significant when the species are considered independently. Leaf PC1 is 
# significant when species are considered independently. 


#################################################################################

#################################################### -- 
# (4) GROWTH RELATIONSHIPS TO ALL TRAIT PCA VALUES
#################################################### -- 

# Same analyses as above to explore the relationships of RGR with the tree traits, but using the PCA 
# scores for all traits considered together. 


# PC1 - All species together 
all_PC1_rgr_all <- ggplot(all_traits_growth_env, aes(x = PC1, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC1 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_all_PC1_rgr_all <- lm(RGR ~ PC1, data = all_traits_growth_env)

summary(lm_all_PC1_rgr_all)


# Multiple R-squared:  0.02399,	Adjusted R-squared:  0.007158 
# F-statistic: 1.425 on 1 and 58 DF,  p-value: 0.2374

# No significant relationship for the species all together 


# PC1 - For separate species
all_PC1_rgr_sep <- ggplot(all_traits_growth_env, aes(x = PC1, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC1 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_all_PC1_rgr_sep <- lm(RGR ~ PC1 * Species, data = all_traits_growth_env)

emtrends(lm_all_PC1_rgr_sep, ~ Species, var = "PC1") %>% test(adjust = "fdr")

# Species PC1.trend      SE df t.ratio p.value
# ABAM     -0.00473 0.00264 46  -1.792  0.1814
# ABGR     -0.00149 0.00209 46  -0.714  0.4788
# ALRU     -0.00500 0.00273 46  -1.832  0.1814
# CONU      0.00644 0.00378 46   1.702  0.1814
# TABR     -0.00225 0.00240 46  -0.935  0.4136
# THPL      0.00368 0.00255 46   1.440  0.2194
# TSHE      0.00361 0.00218 46   1.660  0.1814

# No significant relationships for any species 

anova(lm_all_PC1_rgr_sep)

# 
# Response: RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# PC1          1 0.0001817 0.00018166  2.3951 0.1285695    
# Species      6 0.0027204 0.00045340  5.9779 0.0001118 ***
#   PC1:Species  6 0.0011827 0.00019711  2.5988 0.0297576 *  
#   Residuals   46 0.0034889 0.00007585                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# SIGNIFICANT: RGR is different between species, and their PC1 values may be part of what is driving this. 


## PC2 - All species together 
all_PC2_rgr_all <- ggplot(all_traits_growth_env, aes(x = PC2, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_all_PC2_rgr_all <- lm(RGR ~ PC2, data = all_traits_growth_env)

summary(lm_all_PC2_rgr_all)


# Multiple R-squared:  0.03666,	Adjusted R-squared:  0.02005 
# F-statistic: 2.207 on 1 and 58 DF,  p-value: 0.1428

# No significant relationship for species all together 


# PC2 - For separate species
all_PC2_rgr_sep <- ggplot(all_traits_growth_env, aes(x = PC2, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Combined Traits PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_all_PC2_rgr_sep <- lm(RGR ~ PC2 * Species, data = all_traits_growth_env)

emtrends(lm_all_PC2_rgr_sep, ~ Species, var = "PC2") %>% test(adjust = "fdr")

# Species PC2.trend      SE df t.ratio p.value
# ABAM     0.001023 0.00236 46   0.435  0.6659
# ABGR     0.001740 0.00253 46   0.689  0.6659
# ALRU     0.003627 0.00159 46   2.281  0.1908
# CONU    -0.002806 0.00213 46  -1.316  0.6659
# TABR     0.001493 0.00261 46   0.572  0.6659
# THPL    -0.000918 0.00190 46  -0.482  0.6659
# TSHE    -0.003399 0.00368 46  -0.924  0.6659

# No significant relationships for any species 

anova(lm_all_PC2_rgr_sep)

# 
# Response: RGR
# Df    Sum Sq    Mean Sq F value   Pr(>F)    
# PC2          1 0.0002777 0.00027767  3.2652 0.077308 .  
# Species      6 0.0026755 0.00044592  5.2437 0.000348 ***
#   PC2:Species  6 0.0007086 0.00011811  1.3888 0.239408    
# Residuals   46 0.0039118 0.00008504                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# SIGNIFICANT: RGR is different between species, but the interaction of these is not what is driving variation in RGR


#################################################################################

############################################## -- 
# (4) GATHER SIGNIFICANT RESULTS OF INTEREST  
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



## -- END -- ## 



