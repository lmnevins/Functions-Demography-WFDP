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


# set colors for hosts 
                # ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# Visualize mean RGR between the focal tree species 

spp_rgr <- ggplot(diams, aes(y = mean_RGR, x = Species, fill = Species)) +
  geom_boxplot() +
  geom_point(alpha = 0.6) +
  theme_bw() +
  scale_fill_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

spp_rgr


#summaries 

rgr_spp_aov <- aov(mean_RGR ~ Species, data = diams)

summary(rgr_spp_aov)

# RGR significantly differs between species 
# 
#           Df   Sum Sq   Mean Sq F value         Pr(>F)    
#     Species      6 0.002929 0.0004882   5.553 0.00016 ***
#   Residuals   53 0.004660 0.0000879                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


tuk_rgr_spp <- TukeyHSD(rgr_spp_aov)
tuk_rgr_spp

# 
# diff           lwr           upr     p adj
# ABGR-ABAM -6.952969e-03 -0.0217905642  0.0078846261 0.7799461
# ALRU-ABAM  1.022609e-02 -0.0055115590  0.0259637335 0.4322334
# CONU-ABAM -4.838083e-03 -0.0176878172  0.0080116516 0.9079049
# TABR-ABAM -1.493780e-02 -0.0281396491 -0.0017359550 0.0170683
# THPL-ABAM  4.200038e-03 -0.0086496967  0.0170497720 0.9514257
# TSHE-ABAM -4.903199e-03 -0.0177529337  0.0079465350 0.9024389
# ALRU-ABGR  1.717906e-02 -0.0002195662  0.0345776789 0.0550804
# CONU-ABGR  2.114886e-03 -0.0127227089  0.0169524814 0.9994224
# TABR-ABGR -7.984833e-03 -0.0231283902  0.0071587242 0.6729593
# THPL-ABGR  1.115301e-02 -0.0036845885  0.0259906019 0.2616727
# TSHE-ABGR  2.049770e-03 -0.0127878255  0.0168873649 0.9995173
# CONU-ALRU -1.506417e-02 -0.0308018163  0.0006734762 0.0690773
# TABR-ALRU -2.516389e-02 -0.0411903238 -0.0091374548 0.0002476
# THPL-ALRU -6.026050e-03 -0.0217636959  0.0097115966 0.9009760
# TSHE-ALRU -1.512929e-02 -0.0308669329  0.0006083596 0.0669937
# TABR-CONU -1.009972e-02 -0.0233015663  0.0031021278 0.2429691
# THPL-CONU  9.038120e-03 -0.0038116139  0.0218878548 0.3367882
# TSHE-CONU -6.511657e-05 -0.0129148509  0.0127846178 1.0000000
# THPL-TABR  1.913784e-02  0.0059359926  0.0323396868 0.0008547
# TSHE-TABR  1.003460e-02 -0.0031672444  0.0232364498 0.2497722
# TSHE-THPL -9.103237e-03 -0.0219529714  0.0037464973 0.3283815


# TABR is different from ABAM, ALRU, THPL

#################################################################################

###################################### -- 
# (2) RGR RELATIONSHIP TO TOPOGRAPHY
###################################### -- 

# Considering relationships of both the diameter difference between the first and last time points, 
# and the RGR for each focal tree to slope, aspect, and elevation of each tree. 


# Just using the leaf_growth_env df here because this is just looking at the environmental data, and 
# nothing related to the PCA results yet 



# Exploring separate slopes per species here, because we would expect that growth will vary between species, 
# and may also vary individually with topography. Also testing the linear regressions with an interaction between 
# the topographic variable and species, to get the individual slopes for each species for the relationship 
# between growth and the topographic variable. 


### Explore RGR relationships 

## Slope
slope2 <- ggplot(leaf_growth_env, aes(x = slope, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Slope", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_slope2 <- lm(mean_RGR ~ slope * Species, data = leaf_growth_env)

emtrends(lm_slope2, ~ Species, var = "slope") %>% test(adjust = "fdr")


# Species slope.trend      SE df t.ratio p.value
# ABAM       0.000417 0.00275 46   0.152  0.9599
# ABGR       0.000810 0.00894 46   0.091  0.9599
# ALRU      -0.013660 0.00571 46  -2.391  0.1467
# CONU      -0.000265 0.00524 46  -0.051  0.9599
# TABR       0.000967 0.00256 46   0.378  0.9599
# THPL       0.002027 0.00195 46   1.039  0.9599
# TSHE       0.001395 0.00366 46   0.381  0.9599


anova(lm_slope2)

# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# slope          1 0.0001564 0.00015639  1.7828 0.1883771    
# Species        6 0.0027952 0.00046587  5.3106 0.0003132 ***
#   slope:Species  6 0.0006022 0.00010036  1.1441 0.3525130    
# Residuals     46 0.0040353 0.00008772                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



## Aspect
aspect2 <- ggplot(leaf_growth_env, aes(x = aspect, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Aspect", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_aspect2 <- lm(mean_RGR ~ aspect * Species, data = leaf_growth_env)

emtrends(lm_aspect2, ~ Species, var = "aspect") %>% test(adjust = "fdr")


# Species aspect.trend       SE df t.ratio p.value
# ABAM       -6.95e-05 9.96e-05 46  -0.697  0.9126
# ABGR       -4.72e-05 1.00e-04 46  -0.471  0.9126
# ALRU        2.47e-05 1.82e-04 46   0.135  0.9126
# CONU        1.37e-05 5.91e-05 46   0.233  0.9126
# TABR        7.05e-05 1.50e-04 46   0.469  0.9126
# THPL        1.77e-05 3.29e-05 46   0.537  0.9126
# TSHE        1.39e-05 1.26e-04 46   0.110  0.9126


anova(lm_aspect2)

# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# aspect          1 0.0000018 0.00000181  0.0184 0.8926854    
# Species         6 0.0029389 0.00048981  4.9722 0.0005358 ***
#   aspect:Species  6 0.0001169 0.00001948  0.1978 0.9758017    
# Residuals      46 0.0045315 0.00009851                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Elevation
elev2 <- ggplot(leaf_growth_env, aes(x = elevation_m, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Elevation (m)", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
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
lm_elev2 <- lm(mean_RGR ~ elevation_m * Species, data = leaf_growth_env)

emtrends(lm_elev2, ~ Species, var = "elevation_m") %>% test(adjust = "fdr")

# Species elevation_m.trend       SE df t.ratio p.value
# ABAM            -6.68e-04 0.000356 46  -1.874  0.4714
# ABGR             1.54e-04 0.000452 46   0.340  0.9988
# ALRU            -3.75e-04 0.000532 46  -0.704  0.9988
# CONU            -8.62e-07 0.000550 46  -0.002  0.9988
# TABR             8.30e-05 0.000427 46   0.195  0.9988
# THPL             4.52e-04 0.000406 46   1.113  0.9507
# TSHE            -5.75e-05 0.000333 46  -0.173  0.9988

anova(lm_elev2)

# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# elevation_m          1 0.0001623 0.00016235  1.7918 0.1872884    
# Species              6 0.0027895 0.00046492  5.1312 0.0004158 ***
#   elevation_m:Species  6 0.0004694 0.00007823  0.8634 0.5288464    
# Residuals           46 0.0041679 0.00009061                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Outcome: Same general trends between diam_diff and RGR, species differ in their growth rates, 
# but these are not being driven by any of the topographic variables. 


#################################################################################

################################################ -- 
# (3) GROWTH RELATIONSHIPS TO TRAIT PCA VALUES
################################################ -- 

### RGR ####

## LEAF TRAITS

# PC1 - All species together 
leaf_PC1_rgr_all <- ggplot(leaf_growth_env, aes(x = PC1, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC1 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
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
leaf_PC2_rgr_all <- ggplot(leaf_growth_env, aes(x = PC2, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Leaf Trait PC2 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
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

# PC1 - All species together 
root_PC1_rgr_all <- ggplot(root_growth_env, aes(x = PC1, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC1 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
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
root_PC2_rgr_all <- ggplot(root_growth_env, aes(x = PC2, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Root Trait PC2 Value", y = expression("Mean RGR ("*yr^{-1}*")")) +
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
# (4) GROWTH RELATIONSHIPS TO ALL TRAIT PCA VALUES
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

