# -----------------------------------------------------------------------------#
# Neighborhood - tree growth relationships in WFDP
# Original Author: L. McKinley Nevins 
# December 19, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     sf v 1.0.19
#                     raster v 3.6.32
#                     cowplot v 1.2.0
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

#################################################################################
#                               Main workflow                                   #
#  Use the calculated neighborhood metrics, and explore the relationships of    #
#  the neighborhoods with the growth of the focal trees.                        #
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)


# Read in tree demographic data 
growth <- read.csv("stems_WFDP_20250206_trimmed.csv")

# Will exclude two of the trees that have growth data because they were on the edge of 
# the forest and we don't have complete data for all of their neighbors. These trees are 
# already removed from the summaries down below, so when they are merged with the growth file 
# they will be automatically removed. 


# Load in neighborhood summary files for combined, conspecific, heterospecific, same and 
# different mycorrhizal association for the 9m and 20m radius neighborhoods 

crowding_summary_09 <- read.csv("./Neighborhood_files/combined_crowd_summary_9m.csv", row.names = 1)

crowding_summary_20 <- read.csv( "./Neighborhood_files/combined_crowd_summary_20m.csv", row.names = 1)


con_summary_09 <- read.csv("./Neighborhood_files/conspecific_crowd_summary_9m.csv", row.names = 1)

con_summary_20 <- read.csv("./Neighborhood_files/conspecific_crowd_summary_20m.csv", row.names = 1)


het_summary_09 <- read.csv("./Neighborhood_files/heterospecific_crowd_summary_9m.csv", row.names = 1)

het_summary_20 <- read.csv("./Neighborhood_files/heterospecific_crowd_summary_20m.csv", row.names = 1)


same_myco_summary_09 <- read.csv("./Neighborhood_files/same_myco_crowd_summary_9m.csv", row.names = 1)

same_myco_summary_20 <- read.csv("./Neighborhood_files/same_myco_crowd_summary_20m.csv", row.names = 1)


diff_myco_summary_09 <- read.csv("./Neighborhood_files/diff_myco_crowd_summary_9m.csv", row.names = 1)

diff_myco_summary_20 <- read.csv("./Neighborhood_files/diff_myco_crowd_summary_20m.csv", row.names = 1)


# Total data for 58 trees for all datasets but the diff_myco, where there are a few fewer trees, 
# because a few of the focal trees did not have any neighbors that had a different mycorrhizal 
# association than their own. This is 9 trees less in the 9m neighborhoods, and 4 for the 20m. 

#################################################################################

##################################### -- 
# (2) RELATIONSHIPS TO DEMOGRAPHY 
##################################### -- 

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


# Add column for stem tag to match the crowding data 
diams$focal_stem_tag <- diams$Stem_Tag


#### Merge diams file to all datafiles, and merge the 9m and 20m datafiles into one for each 
# neighborhood type to reduce some of the bulk 

# Pair diams with each of the summary files according to the focal_stem_tag
crowd_growth_summary_09 <- merge(diams, crowding_summary_09, by = "focal_stem_tag")
crowd_growth_summary_20 <- merge(diams, crowding_summary_20, by = "focal_stem_tag")


# add neighborhood radius identifier and bind the dataframe rows 
crowd_growth_summary_09 <- crowd_growth_summary_09 %>%
  mutate(radius = "9m")

crowd_growth_summary_20 <- crowd_growth_summary_20 %>%
  mutate(radius = "20m")

# merge
crowd_growth_summary_all <- bind_rows(crowd_growth_summary_09, crowd_growth_summary_20)


####
con_growth_summary_09 <- merge(diams, con_summary_09, by = "focal_stem_tag")
con_growth_summary_20 <- merge(diams, con_summary_20, by = "focal_stem_tag")


# add neighborhood radius identifier and bind the dataframe rows 
con_growth_summary_09 <- con_growth_summary_09 %>%
  mutate(radius = "9m")

con_growth_summary_20 <- con_growth_summary_20 %>%
  mutate(radius = "20m")

# merge
con_growth_summary_all <- bind_rows(con_growth_summary_09, con_growth_summary_20)


####
het_growth_summary_09 <- merge(diams, het_summary_09, by = "focal_stem_tag")
het_growth_summary_20 <- merge(diams, het_summary_20, by = "focal_stem_tag")


# add neighborhood radius identifier and bind the dataframe rows 
het_growth_summary_09 <- het_growth_summary_09 %>%
  mutate(radius = "9m")

het_growth_summary_20 <- het_growth_summary_20 %>%
  mutate(radius = "20m")

# merge
het_growth_summary_all <- bind_rows(het_growth_summary_09, het_growth_summary_20)


####
same_myco_growth_summary_09 <- merge(diams, same_myco_summary_09, by = "focal_stem_tag")
same_myco_growth_summary_20 <- merge(diams, same_myco_summary_20, by = "focal_stem_tag")


# add neighborhood radius identifier and bind the dataframe rows 
same_myco_growth_summary_09 <- same_myco_growth_summary_09 %>%
  mutate(radius = "9m")

same_myco_growth_summary_20 <- same_myco_growth_summary_20 %>%
  mutate(radius = "20m")

# merge
same_myco_growth_summary_all <- bind_rows(same_myco_growth_summary_09, same_myco_growth_summary_20)


####
diff_myco_growth_summary_09 <- merge(diams, diff_myco_summary_09, by = "focal_stem_tag")
diff_myco_growth_summary_20 <- merge(diams, diff_myco_summary_20, by = "focal_stem_tag")


# add neighborhood radius identifier and bind the dataframe rows 
diff_myco_growth_summary_09 <- diff_myco_growth_summary_09 %>%
  mutate(radius = "9m")

diff_myco_growth_summary_20 <- diff_myco_growth_summary_20 %>%
  mutate(radius = "20m")

# merge
diff_myco_growth_summary_all <- bind_rows(diff_myco_growth_summary_09, diff_myco_growth_summary_20)


### All datasets are now formatted to have the neighbor data for the 9m and 20m radius neighborhoods for 
# each focal tree


######

# Differences in growth rate by focal species 
# Doing this once with the crowding data for the 9m neighborhoods, but this will be the same for all 
# neighborhoods because each focal tree will be the same 


# Load host colors for future plotting:

# set colors for hosts 
                # ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#9b5fe0", "#16a4d8", "#60dbe8", "#8bd346","#efdf48", "#f9a52F", "#d64e12")



RGR_spp <- ggplot(crowd_growth_summary_09, aes(x = focal_species, y = RGR, fill = focal_species)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values=all_hosts, 
                      name="Host Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(
    title = "",
    x = "",
    y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black")) +
  theme(legend.position = "none")

RGR_spp


# Test for significant differences between species 
aov_RGR <- aov(RGR ~ focal_species, data = crowd_growth_summary_09)
summary(aov_RGR)
# p = 0.00055

tuk_RGR <- TukeyHSD(aov_RGR)
tuk_RGR

# 
# diff          lwr          upr     p adj
# ABGR-ABAM -6.961907e-03 -0.022135185  0.008211371 0.7949489
# ALRU-ABAM  1.039501e-02 -0.006988166  0.027778182 0.5313206
# CONU-ABAM -4.782161e-03 -0.017922605  0.008358283 0.9199973
# TABR-ABAM -1.486710e-02 -0.028367619 -0.001366573 0.0221590
# THPL-ABAM  3.972318e-03 -0.009528205  0.017472841 0.9704468
# TSHE-ABAM -4.856806e-03 -0.017997250  0.008283638 0.9143571
# ALRU-ABGR  1.735691e-02 -0.001609682  0.036323512 0.0933962
# CONU-ABGR  2.179746e-03 -0.012993532  0.017353023 0.9993884
# TABR-ABGR -7.905189e-03 -0.023391351  0.007580973 0.7035343
# THPL-ABGR  1.093423e-02 -0.004551937  0.026420387 0.3311547
# TSHE-ABGR  2.105101e-03 -0.013068177  0.017278379 0.9994992
# CONU-ALRU -1.517717e-02 -0.032560343  0.002206004 0.1245227
# TABR-ALRU -2.526210e-02 -0.042919045 -0.007605163 0.0010593
# THPL-ALRU -6.422690e-03 -0.024079631  0.011234251 0.9201703
# TSHE-ALRU -1.525181e-02 -0.032634988  0.002131359 0.1210023
# TABR-CONU -1.008493e-02 -0.023585457  0.003415589 0.2674960
# THPL-CONU  8.754480e-03 -0.004746043  0.022255003 0.4333284
# TSHE-CONU -7.464482e-05 -0.013215089  0.013065799 1.0000000
# THPL-TABR  1.883941e-02  0.004988170  0.032690658 0.0021164
# TSHE-TABR  1.001029e-02 -0.003490233  0.023510812 0.2756100
# TSHE-THPL -8.829124e-03 -0.022329647  0.004671398 0.4229841

# Significant Differences: 

# TABR-ABAM p = 0.022
# TABR-ALRU p = 0.001
# THPL-TABR p = 0.002

# TABR has a negative RGR on average, because some of its bark sloughs off, probably. 

#################################### -- 
# VISUALIZE COMBINED NEIGHBORHOODS  
#################################### -- 

# Set radius as factor and set order for the facet panels for all graphs 
crowd_growth_summary_all <- crowd_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and # neighbors 
RGR_num_neigh <- ggplot(crowd_growth_summary_all, aes(x = num_neighbors, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Neighbors", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

# Having this plot generate the legend that can be used for other plots down below. Everyone else will 
# have the legend hidden so it doesn't duplicate when the plots are compiled. 

RGR_num_neigh


# The regression models need to be performed separately for the two neighborhood radii 
RGR_mod1_09 <- lm(RGR ~ num_neighbors, data = crowd_growth_summary_09)

summary(RGR_mod1_09)

# Significant relationship, trees appear to grow faster with more close neighbors 
# Multiple R-squared:  0.1802,	Adjusted R-squared:  0.1656, 
#F-statistic: 12.31 on 1 and 56 DF,  p-value: 0.0008963

RGR_mod1_20 <- lm(RGR ~ num_neighbors, data = crowd_growth_summary_20)

summary(RGR_mod1_20)

# Significant relationship, trees appear to grow faster with more close neighbors 
# Multiple R-squared:  0.1049,	Adjusted R-squared:  0.0889 
# F-statistic: 6.562 on 1 and 56 DF,  p-value: 0.01314



# Relationship between focal tree RGR and mean neighbor DBH 
RGR_size_neigh <- ggplot(crowd_growth_summary_all, aes(x = mean_neighbor_DBH, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Neighbor DBH", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_size_neigh


RGR_mod2_09 <- lm(RGR ~ mean_neighbor_DBH, data = crowd_growth_summary_09)

summary(RGR_mod2_09)

# Significant relationship, trees appear to grow faster when their neighbors are smaller on average 
# Multiple R-squared:  0.1878,	Adjusted R-squared:  0.1733 
# F-statistic: 12.95 on 1 and 56 DF,  p-value: 0.0006789


RGR_mod2_20 <- lm(RGR ~ mean_neighbor_DBH, data = crowd_growth_summary_20)

summary(RGR_mod2_20)

# Significant relationship, trees appear to grow faster when their neighbors are smaller on average 
# Multiple R-squared:  0.1331,	Adjusted R-squared:  0.1176 
# F-statistic: 8.598 on 1 and 56 DF,  p-value: 0.004868


# Relationship between focal tree RGR and crowding index
RGR_crowd_neigh <- ggplot(crowd_growth_summary_all, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_crowd_neigh


RGR_mod3_09 <- lm(RGR ~ crowding_index, data = crowd_growth_summary_09)

summary(RGR_mod3_09)

# No Significant relationship, p = 0.067


RGR_mod3_20 <- lm(RGR ~ crowding_index, data = crowd_growth_summary_20)

summary(RGR_mod3_20)

# No Significant relationship, p = 0.065


################################################### -- 
# VISUALIZE CON- AND HETEROSPECIFIC NEIGHBORHOODS  
################################################### -- 

### For conspecific neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
con_growth_summary_all <- con_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))

# Relationship between focal tree RGR and # conspecific neighbors 
RGR_num_con_neigh <- ggplot(con_growth_summary_all, aes(x = num_neighbors, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Conspecific Neighbors", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_num_con_neigh


RGR_mod_con1_09 <- lm(RGR ~ num_neighbors, data = con_growth_summary_09)

summary(RGR_mod_con1_09)

# No significant relationship,  
# Multiple R-squared:  0.05006,	Adjusted R-squared:  0.03309 
# F-statistic: 2.951 on 1 and 56 DF,  p-value: 0.09136

RGR_mod_con1_20 <- lm(RGR ~ num_neighbors, data = con_growth_summary_20)

summary(RGR_mod_con1_20)

# No significant relationship,  
# Multiple R-squared:  0.01118,	Adjusted R-squared:  -0.006482 
# F-statistic: 0.6329 on 1 and 56 DF,  p-value: 0.4297


# Relationship between focal tree RGR and mean conspecific neighbor DBH 
RGR_size_con_neigh <- ggplot(con_growth_summary_all, aes(x = mean_neighbor_DBH, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Conspecific Neighbor DBH", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_size_con_neigh


RGR_mod_con2_09 <- lm(RGR ~ mean_neighbor_DBH, data = con_growth_summary_09)

summary(RGR_mod_con2_09)

# Not significant relationship,  
# Multiple R-squared:  0.04093,	Adjusted R-squared:  0.0238 
# F-statistic:  2.39 on 1 and 56 DF,  p-value: 0.1278

RGR_mod_con_20 <- lm(RGR ~ mean_neighbor_DBH, data = con_growth_summary_20)

summary(RGR_mod_con_20)

# Not significant relationship,  
# Multiple R-squared:  0.0006053,	Adjusted R-squared:  -0.01724 
# F-statistic: 0.03392 on 1 and 56 DF,  p-value: 0.8546


# Relationship between focal tree RGR and crowding index
RGR_crowd_con_neigh <- ggplot(con_growth_summary_all, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Conspecific Crowding Index", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_crowd_con_neigh


RGR_mod_con3_09 <- lm(RGR ~ crowding_index, data = con_growth_summary_09)

summary(RGR_mod_con3_09)

# Not significant relationship,  
# Multiple R-squared:  0.05282,	Adjusted R-squared:  0.03591 
# F-statistic: 3.123 on 1 and 56 DF,  p-value: 0.08264

RGR_mod_con3_20 <- lm(RGR ~ crowding_index, data = con_growth_summary_20)

summary(RGR_mod_con3_20)

# Not significant relationship,  
# Multiple R-squared:  0.05325,	Adjusted R-squared:  0.03634 
# F-statistic:  3.15 on 1 and 56 DF,  p-value: 0.08138


### For heterospecific neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
het_growth_summary_all <- het_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and # heterospecific neighbors 
RGR_num_het_neigh <- ggplot(het_growth_summary_all, aes(x = num_neighbors, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Heterospecific Neighbors", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_num_het_neigh


RGR_mod_het1_09 <- lm(RGR ~ num_neighbors, data = het_growth_summary_09)

summary(RGR_mod_het1_09)

# Significant relationship, focal trees grow faster in neighborhoods with more heterospecifics  
# Multiple R-squared:  0.1037,	Adjusted R-squared:  0.08768 
# F-statistic: 6.478 on 1 and 56 DF,  p-value: 0.0137


RGR_mod_het1_20 <- lm(RGR ~ num_neighbors, data = het_growth_summary_20)

summary(RGR_mod_het1_20)

# Significant relationship, focal trees grow faster in neighborhoods with more heterospecifics  
# Multiple R-squared:  0.06719,	Adjusted R-squared:  0.05054 
# F-statistic: 4.034 on 1 and 56 DF,  p-value: 0.04943


# Relationship between focal tree RGR and mean heterospecific neighbor DBH 
RGR_size_het_neigh <- ggplot(het_growth_summary_all, aes(x = mean_neighbor_DBH, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Heterospecific Neighbor DBH", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_size_het_neigh


RGR_mod_het2_09 <- lm(RGR ~ mean_neighbor_DBH, data = het_growth_summary_09)

summary(RGR_mod_het2_09)

# Significant relationship, focal trees grow faster around smaller heterospecific neighbors
# Multiple R-squared:  0.1333,	Adjusted R-squared:  0.1178 
# F-statistic: 8.614 on 1 and 56 DF,  p-value: 0.00483


RGR_mod_het2_20 <- lm(RGR ~ mean_neighbor_DBH, data = het_growth_summary_20)

summary(RGR_mod_het2_20)

# Significant relationship, focal trees grow faster around smaller heterospecific neighbors
# Multiple R-squared:  0.07408,	Adjusted R-squared:  0.05754 
# F-statistic:  4.48 on 1 and 56 DF,  p-value: 0.03874


# Relationship between focal tree RGR and crowding index
RGR_crowd_het_neigh <- ggplot(het_growth_summary_all, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Heterospecific Crowding Index", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_crowd_het_neigh


RGR_mod_het3_09 <- lm(RGR ~ crowding_index, data = het_growth_summary_09)

summary(RGR_mod_het3_09)

# Not significant relationship,  
# Multiple R-squared:  0.006019,	Adjusted R-squared:  -0.01173 
# F-statistic: 0.3391 on 1 and 56 DF,  p-value: 0.5627


RGR_mod_het3_20 <- lm(RGR ~ crowding_index, data = het_growth_summary_20)

summary(RGR_mod_het3_20)

# Not significant relationship,  
# Multiple R-squared:  0.006777,	Adjusted R-squared:  -0.01096 
# F-statistic: 0.3821 on 1 and 56 DF,  p-value: 0.539


# Has that one high crowding outlier again, so can do this with the outlier removed to compare 
het_growth_summary_no_outlier <- het_growth_summary_all %>% filter(focal_cell != "H30") %>% droplevels()


# Visualize again without outlier 
RGR_crowd_het_neigh_no_outlier <- ggplot(het_growth_summary_no_outlier, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Heterospecific Crowding Index", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_crowd_het_neigh_no_outlier


# Remove from the separate radius datasets to test the relationship again 
het_growth_summary_no_outlier_09 <- het_growth_summary_09 %>% filter(focal_cell != "H30") %>% droplevels()



RGR_mod_het3_no_outlier_09 <- lm(RGR ~ crowding_index, data = het_growth_summary_no_outlier_09)

summary(RGR_mod_het3_no_outlier_09)

# Still not a significant relationship 
# Multiple R-squared:  0.02447,	Adjusted R-squared:  0.006735 
# F-statistic:  1.38 on 1 and 55 DF,  p-value: 0.2452


het_growth_summary_no_outlier_20 <- het_growth_summary_20 %>% filter(focal_cell != "H30") %>% droplevels()


RGR_mod_het3_no_outlier_20 <- lm(RGR ~ crowding_index, data = het_growth_summary_no_outlier_20)

summary(RGR_mod_het3_no_outlier_20)

# Still not a significant relationship 
# Multiple R-squared:  0.03302,	Adjusted R-squared:  0.01543 
# F-statistic: 1.878 on 1 and 55 DF,  p-value: 0.1761



# Take a peak at just the data for THPL to see how this one outlier compares to the rest
# Doing this just in the 9m radius neighborhoods because the growth rate for the focal tree itself doesn't differ 
# between the neighborhood radii 

het_growth_summary_09$focal_species <- as.factor(het_growth_summary_09$focal_species)


THPL_het_09 <- het_growth_summary_09 %>% filter(focal_species == "THPL") %>% droplevels()


# Visualize crowding and growth for just THPL

RGR_crowd_het_THPL_09 <- ggplot(THPL_het_09, aes(x = crowding_index, y = RGR)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  labs(x = "Crowding Index - 9m Neighborhoods only THPL", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black"))


RGR_crowd_het_THPL_09

# This shows the one THPL still drives the trend, but it is definitely the lowest RGR out of all of them 


####################################################################### -- 
# VISUALIZE SAME AND DIFFERENT MYCORRHIZAL ASSOCIATIONS NEIGHBORHOODS  
####################################################################### -- 

### For same mycorrhizal neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
same_myco_growth_summary_all <- same_myco_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and # same mycorrhizal neighbors 
RGR_num_same_myco_neigh <- ggplot(same_myco_growth_summary_all, aes(x = num_neighbors, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Same Mycorrhizal Neighbors", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_num_same_myco_neigh


RGR_mod_same_myco1_09 <- lm(RGR ~ num_neighbors, data = same_myco_growth_summary_09)

summary(RGR_mod_same_myco1_09)

# Not significant relationship,  
# Multiple R-squared:  0.05997,	Adjusted R-squared:  0.04319 
# F-statistic: 3.573 on 1 and 56 DF,  p-value: 0.06392


RGR_mod_same_myco1_20 <- lm(RGR ~ num_neighbors, data = same_myco_growth_summary_20)

summary(RGR_mod_same_myco1_20)

# Not significant relationship,  
# Multiple R-squared:  0.02168,	Adjusted R-squared:  0.004211 
# F-statistic: 1.241 on 1 and 56 DF,  p-value: 0.27


# Relationship between focal tree RGR and mean same mycorrhizal neighbor DBH 
RGR_size_same_myco_neigh <- ggplot(same_myco_growth_summary_all, aes(x = mean_neighbor_DBH, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Same Mycorrhizal Neighbor DBH", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_size_same_myco_neigh


RGR_mod_same_myco2_09 <- lm(RGR ~ mean_neighbor_DBH, data = same_myco_growth_summary_09)

summary(RGR_mod_same_myco2_09)

# Significant relationship,  
# Multiple R-squared:  0.1785,	Adjusted R-squared:  0.1639 
# F-statistic: 12.17 on 1 and 56 DF,  p-value: 0.0009525


RGR_mod_same_myco2_20 <- lm(RGR ~ mean_neighbor_DBH, data = same_myco_growth_summary_20)

summary(RGR_mod_same_myco2_20)

# Significant relationship,  
# Multiple R-squared:  0.07809,	Adjusted R-squared:  0.06163 
# F-statistic: 4.743 on 1 and 56 DF,  p-value: 0.03364


# Relationship between focal tree RGR and crowding index
RGR_crowd_same_myco_neigh <- ggplot(same_myco_growth_summary_all, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index - Same Mycorrhizal Neighborhoods", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_crowd_same_myco_neigh


RGR_mod_same_myco3_09 <- lm(RGR ~ crowding_index, data = same_myco_growth_summary_09)

summary(RGR_mod_same_myco3_09)

# Not significant relationship,  
# Multiple R-squared:  0.05733,	Adjusted R-squared:  0.0405 
# F-statistic: 3.406 on 1 and 56 DF,  p-value: 0.07025


RGR_mod_same_myco3_20 <- lm(RGR ~ crowding_index, data = same_myco_growth_summary_20)

summary(RGR_mod_same_myco3_20)

# Not significant relationship,  
# Multiple R-squared:  0.0592,	Adjusted R-squared:  0.0424 
# F-statistic: 3.523 on 1 and 56 DF,  p-value: 0.06572


### For different mycorrhizal neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
diff_myco_growth_summary_all <- diff_myco_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and different mycorrhizal neighbors 
RGR_num_diff_myco_neigh <- ggplot(diff_myco_growth_summary_all, aes(x = num_neighbors, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Different Mycorrhizal Neighbors", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_num_diff_myco_neigh


RGR_mod_diff_myco1_09 <- lm(RGR ~ num_neighbors, data = diff_myco_growth_summary_09)

summary(RGR_mod_diff_myco1_09)

# Significant relationship,  
# Multiple R-squared:  0.1309,	Adjusted R-squared:  0.1124 
# F-statistic: 7.079 on 1 and 47 DF,  p-value: 0.01064


RGR_mod_diff_myco1_20 <- lm(RGR ~ num_neighbors, data = diff_myco_growth_summary_20)

summary(RGR_mod_diff_myco1_20)

# Not significant relationship,  
# Multiple R-squared:  0.05439,	Adjusted R-squared:  0.0362 
# F-statistic: 2.991 on 1 and 52 DF,  p-value: 0.08967


# Relationship between focal tree RGR and mean different mycorrhizal neighbor DBH 
RGR_size_diff_myco_neigh <- ggplot(diff_myco_growth_summary_all, aes(x = mean_neighbor_DBH, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Different Mycorrhizal Neighbor DBH", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_size_diff_myco_neigh


RGR_mod_diff_myco2_09 <- lm(RGR ~ mean_neighbor_DBH, data = diff_myco_growth_summary_09)

summary(RGR_mod_diff_myco2_09)

# No significant relationship,  
# Multiple R-squared:  2.69e-05,	Adjusted R-squared:  -0.02125 
# F-statistic: 0.001264 on 1 and 47 DF,  p-value: 0.9718


RGR_mod_diff_myco2_20 <- lm(RGR ~ mean_neighbor_DBH, data = diff_myco_growth_summary_20)

summary(RGR_mod_diff_myco2_20)

# No significant relationship,  
# Multiple R-squared:  0.01506,	Adjusted R-squared:  -0.003885 
# F-statistic: 0.7949 on 1 and 52 DF,  p-value: 0.3767


# Relationship between focal tree RGR and crowding index
RGR_crowd_diff_myco_neigh <- ggplot(diff_myco_growth_summary_all, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index - Different Mycorrhizal Neighborhoods", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_crowd_diff_myco_neigh


RGR_mod_diff_myco3_09 <- lm(RGR ~ crowding_index, data = diff_myco_growth_summary_09)

summary(RGR_mod_diff_myco3_09)

# Not significant relationship,  
# Multiple R-squared:  0.003669,	Adjusted R-squared:  -0.01753 
#F-statistic: 0.1731 on 1 and 47 DF,  p-value: 0.6793


RGR_mod_diff_myco3_20 <- lm(RGR ~ crowding_index, data = diff_myco_growth_summary_20)

summary(RGR_mod_diff_myco3_20)

# Not significant relationship,  
# Multiple R-squared:  0.001448,	Adjusted R-squared:  -0.01775 
# F-statistic: 0.07541 on 1 and 52 DF,  p-value: 0.7847


# Being impacted by the outlier THPL again 
diff_myco_growth_summary_no_outlier <- diff_myco_growth_summary_all %>% filter(focal_cell != "H30") %>% droplevels()


# Visualize again without outlier 
RGR_crowd_diff_myco_no_outlier <- ggplot(diff_myco_growth_summary_no_outlier, aes(x = crowding_index, y = RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index - Different Mycorrhizal Neighborhoods", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")


RGR_crowd_diff_myco_no_outlier


# Remove from the separate radius datasets to test the relationship again 
diff_myco_growth_summary_no_outlier_09 <- diff_myco_growth_summary_09 %>% filter(focal_cell != "H30") %>% droplevels()


RGR_mod_diff_myco3_no_outlier_09 <- lm(RGR ~ crowding_index, data = diff_myco_growth_summary_no_outlier_09)

summary(RGR_mod_diff_myco3_no_outlier_09)

# Still not a significant relationship 
# Multiple R-squared:  0.004383,	Adjusted R-squared:  -0.01726 
# F-statistic: 0.2025 on 1 and 46 DF,  p-value: 0.6548


diff_myco_growth_summary_no_outlier_20 <- diff_myco_growth_summary_20 %>% filter(focal_cell != "H30") %>% droplevels()


RGR_mod_diff_myco3_no_outlier_20 <- lm(RGR ~ crowding_index, data = diff_myco_growth_summary_no_outlier_20)

summary(RGR_mod_diff_myco3_no_outlier_20)

# Still not a significant relationship 
# Multiple R-squared:  0.01116,	Adjusted R-squared:  -0.008232 
# F-statistic: 0.5754 on 1 and 51 DF,  p-value: 0.4516


# Take a peak at just the data for THPL to see how this one outlier compares to the rest
# Doing this just in the 9m radius neighborhoods because the growth rate for the focal tree itself doesn't differ 
# between the neighborhood radii 

diff_myco_growth_summary_09$focal_species <- as.factor(diff_myco_growth_summary_09$focal_species)


THPL_diff_myco_09 <- diff_myco_growth_summary_09 %>% filter(focal_species == "THPL") %>% droplevels()


# Visualize crowding and growth for just THPL
RGR_crowd_diff_myco_THPL_09 <- ggplot(THPL_het_09, aes(x = crowding_index, y = RGR)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  labs(x = "Crowding Index - Different Mycorrhizal Neighborhoods only THPL", y = expression("Relative Growth Rate ("*yr^{-1}*")"))
  

RGR_crowd_diff_myco_THPL_09

# This shows the one THPL still drives the trend, but it is definitely the lowest RGR out of all of them 

#################################################################################

############################################## -- 
# (3) GATHER SIGNIFICANT RESULTS OF INTEREST  
############################################## -- 

# Gather up the most interesting relationships and organize the plots to compare between the 9 and 20m neighborhoods, 
# and the different types of neighborhoods - con/heterospecific and same/different mycorrhizal associations 
# Not using the combined neighborhoods because our interest is in parsing out the effects of the different 
# types of neighbors on the focal tree growth. 


# Lists of plots: 
# RGR_num_con_neigh, RGR_size_con_neigh, RGR_crowd_con_neigh

# RGR_num_het_neigh, RGR_size_het_neigh, RGR_crowd_het_neigh

# RGR_num_same_myco_neigh, RGR_size_same_myco_neigh, RGR_crowd_same_myco_neigh, 

# RGR_num_diff_myco_neigh, RGR_size_diff_myco_neigh, RGR_crowd_diff_myco_neigh


# When plotting the grid, need to shuffle the order of the plots a little so the same plots are next to eachother 


# This is using the neighborhoods with the high THPL outlier still in them, because removing it did not 
# change the significance of the results. 


# Conspecific vs heterospecific neighbor plots 
con_het_plots <- plot_grid(RGR_num_con_neigh, RGR_num_het_neigh, RGR_size_con_neigh, RGR_size_het_neigh,
                           RGR_crowd_con_neigh, RGR_crowd_het_neigh,
          ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', 
                                         '(e)', '(f)'))

con_het_plots


# Same vs different mycorrhizal association plots 
same_diff_myco_plots <- plot_grid(RGR_num_same_myco_neigh, RGR_num_diff_myco_neigh, 
                                  RGR_size_same_myco_neigh, RGR_size_diff_myco_neigh,
                                  RGR_crowd_same_myco_neigh, RGR_crowd_diff_myco_neigh,
                                  ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', 
                                                                 '(e)', '(f)'))

same_diff_myco_plots

## -- END -- ## 

