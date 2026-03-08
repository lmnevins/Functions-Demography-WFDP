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


# Add column to match the neighborhood data 
diams$focal_stem_tag <- diams$Stem_Tag


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
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")



RGR_spp <- ggplot(crowd_growth_summary_09, aes(x = focal_species, y = mean_RGR, fill = focal_species)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values=all_hosts, 
                      name="Host Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(
    title = "",
    x = "",
    y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black")) +
  theme(legend.position = "none")

RGR_spp


# Test for significant differences between species 
aov_RGR <- aov(mean_RGR ~ focal_species, data = crowd_growth_summary_09)
summary(aov_RGR)
# p = 0.00048

tuk_RGR <- TukeyHSD(aov_RGR)
tuk_RGR

# 
# $focal_species
# diff          lwr          upr     p adj
# ABGR-ABAM -6.952969e-03 -0.022096994  0.008191056 0.7944637
# ALRU-ABAM  1.045281e-02 -0.006896851  0.027802470 0.5224181
# CONU-ABAM -4.838083e-03 -0.017953193  0.008277028 0.9150790
# TABR-ABAM -1.493780e-02 -0.028412297 -0.001463307 0.0208287
# THPL-ABAM  4.057208e-03 -0.009417287  0.017531703 0.9669103
# TSHE-ABAM -4.903199e-03 -0.018018310  0.008211911 0.9099780
# ALRU-ABGR  1.740578e-02 -0.001524253  0.036335810 0.0905799
# CONU-ABGR  2.114886e-03 -0.013029139  0.017258911 0.9994799
# TABR-ABGR -7.984833e-03 -0.023441139  0.007471473 0.6919172
# THPL-ABGR  1.101018e-02 -0.004446129  0.026466483 0.3208679
# TSHE-ABGR  2.049770e-03 -0.013094256  0.017193795 0.9995655
# CONU-ALRU -1.529089e-02 -0.032640553  0.002058768 0.1178372
# TABR-ALRU -2.539061e-02 -0.043013512 -0.007767712 0.0009584
# THPL-ALRU -6.395601e-03 -0.024018502  0.011227299 0.9209775
# TSHE-ALRU -1.535601e-02 -0.032705670  0.001993651 0.1148899
# TABR-CONU -1.009972e-02 -0.023574214  0.003374776 0.2638202
# THPL-CONU  8.895291e-03 -0.004579204  0.022369786 0.4115527
# TSHE-CONU -6.511657e-05 -0.013180227  0.013049994 1.0000000
# THPL-TABR  1.899501e-02  0.005170470  0.032819551 0.0018498
# TSHE-TABR  1.003460e-02 -0.003439892  0.023509098 0.2708414
# TSHE-THPL -8.960408e-03 -0.022434903  0.004514088 0.4026969

# Significant Differences: 

# TABR-ABAM p = 0.020
# TABR-ALRU p = 0.001
# THPL-TABR p = 0.002

# TABR has a negative RGR on average, because some of its bark sloughs off, probably. 

#################################### -- 
# VISUALIZE COMBINED NEIGHBORHOODS  
#################################### -- 

## All of these have been edited to only plot the 9m radius neighborhoods



# Set radius as factor and set order for the facet panels for all graphs 
crowd_growth_summary_all <- crowd_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and # neighbors 
RGR_num_neigh <- ggplot(crowd_growth_summary_all, aes(x = num_neighbors, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Neighbors", y = expression("Mean RGR ("*yr^{-1}*")")) +
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
RGR_mod1_09 <- lm(mean_RGR ~ num_neighbors, data = crowd_growth_summary_09)

summary(RGR_mod1_09)

# Significant relationship, trees appear to grow faster with more close neighbors 
# Multiple R-squared:  0.1802,	Adjusted R-squared:  0.1656, 
#F-statistic: 12.31 on 1 and 56 DF,  p-value: 0.0008963

RGR_mod1_20 <- lm(mean_RGR ~ num_neighbors, data = crowd_growth_summary_20)

summary(RGR_mod1_20)

# Significant relationship, trees appear to grow faster with more close neighbors 
# Multiple R-squared:  0.1053,	Adjusted R-squared:  0.08929 
# F-statistic: 6.589 on 1 and 56 DF,  p-value: 0.01296



# Relationship between focal tree RGR and mean neighbor DBH 
RGR_size_neigh <- ggplot(crowd_growth_summary_all, aes(x = mean_neighbor_DBH, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Neighbor DBH", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod2_09 <- lm(mean_RGR ~ mean_neighbor_DBH, data = crowd_growth_summary_09)

summary(RGR_mod2_09)

# Significant relationship, trees appear to grow faster when their neighbors are smaller on average 
# Multiple R-squared:  0.1887,	Adjusted R-squared:  0.1742 
# F-statistic: 13.02 on 1 and 56 DF,  p-value: 0.0006567


RGR_mod2_20 <- lm(mean_RGR ~ mean_neighbor_DBH, data = crowd_growth_summary_20)

summary(RGR_mod2_20)

# Significant relationship, trees appear to grow faster when their neighbors are smaller on average 
# Multiple R-squared:  0.1338,	Adjusted R-squared:  0.1183 
# F-statistic: 8.651 on 1 and 56 DF,  p-value: 0.004747


# Relationship between focal tree RGR and crowding index
RGR_crowd_neigh <- ggplot(crowd_growth_summary_all, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod3_09 <- lm(mean_RGR ~ crowding_index, data = crowd_growth_summary_09)

summary(RGR_mod3_09)

# No Significant relationship, p = 0.066


RGR_mod3_20 <- lm(mean_RGR ~ crowding_index, data = crowd_growth_summary_20)

summary(RGR_mod3_20)

# No Significant relationship, p = 0.063


################################################### -- 
# VISUALIZE CON- AND HETEROSPECIFIC NEIGHBORHOODS  
################################################### -- 

### For conspecific neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
con_growth_summary_all <- con_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))

# Relationship between focal tree RGR and # conspecific neighbors 
RGR_num_con_neigh <- ggplot(con_growth_summary_all, aes(x = num_neighbors, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Conspecific Neighbors", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_con1_09 <- lm(mean_RGR ~ num_neighbors, data = con_growth_summary_09)

summary(RGR_mod_con1_09)

# No significant relationship,  
# Multiple R-squared:  0.04841,	Adjusted R-squared:  0.03142 
# F-statistic: 2.849 on 1 and 56 DF,  p-value: 0.09699

RGR_mod_con1_20 <- lm(mean_RGR ~ num_neighbors, data = con_growth_summary_20)

summary(RGR_mod_con1_20)

# No significant relationship,  
# Multiple R-squared:  0.01043,	Adjusted R-squared:  -0.007236 
# F-statistic: 0.5905 on 1 and 56 DF,  p-value: 0.4455


# Relationship between focal tree RGR and mean conspecific neighbor DBH 
RGR_size_con_neigh <- ggplot(con_growth_summary_all, aes(x = mean_neighbor_DBH, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Conspecific Neighbor DBH", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_con2_09 <- lm(mean_RGR ~ mean_neighbor_DBH, data = con_growth_summary_09)

summary(RGR_mod_con2_09)

# Not significant relationship,  
# Multiple R-squared:  0.04082,	Adjusted R-squared:  0.02369 
# F-statistic: 2.383 on 1 and 56 DF,  p-value: 0.1283

RGR_mod_con_20 <- lm(mean_RGR ~ mean_neighbor_DBH, data = con_growth_summary_20)

summary(RGR_mod_con_20)

# Not significant relationship,  
# Multiple R-squared:  0.0005199,	Adjusted R-squared:  -0.01733 
# F-statistic: 0.02913 on 1 and 56 DF,  p-value: 0.8651


# Relationship between focal tree RGR and crowding index
RGR_crowd_con_neigh <- ggplot(con_growth_summary_all, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Conspecific Crowding Index", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_con3_09 <- lm(mean_RGR ~ crowding_index, data = con_growth_summary_09)

summary(RGR_mod_con3_09)

# Not significant relationship,  
# Multiple R-squared:  0.05323,	Adjusted R-squared:  0.03632 
# F-statistic: 3.148 on 1 and 56 DF,  p-value: 0.08145

RGR_mod_con3_20 <- lm(mean_RGR ~ crowding_index, data = con_growth_summary_20)

summary(RGR_mod_con3_20)

# Not significant relationship,  
# Multiple R-squared:  0.05365,	Adjusted R-squared:  0.03675 
# F-statistic: 3.174 on 1 and 56 DF,  p-value: 0.08022


### For heterospecific neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
het_growth_summary_all <- het_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and # heterospecific neighbors 
RGR_num_het_neigh <- ggplot(het_growth_summary_all, aes(x = num_neighbors, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Heterospecific Neighbors", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_het1_09 <- lm(mean_RGR ~ num_neighbors, data = het_growth_summary_09)

summary(RGR_mod_het1_09)

# Significant relationship, focal trees grow faster in neighborhoods with more heterospecifics  
# Multiple R-squared:  0.1049,	Adjusted R-squared:  0.08895 
# F-statistic: 6.565 on 1 and 56 DF,  p-value: 0.01312


RGR_mod_het1_20 <- lm(mean_RGR ~ num_neighbors, data = het_growth_summary_20)

summary(RGR_mod_het1_20)

# Significant relationship, focal trees grow faster in neighborhoods with more heterospecifics  
# Multiple R-squared:  0.06821,	Adjusted R-squared:  0.05157 
# F-statistic: 4.099 on 1 and 56 DF,  p-value: 0.04768


# Relationship between focal tree RGR and mean heterospecific neighbor DBH 
RGR_size_het_neigh <- ggplot(het_growth_summary_all, aes(x = mean_neighbor_DBH, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Heterospecific Neighbor DBH", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_het2_09 <- lm(mean_RGR ~ mean_neighbor_DBH, data = het_growth_summary_09)

summary(RGR_mod_het2_09)

# Significant relationship, focal trees grow faster around smaller heterospecific neighbors
# Multiple R-squared:  0.1338,	Adjusted R-squared:  0.1183 
# F-statistic:  8.65 on 1 and 56 DF,  p-value: 0.004747


RGR_mod_het2_20 <- lm(mean_RGR ~ mean_neighbor_DBH, data = het_growth_summary_20)

summary(RGR_mod_het2_20)

# Significant relationship, focal trees grow faster around smaller heterospecific neighbors
# Multiple R-squared:  0.07458,	Adjusted R-squared:  0.05805 
# F-statistic: 4.513 on 1 and 56 DF,  p-value: 0.03807


# Relationship between focal tree RGR and crowding index
RGR_crowd_het_neigh <- ggplot(het_growth_summary_all, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Heterospecific Crowding Index", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_het3_09 <- lm(mean_RGR ~ crowding_index, data = het_growth_summary_09)

summary(RGR_mod_het3_09)

# Not significant relationship,  
# Multiple R-squared:  0.006318,	Adjusted R-squared:  -0.01143 
# F-statistic: 0.3561 on 1 and 56 DF,  p-value: 0.5531


RGR_mod_het3_20 <- lm(mean_RGR ~ crowding_index, data = het_growth_summary_20)

summary(RGR_mod_het3_20)

# Not significant relationship,  
# Multiple R-squared:  0.007074,	Adjusted R-squared:  -0.01066 
# F-statistic: 0.3989 on 1 and 56 DF,  p-value: 0.5302


# Has that one high crowding outlier again, so can do this with the outlier removed to compare 
het_growth_summary_no_outlier <- het_growth_summary_all %>% filter(focal_cell != "H30") %>% droplevels()


# Visualize again without outlier 
RGR_crowd_het_neigh_no_outlier <- ggplot(het_growth_summary_no_outlier, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Heterospecific Crowding Index", y = expression("Mean RGR ("*yr^{-1}*")")) +
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



RGR_mod_het3_no_outlier_09 <- lm(mean_RGR ~ crowding_index, data = het_growth_summary_no_outlier_09)

summary(RGR_mod_het3_no_outlier_09)

# Still not a significant relationship 
# Multiple R-squared:  0.02443,	Adjusted R-squared:  0.006693 
# F-statistic: 1.377 on 1 and 55 DF,  p-value: 0.2456


het_growth_summary_no_outlier_20 <- het_growth_summary_20 %>% filter(focal_cell != "H30") %>% droplevels()


RGR_mod_het3_no_outlier_20 <- lm(mean_RGR ~ crowding_index, data = het_growth_summary_no_outlier_20)

summary(RGR_mod_het3_no_outlier_20)

# Still not a significant relationship 
# Multiple R-squared:  0.03272,	Adjusted R-squared:  0.01514 
# F-statistic: 1.861 on 1 and 55 DF,  p-value: 0.1781



# Take a peak at just the data for THPL to see how this one outlier compares to the rest
# Doing this just in the 9m radius neighborhoods because the growth rate for the focal tree itself doesn't differ 
# between the neighborhood radii 

het_growth_summary_09$focal_species <- as.factor(het_growth_summary_09$focal_species)


THPL_het_09 <- het_growth_summary_09 %>% filter(focal_species == "THPL") %>% droplevels()


# Visualize crowding and growth for just THPL

RGR_crowd_het_THPL_09 <- ggplot(THPL_het_09, aes(x = crowding_index, y = mean_RGR)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  labs(x = "Crowding Index - 9m Neighborhoods only THPL", y = expression("Mean RGR ("*yr^{-1}*")")) +
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
RGR_num_same_myco_neigh <- ggplot(same_myco_growth_summary_all, aes(x = num_neighbors, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Same Mycorrhizal Neighbors", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_same_myco1_09 <- lm(mean_RGR ~ num_neighbors, data = same_myco_growth_summary_09)

summary(RGR_mod_same_myco1_09)

# Not significant relationship,  
# Multiple R-squared:  0.05924,	Adjusted R-squared:  0.04245 
# F-statistic: 3.527 on 1 and 56 DF,  p-value: 0.0656


RGR_mod_same_myco1_20 <- lm(mean_RGR ~ num_neighbors, data = same_myco_growth_summary_20)

summary(RGR_mod_same_myco1_20)

# Not significant relationship,  
# Multiple R-squared:  0.02119,	Adjusted R-squared:  0.003713 
# F-statistic: 1.212 on 1 and 56 DF,  p-value: 0.2756


# Relationship between focal tree RGR and mean same mycorrhizal neighbor DBH 
RGR_size_same_myco_neigh <- ggplot(same_myco_growth_summary_all, aes(x = mean_neighbor_DBH, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Same Mycorrhizal Neighbor DBH", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_same_myco2_09 <- lm(mean_RGR ~ mean_neighbor_DBH, data = same_myco_growth_summary_09)

summary(RGR_mod_same_myco2_09)

# Significant relationship,  
# Multiple R-squared:  0.1792,	Adjusted R-squared:  0.1645 
# F-statistic: 12.22 on 1 and 56 DF,  p-value: 0.0009313


RGR_mod_same_myco2_20 <- lm(mean_RGR ~ mean_neighbor_DBH, data = same_myco_growth_summary_20)

summary(RGR_mod_same_myco2_20)

# Significant relationship,  
# Multiple R-squared:  0.07819,	Adjusted R-squared:  0.06173 
# F-statistic:  4.75 on 1 and 56 DF,  p-value: 0.03352


# Relationship between focal tree RGR and crowding index
RGR_crowd_same_myco_neigh <- ggplot(same_myco_growth_summary_all, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index - Same Mycorrhizal Neighborhoods", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_same_myco3_09 <- lm(mean_RGR ~ crowding_index, data = same_myco_growth_summary_09)

summary(RGR_mod_same_myco3_09)

# Not significant relationship,  
# Multiple R-squared:  0.05779,	Adjusted R-squared:  0.04096 
# F-statistic: 3.434 on 1 and 56 DF,  p-value: 0.06912


RGR_mod_same_myco3_20 <- lm(mean_RGR ~ crowding_index, data = same_myco_growth_summary_20)

summary(RGR_mod_same_myco3_20)

# Not significant relationship,  
# Multiple R-squared:  0.05965,	Adjusted R-squared:  0.04286 
# F-statistic: 3.552 on 1 and 56 DF,  p-value: 0.06466


### For different mycorrhizal neighborhoods: 

# Set radius as factor and set order for the facet panels for all graphs 
diff_myco_growth_summary_all <- diff_myco_growth_summary_all %>%
  mutate(radius = factor(radius, levels = c("9m", "20m")))


# Relationship between focal tree RGR and different mycorrhizal neighbors 
RGR_num_diff_myco_neigh <- ggplot(diff_myco_growth_summary_all, aes(x = num_neighbors, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Number of Different Mycorrhizal Neighbors", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_diff_myco1_09 <- lm(mean_RGR ~ num_neighbors, data = diff_myco_growth_summary_09)

summary(RGR_mod_diff_myco1_09)

# Significant relationship,  
# Multiple R-squared:  0.1324,	Adjusted R-squared:  0.1139 
# F-statistic: 7.173 on 1 and 47 DF,  p-value: 0.01017


RGR_mod_diff_myco1_20 <- lm(mean_RGR ~ num_neighbors, data = diff_myco_growth_summary_20)

summary(RGR_mod_diff_myco1_20)

# Not significant relationship,  
# Multiple R-squared:  0.0552,	Adjusted R-squared:  0.03703 
# F-statistic: 3.038 on 1 and 52 DF,  p-value: 0.08724


# Relationship between focal tree RGR and mean different mycorrhizal neighbor DBH 
RGR_size_diff_myco_neigh <- ggplot(diff_myco_growth_summary_all, aes(x = mean_neighbor_DBH, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Mean Different Mycorrhizal Neighbor DBH", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_diff_myco2_09 <- lm(mean_RGR ~ mean_neighbor_DBH, data = diff_myco_growth_summary_09)

summary(RGR_mod_diff_myco2_09)

# No significant relationship,  
# Multiple R-squared:  4.183e-05,	Adjusted R-squared:  -0.02123 
# F-statistic: 0.001966 on 1 and 47 DF,  p-value: 0.9648


RGR_mod_diff_myco2_20 <- lm(mean_RGR ~ mean_neighbor_DBH, data = diff_myco_growth_summary_20)

summary(RGR_mod_diff_myco2_20)

# No significant relationship,  
# Multiple R-squared:  0.01547,	Adjusted R-squared:  -0.00346 
# F-statistic: 0.8173 on 1 and 52 DF,  p-value: 0.3701


# Relationship between focal tree RGR and crowding index
RGR_crowd_diff_myco_neigh <- ggplot(diff_myco_growth_summary_all, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index - Different Mycorrhizal Neighborhoods", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_diff_myco3_09 <- lm(mean_RGR ~ crowding_index, data = diff_myco_growth_summary_09)

summary(RGR_mod_diff_myco3_09)

# Not significant relationship,  
# Multiple R-squared:  0.003952,	Adjusted R-squared:  -0.01724 
# F-statistic: 0.1865 on 1 and 47 DF,  p-value: 0.6678


RGR_mod_diff_myco3_20 <- lm(mean_RGR ~ crowding_index, data = diff_myco_growth_summary_20)

summary(RGR_mod_diff_myco3_20)

# Not significant relationship,  
# Multiple R-squared:  0.001604,	Adjusted R-squared:  -0.0176 
# F-statistic: 0.08352 on 1 and 52 DF,  p-value: 0.7737


# Being impacted by the outlier THPL again 
diff_myco_growth_summary_no_outlier <- diff_myco_growth_summary_all %>% filter(focal_cell != "H30") %>% droplevels()


# Visualize again without outlier 
RGR_crowd_diff_myco_no_outlier <- ggplot(diff_myco_growth_summary_no_outlier, aes(x = crowding_index, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ radius, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Crowding Index - Different Mycorrhizal Neighborhoods", y = expression("Mean RGR ("*yr^{-1}*")")) +
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


RGR_mod_diff_myco3_no_outlier_09 <- lm(mean_RGR ~ crowding_index, data = diff_myco_growth_summary_no_outlier_09)

summary(RGR_mod_diff_myco3_no_outlier_09)

# Still not a significant relationship 
# Multiple R-squared:  0.004705,	Adjusted R-squared:  -0.01693 
# F-statistic: 0.2174 on 1 and 46 DF,  p-value: 0.6432


diff_myco_growth_summary_no_outlier_20 <- diff_myco_growth_summary_20 %>% filter(focal_cell != "H30") %>% droplevels()


RGR_mod_diff_myco3_no_outlier_20 <- lm(mean_RGR ~ crowding_index, data = diff_myco_growth_summary_no_outlier_20)

summary(RGR_mod_diff_myco3_no_outlier_20)

# Still not a significant relationship 
# Multiple R-squared:  0.01164,	Adjusted R-squared:  -0.007742 
# F-statistic: 0.6005 on 1 and 51 DF,  p-value: 0.442


# Take a peak at just the data for THPL to see how this one outlier compares to the rest
# Doing this just in the 9m radius neighborhoods because the growth rate for the focal tree itself doesn't differ 
# between the neighborhood radii 

diff_myco_growth_summary_09$focal_species <- as.factor(diff_myco_growth_summary_09$focal_species)


THPL_diff_myco_09 <- diff_myco_growth_summary_09 %>% filter(focal_species == "THPL") %>% droplevels()


# Visualize crowding and growth for just THPL
RGR_crowd_diff_myco_THPL_09 <- ggplot(THPL_het_09, aes(x = crowding_index, y = mean_RGR)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  labs(x = "Crowding Index - Different Mycorrhizal Neighborhoods only THPL", y = expression("Mean RGR ("*yr^{-1}*")"))
  

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



# Gather the plots for the supplement to visualize the effects of removing the THPL with the high outlier CI value 
outlier_plots <- plot_grid(RGR_crowd_het_neigh, RGR_crowd_het_neigh_no_outlier, RGR_crowd_diff_myco_neigh, 
                           RGR_crowd_diff_myco_no_outlier, nrow = 2, ncol = 2, labels = c('(a)', '(b)', '(c)', '(d)'))

outlier_plots


## -- END -- ## 

