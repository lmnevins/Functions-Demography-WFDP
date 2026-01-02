# -----------------------------------------------------------------------------#
# Exploring WFDP tree growth data 
# Original Author: L. McKinley Nevins 
# June 25, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 3.5.1
#                    
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")


#################################################################################
#                               Main workflow                                   #
#  Explore the growth data for all of the trait trees sampled at WFDP. Assess   #
#  for differences between species, or any relationships to environmental       #
#  variables so far. Then relate to trait information for the trees themselves, #
#  and for their fungal communities. Model how tree and fungal traits impact    #
#  tree growth.                                                                 #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
############### 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)


# Load in tree demography file 
growth <- read.csv("stems_WFDP_20250206_trimmed.csv")
# using the trimmed datafile, which has had two trees removed from what I originally 
# samples as trait trees - one PSME that was the only one, and one TABR that did not have 
# any associated fungal community data 


# Load in environmental data so far 

# This can be updated once the krieged environmental data is in hand 

env_data <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# make STEM_TAG column for merging 
env_data$Stem_Tag <- env_data$WFDP_Code

# Trim to just relevant columns 
env_data <- select(env_data, Cell, EM_Sample_Name, slope, aspect, elevation_m, Association, Stem_Tag)

# Make association a factor 
env_data$Association <- as.factor(env_data$Association)

# Load in the tree PC values from the trait PCA
tree_pc <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/tree_PC_scores.csv")

# make STEM_TAG column for merging 
tree_pc$Stem_Tag <- tree_pc$WFDP_Code

##########################################
# each tree is identified by its WFDP stem_tag. There are three census years - 2011, 
# 2016, and 2021. I want to collect the diameters from the first and last census year 
# to calculate how much they have grown. 


# get starting and ending diameters for all trees
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

#### 
# MERGE

# Pair diams with the enviro data according to the tree STEM_TAG
diam_env <- merge(diams, env_data, by = "Stem_Tag")

# Add in the PC values for each tree 
diam_env_pc <- merge(diam_env, tree_pc, by = "Stem_Tag")

#################################################################################
####################### -- 
# (1) DATA EXPLORATION
####################### -- 

# Explore diameter difference relationships with slope, aspect, and elevation 

## Slope
slope <- ggplot(diam_env_pc, aes(x = slope, y = diam_diff)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope

# test relationships 
lm_slope <- lm(diam_diff ~ slope, data = diam_env_pc)
summary(lm_slope) # NOT SIGNIFICANT


## Aspect
aspect <- ggplot(diam_env_pc, aes(x = aspect, y = diam_diff)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect

# test relationships 
lm_aspect <- lm(diam_diff ~ aspect, data = diam_env_pc)
summary(lm_aspect) # NOT SIGNIFICANT


## Elevation
elev <- ggplot(diam_env_pc, aes(x = elevation_m, y = diam_diff)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
elev

# test relationships 
lm_elev <- lm(diam_diff ~ elevation_m, data = diam_env_pc)
summary(lm_elev) # NOT SIGNIFICANT


# Explore RGR relationships with slope, aspect, and elevation 

## Slope
slope2 <- ggplot(diam_env_pc, aes(x = slope, y = RGR)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope2

# test relationships 
lm_slope2 <- lm(RGR ~ slope, data = diam_env_pc)
summary(lm_slope2) # NOT SIGNIFICANT


## Aspect
aspect2 <- ggplot(diam_env_pc, aes(x = aspect, y = RGR)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect2

# test relationships 
lm_aspect2 <- lm(RGR ~ aspect, data = diam_env_pc)
summary(lm_aspect2) # NOT SIGNIFICANT


## Elevation
elev2 <- ggplot(diam_env_pc, aes(x = elevation_m, y = RGR)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
elev2

# test relationships 
lm_elev2 <- lm(RGR ~ elevation_m, data = diam_env_pc)
summary(lm_elev2) # NOT SIGNIFICANT


# Explore diameter difference relationships with the PC1 and PC2 trait values 

## PC1
PC1_diam <- ggplot(diam_env_pc, aes(x = PC1, y = diam_diff)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
PC1_diam

# test relationships 
lm_PC1_diam <- lm(diam_diff ~ PC1, data = diam_env_pc)
summary(lm_PC1_diam) # NOT SIGNIFICANT


## PC2
PC2_diam <- ggplot(diam_env_pc, aes(x = PC2, y = diam_diff)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
PC2_diam

# test relationships 
lm_PC2_diam <- lm(diam_diff ~ PC2, data = diam_env_pc)
summary(lm_PC2_diam) # VERY SIGNIFICANT

# Adjusted R-squared:  0.1433,  p-value: 0.001675


# Explore RGR relationships with the PC1 and PC2 trait values 

## PC1
PC1_RGR <- ggplot(diam_env_pc, aes(x = PC1, y = RGR)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
PC1_RGR

# test relationships 
lm_PC1_RGR <- lm(RGR ~ PC1, data = diam_env_pc)
summary(lm_PC1_RGR) # NOT SIGNIFICANT


## PC2
PC2_RGR <- ggplot(diam_env_pc, aes(x = PC2, y = RGR)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
PC2_RGR

# test relationships 
lm_PC2_RGR <- lm(RGR ~ PC2, data = diam_env_pc)
summary(lm_PC2_RGR) # SIGNIFICANT

# Adjusted R-squared:  0.07245,  p-value: 0.02123


#### Nicer Plotting 

#set colors for hosts
                  # ABAM        ABG         ALRU         CONU        TABR          THPL        TSHE        
all_hosts <- c("#0D0887FF", "#5402A3FF", "#B93289FF", "#DB5C68FF", "#F48849FF", "#ffe24cFF", "#fffd66")


# Plot the Results by host
PC2_RGR_nice <- ggplot(diam_env_pc, aes(x = PC2, y = RGR, color = Host_ID, shape = Association)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  scale_colour_manual(values=all_hosts, 
                      name="Host Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  scale_shape_manual(
    values = c("EM" = 16, "Both" = 17), name = "Mycorrhizal Association") +
  labs(x = "PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")"), color = "Host_ID") +
  theme(legend.title = element_text(colour="black", size=14, face = "bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text = element_text(colour="black", size = 12)) +
  theme(axis.title = element_text(colour="black", size = 14))

PC2_RGR_nice



# Faceted option to see individual species trends 
PC2_RGR_nice2 <- ggplot(diam_env_pc, aes(x = PC2, y = diam_diff)) +
  geom_point(aes(color = Host_ID), size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Host_ID) +
  theme_bw() +
  scale_colour_manual(
    values = all_hosts,
    name = "Host Species",
    breaks = c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
    labels = c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")
  ) +
  labs(x = "PC2 Value", y = expression("Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

PC2_RGR_nice2



