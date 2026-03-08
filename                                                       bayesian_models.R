# -----------------------------------------------------------------------------#
# Bayesian Modeling of tree growth  
# Original Author: L. McKinley Nevins 
# February 22, 2026
# Software versions:  R v 4.5.2
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     phyloseq v 1.54.0
#                     rstatix v 0.7.3
#                     ggpubr v 0.6.2
#                     cowplot v 1.2.0
#                     brms v 2.23.0
#                     bayesplot v 1.15.0
#                     rstan v 2.32.7
#                     StanHeaders v 2.32.10
#                     cmdstanr v 0.9.0
#                     tidybayes v 3.0.7
#                     pdp v 0.8.3
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(phyloseq); packageVersion("phyloseq")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(cowplot); packageVersion("cowplot")
library(brms); packageVersion("brms")
library(bayesplot); packageVersion("bayesplot")
library(rstan); packageVersion("rstan")
library(StanHeaders); packageVersion("StanHeaders")
library(cmdstanr); packageVersion("cmdstanr")
library(tidybayes); packageVersion("tidybayes")
library(pdp); packageVersion("pdp")

# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

cmdstanr::install_cmdstan(overwrite = TRUE)

#################################################################################
#                               Main workflow                                   #
#  Create bayesian regression models for the effects of tree traits, fungal     #
#  functions, and other factors on focal tree growth.                           #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/"
setwd(wd)


# Load in all relevant datasets that will eventually be organized at the individual tree 
# level into one dataframe 


# 1. Growth

# Growth dataset and calculate relative growth rate for each focal tree 
# Read in tree demographic data 
growth <- read.csv("./Demography/stems_WFDP_20250206_trimmed.csv")

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
diams$WFDP_Code <- diams$Stem_Tag

# Subset to just relevant columns - and keep Species column for later 

diams <- dplyr::select(diams, Species, mean_RGR, WFDP_Code)


# 2. PC1 and PC2 axis scores for both leaf and root traits for each focal tree

leaf_PCA <- read.csv("./Trait_Data/PCA/PCA_scores_leaf_traits.csv")

root_PCA <- read.csv("./Trait_Data/PCA/PCA_scores_root_traits.csv")

# Subset to just columns of interest 
leaf_PCs <- dplyr::select(leaf_PCA, leaf_PC1 = PC1, leaf_PC2 = PC2, WFDP_Code)

root_PCs <- dplyr::select(root_PCA, root_PC1 = PC1, root_PC2 = PC2, WFDP_Code)

# Merge together for later 
trait_PCs <- merge(leaf_PCs, root_PCs, by = "WFDP_Code")


# 3. PC axis scores for AM and EMF functions for each focal tree 

AM_trait_PCA <- read.csv("./Fungal_Communities/PCA_scores_AM_traits.csv")

# Subset columns 
AM_traits <- dplyr::select(AM_trait_PCA, AM_PC1 = PC1, AM_PC2 = PC2, WFDP_Code)

EM_trait_PCA <- read.csv("./Fungal_Communities/PCA_scores_EMF_traits.csv")

# Subset columns 
EMF_traits <- dplyr::select(EM_trait_PCA, EM_PC1 = PC1, EM_PC2 = PC2, WFDP_Code)


# 4. Neighborhood data

# Load in summaries for same mycorrhizal/different mycorrhizal neighborhoods for each focal tree in 9m radius neighborhoods 

same_myco_summary_09 <- read.csv("./Demography/Neighborhood_files/same_myco_crowd_summary_9m.csv", row.names = 1)

diff_myco_summary_09 <- read.csv("./Demography/Neighborhood_files/diff_myco_crowd_summary_9m.csv", row.names = 1)

same_myco_neigh <- dplyr::select(same_myco_summary_09, same_myco_09_mean_neigh_DBH = mean_neighbor_DBH, 
                                 same_myco_09_num_neigh = num_neighbors, same_myco_09_CI = crowding_index, WFDP_Code = focal_stem_tag)


diff_myco_neigh <- dplyr::select(diff_myco_summary_09, diff_myco_09_mean_neigh_DBH = mean_neighbor_DBH, 
                                 diff_myco_09_num_neigh = num_neighbors, diff_myco_09_CI = crowding_index, WFDP_Code = focal_stem_tag)


# Note: these neighborhoods have 58 trees in the same myco dataset because only the two trees on the edge of 
# the plot were excluded. But the diff myco neigh dataset only has 49 trees in it because some of the focal 
# trees didn't have any neighbors within the 9m radius that were different mycorrhizal associations than their own

# Canopy data 
canopy_all <- read.csv("./Enviro_Data/WFDP_all_canopy_data.csv")


# Right now just focus on the ones that had significant relationships overall with growth - canopy openness at 
# a height of 10m in the 20m radius, and mean canopy height in the 9m radius neighborhoods 

canopy <- dplyr::select(canopy_all, WFDP_Code = Stem_Tag, mean_can_height_9 = X9m_mean, can_open_10_20 = p_10_20)


# 5. Environmental data? 

# Load in general environmental data right now, can add in tree specific soil microclimate and nutrients if/when 
# we have that data 

enviro <- read.csv("./Enviro_Data/WFDP_enviro_data_all.csv")

# Grab just cell, slope, aspect, elevation for now 

env <- dplyr::select(enviro, WFDP_Code, Association, Cell, slope, aspect, elevation_m)

#################################################################################

####################
# (2) PACKAGE DATA
####################

# All component data pieces now have WFDP_Code as a column for merging 

# Growth and tree leaf and root traits 
mod_df <- merge(diams, trait_PCs, by = "WFDP_Code")

# Left join so as not to lose rows for trees that don't have either AM or EM trait data 
mod_df <- left_join(mod_df, AM_traits, by = "WFDP_Code")
mod_df <- left_join(mod_df, EMF_traits, by = "WFDP_Code")

# Left join with the mycorrhizal neighborhood data 
mod_df <- left_join(mod_df, same_myco_neigh, by = "WFDP_Code")
mod_df <- left_join(mod_df, diff_myco_neigh, by = "WFDP_Code")


# For the diff_myco_neigh data, the column for number of neighbors now has NAs for trees that don't 
# have any diff myco neighbors, but these can be changed to zeros for just this column 

mod_df <- mod_df %>%
  mutate(diff_myco_09_num_neigh = replace_na(diff_myco_09_num_neigh, 0))
  
# Do the same thing to load zeros into the AM and EMF trait PCA columns, becasue the model structure 
# will be built to only refer to columns that are relevant for each host as either AM, EMF, or Dual 

mod_df <- mod_df %>%
  replace_na(list(AM_PC1 = 0, AM_PC2 = 0, EM_PC1 = 0, EM_PC2 = 0))


# Add canopy data 
mod_df <- left_join(mod_df, canopy, by = "WFDP_Code")

# Add environmental data 
mod_df <- left_join(mod_df, env, by = "WFDP_Code")


# Conclusion: One final dataset that has all of the data for each focal tree, with some gaps. 

#################################################################################

###################
# (3) MODEL PREP             
###################

# Goal is to generate a series of Bayesian models in brms to compare the effects of adding 
# different abiotic and biotic components to modeling focal tree relative growth rate (RGR)

# Start very simple with models that are just including tree traits, maybe as leaves and roots 
# separately and then together

# Then adding in the mycorrhizal community 

# Then the tree neighborhoods 

# Then environmental factors like the canopy conditions and slope, aspect, and elevation of each tree? 

#   RGR ~ Leaf_PC1 + Leaf_PC2 + Root_PC1 + Root_PC2 +
#   Tree_type * (AM_PC1 + AM_PC2 + EM_PC1 + EM_PC2) +
#   Neighborhood_terms + Environment_terms

####### -- 

# Check the distribution of the growth data for the trees 
hist(mod_df$mean_RGR,breaks=50)
qqnorm(mod_df$mean_RGR)

# Not a super normal distribution, but I have some negative RGR values that would not allow me to log transform
# so I'm going to run this like it is normal enough right now 


# The PCA-related variables for the leaf and root traits have been centered already, and these are 
# the tree positions on the axes, but they do need to be scaled 
mod_df$leaf_PC1_std <- scale(mod_df$leaf_PC1, scale = TRUE)
mod_df$leaf_PC2_std <- scale(mod_df$leaf_PC2, scale = TRUE)
mod_df$root_PC1_std <- scale(mod_df$root_PC1, scale = TRUE)
mod_df$root_PC2_std <- scale(mod_df$root_PC2, scale = TRUE)


# The AM and EMF traits do need to be centered and scaled because I added in zeros for the trees that didn't 
# have the corresponding community 
mod_df$AM_PC1_std <- scale(mod_df$AM_PC1, center = TRUE, scale = TRUE)
mod_df$AM_PC2_std <- scale(mod_df$AM_PC2, center = TRUE, scale = TRUE)
mod_df$EM_PC1_std <- scale(mod_df$EM_PC1, center = TRUE, scale = TRUE)
mod_df$EM_PC2_std <- scale(mod_df$EM_PC2, center = TRUE, scale = TRUE)


# Centering and scaling other variables prior to modeling 
mod_df$same_myco_09_mean_neigh_DBH_std<-scale(mod_df$same_myco_09_mean_neigh_DBH, center = TRUE, scale = TRUE)
mod_df$same_myco_09_num_neigh_std<-scale(mod_df$same_myco_09_num_neigh, center = TRUE, scale = TRUE)
mod_df$same_myco_09_CI_std<-scale(mod_df$same_myco_09_CI, center = TRUE, scale = TRUE)
mod_df$diff_myco_09_mean_neigh_DBH_std<-scale(mod_df$diff_myco_09_mean_neigh_DBH, center = TRUE, scale = TRUE)
mod_df$diff_myco_09_num_neigh_std<-scale(mod_df$diff_myco_09_num_neigh, center = TRUE, scale = TRUE)
mod_df$diff_myco_09_CI_std<-scale(mod_df$diff_myco_09_CI, center = TRUE, scale = TRUE)
mod_df$mean_can_height_9_std <- scale(mod_df$mean_can_height_9, center = TRUE, scale = TRUE)
mod_df$can_open_10_20_std <- scale(mod_df$can_open_10_20, center = TRUE, scale = TRUE)
mod_df$slope_std <- scale(mod_df$slope, center = TRUE, scale = TRUE)
mod_df$aspect_std <- scale(mod_df$aspect, center = TRUE, scale = TRUE)
mod_df$elevation_std <- scale(mod_df$elevation_m, center = TRUE, scale = TRUE)

# Set some of the variables as factors 
mod_df$WFDP_Code <- as.factor(mod_df$WFDP_Code)
mod_df$Species <- as.factor(mod_df$Species)
mod_df$Cell <- as.factor(mod_df$Cell)
mod_df$Association <- as.factor(mod_df$Association)


# Do some checks of the variables 

summary(mod_df$leaf_PC1)
sd(mod_df$leaf_PC1)

# mean = 0, SD = 2.30

#################################################################################

###################
# (4) RUN MODELS         
###################

options(mc.cores=parallel::detectCores()) 

options(brms.backend = "cmdstanr")


# potentially use family = student() to account for slight right-tailing of the RGR data 



## Add species as a random effect? (1| Species) 


# MODEL 1 - Relationship of growth to Leaf PC1 values 

priors <- c(
  prior(normal(0, 0.1), class = "Intercept"),
  prior(normal(0, 0.1), class = "b"), # temporarily loosening this prior to see if it changes the estimates 
  prior(student_t(3, 0, 0.1), class = "sigma"),
  prior(gamma(2, 0.1), class = "nu"))


# Simulate from the priors to check how they behave 
mod_prior <- brm(mean_RGR ~ leaf_PC1,
  data = mod_df, family = student(),
  prior = priors, sample_prior = "only", backend = "cmdstanr")

pp_check(mod_prior)

# Run the model 
mod_01 <- brm(mean_RGR ~ leaf_PC1, 
  data = mod_df, family = student(),
  prior = priors, backend = "cmdstanr", chains = 4, iter = 1500, warmup = 500)

summary(mod_01)
pp_check(mod_01)




# 2. Generate and plot conditional effects
# This plots the marginal effect of leaf_PC1
plot(conditional_effects(mod_01, effects = "leaf_PC1"), points = TRUE)

# 3. For interactions, none in this model 
plot(conditional_effects(mod_01, effects = "zBase:Trt"), points = TRUE)



## Need some edits here 

# 1. Compute partial dependence data
pd <- pdp::partial(mod_01, pred.var = "leaf_PC1", grid.resolution = 30)

# 2. Plot using autoplot
autoplot(pd, rug = TRUE, train = epilepsy) + theme_minimal()

# 3. For 2D interaction (3D surface)
pd2d <- partial(fit, pred.var = c("zAge", "zBase"))
plotPartial(pd2d, levelplot = TRUE, contour = TRUE)









# MODEL 2 - Relationship of growth to Leaf PC1 and PC2 values 

# Simulate from the priors to check how they behave 
mod_prior2 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2, 
                 data = mod_df, family = student(),
                 prior = priors, sample_prior = "only", backend = "cmdstanr")

pp_check(mod_prior2)

# Run the model 
mod_02 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2, 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 4, iter = 1500, warmup = 500)

summary(mod_02)
pp_check(mod_02)


# This plots the marginal effect of leaf_PC1 and leaf_PC2
plot(conditional_effects(mod_02, effects = "leaf_PC2"), points = TRUE)








# MODEL 3 - Relationship of growth to Root PC1 values 

# Simulate from the priors to check how they behave 
mod_prior3 <- brm(mean_RGR ~ root_PC1, 
                  data = mod_df, family = student(),
                  prior = priors, sample_prior = "only", backend = "cmdstanr")

pp_check(mod_prior3)

# Run the model 
mod_03 <- brm(mean_RGR ~ root_PC1, 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 4, iter = 1500, warmup = 500)

summary(mod_03)
pp_check(mod_03)


# This plots the marginal effect of root_PC1
plot(conditional_effects(mod_03, effects = "root_PC1"), points = TRUE)






# MODEL 4 - Relationship of growth to Root PC1 and Root PC2 values 

# Simulate from the priors to check how they behave 
mod_prior4 <- brm(mean_RGR ~ root_PC1 + root_PC2, 
                  data = mod_df, family = student(),
                  prior = priors, sample_prior = "only", backend = "cmdstanr")

pp_check(mod_prior4)

# Run the model 
mod_04 <- brm(mean_RGR ~ root_PC1 + root_PC2,  
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 4, iter = 1500, warmup = 500)

summary(mod_04)
pp_check(mod_04)


# This plots the marginal effect of root_PC2
plot(conditional_effects(mod_04, effects = "root_PC2"), points = TRUE)






# MODEL 5 - 

mod_05 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + root_PC1 + root_PC2 +
    Association * (AM_PC1 + AM_PC2 + EM_PC1 + EM_PC2), 
  data = mod_df, family = student(),
  prior = priors, backend = "cmdstanr", chains = 4, iter = 1500, warmup = 500)


summary(mod_05)
pp_check(mod_05)




# This plots the marginal effect of different terms 
plot(conditional_effects(mod_05, effects = "root_PC2"), points = TRUE)
plot(conditional_effects(mod_05, effects = "AM_PC1"), points = TRUE) # These have fewer points and maybe some weird 
plot(conditional_effects(mod_05, effects = "AM_PC2"), points = TRUE) # behavior because there are not that many AM trees
plot(conditional_effects(mod_05, effects = "EM_PC1"), points = TRUE)
plot(conditional_effects(mod_05, effects = "EM_PC2"), points = TRUE) # These two look very weird 







