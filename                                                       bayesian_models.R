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
#                     forcats v 1.0.1
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
library(forcats); packageVersion("forcats")

# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# cmdstanr::install_cmdstan(overwrite = TRUE)

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

## EDIT: using updated values that do not have 13C

leaf_PCA <- read.csv("./Trait_Data/PCA/PCA_scores_leaf_traits_no13C.csv")

root_PCA <- read.csv("./Trait_Data/PCA/PCA_scores_root_traits_no13C.csv")

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
# a height of 10m in the 9m radius, and mean canopy height in the 9m radius neighborhoods 

canopy <- dplyr::select(canopy_all, WFDP_Code = Stem_Tag, mean_can_height_9 = X9m_mean, can_open_10_09 = prop_10m_9)


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


# Exclude the two focal trees that couldn't have their neighborhoods modeled, so all models have the same 
# number of replicates 

mod_df$WFDP_Code <- as.factor(mod_df$WFDP_Code)

mod_df <- mod_df %>%
  filter(WFDP_Code != "29-0766") %>%
  filter(WFDP_Code != "38-0935") %>%
  droplevels()


# Conclusion: One final dataset that has all of the data for each of the 58 suitable focal tree, with some gaps
# that can be managed in the next step. 

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
#   AM_PC1 + EM_PC1 + Neighborhood_terms + Environment_terms

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
mod_df$can_open_10_09_std <- scale(mod_df$can_open_10_09, center = TRUE, scale = TRUE)
mod_df$slope_std <- scale(mod_df$slope, center = TRUE, scale = TRUE)
mod_df$aspect_std <- scale(mod_df$aspect, center = TRUE, scale = TRUE)
mod_df$elevation_std <- scale(mod_df$elevation_m, center = TRUE, scale = TRUE)

# Set some of the variables as factors 
mod_df$WFDP_Code <- as.factor(mod_df$WFDP_Code)
mod_df$Species <- as.factor(mod_df$Species)
mod_df$Cell <- as.factor(mod_df$Cell)
mod_df$Association <- as.factor(mod_df$Association)


#################################################################################

###################
# (4) RUN MODELS         
###################

options(mc.cores=parallel::detectCores()) 

options(brms.backend = "cmdstanr")

# Using family = student() to account for slight right-tailing of the RGR data 

# Species as a random effect, which means that trait and other effects are estimated after 
# controlling for species identity, where we would inherently expect that species vary in their 
# growth rates (and we know this is true for the species in this study)


## Define priors:


# Tree RGR ranges from -0.01 to 0.04 
# All predictors are centered and scaled around zero 
# Growth rates between species can vary about 0.02 - 0.03 

# Overall goal to use weakly informative priors 

priors <- c(
  prior(normal(0, 0.05), class = "Intercept"),
  prior(normal(0, 0.01), class = "b"), # prior for term effects on growth, don't expect any will be huge 
  prior(student_t(3, 0, 0.02), class = "sd"), # prior for species 
  prior(student_t(3, 0, 0.05), class = "sigma"), # prior for residual variance 
  prior(gamma(2, 0.1), class = "nu"))


# REPS: 

# 20,000 iterations on 5 chains, with a 3,000 burn-in warmup. 
# was done in Germain & Lutz 2021, Ecology 

# Doing 10,000 iterations on 5 chains with a 1,000 burn-in warmup because my sample size is much 
# smaller and that many iterations seems like overkill 



# MODEL 1 - Relationship of growth to Leaf PC1 and PC2 values 

# Simulate from the priors to check how they behave 
mod_prior <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + (1| Species), 
  data = mod_df, family = student(),
  prior = priors, sample_prior = "only", backend = "cmdstanr")

pp_check(mod_prior)

# Run the model 
mod_01 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + (1| Species), 
  data = mod_df, family = student(),
  prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.85)

summary(mod_01)
pp_check(mod_01)
pairs(mod_01)
mcmc_trace(mod_01)
prior_summary(mod_01)


# Generate and plot conditional effects
# This plots the marginal effect of leaf_PC1
plot(conditional_effects(mod_01, effects = "leaf_PC1"), points = TRUE)
plot(conditional_effects(mod_01, effects = "leaf_PC2"), points = TRUE)


# MODEL 2 - Relationship of growth to Root PC1 and PC2 values 

# Run the model 
mod_02 <- brm(mean_RGR ~ root_PC1 + root_PC2 + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.95)

summary(mod_02)
pp_check(mod_02)
mcmc_trace(mod_02)


# This plots the marginal effect of root_PC1
plot(conditional_effects(mod_02, effects = "root_PC1"), points = TRUE)
plot(conditional_effects(mod_02, effects = "root_PC2"), points = TRUE)



## From these first two models, the posterior distributions reflect the relationships of the PC trait 
# values to tree growth. This is assessing the same thing as Hyp. 2, so these plots could be 
# used instead 















# MODEL 3 - Relationship of growth to all trait PCs

# Run the model 
mod_03 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + root_PC1 + root_PC2 + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.95)

summary(mod_03)
pp_check(mod_03)
mcmc_trace(mod_03)



# MODEL 4 - Relationship of growth to myco trait PC1 values 

# Run the model 
mod_04 <- brm(mean_RGR ~ AM_PC1 + EM_PC1 + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", save_pars = save_pars(all = TRUE), 
              chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.975)

summary(mod_04)
pp_check(mod_04)
mcmc_trace(mod_04)



# MODEL 5 - Relationship of growth to trait PCs and mycorrhizal PCs

# Run the model 
mod_05 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + root_PC1 + root_PC2 +
                AM_PC1 + EM_PC1 + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.95)

summary(mod_05)
pp_check(mod_05)
mcmc_trace(mod_05)


# This plots the marginal effect of different terms 
plot(conditional_effects(mod_05, effects = "AM_PC1"), points = TRUE) 
plot(conditional_effects(mod_05, effects = "EM_PC1"), points = TRUE)

# Check the data structures for the mycorrhizal PC values 
ggplot(mod_df, aes(Association, EM_PC1)) + geom_jitter()
ggplot(mod_df, aes(Association, AM_PC1)) + geom_jitter()




# MODEL 6 - Relationship of growth to neighborhood context

# Including two neighborhood factors that had significant relationships with tree RGR

mod_06 <- brm(mean_RGR ~ same_myco_09_mean_neigh_DBH_std + diff_myco_09_num_neigh_std + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.95)

summary(mod_06)
pp_check(mod_06)
mcmc_trace(mod_06)


# This plots the marginal effect of different terms 
plot(conditional_effects(mod_06, effects = "same_myco_09_mean_neigh_DBH_std"), points = TRUE) 
plot(conditional_effects(mod_06, effects = "diff_myco_09_num_neigh_std"), points = TRUE)




# MODEL 7 - Relationship of growth to trait PCs, mycorrhizal PCs, and neighborhood context

# Including two neighborhood factors that had significant relationships with tree RGR

mod_07 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + root_PC1 + root_PC2 +
    AM_PC1 + EM_PC1 + same_myco_09_mean_neigh_DBH_std + diff_myco_09_num_neigh_std + (1| Species), 
  data = mod_df, family = student(),
  prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.85)

summary(mod_07)
pp_check(mod_07)
mcmc_trace(mod_07)




# MODEL 8 - Relationship of growth to canopy context

mod_08 <- brm(mean_RGR ~ mean_can_height_9_std + can_open_10_09_std + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.85)

summary(mod_08)
pp_check(mod_08)
mcmc_trace(mod_08)


# This plots the marginal effect of different terms 
plot(conditional_effects(mod_08, effects = "mean_can_height_9_std"), points = TRUE) 
plot(conditional_effects(mod_08, effects = "can_open_10_09_std"), points = TRUE)



# MODEL 9 - Relationship of growth to trait PCs, mycorrhizal PCs, neighborhood context, and canopy context

mod_09 <- brm(mean_RGR ~ leaf_PC1 + leaf_PC2 + root_PC1 + root_PC2 +
                AM_PC1 + EM_PC1 + same_myco_09_mean_neigh_DBH_std + diff_myco_09_num_neigh_std + 
                mean_can_height_9_std + can_open_10_09_std + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", save_pars = save_pars(all = TRUE), 
              chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.85)

summary(mod_09)
pp_check(mod_09)
mcmc_trace(mod_09)




# MODEL 10 - Relationship of growth to the interactive effects of the three model terms that had significant 
# effects on growth in the most comprehensive model 

mod_10 <- brm(mean_RGR ~ leaf_PC1*diff_myco_09_num_neigh_std*can_open_10_09_std + (1| Species), 
              data = mod_df, family = student(),
              prior = priors, backend = "cmdstanr", save_pars = save_pars(all = TRUE), 
              chains = 5, iter = 10000, warmup = 1000, adapt_delta = 0.85)

summary(mod_10)
pp_check(mod_10)
mcmc_trace(mod_10)


#################################################################################

#######################
# (5) COMPARE MODELS         
#######################

### Model comparisons using leave-one-out cross validation 

loo_01 <- loo(mod_01)
loo_02 <- loo(mod_02)
loo_03 <- loo(mod_03)
loo_04 <- loo(mod_04, moment_match = TRUE)
loo_05 <- loo(mod_05)
loo_06 <- loo(mod_06) 
loo_07 <- loo(mod_07) 
loo_08 <- loo(mod_08)
loo_09 <- loo(mod_09, moment_match = TRUE) # received message of 2 observations with pareto_k > 0.7, so ran 
# moment match and recompiled the model. 
loo_10 <- loo(mod_10, moment_match = TRUE, reloo = TRUE)

# This requires having save_pars = save_pars(all = TRUE),  in the original model call 


## Some remaining issues with the refitting of the model 10. Going forward just for exploration, but this model 
# structure isn't working very well. 


model_comp <- loo_compare(loo_01, loo_02, loo_03, loo_04, loo_05, loo_06, loo_07, loo_08, loo_09, loo_10)
print(model_comp)

plot(model_comp)

comp_df <- as.data.frame(model_comp) %>%
  rownames_to_column("mod_num")


# create dataframe of model number and the model terms for matching 

mod_num <- c("mod_01", "mod_02", "mod_03", "mod_04", "mod_05", "mod_06", "mod_07", "mod_08", "mod_09", "mod_10")

mod_terms <- c("Leaf Traits", "Root Traits", "Leaf + Root Traits", "Myco Traits", 
                   "Leaf + Root + Myco Traits", "Neighborhood", "Leaf + Root + Myco + Neigh", 
                   "Canopy", "Leaf + Root + Myco + Neigh + Canopy", "Leaf PC1 x Neigh x Canopy")

mod_details <- data.frame(mod_num, mod_terms)


# Merge comparisons with model details 
comp_df <- merge(mod_details, comp_df, by = "mod_num")



# compute model weights, which show the probability each model is the best predictive model
weights <- loo_model_weights(list(loo_01, loo_02, loo_03, loo_04, loo_05, loo_06, loo_07, loo_08, 
                                  loo_09, loo_10),method = "stacking")


weight_df <- data.frame(model = c("Leaf Traits", "Root Traits", "Leaf + Root Traits", "Myco Traits", 
                                  "Leaf + Root + Myco Traits", "Neighborhood", "Leaf + Root + Myco + Neigh", 
                                  "Canopy", "Leaf + Root + Myco + Neigh + Canopy", "Leaf PC1 x Neigh x Canopy"), weight = weights)


# Set order for plotting - putting them in the same order as the overall model performance 
# values down below 

## !! Not changing this to reflect Mod 10 differences, just putting it at the top right now !! 
new_mod_order <- c("Leaf + Root + Myco Traits", "Leaf + Root Traits", "Root Traits", "Myco Traits", "Leaf Traits",
                   "Leaf + Root + Myco + Neigh", "Neighborhood", "Canopy", "Leaf + Root + Myco + Neigh + Canopy", 
                   "Leaf PC1 x Neigh x Canopy")

weight_df$model <- factor(weight_df$model, levels = new_mod_order)



# Plot weights 
weight_plot <- ggplot(weight_df, aes(y = model, x = weight)) +
  geom_col() +
  labs(x = "Model Predictive Weight",
       title = "", 
       y = "Model") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 14, colour="black"), 
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.text = element_text(size = 14, colour="black"), 
        legend.title = element_text(size = 14, face = "bold", colour="black"))

weight_plot

# Weight ≈ 1 clearly best model; Split weights = multiple models explain data similarly




# plot predictive performance - ordered relative to the best model 

p1 <- ggplot(comp_df,
             aes(x = reorder(mod_terms, elpd_diff),
                 y = elpd_diff)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = elpd_diff - se_diff,
                    ymax = elpd_diff + se_diff),
                width = 0.25) +
  coord_flip() +
  labs(x = "Model",
       y = "Δ in Expected Log Predictive Density",
       title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 14, colour="black"), 
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.text = element_text(size = 14, colour="black"), 
        legend.title = element_text(size = 14, face = "bold", colour="black"))

p1

# Interpretation: The model at zero is considered the best one, and then each model is compared relative to that 
# ELDP less than 4 suggests models are similar in performance, but greater than 4 shows they are very different

# So I believe the other models being different but having a lot of overlap shows they are different but 
# not super different, but the Mod_01 of just the leaf traits performs much worse. 



# Gather posterior trait effects

leaf_draws <- mod_01 %>%
  gather_draws(b_leaf_PC1, b_leaf_PC2) %>%
  mutate(model = "Leaf Traits")


root_draws <- mod_02 %>%
  gather_draws(b_root_PC1, b_root_PC2) %>%
  mutate(model = "Root Traits")


leaf_root_draws <- mod_03 %>%
  gather_draws(b_leaf_PC1, b_leaf_PC2, b_root_PC1, b_root_PC2) %>%
  mutate(model = "Leaf + Root Traits")


myco_draws <- mod_04 %>%
  gather_draws(b_AM_PC1, b_EM_PC1) %>%
  mutate(model = "Myco Traits")


leaf_root_myco_draws <- mod_05 %>%
  gather_draws(b_leaf_PC1, b_leaf_PC2, b_root_PC1, b_root_PC2, b_AM_PC1, b_EM_PC1) %>%
  mutate(model = "Leaf + Root + Myco Traits")


neigh_draws <- mod_06 %>%
  gather_draws(b_same_myco_09_mean_neigh_DBH_std, b_diff_myco_09_num_neigh_std) %>%
  mutate(model = "Neighborhood")


leaf_root_myco_neigh_draws <- mod_07 %>%
  gather_draws(b_leaf_PC1, b_leaf_PC2, b_root_PC1, b_root_PC2, b_AM_PC1, b_EM_PC1, 
               b_same_myco_09_mean_neigh_DBH_std, b_diff_myco_09_num_neigh_std) %>%
  mutate(model = "Leaf + Root + Myco Traits + Neigh")


canopy_draws <- mod_08 %>%
  gather_draws(b_mean_can_height_9_std, b_can_open_10_09_std) %>%
  mutate(model = "Canopy")


leaf_root_myco_neigh_canopy_draws <- mod_09 %>%
  gather_draws(b_leaf_PC1, b_leaf_PC2, b_root_PC1, b_root_PC2, b_AM_PC1, b_EM_PC1, 
               b_same_myco_09_mean_neigh_DBH_std, b_diff_myco_09_num_neigh_std, 
               b_mean_can_height_9_std, b_can_open_10_09_std) %>%
  mutate(model = "Leaf + Root + Myco Traits + Neigh + Canopy")


leaf_neigh_canopy_interaction_draws <- mod_10 %>%
  gather_draws(b_leaf_PC1, b_diff_myco_09_num_neigh_std, b_can_open_10_09_std) %>%
  mutate(model = "Leaf PC1 x Neigh x Canopy")




# combine all 
draws_df <- bind_rows(
  leaf_draws,
  root_draws,
  leaf_root_draws, 
  myco_draws,
  leaf_root_myco_draws, 
  neigh_draws,
  leaf_root_myco_neigh_draws, 
  canopy_draws,
  leaf_root_myco_neigh_canopy_draws, 
  leaf_neigh_canopy_interaction_draws)



# Clean up labels of coefficients
draws_df <- draws_df %>%
  mutate(term = dplyr::recode(.variable,
                       "b_leaf_PC1" = "Leaf Trait PC1",
                       "b_leaf_PC2" = "Leaf Trait PC2",
                       "b_root_PC1" = "Root Trait PC1",
                       "b_root_PC2" = "Root Trait PC2", 
                       "b_AM_PC1" = "AM Trait PC1", 
                       "b_EM_PC1" = "ECM Trait PC1", 
                       "b_same_myco_09_mean_neigh_DBH_std" = "Con-Myco Neighbor Size", 
                       "b_diff_myco_09_num_neigh_std" = "Num Hetero-Myco Neighbors", 
                       "b_mean_can_height_9_std" = "Mean Neighbor Canopy Height", 
                       "b_can_open_10_09_std" = "Prop. Canopy Openness - 10 m"))



# Plot effects 
p2 <- ggplot(draws_df,
       aes(x = .value, y = term, fill = model)) +
  stat_halfeye(.width = c(.5, .8, .95), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~model) +
  labs(
    x = "Effect on RGR",
    y = "Model Term") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "bottom", legend.justification = "left")

p2

# Format is not as effective because it splits up the terms, and the facetting makes 
# things very small 


# library(patchwork)

# p1 + p2 + plot_layout(ncol = 2)



# Create summaries and calculate 95% confidence intervals around the median 
coef_summary <- draws_df %>%
  group_by(model, term) %>%
  summarise(
    median = median(.value),
    lower = quantile(.value, 0.025),
    upper = quantile(.value, 0.975),
    .groups = "drop"
  )


new_term_order <- c("Prop. Canopy Openness - 10 m", "Mean Neighbor Canopy Height", 
                    "Num Hetero-Myco Neighbors", "Con-Myco Neighbor Size", 
                    "ECM Trait PC1", "AM Trait PC1", "Root Trait PC2", 
                    "Root Trait PC1", "Leaf Trait PC2","Leaf Trait PC1")

coef_summary$term <- factor(coef_summary$term, levels = new_term_order)

# Code for if each term is significant or not according to the confidence intervals 
coef_summary <- coef_summary %>%
  mutate(sig = ifelse(lower > 0 | upper < 0, "significant", "nonsignificant"))


# set model colors

                # Leaf      Root   Leaf + Root    Myco      Neigh       Canopy                                      #interaction
my_cols <- c( "darkgreen", "brown",  "gray60", "#0072C2", "purple3",  "#F5C710", "#E69F00",  "#D55E00",  "#CC79A7", "navy")



## Added in the tenth model for the interactions here. Would need to take it out just by removing the last color and the last 
# items in the breaks and labels lists, along with it from the dataset generated up above. 


# Create a forest plot 
p3 <- ggplot(coef_summary, aes(x = median, y = term, color = model, group = model, shape = sig)) +
 geom_point(position = position_dodge(width = 0.85), size = 3, stroke = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, linetype = sig),
                 position = position_dodge(width = 0.85),
                 height = 0.25) +
   scale_color_manual(values = my_cols, 
                      name="Model",
                      breaks = c("Leaf Traits", "Root Traits", "Leaf + Root Traits", "Myco Traits", 
                                "Neighborhood", "Canopy", "Leaf + Root + Myco Traits",
                                "Leaf + Root + Myco Traits + Neigh", "Leaf + Root + Myco Traits + Neigh + Canopy", 
                                "Leaf PC1 x Neigh x Canopy"),
                      labels = c("Leaf Traits", "Root Traits", "Leaf + Root Traits", "Myco Traits", 
                                 "Neighborhood", "Canopy", "Leaf + Root + Myco Traits",
                                 "Leaf + Root + Myco Traits + Neigh", "Leaf + Root + Myco Traits + Neigh + Canopy", 
                                 "Leaf PC1 x Neigh x Canopy")) +
  scale_shape_manual(values = c(19,1),
                     breaks = c("significant", "nonsignificant")) +              
  guides(linetype = "none", shape = "none") +      
  guides(color = guide_legend(nrow = 3, ncol = 4)) +
  scale_linetype_manual(values = c("significant" = "solid", "nonsignificant" = "dashed")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Posterior Distribution of Effect on RGR", y = "Model Term", color = "Model") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 14, colour="black"), 
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.text = element_text(size = 14, colour="black"), 
        legend.title = element_text(size = 14, face = "bold", colour="black")) +
  theme(legend.position = "none", legend.justification = "left") +
  theme(legend.key.spacing.x = unit(1, 'cm'))

p3 




# Interpretation: 

# Showing the posterior estimate of the regression coefficient (β) for each predictor in 
# the model.Essentially the change in RGR associated with a one-unit increase in the predictor, 
# holding other predictors constant.

# 1 PCA unit would be one standard deviation along the trait axis 

# Positive = factor increases growth 
# Negative = factor decreases growth 

# Crosses the zero shows that it had a mix of positive and negative effects on growth 


# Leaf traits became more important in the full model, likely because they relate to competition for light 
# that is partially reflected by both the neighborhood and condition of the canopy. 


# This visualization isn't really fully getting at the interaction effects either. 

############################################## --

################################
# (6) SAVE RESULTS OF INTEREST         
################################

# Save weight plot 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/bayesian_weights_plot.png", 
       plot = weight_plot, width = 8, height = 10, units = "in", dpi = 300)


# Save ELPD plot 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/bayesian_ELPD_plot.png", 
       plot = p1, width = 8, height = 10, units = "in", dpi = 300)


# Save forest plot 

# save without the legend and then add it separately 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/bayesian_forest_plot.png", 
       plot = p3, width = 10, height = 10, units = "in", dpi = 300)


