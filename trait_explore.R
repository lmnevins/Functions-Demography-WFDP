# -----------------------------------------------------------------------------#
# Exploring WFDP tree root and leaf trait data 
# Original Author: L. McKinley Nevins 
# January 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 3.5.1
#                     rstatix v 0.7.2
#                     ggpubr v 0.6.0
#                     fundiversity v 1.1.1
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(ggpubr); packageVersion("ggpubr")
library(fundiversity); packageVersion("fundiversity")
library(summarytools)
require(summarytools)

#################################################################################
#                               Main workflow                                   #
#  Explore the leaf and root trait data. Determine if there are differences     #
#  between species, and maybe try to map some of the variation across the plot? #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
############### 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/"
setwd(wd)

# Load in leaf and root trait datasets 

leaf <- read.csv("WFDP_leaf_traits.csv")

root <- read.csv("WFDP_root_traits.csv")

# combine into one dataset for the trees 
traits <- merge(leaf, root, by = 'code')

# subset just variables of interest, excluding petiole data since it is absent for
# needle-leaf species 

traits <- select(traits, code, species = spp.x, SLA_leaf, LDMC_leaf, LMA_leaf, 
                 leaf_pct_N, leaf_pct_C, specific_root_length, specific_root_area,
                 avg_root_dia, root_pct_N, root_pct_C)

# bit more formatting to make it suitable for the fundiversity package 
# need species as the rownames and columns to only be the traits 

# remove code column 
traits_tree <- subset(traits, select = -(code))

#need to get average traits for each species 
traits_tree <- traits_tree %>% group_by(species) %>%
  summarise(across(everything(), mean))

# move species into the row names 
traits_tree <- traits_tree %>%
  column_to_rownames(var = "species")

#check the spread of the raw trait data 
summary(traits_tree)

# traits do not have the same sample sizes, and should also probably be standardized so they 
# are comparable 

#scaling is one option 
traits_tree_sc <- scale(traits_tree)
summary(traits_tree_sc)


#could also scale all traits between zero and 1 

min_values <- as.numeric(lapply(as.data.frame(traits_tree), min))
max_values <- as.numeric(lapply(as.data.frame(traits_tree), max))

traits_tree_minmax <- apply(traits_tree, 1, function(x) {
  (x - min_values)/(max_values - min_values)
})
traits_tree_minmax <- t(traits_tree_minmax)
summary(traits_tree_minmax)


#################################################################################

###############
# (2) FUNCTIONAL ANALYSES
###############


## Functional analyses using the fundiversity package 
# following this vignette: https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html
###############

# not using any species x site matrix, as I'm not making comparisons between sites (all data are 
# from WFDP). If I want to in the future I can make comparisons across some category of 
# environmental variable and perhaps use that as the 'site' 


## Richness

fd_fric(traits_EM)





