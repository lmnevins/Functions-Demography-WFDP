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
growth <- read.csv("stems_for_demo_data_WFDP_20250206.csv")






