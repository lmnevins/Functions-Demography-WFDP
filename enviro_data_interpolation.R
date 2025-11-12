# -----------------------------------------------------------------------------#
# Environmental data interpolation with kriging  
# Original Author: L. McKinley Nevins 
# November 12, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     sf v 1.0.19
#                     raster v 3.6.32
#                     gstat v 2.1.4
#                     terra v 1.8.60
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(sf); packageVersion("sf")
library(raster); packageVersion("raster")
library(gstat); packageVersion("gstat")
library(terra); packageVersion("terra")

#################################################################################
#                               Main workflow                                   #
#  Use the prepared spatial data for the WFDP plot and the soil nutrient and    #
#  microclimate data collected across the plot. Interpolate the variation in    #
#  key variables using kriging. Generate values for environmental variables     #
#  and match these spatially with the focal trees and their neighborhoods.      # 
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/"
setwd(wd)

# Load WFDP plot boundary (polygon shapefile)
wfdp_poly <- st_read("./GIS/WFDP_boundary_utm.shp")

# Load dataloggers data of soil temperature and moisture 
dataloggers <- read.csv("./Enviro_Data/dataloggers_utm_mean_data.csv")

# Load soil samples data of soil nutrients, pH, etc. 
soil_samples <- read.csv("./Enviro_Data/all_soils_utm.csv")


# Convert to sf objects (using UTM Zone 10N)
dataloggers_sf <- st_as_sf(dataloggers, coords = c("UTM_X", "UTM_Y"), crs = 32610)
soil_sf        <- st_as_sf(soil_samples, coords = c("UTM_X", "UTM_Y"), crs = 32610)


# Plot all together to check alignment 
check_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black") +
  geom_sf(data = dataloggers_sf, aes(color = "Dataloggers"), size = 2) +
  geom_sf(data = soil_sf, aes(color = "Soil samples"), shape = 17, size = 2) +
  scale_color_manual(values = c("Dataloggers" = "blue", "Soil samples" = "brown")) +
  theme_minimal() +
  labs(title = "WFDP sensors and sampling points", color = "Point type")

check_plot

# Three datalogger gaps for technical difficulties, and two pulled up by wild animals. 




