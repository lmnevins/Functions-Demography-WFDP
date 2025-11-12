# -----------------------------------------------------------------------------#
# Soil sample data prep
# Original Author: L. McKinley Nevins 
# November 10, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 3.5.1
#                     gstat v 2.1.4
#                     terra v 1.8.60
#                     sf v 1.0.19
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(myClim); packageVersion("myClim")
library(gstat); packageVersion("gstat")
library(terra); packageVersion("terra")
library(sf); packageVersion("sf")

#################################################################################
#                               Main workflow                                   #
#  Process the soil sample data from WFDP to get the coordinates in the right   #
#  format for spatial mapping and kriging interpolation in QGIS.               #
#                                                                               #
#################################################################################

############### -- 
# (1) DATA PREP
############### -- 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/"
setwd(wd)


# Load in all soil data, which has data from the two sampling depths: 0-10 cm and 10-20 cm organized 
# for each of the 40 sampling locations 
soil_data <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/all_soil_data.csv")

# Check that data is intitially in EPSG: 4326, which is standard for GPS units 
soil <- st_as_sf(soil_data, coords = c("Lon","Lat"), crs = 4326)

# check
st_crs(soil)
st_bbox(soil)

# In CRS 4236 currently, but I want to reproject to ESG 32610 so I can perform the krieging 
# analyses on the same map system as my other points, and the WFDP boundary


#################################

#################### -- 
# (2) SPATIAL PREP
#################### -- 

soil_utm <- st_transform(soil, 32610)

# check UTM values (should be ~580000, ~5074000)
st_crs(soil_utm)
st_bbox(soil_utm)
head(st_coordinates(soil_utm))


# Export the UTM coordinates to read into QGIS
soil_utm <- soil_utm %>%
  mutate(UTM_X = st_coordinates(.)[,1],
         UTM_Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()


# Save 
write.csv(soil_utm, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/all_soils_utm.csv", row.names = FALSE)

#################################

# -- END -- # 




