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
#                     automap v 1.1.20
#                     stars v 0.6.8
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
library(automap); packageVersion("automap")
library(stars); packageVersion("stars")

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

# Three datalogger gaps: one for technical difficulties, and two pulled up by wild animals. 


#################################################################################

######################## --
# (2) DATA EXPLORATION
######################## --

# Explore the key soil variables of interest to check normality and assumptions 
# before kriging

### Soil Microclimate Data ### ---

## Soil Moisture - Mean

# 'TMS_moist_mean'

hist(dataloggers_sf$TMS_moist_mean)
qqnorm(dataloggers_sf$TMS_moist_mean)
qqline(dataloggers_sf$TMS_moist_mean)
shapiro.test(dataloggers_sf$TMS_moist_mean) 
# p-value = 0.239
# as p-value > 0.05, we can assume normal distribution

#NORMAL 

## Soil Temperature - Mean

# ' TMS_T1_mean' 

hist(dataloggers_sf$TMS_T1_mean)
qqnorm(dataloggers_sf$TMS_T1_mean)
qqline(dataloggers_sf$TMS_T1_mean)
shapiro.test(dataloggers_sf$TMS_T1_mean) 

#NOT NORMAL - a few outliers 


## Air Temperature - Mean

# ' TMS_T3_mean' 
hist(dataloggers_sf$TMS_T3_mean)
qqnorm(dataloggers_sf$TMS_T3_mean)
qqline(dataloggers_sf$TMS_T3_mean)
shapiro.test(dataloggers_sf$TMS_T3_mean) 

#NORMAL

### Soil Nutrient and Properties Data ### ---

## Soil Organic Matter - 10 cm 

# 'OM_10cm'

hist(soil_sf$OM_10cm)
qqnorm(soil_sf$OM_10cm)
qqline(soil_sf$OM_10cm)
shapiro.test(soil_sf$OM_10cm) 

##VERY NOT NORMAL 


## Soil Organic Matter - 20 cm 

# 'OM_20cm'
hist(soil_sf$OM_20cm)
qqnorm(soil_sf$OM_20cm)
qqline(soil_sf$OM_20cm)
shapiro.test(soil_sf$OM_20cm) 

##VERY NOT NORMAL 


## Soil pH - 10 cm 

# 'pH_10cm'
hist(soil_sf$pH_10cm)
qqnorm(soil_sf$pH_10cm)
qqline(soil_sf$pH_10cm)
shapiro.test(soil_sf$pH_10cm) 

##VERY NOT NORMAL 


## Soil pH - 20 cm 

# 'pH_20cm'
hist(soil_sf$pH_20cm)
qqnorm(soil_sf$pH_20cm)
qqline(soil_sf$pH_20cm)
shapiro.test(soil_sf$pH_20cm) 

##VERY NOT NORMAL 


## Soil Cation Exchange Capacity - 10 cm 

# 'EC_10cm'
hist(soil_sf$EC_10cm)
qqnorm(soil_sf$EC_10cm)
qqline(soil_sf$EC_10cm)
shapiro.test(soil_sf$EC_10cm) 

##VERY NOT NORMAL 

## Soil Cation Exchange Capacity - 20 cm 

# 'EC_20cm'
hist(soil_sf$EC_20cm)
qqnorm(soil_sf$EC_20cm)
qqline(soil_sf$EC_20cm)
shapiro.test(soil_sf$EC_20cm) 

##VERY NOT NORMAL 


## Soil Percent Nitrogen - 10 cm 

# 'PctN_10cm'
hist(soil_sf$PctN_10cm)
qqnorm(soil_sf$PctN_10cm)
qqline(soil_sf$PctN_10cm)
shapiro.test(soil_sf$PctN_10cm) 

##VERY NOT NORMAL 


## Soil Percent Nitrogen - 20 cm 

# 'PctN_20cm'
hist(soil_sf$PctN_20cm)
qqnorm(soil_sf$PctN_20cm)
qqline(soil_sf$PctN_20cm)
shapiro.test(soil_sf$PctN_20cm) 

##VERY NOT NORMAL 

## Soil Percent Carbon - 10 cm 

# 'PctC_10cm'
hist(soil_sf$PctC_10cm)
qqnorm(soil_sf$PctC_10cm)
qqline(soil_sf$PctC_10cm)
shapiro.test(soil_sf$PctC_10cm) 

##VERY NOT NORMAL 

## Soil Percent Carbon - 20 cm 

# 'PctC_20cm'
hist(soil_sf$PctC_20cm)
qqnorm(soil_sf$PctC_20cm)
qqline(soil_sf$PctC_20cm)
shapiro.test(soil_sf$PctC_20cm) 

##VERY NOT NORMAL 


### All of the soil chemical and physical data is not normal, and 
# has high outliers 









####################################################################

# Use mean soil moisture as a test variable 

# Create a prediction grid at 10 m resolution within the WFDP map 
cell_res <- 10

grid <- st_make_grid(
  wfdp_poly,
  cellsize = cell_res,
  what = "centers",
  square = TRUE
)

# Keep only grid cells inside polygon
grid <- st_sf(geometry = grid)
grid <- grid[wfdp_poly, ]



# use automap to krige the mean soil moisture variable 

var <- "TMS_moist_max"
form <- as.formula(paste(var, "~ 1"))

kr <- autoKrige(formula = form, input_data = dataloggers_sf,
                new_data = grid, model = c("Sph", "Exp", "Gau"), verbose = TRUE)


kr_stars <- st_as_stars(kr$krige_output)
plot(kr_stars["var1.pred"], breaks = "equal")





# plot a quick variogram to inspect 
vg <- variogram(TMS_moist_max ~ 1, data = dataloggers_sf)
plot(vg)

# fit variogram models for comparison 
vg.fit <- fit.variogram(vg, model=vgm("Sph", nugget = 2500, psill = 75, range = 150))

# Plot the sample values along with the fitted model 
plot(vg, vg.fit)


# load spatial domain to interpolate over
                                                           




