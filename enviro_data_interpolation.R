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
#                     car v 3.1.3
#                     FNN v 1.1.4.1
#                     spdep v 1.4.1
#                     ade4 v 1.7.23
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
library(gstat); packageVersion("gstat")
library(car); packageVersion("car")
library(FNN); packageVersion("FNN")
library(spdep); packageVersion("spdep")
library(ade4); packageVersion("ade4")
library(cowplot); packageVersion("cowplot")

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
wfdp_poly <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_boundary_UTM.shp")

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

## Can skip these first two sections and go to section 4 if just interested in the plots of 
# environmental variation across WFDP, and the tests for spatial autocorrelation. 


#################################################################################

######################## --
# (2) DATA EXPLORATION
######################## --

# Explore the key soil variables of interest to check normality and assumptions 
# before kriging

### Soil Microclimate Data ### ---

## Soil volumetric water content - Mean

# 'VWC_moisture_mean'

hist(dataloggers_sf$VWC_moisture_mean)
qqnorm(dataloggers_sf$VWC_moisture_mean)
qqline(dataloggers_sf$VWC_moisture_mean)
shapiro.test(dataloggers_sf$VWC_moisture_mean) 
# p-value = 0.206
# as p-value > 0.05, we can assume normal distribution

#NORMAL 

## Soil Temperature - Mean

# ' TMS_T1_mean' 

hist(dataloggers_sf$TMS_T1_mean)
qqnorm(dataloggers_sf$TMS_T1_mean)
qqline(dataloggers_sf$TMS_T1_mean)
shapiro.test(dataloggers_sf$TMS_T1_mean) 
# p = 0.020

#NOT NORMAL - a few outliers 


## Surface Temperature - Mean

# ' TMS_T2_mean' 

hist(dataloggers_sf$TMS_T2_mean)
qqnorm(dataloggers_sf$TMS_T2_mean)
qqline(dataloggers_sf$TMS_T2_mean)
shapiro.test(dataloggers_sf$TMS_T2_mean) 
# p = 0.793

#NORMAL


## Air Temperature - Mean

# ' TMS_T3_mean' 
hist(dataloggers_sf$TMS_T3_mean)
qqnorm(dataloggers_sf$TMS_T3_mean)
qqline(dataloggers_sf$TMS_T3_mean)
shapiro.test(dataloggers_sf$TMS_T3_mean) 
# p = 1167

#NORMAL

### Soil Nutrient and Properties Data ### ---

## Soil Organic Matter - 10 cm 

# 'OM_10cm'

hist(soil_sf$OM_10cm)
qqnorm(soil_sf$OM_10cm)
qqline(soil_sf$OM_10cm)
shapiro.test(soil_sf$OM_10cm) 
# p < 0.001

##VERY NOT NORMAL 


## Soil Organic Matter - 20 cm 

# 'OM_20cm'
hist(soil_sf$OM_20cm)
qqnorm(soil_sf$OM_20cm)
qqline(soil_sf$OM_20cm)
shapiro.test(soil_sf$OM_20cm) 
# p < 0.001

##VERY NOT NORMAL 


## Soil pH - 10 cm 

# 'pH_10cm'
hist(soil_sf$pH_10cm)
qqnorm(soil_sf$pH_10cm)
qqline(soil_sf$pH_10cm)
shapiro.test(soil_sf$pH_10cm)
# p < 0.001

##VERY NOT NORMAL 


## Soil pH - 20 cm 

# 'pH_20cm'
hist(soil_sf$pH_20cm)
qqnorm(soil_sf$pH_20cm)
qqline(soil_sf$pH_20cm)
shapiro.test(soil_sf$pH_20cm) 
# p < 0.001

##VERY NOT NORMAL 


## Soil Percent Nitrogen - 10 cm 

# 'PctN_10cm'
hist(soil_sf$PctN_10cm)
qqnorm(soil_sf$PctN_10cm)
qqline(soil_sf$PctN_10cm)
shapiro.test(soil_sf$PctN_10cm) 
# p < 0.001

##VERY NOT NORMAL 


## Soil Percent Nitrogen - 20 cm 

# 'PctN_20cm'
hist(soil_sf$PctN_20cm)
qqnorm(soil_sf$PctN_20cm)
qqline(soil_sf$PctN_20cm)
shapiro.test(soil_sf$PctN_20cm) 
# p < 0.001

##VERY NOT NORMAL 

## Soil Percent Carbon - 10 cm 

# 'PctC_10cm'
hist(soil_sf$PctC_10cm)
qqnorm(soil_sf$PctC_10cm)
qqline(soil_sf$PctC_10cm)
shapiro.test(soil_sf$PctC_10cm) 
# p < 0.001

##VERY NOT NORMAL 

## Soil Percent Carbon - 20 cm 

# 'PctC_20cm'
hist(soil_sf$PctC_20cm)
qqnorm(soil_sf$PctC_20cm)
qqline(soil_sf$PctC_20cm)
shapiro.test(soil_sf$PctC_20cm) 
# p < 0.001

##VERY NOT NORMAL 


### All of the soil chemical and physical data is not normal, and 
# has high outliers 



#################################################################################

########################### --
# (3) DATA TRANSFORMATION
########################### --

# powerTransform function of the car package 
# "Uses the maximum likelihood-like approach of Box and Cox (1964) to select a transformation of a 
# univariate or multivariate response for normality, linearity and/or constant variance. Available 
# families of transformations are the default Box-Cox power family and two additioal families that are 
# modifications of the Box-Cox family that allow for (a few) negative responses. The summary method 
# automatically computes two or three likelihood ratio type tests concerning the transformation powers."


# create a vector of variables I want to transform

datalogger_vars <- c("TMS_moist_mean", "TMS_T1_mean", "TMS_T2_mean", "TMS_T3_mean")
# None of these are negative 


soil_vars <- c("OM_10cm", "OM_20cm", "pH_10cm", "pH_20cm", "PCtN_10cm", "PCtN_20cm", 
               "PCtC_10cm", "PCtC_20cm")
# None of these are negative 


# Transform datalogger variables 
datalogger_trans <- dataloggers %>%
  mutate(
    TMS_moist_mean_bc    = bcPower(TMS_moist_mean, powerTransform(TMS_moist_mean)$lambda),
    TMS_T1_mean_bc    = bcPower(TMS_T1_mean, powerTransform(TMS_T1_mean)$lambda),
    TMS_T2_mean_bc    = bcPower(TMS_T2_mean, powerTransform(TMS_T2_mean)$lambda), 
    TMS_T3_mean_bc    = bcPower(TMS_T3_mean, powerTransform(TMS_T3_mean)$lambda))



soil_trans <- soil_samples %>%
  mutate(
    OM_10cm_bc    = bcPower(OM_10cm, powerTransform(OM_10cm)$lambda),
    OM_20cm_bc    = bcPower(OM_20cm, powerTransform(OM_20cm)$lambda),
    pH_10cm_bc    = bcPower(pH_10cm, powerTransform(pH_10cm)$lambda),
    pH_20cm_bc    = bcPower(pH_20cm, powerTransform(pH_20cm)$lambda),
    PctN_10cm_bc    = bcPower(PctN_10cm, powerTransform(PctN_10cm)$lambda),
    PctN_20cm_bc    = bcPower(PctN_20cm, powerTransform(PctN_20cm)$lambda),
    PctC_10cm_bc    = bcPower(PctC_10cm, powerTransform(PctC_10cm)$lambda),
    PctC_20cm_bc    = bcPower(PctC_20cm, powerTransform(PctC_20cm)$lambda))

# Check if transformations worked 

# Dataloggers

## Soil Moisture - Mean

# 'TMS_moist_mean'

hist(dataloggers_sf$TMS_moist_mean, main = "Original")
hist(datalogger_trans$TMS_moist_mean_bc, main = "Box-Cox")


qqPlot(datalogger_trans$TMS_moist_mean_bc)
shapiro.test(datalogger_trans$TMS_moist_mean_bc) 
# p-value = 0.426

# Was already NORMAL, but now even more normal, and this will stay consistent 
# with the other variables 


## Soil Temperature - Mean

# ' TMS_T1_mean' 

hist(dataloggers_sf$TMS_T1_mean, main = "Original")
hist(datalogger_trans$TMS_T1_mean_bc, main = "Box-Cox")


qqPlot(datalogger_trans$TMS_T1_mean_bc)
shapiro.test(datalogger_trans$TMS_T1_mean_bc) 
# p = 0.081

#NORMAL


## Soil Temperature - Mean

# ' TMS_T2_mean' 

hist(dataloggers_sf$TMS_T2_mean, main = "Original")
hist(datalogger_trans$TMS_T2_mean_bc, main = "Box-Cox")


qqPlot(datalogger_trans$TMS_T2_mean_bc)
shapiro.test(datalogger_trans$TMS_T2_mean_bc) 
# p = 0.905

# Was already NORMAL, but now even more normal, and this will stay consistent 
# with the other variables 


## Air Temperature - Mean

# ' TMS_T3_mean' 
hist(dataloggers_sf$TMS_T3_mean, main = "Original")
hist(datalogger_trans$TMS_T3_mean_bc, main = "Box-Cox")


qqPlot(datalogger_trans$TMS_T3_mean_bc)
shapiro.test(datalogger_trans$TMS_T3_mean_bc) 
# p = 0.128

# Was already NORMAL, but now even more normal, and this will stay consistent 
# with the other variables 



### Soil Nutrient and Properties Data ### ---

## Soil Organic Matter - 10 cm 

# 'OM_10cm'
hist(soil_samples$OM_10cm, main = "Original")
hist(soil_trans$OM_10cm_bc, main = "Box-Cox")


qqPlot(soil_trans$OM_10cm_bc)
shapiro.test(soil_trans$OM_10cm_bc) 
# p = 0.09

# NOW NORMAL 


## Soil Organic Matter - 20 cm 

# 'OM_20cm'
hist(soil_samples$OM_20cm, main = "Original")
hist(soil_trans$OM_20cm_bc, main = "Box-Cox")


qqPlot(soil_trans$OM_20cm_bc)
shapiro.test(soil_trans$OM_20cm_bc) 
# p = 0.673

# NOW NORMAL 


## Soil pH - 10 cm 

# 'pH_10cm'
hist(soil_samples$pH_10cm, main = "Original")
hist(soil_trans$pH_10cm_bc, main = "Box-Cox")


qqPlot(soil_trans$pH_10cm_bc)
shapiro.test(soil_trans$pH_10cm_bc) 
# p = 0.304

# NOW NORMAL 


## Soil pH - 20 cm 

# 'pH_20cm'
hist(soil_samples$pH_20cm, main = "Original")
hist(soil_trans$pH_20cm_bc, main = "Box-Cox")


qqPlot(soil_trans$pH_20cm_bc)
shapiro.test(soil_trans$pH_20cm_bc) 
# p = 0.03

# Still not normal, but better so can give it a try with kriging 


## Soil Percent Nitrogen - 10 cm 

# 'PctN_10cm'
hist(soil_samples$PctN_10cm, main = "Original")
hist(soil_trans$PctN_10cm_bc, main = "Box-Cox")


qqPlot(soil_trans$PctN_10cm_bc)
shapiro.test(soil_trans$PctN_10cm_bc) 
# p = 0.165

# NOW NORMAL 


## Soil Percent Nitrogen - 20 cm 

# 'PctN_20cm'
hist(soil_samples$PctN_20cm, main = "Original")
hist(soil_trans$PctN_20cm_bc, main = "Box-Cox")


qqPlot(soil_trans$PctN_20cm_bc)
shapiro.test(soil_trans$PctN_20cm_bc) 
# p = 0.849

# NOW NORMAL 

## Soil Percent Carbon - 10 cm 

# 'PctC_10cm'
hist(soil_samples$PctC_10cm, main = "Original")
hist(soil_trans$PctC_10cm_bc, main = "Box-Cox")


qqPlot(soil_trans$PctC_10cm_bc)
shapiro.test(soil_trans$PctC_10cm_bc) 
# p = 0.269

# NOW NORMAL 

## Soil Percent Carbon - 20 cm 

# 'PctC_20cm'
hist(soil_samples$PctC_20cm, main = "Original")
hist(soil_trans$PctC_20cm_bc, main = "Box-Cox")


qqPlot(soil_trans$PctC_20cm_bc)
shapiro.test(soil_trans$PctC_20cm_bc) 
# p = 0.742

# NOW NORMAL 


## Save transformed datasets for QGIS

write.csv(datalogger_trans, "./Enviro_Data/dataloggers_transformed_data.csv")

write.csv(soil_trans, "./Enviro_Data/soil_transformed_data.csv")


# make the transformed datasets sf objects again 
datalogger_trans_sf <- st_as_sf(
  datalogger_trans,
  coords = c("UTM_X", "UTM_Y"),
  crs = 32610
)


soil_trans_sf <- st_as_sf(
  soil_trans,
  coords = c("UTM_X", "UTM_Y"),
  crs = 32610
)


####################################################################

######################################## --
# (4) TEST FOR SPATIAL AUTOCORRELATION 
######################################## --

# Need to assess if there is spatial structuring to the data, as that influences 
# the type of interpolation technique that can be used, if one can be used at all. 


# For each variable of interest, assess variation in the raw data across the plot, 
# then perform a Moran's I test, and visualize with a mantel correlogram. 


# Interested in the mean soil temperature and moisture values across the plot, as well 
# as the range, as this could reflect areas that are more variable for the trees and 
# fungal communities. 

# Excluding surface temp because I'm not really that interested in it 

# TMS_T1_mean - soil 
# TMS_T1_range
# TMS_T3_mean - air
# TMS_T3_range
# VWC_moisture_mean
# VWC_moisture_range


# convert VWC values to percentages for plotting 

dataloggers <- dataloggers %>%
  dplyr::mutate(
    VWC_mean_percent = VWC_moisture_mean * 100, 
    VWC_range_percent = VWC_moisture_range * 100)

dataloggers_sf <- dataloggers_sf %>%
  dplyr::mutate(
    VWC_mean_percent = VWC_moisture_mean * 100, 
    VWC_range_percent = VWC_moisture_range * 100)



# And for soil data I am interested in a subset of variables - can split into chemical and physical for 
# thinking about them, to make plotting dimensions easier 

# Not considering PctC because it is so strongly correlated to soil organic matter 

# pH_10cm
# pH_20cm

# OM_10cm
# OM_20cm 
# PctN_10cm
# PctN_20cm

######################## -- 

## PREP
# get spatial dist for all dataloggers 
logger.dists <- dist(cbind(dataloggers$UTM_X, dataloggers$UTM_Y))


# get spatial dist for all soil sites 
soil.dists <- dist(cbind(soil_samples$UTM_X, soil_samples$UTM_Y))


# Function to plot the results of the mantel correlogram and make it customizable with ggplot
ggplot_mantel <- function(mant, fill = "gray50") {
  df <- data.frame(x = mant$plot$hist$mids,
                   y = mant$plot$hist$counts)
  ggplot(df, aes(x, y)) + 
    geom_col(orientation = "x", 
             width = diff(mant$plot$hist$breaks)[1],
             fill = fill, color = "gray30") +
    labs(x = mant$plot$hist$xname, y = "Frequency") +
    scale_x_continuous(limits = mant$plot$xlim) +
    geom_segment(aes(x = mant$obs, xend = mant$obs, y = 0,
                     yend = 0.75 * max(y))) +
    geom_point(aes(x = mant$obs, y = 0.75 * max(y)), size = 5,
               shape = 18) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 11, colour="black"),
      axis.text.y = element_text(size = 11, colour="black"),
      axis.title.y = element_text(size = 12, colour="black"),
      axis.title.x = element_text(size = 12, colour="black"))
}



## Soil datalogger variables ## ---

# 1. Mean Soil Temp

## Exploring soil temperature 
T1_mean_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = dataloggers_sf, aes(color = TMS_T1_mean), size = 3) +
  scale_color_gradient(low = "tan", high = "red4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Mean Soil Temp (°C)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

T1_mean_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(dataloggers[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

dataloggers$nn_dist <- nn$nn.dist[,1]
dataloggers$nn_diff <- abs(
  dataloggers$TMS_T1_mean -
    dataloggers$TMS_T1_mean[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(dataloggers, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in mean soil temp"
  )

nn_plot 

# No pattern whatsoever

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(dataloggers$TMS_T1_mean, lw)

# NO SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = 0.64327, p-value = 0.26
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#         0.035969643      -0.027777778       0.009820503 


# get euclidian distance for temp data for loggers 
env_dist_T1_mean <- dist(dataloggers$TMS_T1_mean)


# Run mantel correlogram 
mc_T1_mean <- ade4::mantel.randtest(env_dist_T1_mean, logger.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_T1_mean, m2 = logger.dists, 
#                             nrepet = 9999)
# 
# Observation: -0.03244027 
# 
# Based on 9999 replicates
# Simulated p-value: 0.6624 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# -0.483570334  0.001657415  0.004971993 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_T1_mean <- ggplot_mantel(mc_T1_mean)

mantel_T1_mean

# Throws some errors, but this plots the same as if you just did plot() so it's okay 

# Vertical line reflects observed value, and the frequency is the distribution of the 
# simulated values. If the distribution is skewed to one side or the other of the vertical 
# line, then this reflects a significant trend in spatial autocorrelation, either positively 
# or negatively 



# 2. Soil Temp Range

## Exploring the range in soil temperature 
T1_range_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = dataloggers_sf, aes(color = TMS_T1_range), size = 3) +
  scale_color_gradient(low = "tan", high = "red4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Range of Mean \nSoil Temp (°C)") + #using manual linebreak to wrap text 
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

T1_range_plot

# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(dataloggers[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

dataloggers$nn_dist <- nn$nn.dist[,1]
dataloggers$nn_diff <- abs(
  dataloggers$TMS_T1_range -
    dataloggers$TMS_T1_range[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(dataloggers, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil temp range"
  )

nn_plot 


# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(dataloggers$TMS_T1_range, lw)

# Weak but significant spatial autocorrelation 


# Moran I statistic standard deviate = 1.9268, p-value = 0.027
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#          0.155328016      -0.027777778       0.009030916 


# get euclidian distance for temp data for loggers 
env_dist_T1_range <- dist(dataloggers$TMS_T1_range)


# Run mantel correlogram 
mc_T1_range <- ade4::mantel.randtest(env_dist_T1_range, logger.dists, nrepet = 9999)

# Highlight just the model portion to get the output

# 
# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_T1_range, m2 = logger.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.0556388 
# 
# Based on 9999 replicates
# Simulated p-value: 0.2423 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 0.7633177978 0.0005204531 0.0052141267 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 

mantel_T1_range <- ggplot_mantel(mc_T1_range)

mantel_T1_range

# Vertical line reflects observed value, and the frequency is the distribution of the 
# simulated values. If the distribution is skewed to one side or the other of the vertical 
# line, then this reflects a significant trend in spatial autocorrelation, either positively 
# or negatively 



# 3. Mean Air Temp

## Exploring air temperature 
T3_mean_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = dataloggers_sf, aes(color = TMS_T3_mean), size = 3) +
  scale_color_gradient(low = "tan", high = "lightblue3") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Mean Air Temp (°C)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

T3_mean_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(dataloggers[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

dataloggers$nn_dist <- nn$nn.dist[,1]
dataloggers$nn_diff <- abs(
  dataloggers$TMS_T3_mean -
    dataloggers$TMS_T3_mean[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(dataloggers, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in air temp"
  )

nn_plot 

# No pattern whatsoever

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(dataloggers$TMS_T3_mean, lw)

# SIGNIFICANT spatial autocorrelation, but relatively weak 


# Moran I statistic standard deviate = 2.0466, p-value = 0.02035
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.18188583       -0.02777778        0.01049459 


# get euclidian distance for temp data for loggers 
env_dist_T3_mean <- dist(dataloggers$TMS_T3_mean)


# Run mantel correlogram 
mc_T3_mean <- ade4::mantel.randtest(env_dist_T3_mean, logger.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_T3_mean, m2 = logger.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.07837364 
# 
# Based on 9999 replicates
# Simulated p-value: 0.1269 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 1.1771697154 0.0002069525 0.0044092550 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_T3_mean <- ggplot_mantel(mc_T3_mean)

mantel_T3_mean

# Throws some errors, but this plots the same as if you just did plot() so it's okay 



# 4. Air Temp Range

## Exploring the range in air temperature 
T3_range_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = dataloggers_sf, aes(color = TMS_T3_range), size = 3) +
  scale_color_gradient(low = "tan", high = "lightblue3") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Range of Mean \nAir Temp (°C)") + #using manual linebreak to wrap text 
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

T3_range_plot

# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(dataloggers[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

dataloggers$nn_dist <- nn$nn.dist[,1]
dataloggers$nn_diff <- abs(
  dataloggers$TMS_T3_range -
    dataloggers$TMS_T3_range[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(dataloggers, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in air temp range"
  )

nn_plot 

# really no clear trend 


# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(dataloggers$TMS_T3_range, lw)

# Weak but significant spatial autocorrelation 


# Moran I statistic standard deviate = 2.0395, p-value = 0.0207
# alternative hypothesis: greater
# sample estimates:
# Moran I statistic       Expectation          Variance 
#        0.18130695       -0.02777778        0.01050935 


# get euclidian distance for temp data for loggers 
env_dist_T3_range <- dist(dataloggers$TMS_T3_range)


# Run mantel correlogram 
mc_T3_range <- ade4::mantel.randtest(env_dist_T3_range, logger.dists, nrepet = 9999)

# Highlight just the model portion to get the output

# 
# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_T3_range, m2 = logger.dists, 
#                             nrepet = 9999)
# 
# Observation: -0.06749089 
# 
# Based on 9999 replicates
# Simulated p-value: 0.8499 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# -1.050350237 -0.000936913  0.004014946 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 

mantel_T3_range <- ggplot_mantel(mc_T3_range)

mantel_T3_range



# 5. Mean soil volumetric water content 

## Exploring soil moisture 
VWC_mean_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = dataloggers_sf, aes(color = VWC_mean_percent), size = 3) +
  scale_color_gradient(low = "tan", high = "blue4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Mean Soil Volumetric\n Water Content(%)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

VWC_mean_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(dataloggers[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

dataloggers$nn_dist <- nn$nn.dist[,1]
dataloggers$nn_diff <- abs(
  dataloggers$VWC_mean_percent -
    dataloggers$VWC_mean_percent[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(dataloggers, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in mean VWC"
  )

nn_plot 

# Maybe a bit of a trend 

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(dataloggers$VWC_mean_percent, lw)

# NO SIGNIFICANT spatial autocorrelation


# Moran I statistic standard deviate = -1.6783, p-value = 0.9534
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.20307871       -0.02777778        0.01090966 


# get euclidian distance for temp data for loggers 
env_dist_VWC_mean <- dist(dataloggers$VWC_mean_percent)


# Run mantel correlogram 
mc_VWC_mean <- ade4::mantel.randtest(env_dist_VWC_mean, logger.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_VWC_mean, m2 = logger.dists, 
#                             nrepet = 9999)
# 
# Observation: -0.01925831 
# 
# Based on 9999 replicates
# Simulated p-value: 0.618 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# -0.3546149938  0.0007874281  0.0031954350 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_VWC_mean <- ggplot_mantel(mc_VWC_mean)

mantel_VWC_mean

# Throws some errors, but this plots the same as if you just did plot() so it's okay 



# 6. Soil VWC Range

## Exploring the range in soil VWC
VWC_range_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = dataloggers_sf, aes(color = VWC_range_percent), size = 3) +
  scale_color_gradient(low = "tan", high = "blue4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Range of Mean Soil\n Volumetric Water\n Content(%)") + #using manual linebreak to wrap text 
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

VWC_range_plot

# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(dataloggers[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

dataloggers$nn_dist <- nn$nn.dist[,1]
dataloggers$nn_diff <- abs(
  dataloggers$VWC_range_percent -
    dataloggers$VWC_range_percent[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(dataloggers, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil VWC range"
  )

nn_plot 

# really no clear trend 


# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(dataloggers$VWC_range_percent, lw)

# NO significant spatial autocorrelation 


# Moran I statistic standard deviate = 0.31949, p-value = 0.3747
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#         0.004910257      -0.027777778       0.010468059 


# get euclidian distance for temp data for loggers 
env_dist_VWC_range <- dist(dataloggers$VWC_range_percent)


# Run mantel correlogram 
mc_VWC_range <- ade4::mantel.randtest(env_dist_VWC_range, logger.dists, nrepet = 9999)

# Highlight just the model portion to get the output

# 
# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_VWC_range, m2 = logger.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.00404931 
# 
# Based on 9999 replicates
# Simulated p-value: 0.4611 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 0.0494212817 0.0008518961 0.0041857164 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 

mantel_VWC_range <- ggplot_mantel(mc_VWC_range)

mantel_VWC_range





## Soil sample variables ## ---

# 1. Soil pH  - 10 cm 

## Exploring soil pH 10 cm 
pH_10cm_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = soil_sf, aes(color = pH_10cm), size = 3) +
  scale_color_gradient(low = "tan", high = "purple4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Soil pH (0-10 cm)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

pH_10cm_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(soil_samples[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

soil_samples$nn_dist <- nn$nn.dist[,1]
soil_samples$nn_diff <- abs(
  soil_samples$pH_10cm -
    soil_samples$pH_10cm[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(soil_samples, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil pH (0-10 cm)"
  )

nn_plot 

# No pattern whatsoever

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(soil_samples$pH_10cm, lw)

# NO SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = 0.95492, p-value = 0.1698
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#        0.061590087      -0.025641026       0.008344714 


# get euclidian distance for temp data for loggers 
env_dist_pH_10cm <- dist(soil_samples$pH_10cm)


# Run mantel correlogram 
mc_pH_10cm <- ade4::mantel.randtest(env_dist_pH_10cm, soil.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_pH_10cm, m2 = soil.dists, 
#                             nrepet = 9999)
# 
# Observation: -0.02537348 
# 
# Based on 9999 replicates
# Simulated p-value: 0.5951 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# -0.3452528244 -0.0005691883  0.0051615384 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_pH_10cm <- ggplot_mantel(mc_pH_10cm)

mantel_pH_10cm


# 2. Soil pH  - 20 cm 

## Exploring soil pH 20 cm 
pH_20cm_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = soil_sf, aes(color = pH_20cm), size = 3) +
  scale_color_gradient(low = "tan", high = "purple4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Soil pH (10-20 cm)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

pH_20cm_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(soil_samples[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

soil_samples$nn_dist <- nn$nn.dist[,1]
soil_samples$nn_diff <- abs(
  soil_samples$pH_20cm -
    soil_samples$pH_20cm[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(soil_samples, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil pH (10-20 cm)"
  )

nn_plot 

# No pattern whatsoever

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(soil_samples$pH_20cm, lw)

# SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = 2.432, p-value = 0.007509
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
#    0.198558810      -0.025641026       0.008498749 


# get euclidian distance for temp data for loggers 
env_dist_pH_20cm <- dist(soil_samples$pH_20cm)


# Run mantel correlogram 
mc_pH_20cm <- ade4::mantel.randtest(env_dist_pH_20cm, soil.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_pH_20cm, m2 = soil.dists, 
#                             nrepet = 9999)
# 
# Observation: -0.03286255 
# 
# Based on 9999 replicates
# Simulated p-value: 0.638 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# -0.4484984857 -0.0004711068  0.0052160126 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_pH_20cm <- ggplot_mantel(mc_pH_20cm)

mantel_pH_20cm



# 3. Soil OM  - 10 cm 

## Exploring soil EC 10 cm 
OM_10cm_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = soil_sf, aes(color = OM_10cm), size = 3) +
  scale_color_gradient(low = "tan", high = "green4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Soil Organic Matter (%)\n (0-10 cm)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

OM_10cm_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(soil_samples[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

soil_samples$nn_dist <- nn$nn.dist[,1]
soil_samples$nn_diff <- abs(
  soil_samples$OM_10cm -
    soil_samples$OM_10cm[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(soil_samples, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil Organic Matter (%) (0-10 cm)"
  )

nn_plot 

# Maybe a tiny bit of a trend?

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(soil_samples$OM_10cm, lw)

# NO SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = -0.68868, p-value = 0.7545
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.09497620       -0.02564103        0.01013604 


# get euclidian distance for temp data for loggers 
env_dist_OM_10cm <- dist(soil_samples$OM_10cm)


# Run mantel correlogram 
mc_OM_10cm <- ade4::mantel.randtest(env_dist_OM_10cm, soil.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_OM_10cm, m2 = soil.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.01511666 
# 
# Based on 9999 replicates
# Simulated p-value: 0.3842 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 0.2379176068 0.0008608812 0.0035902857 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_OM_10cm <- ggplot_mantel(mc_OM_10cm)

mantel_OM_10cm


# 4. Soil OM  - 20 cm 

## Exploring soil OM 20 cm 
OM_20cm_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = soil_sf, aes(color = OM_20cm), size = 3) +
  scale_color_gradient(low = "tan", high = "green4") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Soil Organic Matter (%)\n (10-20 cm)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

OM_20cm_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(soil_samples[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

soil_samples$nn_dist <- nn$nn.dist[,1]
soil_samples$nn_diff <- abs(
  soil_samples$OM_20cm -
    soil_samples$OM_20cm[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(soil_samples, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil Organic Matter (%) (10-20 cm)"
  )

nn_plot 

# No pattern whatsoever

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(soil_samples$OM_20cm, lw)

# SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = 0.96106, p-value = 0.1683
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.061097441      -0.025641026       0.008145616 


# get euclidian distance for temp data for loggers 
env_dist_OM_20cm <- dist(soil_samples$OM_20cm)


# Run mantel correlogram 
mc_OM_20cm <- ade4::mantel.randtest(env_dist_OM_20cm, soil.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_OM_20cm, m2 = soil.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.1878742 
# 
# Based on 9999 replicates
# Simulated p-value: 0.0097 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# 2.5414160578 -0.0005974297  0.0054997154 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_OM_20cm <- ggplot_mantel(mc_OM_20cm)

mantel_OM_20cm


# 3. Soil Pct N  - 10 cm 

## Exploring soil percent N 10 cm 
N_10cm_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = soil_sf, aes(color = PctN_10cm), size = 3) +
  scale_color_gradient(low = "tan", high = "brown3") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Soil N Content (%)\n (0-10 cm)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

N_10cm_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(soil_samples[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

soil_samples$nn_dist <- nn$nn.dist[,1]
soil_samples$nn_diff <- abs(
  soil_samples$PctN_10cm -
    soil_samples$PctN_10cm[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(soil_samples, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil N content (%) (0-10 cm)"
  )

nn_plot 

# Maybe a tiny bit of a trend?

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(soil_samples$PctN_10cm, lw)

# NO SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = -0.14233, p-value = 0.5566
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.04002031       -0.02564103        0.01020611 


# get euclidian distance for temp data for loggers 
env_dist_N_10cm <- dist(soil_samples$PctN_10cm)


# Run mantel correlogram 
mc_N_10cm <- ade4::mantel.randtest(env_dist_N_10cm, soil.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_N_10cm, m2 = soil.dists, 
#                             nrepet = 9999)
# 
# Observation: -0.01251803 
# 
# Based on 9999 replicates
# Simulated p-value: 0.5575 
# Alternative hypothesis: greater 
# 
# Std.Obs   Expectation      Variance 
# -0.2107945853 -0.0001926254  0.0034188756 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_N_10cm <- ggplot_mantel(mc_N_10cm)

mantel_N_10cm


# 6. Soil N content  - 20 cm 

## Exploring soil N content 20 cm 
N_20cm_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = soil_sf, aes(color = PctN_20cm), size = 3) +
  scale_color_gradient(low = "tan", high = "brown3") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal() +
  labs(color = "Soil N Content (%)\n (10-20 cm)") +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=11),
    axis.text    = element_text(size=11, colour="black"),
    axis.title   = element_text(size=12, colour="black")
  )

N_20cm_plot


# Perform nearest-neighbor comparison of the data 
coords <- as.matrix(soil_samples[, c("UTM_X", "UTM_Y")])
nn <- get.knn(coords, k = 2)

soil_samples$nn_dist <- nn$nn.dist[,1]
soil_samples$nn_diff <- abs(
  soil_samples$PctN_20cm -
    soil_samples$PctN_20cm[nn$nn.index[,1]]
)


nn_plot <- ggplot2::ggplot(soil_samples, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest logger (m)",
    y = "Difference in soil N Content (%) (10-20 cm)"
  )

nn_plot 

# No pattern whatsoever

# Moran's I test for spatial structure 
coords_nb <- knearneigh(coords, k = 4)
nb <- knn2nb(coords_nb)
lw <- nb2listw(nb)

moran.test(soil_samples$PctN_20cm, lw)

# SIGNIFICANT spatial autocorrelation 


# Moran I statistic standard deviate = 1.3994, p-value = 0.08084
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.092657946      -0.025641026       0.007145957 


# get euclidian distance for temp data for loggers 
env_dist_N_20cm <- dist(soil_samples$PctN_20cm)


# Run mantel correlogram 
mc_N_20cm <- ade4::mantel.randtest(env_dist_N_20cm, soil.dists, nrepet = 9999)

# Highlight just the model portion to get the output


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_N_20cm, m2 = soil.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.1363918 
# 
# Based on 9999 replicates
# Simulated p-value: 0.0439 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 1.864637803 -0.001325592  0.005454927 

# Observation is the Mantel statistic, which indicates the strength and direction 
# of the correlation, if one exists 


mantel_N_20cm <- ggplot_mantel(mc_N_20cm)

mantel_N_20cm


############################################## -- 
# (4) GATHER SIGNIFICANT RESULTS OF INTEREST  
############################################## -- 

# Gather up the most interesting plots and organize the plots to show the environmental variation 
# across the plot, and then to summarize the Mantel results for spatial autocorrelation


# Gather summary plots for the dataloggers - the mean data:
# T1_mean_plot, T3_mean_plot, VWC_mean_plot

# And the ranges of the variables to show the amount of variability across the plot: 
# T1_range_plot, T3_range_plot, VWC_range_plot


# For the soil data there are a lot more plots, and they need to be organized so the 10cm and 10cm 
# depths are next to eachother 
# pH_10cm_plot, pH_20cm_plot, OM_10cm_plot, OM_20cm_plot, N_10cm_plot, N_20cm_plot


# Gather mantel results - for the datalogger mean data 
# mantel_T1_mean, mantel_T3_mean, mantel_VWC_mean


# For the range data 
# mantel_T1_range, mantel_T3_range, mantel_VWC_range


# For the soil variables - pH, OM, PctN
# mantel_pH_10cm, mantel_pH_20cm, mantel_OM_10cm, mantel_OM_20cm, mantel_N_10cm, mantel_N_20cm


# datalogger mean data plots 
datalogger_mean_plots <- plot_grid(T1_mean_plot, T3_mean_plot, VWC_mean_plot,
                       ncol = 1, nrow = 3, labels = c('(a)', '(b)', '(c)'))

datalogger_mean_plots


# datalogger data range plots 
datalogger_range_plots <- plot_grid(T1_range_plot, T3_range_plot, VWC_range_plot,
                                   ncol = 1, nrow = 3, labels = c('(a)', '(b)', '(c)'))

datalogger_range_plots


# datalogger mean and range Mantel's correlation plots 
datalogger_mantel_plots <- plot_grid(mantel_T1_mean, mantel_T1_range, mantel_T3_mean, mantel_T3_range,
                                          mantel_VWC_mean, mantel_VWC_range,
                                   ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'), 
                                   align = "hv", hjust = -0.1)

datalogger_mantel_plots


# soil data plots - do 10 and 20 cm depths separately 

# 10 cm 
soil_10cm_plots <- plot_grid(pH_10cm_plot, OM_10cm_plot, N_10cm_plot, 
                        ncol = 1, nrow = 3, labels = c('(a)', '(b)', '(c)'), 
                        align = "hv", hjust = -0.1)

soil_10cm_plots


soil_20cm_plots <- plot_grid(pH_20cm_plot, OM_20cm_plot, N_20cm_plot, 
                             ncol = 1, nrow = 3, labels = c('(a)', '(b)', '(c)'), 
                             align = "hv", hjust = -0.1)

soil_20cm_plots


# soil Mantel's correlation plots 
soil_mantel_plots <- plot_grid(mantel_pH_10cm, mantel_pH_20cm, mantel_OM_10cm, mantel_OM_20cm, mantel_N_10cm, mantel_N_20cm,
                                     ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'), 
                                     align = "hv", hjust = -0.1)

soil_mantel_plots

