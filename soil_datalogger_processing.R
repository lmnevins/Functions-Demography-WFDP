# -----------------------------------------------------------------------------#
# Soil Datalogger Data Processing
# Original Author: L. McKinley Nevins 
# February 2, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 4.0.1
#                     myClim v 1.5.0
#                     sf v 1.0.23
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(myClim); packageVersion("myClim")
library(sf); packageVersion("sf")

#################################################################################
#                               Main workflow                                   #
#  Process the soil data from all of the dataloggers retrieved from WFDP using  #
#  the myClim package, specifically designed to process microlimate datalogger  #
#  data. Prep the datalogger points for spatial analyses.                       #
#                                                                               #
#################################################################################

############### -- 
# (1) DATA PREP
############### --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Dataloggers/"
setwd(wd)


# using the myClim data processing package: https://cran.r-project.org/web/packages/myClim/index.html

# doing myself a favor and reading in tables that organize the datafiles and give some metadata

## Read pre-defined logger with metadata
# load in two tables - path should work with just the csv of each datafile if the directory is set 
# to the appropriate file already 

# Robust read for files_table
ft <- read.csv("files_table.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Force important columns to character
cols_to_char_ft <- c("serial_number", "locality_id")
existing_cols_ft <- intersect(cols_to_char_ft, names(ft))
ft[existing_cols_ft] <- lapply(ft[existing_cols_ft], as.character)

# Robust read for localities_table
lt <- read.csv("localities_table.csv", sep = ",", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

# Force important columns to character
cols_to_char_lt <- c("locality_id")
existing_cols_lt <- intersect(cols_to_char_lt, names(lt))
lt[existing_cols_lt] <- lapply(lt[existing_cols_lt], as.character)

# Now run myClim
tms.m <- myClim::mc_read_data(
  files_table = ft,
  localities_table = lt,
  silent = TRUE
)

# five loggers have errors in the cleaning process, they say they have different values for 
# moisture, and T1, T2, and T3 at the same time. 

# these loggers are: 95132139, 95132132, 95132130, 95132114, 95132113

# not sure what this means yet 

###########################################################################################

##################### --
# (2) DATA CLEANING
##################### --

# This was performed automatically while reading in, but I want the output file 

tms <- mc_prep_clean(tms.m, silent = T)

# calculate time zones 
# calculate solar time - default time is UTC
# needs longitude to be able to do this, derived from the localities table 
tms <- mc_prep_solar_tz(tms)

mc_info_count(tms) #which returns the number of localities, loggers and sensors in myClim object
mc_info_clean(tms) #returning the data frame with cleaning log
info <- mc_info(tms) #returning data frame with summary per sensor


## crop the time-series to match when they were deployed and retrieved 
start <- as.POSIXct("2024-04-01", tz = "UTC")
end   <- as.POSIXct("2024-11-09", tz = "UTC")

tms   <- mc_prep_crop(tms, start, end)


### visualizing ###

## plot for B03 in the dry southwest corner
tms.plot <- mc_filter(tms, localities = "B03")

p <- mc_plot_line(tms.plot, sensors = c("TMS_T1", "TMS_moist"))
p <- p+ggplot2::scale_x_datetime(date_breaks = "2 weeks", date_labels = "%W")
p <- p+ggplot2::xlab("week")
p <- p+ggplot2::aes(size = sensor_name)
p <- p+ggplot2::scale_size_manual(values = c(0.4 ,0.4))
p <- p+ggplot2::guides(size = "none")
p <- p+ggplot2::scale_color_manual(values = c("darkblue", "red"), name = NULL)
P <- p+ggplot2::geom_line(linewidth = 1)

p


## plot for N39 in the wet northeast corner
tms.plot2 <- mc_filter(tms, localities = "N39")

# running into a repeated issue with the new scale color

q <- mc_plot_line(tms.plot2, sensors = c("TMS_T1", "TMS_moist"))
q <- q+ggplot2::scale_color_manual(values = c("darkblue", "red"), name = NULL)
q <- q+ggplot2::scale_x_datetime(date_breaks = "2 weeks", date_labels = "%W")
q <- q+ggplot2::xlab("week")
q <- q+ggplot2::scale_size_manual(values = c(0.3 ,0.3))
q <- q+ggplot2::guides(size = "none")
q <- q+ggplot2::geom_line(aes(color = sensor_name))

q


# Visualize two that were pulled up by animals 

## plot for J19
tms.plot3 <- mc_filter(tms, localities = "J19")

# running into a repeated issue with the new scale color

r <- mc_plot_line(tms.plot3, sensors = c("TMS_T1", "TMS_moist"))
r <- r+ggplot2::scale_color_manual(values = c("darkblue", "red"), name = NULL)
r <- r+ggplot2::scale_x_datetime(date_breaks = "2 weeks", date_labels = "%W")
r <- r+ggplot2::xlab("week")
r <- r+ggplot2::scale_size_manual(values = c(0.3 ,0.3))
r <- r+ggplot2::guides(size = "none")
r <- r+ggplot2::geom_line(aes(color = sensor_name))

r

# Yep, major variation in temperature after about week 37, and the moisture drops a ton once 
# it's no longer in the ground 


## plot for N11
tms.plot4 <- mc_filter(tms, localities = "N11")

# running into a repeated issue with the new scale color

s <- mc_plot_line(tms.plot4, sensors = c("TMS_T1", "TMS_moist"))
s <- s+ggplot2::scale_color_manual(values = c("darkblue", "red"), name = NULL)
s <- s+ggplot2::scale_x_datetime(date_breaks = "2 weeks", date_labels = "%W")
s <- s+ggplot2::xlab("week")
s <- s+ggplot2::scale_size_manual(values = c(0.3 ,0.3))
s <- s+ggplot2::guides(size = "none")
s <- s+ggplot2::geom_line(aes(color = sensor_name))

s

# Same as the other one, temp and moisture get crazy. Looks like they were pulled up right around the same time too. 


## raster
# this is showing temperature in the T3 sensor (air) over the whole growing season
mc_plot_raster(tms, sensors = c("TMS_T3"))

# this is soil temp
mc_plot_raster(tms, sensors = c("TMS_T1"))

# this is soil moisture 
mc_plot_raster(tms, sensors = c("TMS_moist"))

# this is the plot that shows the two that were pulled out of the ground (J19 and N11)
# Their moisture sensors were picking up almost no moisture, while everyone else was 
# clocking the soil inundation with the fall precipitation 


#### -- 
## Need to calculate volumetric water content from the soil moisture values, which right 
# now just reflect electrical impulses 

# From NEON data, WFDP has sand:silt:clay 59:33.5:7.5, which almost exactly matches to 
# 'Sandy Loam A' tested by Wild et al. from TOMST: https://www.sciencedirect.com/science/article/pii/S0168192318304118

localities <- (lt$locality_id)


VWC <- myClim::mc_calc_vwc(tms, moist_sensor = "TMS_moist", temp_sensor = "TMS_T1",
  output_sensor = "VWC_moisture", soiltype = "sandy loam A",
  localities = localities,
 frozen2NA = TRUE # if TRUE then VWC values are set to NA when the soil temperature is below 0 °C (default TRUE)
)

# This created a new VWC_moisture 'logger' that contains the volumetric water content data for 


## plot for B03 in the dry southwest corner 
tms.plot.test <- mc_filter(VWC, localities = "B03")

t <- mc_plot_line(tms.plot.test, sensors = c("TMS_moist", "VWC_moisture"))
t <- t+ggplot2::scale_x_datetime(date_breaks = "4 weeks", date_labels = "%W")
t <- t+ggplot2::xlab("week")
t <- t+ggplot2::aes(size = sensor_name)
t <- t+ggplot2::scale_size_manual(values = c(0.4 ,0.4))
t <- t+ggplot2::guides(size = "none")
t <- t+ggplot2::scale_color_manual(values = c("darkblue", "lightblue"), name = NULL)

t

# Just VWC so I can view the scale better 
t2 <- mc_plot_line(tms.plot.test, sensors = "VWC_moisture")
t2 <- t2+ggplot2::scale_x_datetime(date_breaks = "4 weeks", date_labels = "%W")
t2 <- t2+ggplot2::xlab("week")
t2 <- t2+ggplot2::aes(size = sensor_name)
t2 <- t2+ggplot2::scale_size_manual(values = c(0.4 ,0.4))
t2 <- t2+ggplot2::guides(size = "none")
t2 <- t2+ggplot2::scale_color_manual(values = c("darkblue"), name = NULL)

t2


## plot for N39 in the wettest corner 
tms.plot.test2 <- mc_filter(VWC, localities = "N39")

u <- mc_plot_line(tms.plot.test2, sensors = c("TMS_moist", "VWC_moisture"))
u <- u+ggplot2::scale_x_datetime(date_breaks = "4 weeks", date_labels = "%W")
u <- u+ggplot2::xlab("week")
u <- u+ggplot2::aes(size = sensor_name)
u <- u+ggplot2::scale_size_manual(values = c(0.4 ,0.4))
u <- u+ggplot2::guides(size = "none")
u <- u+ggplot2::scale_color_manual(values = c("darkblue", "lightblue"), name = NULL)

u

# Just for VWC
u2 <- mc_plot_line(tms.plot.test2, sensors = "VWC_moisture")
u2 <- u2+ggplot2::scale_x_datetime(date_breaks = "4 weeks", date_labels = "%W")
u2 <- u2+ggplot2::xlab("week")
u2 <- u2+ggplot2::aes(size = sensor_name)
u2 <- u2+ggplot2::scale_size_manual(values = c(0.4 ,0.4))
u2 <- u2+ggplot2::guides(size = "none")
u2 <- u2+ggplot2::scale_color_manual(values = c("darkblue"), name = NULL)

u2

# These checks look good, and the VWC values fall within the expected range for temperate 
# forest soils, and exhibit expected variation between the dry and wet sides of the plot. 


### data aggregation to get some summaries ###

# Using the VWC myClim object now


# with defaults only convert Raw-format  to Agg-format
tms.ag <- mc_agg(VWC,fun = NULL, period = NULL)

# aggregate to daily mean, range, coverage, and 95 percentile. 
tms.day <- mc_agg(VWC, fun = c("mean", "range", "coverage", "percentile"),
                  percentiles = 95, period = "day", min_coverage = 0.95)

# aggregate all time-series, return one value per sensor.
tms.all <- mc_agg(VWC, fun = c("mean", "max", "min", "range", "coverage", "percentile"),
                  percentiles = 95, period = "all", min_coverage = 0.95)

# B03, driest and hottest corner 
B03 <- tms.all$localities$B03

# Soil temperature 
# TMS1 mean = 11.88
# TMS1 Max = 20.75
# TMS1 Min = 3.2
# TMS1 range = 17.5

# Air Temperature
# TMS3 mean = 13.22
# TMS3 Max = 39.88
# TMS3 Min = -1.5
# TMS3 range = 41.38

# Soil Moisture 
# VWC_moist mean = 0.23
# VWC_moist Max = 0.39
# VWC_moist Min = 0.10
# VWC_moist range = 0.28


# N39, wettest and coolest corner 
N39 <- tms.all$localities$N39

# Soil temperature 
# TMS1 mean = 11.89
# TMS1 Max = 18.13
# TMS1 Min = 3.88
# TMS1 range = 14.25

# Air Temperature
# TMS3 mean = 12.94
# TMS3 Max = 35.63
# TMS3 Min = -0.81
# TMS3 range = 36.44

# Soil Moisture 
# VWC_moist mean = 0.23
# VWC_moist Max = 0.43
# VWC_moist Min = 0.13
# VWC_moist range = 0.30

# This generates output plots to the specified directory that shows variation in temp
# and moisture over the whole study period 


# Can subset this to just certain localities or certain sensors of interest 
mc_plot_loggers(VWC, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Dataloggers/")


################################################################

######################### --
# (3) SPATIAL DATA PREP 
######################### --

# Load in datalogger coordinates 
logger_data <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/datalogger_coords_GIS.csv")


# Check that data is intitially in EPSG: 4326, which is standard for GPS units 
loggers <- st_as_sf(logger_data, coords = c("Lon_DD","Lat_DD"), crs = 4326)

# check
st_crs(loggers)
st_bbox(loggers)

# In CRS 4236 currently, but I want to reproject to ESG 32610 so I can perform
# analyses on the same map system as my other points, and the WFDP boundary


loggers_utm <- st_transform(loggers, 32610)

# check UTM values (should be ~580000, ~5074000)
st_crs(loggers_utm)
st_bbox(loggers_utm)
head(st_coordinates(loggers_utm))


# TMS_T1 - soil temperature sensor in Tomst TMS (°C)
# TMS_T2 - surface temperature sensor in Tomst TMS (°C)
# TMS_T3 - air temperature sensor in Tomst TMS (°C)
# VWC_moisture - soil volumetric water content calculated from soil moisture sensor ()


# Export the UTM coordinates to read into QGIS
loggers_utm <- loggers_utm %>%
  mutate(UTM_X = st_coordinates(.)[,1],
         UTM_Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()


## Final step - the two dataloggers that were pulled up by animals need to be removed 
# so they are not wrongly used in the analyses - J19 and N11

loggers_utm$Cell <- as.factor(loggers_utm$Cell)

loggers_utm_sub <- loggers_utm %>%
  filter(
    !(Cell %in% c("J19", "N11")))


# Save 
write.csv(loggers_utm_sub, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/dataloggers_utm.csv", row.names = FALSE)

#################################

######################### --
# (4) MERGE ENVIRO DATA  
######################### --

# Merge summary data of the envrionmental variation across the plot with the spatial coordinates 
# This will give a set of variables that can be used for interpolation, or just to visualize 
# variation in soil temp and moisture over space and time in the plot. 


# Create compatible id column 
loggers_utm_sub$locality_id <- loggers_utm_sub$Cell


# Pull out aggregated data from myClim 

# All time series aggregated with one value per sensor 
tms_all_df <- myClim::mc_reshape_long(tms.all)

# Daily aggregate of mean, range, coverage, and 95 percentile 
tms_day_df <- myClim::mc_reshape_long(tms.day)



# Switch to wide format for just the all_data
# The daily data is giant and I need to decide if I need it 
# This is excluding height data, but that is just where the different sensors 
# are located, so that info is readily available 
tms_all_means <- tms_all_df %>%
  dplyr::select(locality_id, sensor_name, value) %>%
  pivot_wider(
    names_from = sensor_name,
    values_from = value)


# This is now organized to have the aggregated data for each sensor locality 

# Merge all mean data with the spatial dataset
loggers_merged <- merge(loggers_utm_sub, tms_all_means, by = "locality_id")

# Save
write.csv(loggers_merged, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_data/dataloggers_utm_mean_data.csv")


## Note: If I want to look at variation over time I can use the tms_day_df, but 
# it'll be a giant file


# -- END -- # 
