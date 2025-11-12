# -----------------------------------------------------------------------------#
# Soil Datalogger Data Processing
# Original Author: L. McKinley Nevins 
# February 2, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 3.5.1
#                     climR v 1.3.0
#                     geosphere v 1.5.20
#                     sf v 1.0.19
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
#  data. Create a vector of WFDP boundary.                                      #
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
ft <- read.table("files_table.csv", sep=",", header = T)
lt <- read.table("localities_table.csv", sep=",", header = T, fill = T)


tms.m <- mc_read_data(files_table = "files_table.csv",
                      localities_table = lt,
                      silent = T)

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


### data aggregation to get some summaries ###

# with defaults only convert Raw-format  to Agg-format
tms.ag <- mc_agg(tms.m,fun = NULL, period = NULL)

# aggregate to daily mean, range, coverage, and 95 percentile. 
tms.day <- mc_agg(tms, fun = c("mean", "range", "coverage", "percentile"),
                  percentiles = 95, period = "day", min_coverage = 0.95)

# aggregate all time-series, return one value per sensor.
tms.all <- mc_agg(tms, fun = c("mean", "max", "min", "range", "coverage", "percentile"),
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
# TMS_moist mean = 1632
# TMS_moist Max = 2434
# TMS_moist Min = 1051
# TMS_moist range = 1383


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
# TMS_moist mean = 1574
# TMS_moist Max = 2706
# TMS_moist Min = 1144
# TMS_moist range = 1562

# This generates output plots to the specified directory that shows variation in temp
# and moisture over the whole study period 


# Can subset this to just certain localities or certain sensors of interest 
mc_plot_loggers(tms, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Dataloggers/")


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

# In CRS 4236 currently, but I want to reproject to ESG 32610 so I can perform the kriging 
# analyses on the same map system as my other points, and the WFDP boundary


loggers_utm <- st_transform(loggers, 32610)

# check UTM values (should be ~580000, ~5074000)
st_crs(loggers_utm)
st_bbox(loggers_utm)
head(st_coordinates(loggers_utm))


# TMS_T1 - soil temperature sensor in Tomst TMS (°C)
# TMS_T2 - surface temperature sensor in Tomst TMS (°C)
# TMS_T3 - air temperature sensor in Tomst TMS (°C)
# TMS_moist - soil moisture sensor in Tomst TMS (raw TMS units)


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
  select(locality_id, sensor_name, value) %>%
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
