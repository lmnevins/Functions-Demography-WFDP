# -----------------------------------------------------------------------------#
# Datalogger soil data 
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
library(geosphere); packageVersion("geosphere")
library(sf); packageVersion("sf")

#################################################################################
#                               Main workflow                                   #
#  Process the soil data from all of the dataloggers retrieved from WFDP using  #
#  the myClim package, specifically designed to process microlimate datalogger  #
#  data. Create a vector of WFDP and interpolate environmental variation        #
#  across the plot.                                                             #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
############### 

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

#####################
# (2) DATA CLEANING
##################### 

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

# This generates output plots to the specified directory that shows varation in temp
# and moisture over the whole study period 

# Can subset this to just certain localities or certain sensors of interest 
mc_plot_loggers(tms, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Dataloggers/")


################################################################

###########################
# (3) GENERATE WFDP VECTOR
########################### 

# Step 1: Define center in UTM (EPSG:32610)
# Approximate UTM coordinates for 45.8197N, -121.9558W
# Convert lat/lon to UTM
center_latlon <- st_sfc(st_point(c(-121.9558, 45.8197)), crs = 4326)
center_utm <- st_transform(center_latlon, 32610)  # UTM Zone 10N

# Step 2: Define half-width/height in meters
half_width <- 400   # E-W
half_height <- 170  # N-S

# Step 3: Build rectangle in UTM coordinates
center_coords <- st_coordinates(center_utm)
xmin <- center_coords[1] - half_width
xmax <- center_coords[1] + half_width
ymin <- center_coords[2] - half_height
ymax <- center_coords[2] + half_height

# Rectangle corners (clockwise)
rect_coords <- matrix(
  c(xmin, ymax,
    xmax, ymax,
    xmax, ymin,
    xmin, ymin,
    xmin, ymax), # close polygon
  ncol = 2, byrow = TRUE
)

# Step 4: Create sf polygon in UTM
plot_boundary_utm <- st_polygon(list(rect_coords)) |> 
  st_sfc(crs = 32610) |> 
  st_sf(plot_name = "WFDP")

# Step 5: Export shapefile
st_write(plot_boundary_utm, "WFDP_boundary_utm.shp", delete_layer = TRUE)

# Step 6. Visualize
plot(st_geometry(plot_boundary_utm), col = NA, border = "darkgreen", lwd = 2)
points(center_lon, center_lat, pch = 19, col = "red")
text(center_lon, center_lat, labels = "Center", pos = 3)

# Right now this is oriented perfectly N-S E-W, but the actual plot is offset a bit 
# Rotate plot to the -5.2 degree offset (switch to positive in relation to north)

# center coordinates
center_coords <- st_coordinates(center_utm)

# Rotation angle in degrees (clockwise from north)
theta <- 5.26
theta_rad <- theta * pi / 180

# Function to rotate points around center
rotate_points <- function(xy, center, angle) {
  x_shift <- xy[,1] - center[1]
  y_shift <- xy[,2] - center[2]
  
  x_rot <- x_shift * cos(angle) - y_shift * sin(angle) + center[1]
  y_rot <- x_shift * sin(angle) + y_shift * cos(angle) + center[2]
  
  cbind(x_rot, y_rot)
}

# Rotate rectangle coordinates
rect_coords_rot <- rotate_points(rect_coords, center_coords, theta_rad)

# Recreate polygon with rotation
plot_boundary_utm_rot <- st_polygon(list(rect_coords_rot)) |>
  st_sfc(crs = 32610) |>
  st_sf(plot_name = "WFDP")

# Optional: export
st_write(plot_boundary_utm_rot, "WFDP_boundary_utm_rotated.shp", delete_layer = TRUE)

# Visual check
plot(st_geometry(plot_boundary_utm_rot), col = NA, border = "blue", lwd = 2)
points(center_coords[1], center_coords[2], pch = 19, col = "red")

# Load into shorter name 
wfdp <- plot_boundary_utm_rot


################################################################

##################################
# (4) DATA PREP FOR INTERPOLATION
################################## 

# Load in datalogger coordinates 
logger_data <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/datalogger_coords_GIS.csv")

# Create compatible id column 
logger_data$locality_id <- logger_data$Cell


# Get a summary of the mean and SD soil temperature and soil moisture across the time period
# tms.all is the aggregated data for each sensor with mean, min, max, and percentiles
# calculated for each variable 
all_data <- mc_reshape_long(tms.all)

all_data$sensor_name <- as.factor(all_data$sensor_name)
all_data$locality_id <- as.factor(all_data$locality_id)


# Join logger data to enviro data file 
env <- merge(all_data, logger_data, by = "locality_id")


# create a 'geometry' column that contains the point coordinates 
env_sf <- st_as_sf(env, coords = c("Lon", "Lat"), crs = 4326)

# Check coordinates read in correctly
coords <- st_coordinates(env_sf) # Looks good

view(coords)

# TMS_T1 - soil temperature sensor in Tomst TMS (°C)
# TMS_T2 - surface temperature sensor in Tomst TMS (°C)
# TMS_T3 - air temperature sensor in Tomst TMS (°C)
# TMS_moist - soil moisture sensor in Tomst TMS (raw TMS units)

# Double check that the logger coordinates and WFDP polygon are using the same UTM
# This allows the interpolation to be performed in meters 

utm_crs <- 32610  # UTM zone 10N covers Washington
loggers_utm <- st_transform(env_sf, crs = utm_crs)
wfdp_utm <- st_transform(wfdp, crs = utm_crs)

# create a grid across WFDP 
# Right now this is at a 5 m resolution, could change this. I don't think I could 
# confidently go any finer resolution. 
grid <- st_make_grid(wfdp_utm, cellsize = 5, what = "centers")
grid <- st_intersection(st_sf(geometry = grid), wfdp_utm)


### Logger coordinates look funky, need to check these 



#################################
