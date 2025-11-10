# -----------------------------------------------------------------------------#
# Neighborhood Modeling 
# Original Author: L. McKinley Nevins 
# November 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     sf v 1.0.19
#                     raster v 3.6.32
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(sf); packageVersion("sf")
library(raster); packageVersion("raster")

#################################################################################
#                               Main workflow                                   #
#  Use the X-Y plot coordinates for each tree in WFDP to calculate the          #
#  neighborhood for each studied focal tree. Incorporate DBH and species data   #
#  into the neighborhoods and use distance to generate an overall crowding      #
#  metric for each focal tree.                                                  #
#                                                                               #
#################################################################################

###############
# (1) DATA PREP
############### 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)

# Load in dataset of all trees in WFDP from the 2021 census 
trees <- read.csv("all_stems_2021.csv") 
# 35,929 trees, including dead trees

# Make a few columns factors 
trees$Cell <- as.factor(trees$Cell)
trees$Species <- as.factor(trees$Species)
trees$Status <- as.factor(trees$Status)

# Status codes: 
# CWD DALB DEAD DOG FIND GONE IN_17 IN_18 IN_19 IN_20 INDY LIVE

# CWD = coarse woody debris; DALB = dead above-live below; DEAD = dead; DOG = dead on ground;
# FIND = need to locate; GONE = missing; IN_# = expected to reach a trackable size class in 
# that year; INDY = expected to reach a trackable size, but year isn't given; LIVE = living condition


# Want to focus only on living and healthy neighbors for which we have accurate location 
# information, so filtering to remove trees with status DALB, DEAD, DOG, FIND, GONE, and INDY

# Also removing any trees with missing status, and any trees with 'NULL' for DBH 

trees_filt <- trees %>%
  filter(
    # keep only trees with acceptable statuses
    !(Status %in% c("DALB", "DEAD", "DOG", "FIND", "GONE", "INDY")),
    # remove missing status
    !is.na(Status),
    # remove NULL or missing DBH values
    !is.na(DBH),
    DBH != "NULL")

# convert DBH to numeric as the Null codes had it read in as a character 
trees_filt <- trees_filt %>%
  mutate(DBH = as.numeric(DBH))

# 30,792 stems remaining, with complete data  


# Load in the dataset of focal trees 
focal <- read.csv("stems_WFDP_20250206_trimmed.csv")

# Here all trees have DBH values for either two or three census times: 2010/2011, 2016, and 2021
# These trees have UTM coordinates, as well as the X-Y coordinates same as all of the other stems 

# Make a few columns factors 
focal$Cell <- as.factor(focal$Cell)
focal$Species <- as.factor(focal$Species)

#################################################################################

############################
# (2) SPATIAL CALCULATIONS
########################### 

# Create universal system of locations for each focal tree and potential neighbors, 
# so the neighbors within a given radius can be identified 

# Convert letters to row index (A=1, B=2, ...)
letter_to_index <- function(letter) match(letter, LETTERS)


# First compute offsets for all of the eligible stems in the plot 
# This calculates the distance and offset of each tree from the A-1 marker in 
# the southwest corner of the plot. 
trees_filt <- trees_filt %>%
  mutate(
    row_letter = str_extract(Cell, "^[A-Z]"),
    col_number = as.numeric(str_extract(Cell, "[0-9]+")),
    row_index = letter_to_index(row_letter),
    Offset_X = (col_number - 1) * 20,
    Offset_Y = (row_index - 1) * 20,
    Plot_X = Offset_X + Cell_X,
    Plot_Y = Offset_Y + Cell_Y)

# Merge the UTMs of the focal trees to the full tree dataset using the stem tag, and 
# remove duplicate rows from multiple DBHs
focal_filt <- dplyr::select(focal, Stem_Tag, UTM_X, UTM_Y) %>% unique()

# combine
all_trees <- left_join(trees_filt, focal_filt, by = "Stem_Tag")

# Create a spatial object of the coordinates of all trees 
trees_sf <- st_as_sf(all_trees, coords = c("Plot_X", "Plot_Y"), crs = 32610) # UTM Zone 10N

### Perform a check that the calculated coordinates are accurate ###

# Identify focal trees with real UTM coords
focal_sf <- all_trees %>%
  filter(!is.na(UTM_X), !is.na(UTM_Y)) %>%
  st_as_sf(coords = c("UTM_X", "UTM_Y"), crs = 32610)

# Check against the plot coordinate calculation 
focal_local_sf <- all_trees %>%
  filter(!is.na(UTM_X), !is.na(UTM_Y)) %>%
  st_as_sf(coords = c("Plot_X", "Plot_Y"), crs = 32610)

# Extract coordinate matrices for both the accurate UTM coordinates and the calculated ones 
utm_coords <- st_coordinates(focal_sf)
local_coords <- st_coordinates(focal_local_sf)

# Run Procrustes analysis to align the two coordinate systems
geo_align <- procrustes(utm_coords, local_coords, scale = TRUE)

# Check alignment summary
summary(geo_align)

# Things generally look quite close, as this is showing that the location of each tree was 
# calculated with an accuracy of within ~ 40 cm. But, WFDP has a small rotation from north that 
# needs to be accounted for. Based on this model, the calculated values can be adjusted to 
# reflect the rotation and slight differences from the UTM coordinates to make sure everything 
# is accurate. 

# Extract components from the procrustes output
R <- matrix(c(0.99578848, 0.09168048,
              -0.09168048, 0.99578848), ncol = 2, byrow = TRUE)
translation <- c(580731, 5074350)
scale_factor <- 1.000421

# Apply transformation to all trees
local_coords <- all_trees %>% dplyr::select(Plot_X, Plot_Y) %>% as.matrix()

utm_pred <- scale_factor * (local_coords %*% R) + 
  matrix(rep(translation, each = nrow(local_coords)), ncol = 2, byrow = FALSE)

# Add new columns
trees_geo <- all_trees %>%
  mutate(UTM_X_pred = utm_pred[,1],
         UTM_Y_pred = utm_pred[,2])


### Plotting to check ###

# Extract from your procrustes result
R <- matrix(c(0.99578848, 0.09168048,
              -0.09168048, 0.99578848), ncol = 2, byrow = TRUE)
translation <- c(580731, 5074350)
scale_factor <- 1.000421

# Extract coordinates
utm_coords <- st_coordinates(focal_sf)
local_coords <- st_coordinates(focal_local_sf)

# Apply transformation
utm_pred <- scale_factor * (local_coords %*% R) +
  matrix(rep(translation, each = nrow(local_coords)), ncol = 2, byrow = FALSE)

# Make data frame for plotting
align_df <- data.frame(
  UTM_X_true = utm_coords[,1],
  UTM_Y_true = utm_coords[,2],
  UTM_X_pred = utm_pred[,1],
  UTM_Y_pred = utm_pred[,2]
)


# Plot real and predicted UTMs for focal trees 
pred_plot <- ggplot(align_df) +
  geom_point(aes(UTM_X_true, UTM_Y_true), color = "steelblue", size = 2) +
  geom_point(aes(UTM_X_pred, UTM_Y_pred), color = "firebrick", size = 2) +
  geom_segment(aes(x = UTM_X_pred, y = UTM_Y_pred,
                   xend = UTM_X_true, yend = UTM_Y_true),
               arrow = arrow(length = unit(0.15, "cm")),
               alpha = 0.4) +
  coord_fixed() +
  labs(x = "UTM Easting", y = "UTM Northing",
       title = "WFDP Plot Alignment (Predicted vs True Coordinates)") +
  theme_minimal()

pred_plot

# Looks almost identical 

# Save to be able to load into QGIS
st_write(st_as_sf(trees_geo, coords = c("UTM_X_pred", "UTM_Y_pred"), crs = 32610),
         "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_trees_georef.shp")


# This is a file of the location of every tree in WFDP, this will maybe blow up 
# my QGIS. I only need to really focus on the subset of trees that are within 
# the neighborhood radius of each focal tree. 


## ! Getting an error down below, so this isn't totally working 

# Define WFDP boundary in local coordinates (m)
wfdp_local <- data.frame(
  x = c(0, 800, 800, 0, 0),
  y = c(0, 0, 340, 340, 0)
)

wfdp_local_sf <- st_as_sf(wfdp_local, coords = c("x", "y"), crs = NA) |>
  summarise(geometry = st_combine(geometry)) |>
  st_cast("POLYGON")

R <- matrix(c(0.99578848, 0.09168048,
              -0.09168048, 0.99578848), ncol = 2, byrow = TRUE)
translation <- c(580731, 5074350)
scale_factor <- 1.000421

# Extract coordinates, transform, and rebuild polygon
local_coords <- st_coordinates(wfdp_local_sf[[1]][[1]])

utm_coords <- scale_factor * (local_coords %*% R) +
  matrix(rep(translation, each = nrow(local_coords)), ncol = 2, byrow = FALSE)

wfdp_utm_sf <- st_sf(
  geometry = st_sfc(st_polygon(list(utm_coords))),
  crs = 32610 # UTM zone 10N (assuming WFDP is in WA/OR)
)


st_write(wfdp_utm_sf, "WFDP_boundary_UTM.shp", delete_dsn = TRUE)


#################################################################################

#### In QGIS, I found the neighbors of each focal tree within a 9 m radius using 
# the Buffer tool to set the radius, and the Vector Tool - Select by Location to 
# pick just trees that intersected with the radius polygon of each focal tree 

# Import neighbor trees shape file 



