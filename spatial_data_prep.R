# -----------------------------------------------------------------------------#
# Spatial data preparation for neighborhood modeling in WFDP 
# Original Author: L. McKinley Nevins 
# November 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 4.0.1
#                     sf v 1.0.23
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
#  Use the X-Y plot coordinates for each tree in WFDP to find the locations of  #
#  all living trees in the plot. Generate vector of WFDP perimeter. Save files  #
#  for use in spatial calculations in QGIS to find focal tree neighbors.        #
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

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

############################ -- 
# (2) SPATIAL CALCULATIONS
############################ --

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
# This takes the UTM_X_pred and UTM_Y_pred columns and saves them into the geometry column 
st_write(st_as_sf(trees_geo, coords = c("UTM_X_pred", "UTM_Y_pred"), crs = 32610),
         "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_trees_georef.shp")

#################################################################################

########################### -- 
# (3) GENERATE WFDP VECTOR
########################### --

# Generate local coordinates based off of plot dimensions 
wfdp_local_coords <- matrix(
  c(0,   0,
    800, 0,
    800, 340,
    0,   340,
    0,   0), 
  ncol = 2, byrow = TRUE
)
colnames(wfdp_local_coords) <- c("X", "Y")

# --- ensure numeric n x 2 matrix ---
local_coords <- as.matrix(wfdp_local_coords)
if(ncol(local_coords) < 2) stop("local_coords must have 2 columns (X,Y).")
local_coords <- local_coords[, 1:2, drop = FALSE]  # enforce n x 2

# --- apply transformation ---
n <- nrow(local_coords)
translation_mat <- matrix(rep(translation, each = n), nrow = n, ncol = 2, byrow = FALSE)

utm_coords <- scale_factor * (local_coords %*% R) + translation_mat

# --- build sf polygon in UTM coordinates ---
wfdp_utm_poly <- st_polygon(list(as.matrix(utm_coords)))
wfdp_utm_sf <- st_sfc(wfdp_utm_poly, crs = 32610) %>% st_sf(plot = "WFDP_boundary")

# --- quick plot to check ---
plot(st_geometry(wfdp_utm_sf), border = "darkgreen", lwd = 2, asp = 1)
# if you have focal points in 'focal_sf' (EPSG:32610), plot them:
if(exists("focal_sf")) {
  plot(st_geometry(focal_sf), add = TRUE, col = "red", pch = 16)
}

# --- write shapefile  ---
st_write(wfdp_utm_sf, "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_boundary_UTM.shp", delete_layer = TRUE)


#################################################################################

## -- END -- ##

