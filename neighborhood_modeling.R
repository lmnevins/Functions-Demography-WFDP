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
st_write(wfdp_utm_sf, "WFDP_boundary_UTM.shp", delete_layer = TRUE)


#################################################################################

############################ --
# (4) PROCESS NEIGHBOR DATA
############################ -- 

#### In QGIS, I found the neighbors of each focal tree within a 9 m radius using 
# the Buffer tool to set the radius, and Join Attributes by Location tool to 
# pick just trees that intersected with the radius polygon of each focal tree 

# Import neighbor trees file 
neighbors <- read.csv("./neighbor_trees.csv") #2,147 neighbors in total 

# columns related to the focal tree are generally identified with a 'focal_' prefix, 
# while columns related to the neighbors have a 'neigh_' prefix 


# Explore the neighborhood data a bit 

# summarize number of neighbors per focal tree
neighbor_summary <- neighbors %>%
  group_by(focal_stem_tag) %>%
  summarise(
    focal_species = first(focal_species),
    focal_DBH = first(focal_DBH),
    n_neighbors = n(),
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    sd_neighbor_DBH = sd(neigh_DBH, na.rm = TRUE))

# overall mean and range
summary_stats <- neighbor_summary %>%
  summarise(
    mean_neighbors = mean(n_neighbors),
    sd_neighbors = sd(n_neighbors),
    min_neighbors = min(n_neighbors),
    max_neighbors = max(n_neighbors))

neighbor_summary
summary_stats

# mean_neighbors     sd_neighbors     min_neighbors     max_neighbors
#       35.8              17.4               7                86


# plot count of neighbors 
neigh_count <- ggplot(neighbor_summary, aes(x = n_neighbors)) +
  geom_histogram(binwidth = 2, color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Neighbor Counts per Focal Tree",
    x = "Number of Neighbors (within 10 m)",
    y = "Number of Focal Trees"
  )

neigh_count


# Check neighbors by focal species 
neigh_spp <- ggplot(neighbor_summary, aes(x = focal_species, y = n_neighbors, fill = focal_species)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Neighbor Density by Focal Species",
    x = "Focal Species",
    y = "Number of Neighbors"
  ) +
  theme(legend.position = "none")

neigh_spp 


# Check for relationship between tree size and number of neighbors 
tree_size <- ggplot(neighbor_summary, aes(x = focal_DBH, y = n_neighbors)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    title = "Relationship Between Focal Tree DBH and Number of Neighbors",
    x = "Focal Tree DBH (cm)",
    y = "Number of Neighbors"
  )

tree_size 

size_mod <- lm(focal_DBH ~ n_neighbors, data = neighbor_summary)

summary(size_mod)

# No relationship here 


# Check for relationship between focal tree size and neighbor size 
neigh_size <- ggplot(neighbor_summary, aes(x = focal_DBH, y = mean_neighbor_DBH)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    title = "Relationship Between Focal Tree Size and Average Neighbor Size",
    x = "Focal Tree DBH (cm)",
    y = "Mean Neighbor DBH (cm)")

neigh_size

size_mod2 <- lm(focal_DBH ~ mean_neighbor_DBH, data = neighbor_summary)

summary(size_mod2)

# No relationship here 


#################################################################################

##################################### -- 
# (5) CALCULATE NEIGHBORHOOD METRICS
##################################### -- 

# Need to calculate distances between the focal tree to each of its neighbors

# The trees_geo df has the predicted and corrected UTMs for each of the trees, so I can 
# subset this, then merge it with the neighbors data 

trees_geo$neigh_stem_tag <- trees_geo$Stem_Tag

geo_sub <- dplyr::select(trees_geo, neigh_stem_tag, UTM_X_pred, UTM_Y_pred)

neighbors <- merge(neighbors, geo_sub, by = "neigh_stem_tag")



# Each neighbor tree now has UTM coordinates, so can calculate it's distance to the focal tree


# Make sure coordinates are numeric
neighbors <- neighbors %>%
  mutate(UTM_X_pred = as.numeric(UTM_X_pred),
         UTM_Y_pred = as.numeric(UTM_Y_pred))

# Make a few things factors 
neighbors$focal_cell <- as.factor(neighbors$focal_cell)
neighbors$focal_species <- as.factor(neighbors$focal_species)
neighbors$neigh_species <- as.factor(neighbors$neigh_species)


# Calculate distance and then calculate the crowding index per neighbor
neighbors <- neighbors %>%
  mutate(
    # Step 1: Euclidean distance between focal and neighbor trees
    distance = sqrt((UTM_X_pred - UTM_X)^2 + (UTM_Y_pred - UTM_Y)^2),
    
    # Step 2: Neighborhood crowding index contribution per neighbor
    size_dist2 = (neigh_DBH^2) / (distance^2)
  )


# Summarize crowding 
crowding_summary <- neighbors %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(size_dist2, na.rm = TRUE),
    .groups = "drop")

# A high crowding index means the focal tree is surrounded by many large, nearby neighbors 
# A low crowding index means few, small, or distant neighbors 

### Visualize ### 

# Relationship between focal tree size and crowding 
size_crowd <- ggplot(crowding_summary, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding")

size_crowd 


crowd_mod1 <- lm(focal_DBH ~ crowding_index, data = crowding_summary)

summary(crowd_mod1)

# Significant relationship - bigger trees are more crowded (though this study was focused on 
# trees 10 - 20 cm DBH so this is only a small snapshot)

# Crowding across different species 
spp_crowd <- ggplot(crowding_summary, aes(x = focal_species, y = crowding_index, fill = focal_species)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    x = "Focal Species",
    y = "Neighborhood Crowding Index",
    title = "Variation in Crowding by Focal Species"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

spp_crowd 


# Merge back in UTM info to be able to plot spatially 
geo_sub_focal <- dplyr::select(neighbors, focal_stem_tag, UTM_X, UTM_Y) %>% unique()

crowding_summary <- merge(crowding_summary, geo_sub_focal, by = "focal_stem_tag")


# Spatial map of crowding across the plot 
space_crowd <- ggplot(crowding_summary, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Neighborhood Crowding")

space_crowd 

# Doesn't seem like a consistent trend with crowding across the plot 


#################################################################################

##################################### -- 
# (6) RELATIONSHIPS TO DEMOGRAPHY 
##################################### -- 

# Read in tree demographic data 
growth <- read.csv("stems_WFDP_20250206_trimmed.csv")


# get starting and ending diameters for all trees
diams <- growth %>%
  group_by(Stem_Tag, Species) %>%
  summarise(
    dia_first = DBH[which.min(DBH_DATE)],
    dia_last = DBH[which.max(DBH_DATE)],
    year_first = min(DBH_DATE),
    year_last = max(DBH_DATE),
    .groups = "drop"
  )

# get diameter difference between the two time points 
diams <- diams %>%
  mutate(diam_diff = dia_last - dia_first)

# Calculate relative growth rate for each tree 
diams <- diams %>%
  mutate(RGR = (log(dia_last) - log(dia_first)) / (year_last - year_first))


# Add column for stem tag to match the crowding data 
diams$focal_stem_tag <- diams$Stem_Tag

# Pair diams with the crowding summary file according to the focal_stem_tag
crowd_growth_summary <- merge(diams, crowding_summary, by = "focal_stem_tag")

### Visualize ### 

#Explore relationships between crowding/neighbors and RGR


# Relationship between focal tree RGR and # neighbors 
RGR_num_neigh <- ggplot(crowd_growth_summary, aes(x = num_neighbors, y = RGR)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(x = "Number of Neighbors", y = expression("Relative Growth Rate ("*yr^{-1}*")"))

RGR_num_neigh


RGR_mod1 <- lm(RGR ~ num_neighbors, data = crowd_growth_summary)

summary(RGR_mod1)

# Significant relationship, trees appear to grow faster with more close neighbors 
# Adjusted R-squared:  0.1455, p-value: 0.001544


# Relationship between focal tree RGR and mean neighbor DBH 
RGR_size_neigh <- ggplot(crowd_growth_summary, aes(x = mean_neighbor_DBH, y = RGR)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(x = "Mean Neighbor DBH", y = expression("Relative Growth Rate ("*yr^{-1}*")"))

RGR_size_neigh


RGR_mod2 <- lm(RGR ~ mean_neighbor_DBH, data = crowd_growth_summary)

summary(RGR_mod2)

# Significant relationship, trees appear to grow faster when their neighbors are smaller on average 
# Adjusted R-squared:  0.165, p-value: 0.0007533


# Relationship between focal tree RGR and crowding index
RGR_crowd_neigh <- ggplot(crowd_growth_summary, aes(x = crowding_index, y = RGR)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(x = "Crowding Index", y = expression("Relative Growth Rate ("*yr^{-1}*")"))

RGR_crowd_neigh


RGR_mod3 <- lm(RGR ~ crowding_index, data = crowd_growth_summary)

summary(RGR_mod3)

# Significant relationship, trees appear to grow faster when their neighbors are smaller on average 
# Adjusted R-squared:  0.165, p-value: 0.0007533



