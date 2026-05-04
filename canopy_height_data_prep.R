# -----------------------------------------------------------------------------#
# Load NEON Canopy Height Model data and prepare for QGIS  
# Original Author: L. McKinley Nevins 
# February 5, 2026
# Software versions:  R v 4.5.2
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     tibble v 3.3.0
#                     vegan v 2.7.2
#                     sf v 1.0.23
#                     terra v 1.8.86
#                     neonUtilities v 3.0.3
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(vegan); packageVersion("vegan")
library(sf); packageVersion("sf")
library(terra); packageVersion("terra")
library(neonUtilities); packageVersion("neonUtilities")

#################################################################################
#                               Main workflow                                   #
#  Pull data from the Canopy Height model from NEON and prepare it to be able   #
#  to extract canopy info for 9 and 20 m radii around each focal tree in QGIS.  #
#  Also pull 1 m Digital Elevation Model data for a basemap layer for WFDP.     # 
#                                                                               # 
#################################################################################

###############
# (1) DATA PREP
###############

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/"
setwd(wd)

# Load WFDP plot boundary (polygon shapefile)
wfdp_poly <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_boundary_UTM.shp")

# Load focal tree points 
focal_trees <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_focal_trees.shp")

# View plot 
WFDP_plot <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black") +
  geom_sf(data = focal_trees) +
  theme_minimal() +
  labs(title = "WFDP boundary")

WFDP_plot


## PREP SPATIAL DATA


# Right now this is reflecting the 5.2 rotation of the plot relative to true north. For generating figures 
# (only figures, not for doing any spatial analyses), I want to rotate the plot and point coordinates 
# to be aligned with true north 

# Rotation helper function 
rotation <- function(a) {
  r <- a * pi / 180 # Convert degrees to radians
  matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
}

# 1. Get the geometry and its centroid for both the plot boundary and the points 
wfdp_poly_geom <- st_geometry(wfdp_poly)
wfdp_poly_centroid <- st_centroid(st_union(wfdp_poly_geom))

focal_trees_geom <- st_geometry(focal_trees)
focal_trees_centroid <- st_centroid(st_union(focal_trees_geom))


# 2. Apply rotation (5.2 degrees)
wfdp_poly_rotated_geom <- (wfdp_poly_geom - wfdp_poly_centroid) * rotation(5.2) + wfdp_poly_centroid

focal_trees_rotated_geom <- (focal_trees_geom - focal_trees_centroid) * rotation(5.2) + focal_trees_centroid


# 3. Update the sf object with the new geometry
wfdp_poly_rotated <- st_set_geometry(wfdp_poly, wfdp_poly_rotated_geom)


focal_trees_rotated <- st_set_geometry(focal_trees, focal_trees_rotated_geom)


# Replot to check 
WFDP_plot_rotated <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black") +
  geom_sf(data = focal_trees_rotated) +
  theme_minimal() +
  labs(title = "WFDP boundary rotated")

WFDP_plot_rotated


# NOTE: This essentially breaks the CRS and had converted this back to UTM coordinates. That's okay 
# because I'm not plotting with the coordinates showing anyways. This is just to get the relative positions 
# of the trees in the plot and to eachother. 

####################################################

# Shooting to get the canopy height data to be able to assess this for each focal tree, and 
# the surrounding neighborhood within defined radii 

# Convert boundary into easting and northing values to give to the neonUtilities function to 
# define the range of the canopy data 


bb <- st_bbox(wfdp_poly)
bb

eastings  <- c(bb["xmin"], bb["xmax"])
northings <- c(bb["ymin"], bb["ymax"])

#################################################################################

#############################
# (2) READ NEON CANOPY DATA 
#############################

#Read in LiDAR data tiles from the canopy height model data product 
chm_files <- byTileAOP(
  dpID = "DP3.30015.001",   # Canopy Height Model
  site = "WREF",
  year = 2022,             # Doing 2022 as this was when I sampled the fungal communities, and 
  check.size = FALSE,       # was the stem census data 
  easting = eastings,
  northing = northings,
  buffer = 20             # Add a little boundary around the defined area, just in case 
)

# This saves to a file in my working directory, so need to read in the two gtif files that were saved 

# WREF has been flown 7 times since the beginning of NEON 

chm1 <- rast("/Users/mckinleynevins/Dropbox/WSU/WFDP_Chapter_3_Project/DP3.30015.001/neon-aop-products/2022/FullSite/D16/2022_WREF_5/L3/DiscreteLidar/CanopyHeightModelGtif/NEON_D16_WREF_DP3_580000_5074000_CHM.tif")

chm2 <- rast("/Users/mckinleynevins/Dropbox/WSU/WFDP_Chapter_3_Project/DP3.30015.001/neon-aop-products/2022/FullSite/D16/2022_WREF_5/L3/DiscreteLidar/CanopyHeightModelGtif/NEON_D16_WREF_DP3_581000_5074000_CHM.tif")

# plot the CHM file - western side of the plot 
plot(chm1, col = topo.colors(10))

# plot the CHM file - eastern side of the plot 
plot(chm2, col = topo.colors(10))

# This is showing canopy height across WFDP. 


plot(st_geometry(wfdp_poly), border = "black")
plot(ext(chm1), add = TRUE, border = "red")
plot(ext(chm2), add = TRUE, border = "red")


# Save western raster to load into QGIS 
writeRaster(
  chm1,
  filename = "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WREF_CHM1.tif",
  filetype = "GTiff",
  overwrite = TRUE
)


# Save eastern raster to load into QGIS 
writeRaster(
  chm2,
  filename = "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WREF_CHM2.tif",
  filetype = "GTiff",
  overwrite = TRUE
)


# Merge the two CHM files 

# This averages any areas they may overlap 
chm_merged <- mosaic(chm1, chm2, fun = mean)

# Plot merged 
plot(chm_merged, col = topo.colors(10))

# Crop to the WFDP boundary 
chm_crop <- crop(chm_merged, wfdp_poly)

# check again 
plot(chm_crop, col = topo.colors(10))

# Mask anything outside of the boundary 
chm_plot <- mask(chm_crop, wfdp_poly)

# check again 
plot(chm_plot, col = topo.colors(10))


#################################################################################

#############################
# (3) READ NEON DSM DATA 
#############################

#Read in LiDAR data tiles from the canopy height model data product 
chm_files <- byTileAOP(
  dpID = "DP3.30024.001",   # Digital Elevation Model
  site = "WREF",
  year = 2022,             # Doing 2022, timing doesn't matter as much here 
  check.size = FALSE,       
  easting = eastings,
  northing = northings,
  buffer = 20             # Add a little boundary around the defined area, just in case 
)

# This saves to a file in my working directory, so need to read in the two gtif files that were saved 

# This gives files for both the Digital Surface Model, and the Digital Terrain Model. The DSM is the raw heights 
# from the LiDAR point cloud (so actually could be interesting), the Digital Terrain Model is elevations from the 
# physical terrain surface 

# In that file, there are two files for the two tiles that cover the requested WREF area 

dtm1 <- rast("/Users/mckinleynevins/Dropbox/WSU/WFDP_Chapter_3_Project/DP3.30024.001/neon-aop-products/2022/FullSite/D16/2022_WREF_5/L3/DiscreteLidar/DTMGtif/NEON_D16_WREF_DP3_580000_5074000_DTM.tif")

dtm2 <- rast("/Users/mckinleynevins/Dropbox/WSU/WFDP_Chapter_3_Project/DP3.30024.001/neon-aop-products/2022/FullSite/D16/2022_WREF_5/L3/DiscreteLidar/DTMGtif/NEON_D16_WREF_DP3_581000_5074000_DTM.tif")

# plot the DTM file - western side of the plot 
plot(dtm1, col = topo.colors(10))

# plot the DTM file - eastern side of the plot 
plot(dtm2, col = topo.colors(10))

# This is showing bare terrain elevation across WFDP. 


plot(st_geometry(wfdp_poly), border = "black")
plot(ext(dtm1), add = TRUE, border = "red")
plot(ext(dtm2), add = TRUE, border = "red")


# Save western raster to load into QGIS 
writeRaster(
  dtm1,
  filename = "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WREF_DTM1.tif",
  filetype = "GTiff",
  overwrite = TRUE
)


# Save eastern raster to load into QGIS 
writeRaster(
  dtm2,
  filename = "~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WREF_DTM2.tif",
  filetype = "GTiff",
  overwrite = TRUE
)


#################################################################################

########################
# (4) ANALYZE CHM DATA 
########################

# sample 100k pixels for speed
vals <- values(chm_plot, mat = FALSE)
vals <- vals[!is.na(vals)]

hist(vals, breaks = 100, main = "CHM Height Distribution", xlab = "Height (m)")

# Distribution of heights overall 
quantile(vals, probs = seq(0,1,0.1))

# Distribution at lower heights 
quantile(vals, probs = seq(0,0.3,0.05))


# Based on these distributions, it looks like <2m would be essentially down to the bare ground, 
# <4 m is the threshold for gaps used by Francis et al. 2023
# <5 m would include the shrub layer and very young trees 
# <10 m would be a canopy opening surrounded by mature forest that is on average about 30 m tall 
# <12 m is the cut-off for the 'dim zone' defined by Parker 1997 in WFDP 
# <40 m is the cut-off for the bottom of the 'bright zone' defined by Parker 1997 in WFDP 
  # Between these two zones is a transition area where the light continues to decrease
# 60 m is the tallest tree at WFDP 


# Set buffers for 9 and 20 m neighborhoods around each of the focal trees 
buffers9  <- st_buffer(focal_trees, 9)
buffers20 <- st_buffer(focal_trees, 20)


# Create binary rasters for each of the 3 height levels, so it will read as true or false for if 
# the height is less than each threshold

# 2 m
gap2 <- chm_plot < 2

plot(gap2, col = topo.colors(10))

# 4 m
gap4 <- chm_plot < 4

plot(gap4, col = topo.colors(10))


# 5 m 
gap5 <- chm_plot < 5

plot(gap5, col = topo.colors(10))


# 10 m
gap10 <- chm_plot < 10

plot(gap10, col = topo.colors(10))


# 12 m
gap12 <- chm_plot < 12

plot(gap12, col = topo.colors(10))


# 40 m
gap40 <- chm_plot < 40

plot(gap40, col = topo.colors(10))



# For the 9 m and 20 m buffers, extract the proportions of open to closed canopy from the binary rasters for 
# each of the three thresholds 
# 9m
gap2_9  <- terra::extract(gap2, vect(buffers9),  fun = mean, na.rm = TRUE)
gap4_9  <- terra::extract(gap4, vect(buffers9),  fun = mean, na.rm = TRUE)
gap5_9  <- terra::extract(gap5, vect(buffers9),  fun = mean, na.rm = TRUE)
gap10_9 <- terra::extract(gap10, vect(buffers9), fun = mean, na.rm = TRUE)
gap12_9 <- terra::extract(gap12, vect(buffers9), fun = mean, na.rm = TRUE)
gap40_9  <- terra::extract(gap40, vect(buffers9),  fun = mean, na.rm = TRUE)

# 20m
gap2_20  <- terra::extract(gap2, vect(buffers20),  fun = mean, na.rm = TRUE)
gap4_20  <- terra::extract(gap4, vect(buffers20),  fun = mean, na.rm = TRUE)
gap5_20  <- terra::extract(gap5, vect(buffers20),  fun = mean, na.rm = TRUE)
gap10_20 <- terra::extract(gap10, vect(buffers20), fun = mean, na.rm = TRUE)
gap12_20 <- terra::extract(gap12, vect(buffers20), fun = mean, na.rm = TRUE)
gap40_20  <- terra::extract(gap40, vect(buffers20),  fun = mean, na.rm = TRUE)

# These dataframes don't have the trees labeled anymore, but the ID number is the row number from the buffers file
# Can merge all of these back together again with identifying info 

gaps_9 <- cbind(buffers9, gap2_9, gap4_9, gap5_9, gap10_9, gap12_9, gap40_9)

gaps_9 <- gaps_9 %>%
  rename(prop_2m_9 = NEON_D16_WREF_DP3_580000_5074000_CHM, prop_4m_9 = NEON_D16_WREF_DP3_580000_5074000_CHM.1,
         prop_5m_9 = NEON_D16_WREF_DP3_580000_5074000_CHM.2, prop_10m_9 = NEON_D16_WREF_DP3_580000_5074000_CHM.3, 
          prop_12m_9 = NEON_D16_WREF_DP3_580000_5074000_CHM.4, prop_40m_9 = NEON_D16_WREF_DP3_580000_5074000_CHM.5)


gaps_20 <- cbind(buffers20, gap2_20, gap4_20, gap5_20, gap10_20, gap12_20, gap40_20)

gaps_20 <- gaps_20 %>%
  rename(prop_2m_20 = NEON_D16_WREF_DP3_580000_5074000_CHM,  prop_4m_20 = NEON_D16_WREF_DP3_580000_5074000_CHM.1, 
         prop_5m_20 = NEON_D16_WREF_DP3_580000_5074000_CHM.2, prop_10m_20 = NEON_D16_WREF_DP3_580000_5074000_CHM.3, 
         prop_12m_20 = NEON_D16_WREF_DP3_580000_5074000_CHM.4, prop_40m_20 = NEON_D16_WREF_DP3_580000_5074000_CHM.5)



# These dataframes now have proportions of canopy openness for multiple height categories for each focal tree. 

# Save to analyze 
gaps_9 <- dplyr::select(gaps_9, Stem_Tag, prop_2m_9, prop_4m_9, prop_5m_9, prop_10m_9, prop_12m_9, prop_40m_9) 

gaps_20 <- dplyr::select(gaps_20, Stem_Tag, prop_2m_20, prop_4m_20, prop_5m_20, prop_10m_20, prop_12m_20, prop_40m_20)

# save as sf objects 
st_write(st_as_sf(gaps_9, crs = 32610),
         "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_canopy_openness_9m.shp", append = FALSE)

st_write(st_as_sf(gaps_20, crs = 32610),
         "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_canopy_openness_20m.shp", append = FALSE)


#################################################################################

############################ -- 
# (5) PLOT CANOPY OPENNESS 
############################ -- 

# Right now this is reflecting the 5.2 rotation of the plot relative to true north. For generating figures 
# (only figures, not for doing any spatial analyses), I want to rotate the plot and point coordinates 
# to be aligned with true north 

# Rotation helper function 
rotation <- function(a) {
  r <- a * pi / 180 # Convert degrees to radians
  matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
}

# 1. Get the geometry and its centroid for both the plot boundary and the points 
wfdp_poly_geom <- st_geometry(wfdp_poly)
wfdp_poly_centroid <- st_centroid(st_union(wfdp_poly_geom))

gaps_9_geom <- st_geometry(gaps_9)
gaps_9_centroid <- st_centroid(st_union(gaps_9_geom))


# 2. Apply rotation (5.2 degrees)
wfdp_poly_rotated_geom <- (wfdp_poly_geom - wfdp_poly_centroid) * rotation(5.2) + wfdp_poly_centroid

gaps_9_rotated_geom <- (gaps_9_geom - gaps_9_centroid) * rotation(5.2) + gaps_9_centroid


# 3. Update the sf object with the new geometry
wfdp_poly_rotated <- st_set_geometry(wfdp_poly, wfdp_poly_rotated_geom)


gaps_9_rotated <- st_set_geometry(gaps_9, gaps_9_rotated_geom)


# Replot to check 
WFDP_plot_rotated <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black") +
  geom_sf(data = gaps_9_rotated) +
  theme_minimal() +
  labs(title = "9 m Focal Tree Radii")

WFDP_plot_rotated


# NOTE: This essentially breaks the CRS and had converted this back to UTM coordinates. That's okay 
# because I'm not plotting with the coordinates showing anyways. This is just to get the relative positions 
# of the trees in the plot and to eachother. 



# -- #### END ##### -- 
