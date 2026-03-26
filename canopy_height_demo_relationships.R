# -----------------------------------------------------------------------------#
# Canopy Height analyses 
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
#                     cowplot v 1.2.0
#                     ggspatial v 1.1.10
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(vegan); packageVersion("vegan")
library(sf); packageVersion("sf")
library(terra); packageVersion("terra")
library(cowplot); packageVersion("cowplot")
library(ggspatial); packageVersion("ggspatial")

#################################################################################
#                               Main workflow                                   #
#  Use the Canopy Height data derived from NEON and analyze the canopy          #
#  characteristics in 9 and 20 m radii around each focal tree. Model the        #
#  relationships of canopy height to the growth of each focal tree to see if    #
#  canopy height of neighbors has an affect on the focal tree.                  # 
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


# Read in canopy height dataframe
canopy <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/canopy_height_summaries.csv")           

# Read in 9 openness 
gaps_9m <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_canopy_openness_9m.shp")

gaps_9m <- as.data.frame(gaps_9m)

# Read in 20 openness 
gaps_20m <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_canopy_openness_20m.shp")

# Fix column naming error 
gaps_20m$Stem_Tag <- gaps_20m$Stem_Tg

gaps_20m <- as.data.frame(gaps_20m)


# Read in tree demographic data 
growth <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/stems_WFDP_20250206_trimmed.csv")

# get diameters for each of the three timepoints for all 60 trees
# Some trees have 2010 or 2011 for their first census date 
rgr_intervals <- growth %>%
  arrange(Stem_Tag, DBH_DATE) %>%
  group_by(Stem_Tag, Species) %>%
  mutate(
    DBH_prev  = lag(DBH),
    year_prev = lag(DBH_DATE),
    RGR_interval = (log(DBH) - log(DBH_prev)) / (DBH_DATE - year_prev)
  ) %>%
  filter(!is.na(RGR_interval))

# calc mean RGR across both intervals 
diams <- rgr_intervals %>%
  summarise(
    mean_RGR = mean(RGR_interval, na.rm = TRUE),
    n_intervals = n(),
    .groups = "drop"
  )

#### Merge diams file to the canopy height file 

# Pair diams with the summary files according to the focal_stem_tag
all_canopy_growth <- merge(diams, canopy, by = "Stem_Tag")

all_canopy_growth <- merge(all_canopy_growth, gaps_9m, by = "Stem_Tag")

all_canopy_growth <- merge(all_canopy_growth, gaps_20m, by = "Stem_Tag")


# Convert to sf objects
all_canopy_growth_sf <- st_as_sf(all_canopy_growth, coords = c("UTM_X", "UTM_Y"), crs = 32610)



# Create object to save for later analyses 

all_canopy <- merge(canopy, gaps_9m, by = "Stem_Tag")

all_canopy <- merge(all_canopy, gaps_20m, by = "Stem_Tag")

# subset 

all_canopy <- dplyr::select(all_canopy, Cell, Species, Stem_Tag, X9m_mean, X20m_mean, prop_2m_9, 
                            prop_5m_9, prop_10m_9, pr_2_20, pr_5_20, p_10_20)

write.csv(all_canopy, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_all_canopy_data.csv", 
          row.names = FALSE)


# Get summary of variables for results section text
summary(all_canopy_growth)


# make dataframe of full species names 
sci_name <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
              "T. plicata", "T. heterophylla")

Species <- c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")


taxa <- data.frame(sci_name, Species)

# Merge to all of the files 
all_canopy_growth$Species <- all_canopy_growth$Species.x

all_canopy_growth <- merge(all_canopy_growth, taxa, by = "Species")


# This has a lot of different geometries in it that I'm going to ignore for the next step and just go back 
# to the original files 

##################

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


# Plotting down below uses the spatial coordinates of the focal trees only, so I can merge the rotated 
# focal tree coordinates to a subsetted version of the all_canopy_growth 

canopy_rgr_sub <- dplyr::select(all_canopy_growth, Species, Stem_Tag, mean_RGR, Cell, X9m_mean, X9m_min, X9m_max, 
                               X20m_mean, X20m_min, X20m_max, prop_2m_9, prop_5m_9, prop_10m_9, pr_2_20, 
                               pr_5_20, p_10_20, sci_name)


growth_spatial_df <- merge(focal_trees_rotated, canopy_rgr_sub, by = "Stem_Tag")


#################################################################################

###################################################
# (2) VISUALIZE CANOPY VARIATION ACROSS THE PLOT 
###################################################


# 1. Mean Canopy Height - 9m radius 
plot_mean_9m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = X9m_mean), size = 4) +
  scale_color_gradient2(low = "tan", mid = "green4", high = "#1D2E28", # adjust so both plots use the same scale 
                        midpoint = 25, limits = c(0, 45), 
                       name = "Mean Canopy\nHeight (m)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

# removing legends for plotting too 

plot_mean_9m


# 2. Mean Canopy Height - 20m radius 
plot_mean_20m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = X20m_mean), size = 6) +
  scale_color_gradient2(low = "tan", mid = "green4", high = "#1D2E28", 
                        midpoint = 25, limits = c(0, 45), 
                        name = "Mean Canopy\nHeight (m)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_mean_20m



# 3. Percent Canopy Openness 2 m height - 20 m radius 
plot_2open_20m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = pr_2_20), size = 6) +
  scale_color_gradient2(high = "tan", mid = "lightblue", low = "blue4", 
                        midpoint = 0.25, limits = c(0, 0.45), 
                        name = "Proportion Canopy\nOpenness 2m (%)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_2open_20m



# 4. Percent Canopy Openness 5 m height - 20 m radius 
plot_5open_20m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = pr_5_20), size = 6) +
  scale_color_gradient2(high = "tan", mid = "lightblue", low = "blue4", 
                        midpoint = 0.25, limits = c(0, 0.45), 
                        name = "Proportion Canopy\nOpenness 5m (%)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_5open_20m



# 6. Percent Canopy Openness 10 m height - 20 m radius 
plot_10open_20m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = p_10_20), size = 6) +
  scale_color_gradient2(high = "tan", mid = "lightblue", low = "blue4", 
                        midpoint = 0.25, limits = c(0, 0.45), 
                        name = "Proportion Canopy\nOpenness 10m (%)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_10open_20m



# 7. Percent Canopy Openness 2 m height - 9m radius 
plot_2open_9m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = prop_2m_9), size = 4) +
  scale_color_gradient2(high = "tan", mid = "lightblue", low = "blue4", 
                        midpoint = 0.375, limits = c(0, 0.75), 
                        name = "Proportion Canopy\nOpenness 2m (%)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_2open_9m



# 8. Percent Canopy Openness 5 m height - 9m radius 
plot_5open_9m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = prop_5m_9), size = 4) +
  scale_color_gradient2(high = "tan", mid = "lightblue", low = "blue4", 
                        midpoint = 0.375, limits = c(0, 0.75), 
                        name = "Proportion Canopy\nOpenness 5m (%)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_5open_9m



# 9. Percent Canopy Openness 10 m height - 9m radius 
plot_10open_9m <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = growth_spatial_df, aes(color = prop_10m_9), size = 4) +
  scale_color_gradient2(high = "tan", mid = "lightblue", low = "blue4", 
                        midpoint = 0.375, limits = c(0, 0.75), 
                        name = "Proportion Canopy\nOpenness 10m (%)") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text  = element_text(colour="black", size=12),
    axis.text    = element_text(size=1, colour="white"),
    axis.title   = element_text(size=12, colour="black"))

plot_10open_9m



# Are there spatial relationships with canopy height? 

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# Define shapes for species 

# ABAM, ABGR, ALRU, CONU, TABR, THPL, TSHE  
species_shapes <- c(15, 16, 17, 18, 7, 8, 9)


spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
all_canopy_growth$sci_name <- factor(all_canopy_growth$sci_name, levels = spp_order)


# with easting in 9m neighborhoods 
east_height_9m <- ggplot(all_canopy_growth, aes(x = UTM_X, y = X9m_mean, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Easting", y = "Mean Canopy Height (m)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

# No legend here, will add later 

east_height_9m

# The regression model
east_mod <- lm(X9m_mean ~ UTM_X, data = all_canopy_growth)

summary(east_mod)

# Not significant 
# Multiple R-squared:  0.03683,	Adjusted R-squared:  0.02023 
# F-statistic: 2.218 on 1 and 58 DF,  p-value: 0.1418


# with easting in 9m neighborhoods 
north_height_9m <- ggplot(all_canopy_growth, aes(x = UTM_Y, y = X9m_mean, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Northing", y = "Mean Canopy Height (m)") +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

# No legend here, will add later 

north_height_9m

# The regression model
north_mod <- lm(X9m_mean ~ UTM_Y, data = all_canopy_growth)

summary(north_mod)

# Not significant 
# Multiple R-squared:  0.003428,	Adjusted R-squared:  -0.01375 
# F-statistic: 0.1995 on 1 and 58 DF,  p-value: 0.6568


# No relationships, though this is a bit of a simplification. 


# Just plot mean canopy height around each focal tree, by species 

height_plot <- ggplot(all_canopy_growth, aes(x = sci_name, y = X9m_mean, fill = sci_name)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  labs(title = "", y = "Mean Canopy Height (m)", x = "") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black", face = "italic"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black"))

height_plot


# Test for significant differences between species 
aov_height <- aov(X9m_mean ~ Species.x, data = all_canopy_growth)
summary(aov_height)

# P = 0.00322

tuk_height <- TukeyHSD(aov_height)
tuk_height


# 
#               diff         lwr       upr     p adj
# ABGR-ABAM -2.692543 -13.4188537  8.033767 0.9869208
# ALRU-ABAM  1.007210 -10.3697604 12.384180 0.9999640
# CONU-ABAM  5.132452  -4.1568056 14.421709 0.6235181
# TABR-ABAM  9.767041   0.2232362 19.310845 0.0416597  *
# THPL-ABAM  8.473882  -0.8153752 17.763139 0.0955090
# TSHE-ABAM  7.021299  -2.2679582 16.310556 0.2557194
# ALRU-ABGR  3.699753  -8.8779606 16.277467 0.9707725
# CONU-ABGR  7.824995  -2.9013154 18.551305 0.2947138
# TABR-ABGR 12.459584   1.5120895 23.407079 0.0161203 **
# THPL-ABGR 11.166425   0.4401151 21.892736 0.0362085 ** 
# TSHE-ABGR  9.713842  -1.0124680 20.440153 0.1000271
# CONU-ALRU  4.125242  -7.2517284 15.502212 0.9219267
# TABR-ALRU  8.759831  -2.8259084 20.345570 0.2553861
# THPL-ALRU  7.466672  -3.9102980 18.843642 0.4200637
# TSHE-ALRU  6.014089  -5.3628810 17.391059 0.6703624
# TABR-CONU  4.634589  -4.9092154 14.178394 0.7504738
# THPL-CONU  3.341430  -5.9478269 12.630688 0.9246544
# TSHE-CONU  1.888847  -7.4004099 11.178105 0.9957559
# THPL-TABR -1.293159 -10.8369633  8.250646 0.9995682
# TSHE-TABR -2.745742 -12.2895464  6.798063 0.9738116
# TSHE-THPL -1.452583 -10.7418403  7.836674 0.9990200


# TABR has a taller surrounding canopy that ABAM and ABGR, and THPL has a taller 
# surrounding canopy than ABGR 


#################################################################################

#########################################################################
# (3) MODEL RELATIONSHIPS OF NEIGHBORHOOD CANOPY HEIGHT AND TREE GROWTH 
#########################################################################

# 1. Mean canopy height 

# Relationship between focal tree RGR and mean canopy height in 9m radius 
RGR_height_9m <- ggplot(all_canopy_growth, aes(x = X9m_mean, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Mean Canopy Height (m)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

# No legend here, will add later 

RGR_height_9m


# The regression model
RGR_mod1_09 <- lm(mean_RGR ~ X9m_mean, data = all_canopy_growth)

summary(RGR_mod1_09)

# Significant relationship, trees appear to grow faster when the surrounding canopy 
# is shorter. This supports the effect of light competition with neighbors on 
# reducing tree growth 
# Multiple R-squared:  0.1529,	Adjusted R-squared:  0.1383 
# F-statistic: 10.47 on 1 and 58 DF,  p-value: 0.00201



# Relationship between focal tree RGR and canopy height in 20m radius 
RGR_height_20m <- ggplot(all_canopy_growth, aes(x = X20m_mean, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Mean Canopy Height (m)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_height_20m



# The regression model
RGR_mod1_20 <- lm(mean_RGR ~ X20m_mean, data = all_canopy_growth)

summary(RGR_mod1_20)

# No significant relationship, seems like the negative effects of canopy height on focal 
# tree growth is only felt in the closer radius neighborhood 
# Multiple R-squared:  0.0375,	Adjusted R-squared:  0.02091 
# F-statistic:  2.26 on 1 and 58 DF,  p-value: 0.1382



# 2. Max canopy height 

# Question: Does the presence of taller neighbors affect growth more? 

# Relationship between focal tree RGR and max canopy height in 9m radius 
RGR_max_height_9m <- ggplot(all_canopy_growth, aes(x = X9m_max, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Maximum Canopy Height (m)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_max_height_9m


# The regression model
RGR_mod2_09 <- lm(mean_RGR ~ X9m_max, data = all_canopy_growth)

summary(RGR_mod2_09)

# No Significant relationship
# Multiple R-squared:  0.02009,	Adjusted R-squared:  0.003196 
# F-statistic: 1.189 on 1 and 58 DF,  p-value: 0.28



# Relationship between focal tree RGR and max canopy height in 20m radius 
RGR_max_height_20m <- ggplot(all_canopy_growth, aes(x = X20m_max, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Maximum Canopy Height (m)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_max_height_20m



# The regression model
RGR_mod2_20 <- lm(mean_RGR ~ X20m_max, data = all_canopy_growth)

summary(RGR_mod2_20)

# No significant relationship
# Multiple R-squared:  0.004034,	Adjusted R-squared:  -0.01314 
# F-statistic: 0.2349 on 1 and 58 DF,  p-value: 0.6297


# 3. Canopy Openness - only doing for 20 m radius because 9m can be the radius of a single 
# tree's crown 


# Relationship between focal tree RGR and 2m canopy openness in 20m radius 
RGR_2open_20m <- ggplot(all_canopy_growth, aes(x = pr_2_20, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Proportion Canopy Openness 2m (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_2open_20m



# The regression model
RGR_mod3_20 <- lm(mean_RGR ~ pr_2_20, data = all_canopy_growth)

summary(RGR_mod3_20)

# No significant relationship
# Multiple R-squared:  0.00114,	Adjusted R-squared:  -0.01608 
# F-statistic: 0.06622 on 1 and 58 DF,  p-value: 0.7978


# Relationship between focal tree RGR and 5m canopy openness in 20m radius 
RGR_5open_20m <- ggplot(all_canopy_growth, aes(x = pr_5_20, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Proportion Canopy Openness 5m (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_5open_20m



# The regression model
RGR_mod4_20 <- lm(mean_RGR ~ pr_5_20, data = all_canopy_growth)

summary(RGR_mod4_20)

# No significant relationship
# Multiple R-squared:  0.0283,	Adjusted R-squared:  0.01154 
# F-statistic: 1.689 on 1 and 58 DF,  p-value: 0.1989



# Relationship between focal tree RGR and 10m canopy openness in 20m radius 
RGR_10open_20m <- ggplot(all_canopy_growth, aes(x = p_10_20, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Proportion Canopy Openness 10m (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_10open_20m



# The regression model
RGR_mod5_20 <- lm(mean_RGR ~ p_10_20, data = all_canopy_growth)

summary(RGR_mod5_20)

# Significant relationship
# Multiple R-squared:  0.07291,	Adjusted R-squared:  0.05692 
# F-statistic: 4.561 on 1 and 58 DF,  p-value: 0.03694




##### Then repeating for the 9m radius for canopy openness 

# Relationship between focal tree RGR and 2m canopy openness in 9m radius 
RGR_2open_9m <- ggplot(all_canopy_growth, aes(x = prop_2m_9, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 2) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Proportion Canopy Openness 2m (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_2open_9m



# The regression model
RGR_mod6_09 <- lm(mean_RGR ~ prop_2m_9, data = all_canopy_growth)

summary(RGR_mod6_09)

# No significant relationship
# Multiple R-squared:  0.008002,	Adjusted R-squared:  -0.009101 
# F-statistic: 0.4679 on 1 and 58 DF,  p-value: 0.4967


# Relationship between focal tree RGR and 5m canopy openness in 9m radius 
RGR_5open_9m <- ggplot(all_canopy_growth, aes(x = prop_5m_9, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Proportion Canopy Openness 5m (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_5open_9m



# The regression model
RGR_mod7_09 <- lm(mean_RGR ~ prop_5m_9, data = all_canopy_growth)

summary(RGR_mod7_09)

# Significant relationship
# Multiple R-squared:  0.1238,	Adjusted R-squared:  0.1087 
# F-statistic: 8.196 on 1 and 58 DF,  p-value: 0.005833



# Relationship between focal tree RGR and 10m canopy openness in 9m radius 
RGR_10open_9m <- ggplot(all_canopy_growth, aes(x = prop_10m_9, y = mean_RGR, colour = sci_name)) +
  geom_point(alpha = 1, cex = 2.5, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = 1) +
  theme_bw() +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(x = "Proportion Canopy Openness 10m (%)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

RGR_10open_9m



# The regression model
RGR_mod8_09 <- lm(mean_RGR ~ prop_10m_9, data = all_canopy_growth)

summary(RGR_mod8_09)

# Significant relationship
# Multiple R-squared:  0.2404,	Adjusted R-squared:  0.2273 
# F-statistic: 18.35 on 1 and 58 DF,  p-value: 7e-05



############################################## -- 
# (4) GATHER SIGNIFICANT RESULTS OF INTEREST  
############################################## -- 

# Gather up the most interesting plots and organize the plots to show the canopy height variation 
# across the plot, and then to summarize the relationships to focal tree growth


# Spatial plots with canopy openness for 9m radius 
spatial_plots <- plot_grid(plot_mean_9m, plot_2open_9m, plot_5open_9m, plot_10open_9m,
                           ncol = 1, nrow = 4, labels = c('(a)', '(b)', '(c)', '(d)'))

spatial_plots


# Save figure
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/canopy_spatial_plots.png", 
       plot = spatial_plots, width = 8, height = 13, units = "in", dpi = 300)


# Organize and number plots for the growth relationships with mean canopy height 
canopy_plots <- plot_grid(RGR_height_9m, RGR_2open_9m, RGR_5open_9m, RGR_10open_9m,
                          ncol = 2, nrow = 2, labels = c('(a)', '(b)', '(c)', '(d)', '(e)'))

canopy_plots


# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/canopy_RGR_plots.png", 
       plot = canopy_plots, width = 8, height = 7, units = "in", dpi = 300)


## -- END -- ## 
