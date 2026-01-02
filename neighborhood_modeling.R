# -----------------------------------------------------------------------------#
# Neighborhood Modeling of focal trees in WFDP
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
#  neighborhood for each studied focal tree. Incorporate DBH, species data, and #
#  mycorrhizal associations into the neighborhoods and use distance to generate #
# an overall crowding metric for each focal tree.                               #
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)


# Import neighbor trees file for 9m radius neighborhoods 
neighbors <- read.csv("./neighbor_trees.csv") #2,147 neighbors in total 

# columns related to the focal tree are generally identified with a 'focal_' prefix, 
# while columns related to the neighbors have a 'neigh_' prefix 


# Import neighbor trees file for 20m radius neighborhoods 
neighbors_20 <- read.csv("./neighbor_trees_20.csv") # 9,288 neighbors in total 


# There are two trees that are on the northern edge of the plot, and part of their neighborhood 
# radii reaches outside of the plot, so it contains trees we can't totally account for. These 
# trees will be retained for other analyses, but need to be removed for these neighborhood 
# analyses:
  # ALRU in Q38 - 38-0935
  # THPL in Q29 - 29-0766

# Make cell a factor in both df's, because there are no other focal trees that share the same cell  

neighbors$focal_cell <- as.factor(neighbors$focal_cell)
neighbors_20$focal_cell <- as.factor(neighbors_20$focal_cell)

# Remove these cells from each 
neighbors <- neighbors %>% filter(focal_cell != "Q38") %>% droplevels()
neighbors <- neighbors %>% filter(focal_cell != "Q29") %>% droplevels()
# removes 58 neighbor records 

neighbors_20 <- neighbors_20 %>% filter(focal_cell != "Q38") %>% droplevels()
neighbors_20 <- neighbors_20 %>% filter(focal_cell != "Q29") %>% droplevels()
# removes 153 neighbor records 


# Load file of all of the coordinate positions of each tree in WFDP
trees_geo <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_trees_georef.shp")

# Mutate the geometry column to extract the UMTM_X and UTM_Y data 
trees_geo <- trees_geo %>%
  mutate(
    UTM_X_pred = st_coordinates(.)[, 1],
    UTM_Y_pred = st_coordinates(.)[, 2])
 

# Load WFDP plot boundary (polygon shapefile)
wfdp_poly <- st_read("~/Dropbox/WSU/WFDP_Chapter_3_Project/GIS/WFDP_boundary_utm.shp")

############################ --
# (2) PROCESS NEIGHBOR DATA
############################ -- 

#### In QGIS, I found the neighbors of each focal tree within a 9 m radius using 
# the Buffer tool to set the radius, and Join Attributes by Location tool to 
# pick just trees that intersected with the radius polygon of each focal tree 

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#9b5fe0", "#16a4d8", "#60dbe8", "#8bd346","#efdf48", "#f9a52F", "#d64e12")


### Explore 9m data ### -- 

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
#       36              17.2               9                86


# plot count of neighbors 
neigh_count <- ggplot(neighbor_summary, aes(x = n_neighbors)) +
  geom_histogram(binwidth = 2, color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Neighbor Counts per Focal Tree",
    x = "Number of Neighbors (within 9 m)",
    y = "Number of Focal Trees"
  )

neigh_count


# Check neighbors by focal species - just crowding, nothing separate 
neigh_spp <- ggplot(neighbor_summary, aes(x = focal_species, y = n_neighbors, fill = focal_species)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values=all_hosts, 
                    name="Host Species",
                    breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                    labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(
    title = "",
    x = "",
    y = "Number of Neighbors"
  ) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

neigh_spp 

# Test for significant differences between species 
aov_neigh <- aov(n_neighbors ~ focal_species, data = neighbor_summary)
summary(aov_neigh)
# p = 0.19

tuk_neigh <- TukeyHSD(aov_neigh)
tuk_neigh

# No significant differences between species


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

# No relationship here, p = 0.548


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

# No relationship here, p = 0.399 



### Explore 20m data ### -- 

# summarize number of neighbors per focal tree
neighbor_summary_20 <- neighbors_20 %>%
  group_by(focal_stem_tag) %>%
  summarise(
    focal_species = first(focal_species),
    focal_DBH = first(focal_DBH),
    n_neighbors = n(),
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    sd_neighbor_DBH = sd(neigh_DBH, na.rm = TRUE))

# overall mean and range
summary_stats_20 <- neighbor_summary_20 %>%
  summarise(
    mean_neighbors = mean(n_neighbors),
    sd_neighbors = sd(n_neighbors),
    min_neighbors = min(n_neighbors),
    max_neighbors = max(n_neighbors))

neighbor_summary_20
summary_stats_20

# mean_neighbors     sd_neighbors     min_neighbors     max_neighbors
#       158              51.1               70                287


# plot count of neighbors 
neigh_count_20 <- ggplot(neighbor_summary_20, aes(x = n_neighbors)) +
  geom_histogram(binwidth = 2, color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Neighbor Counts per Focal Tree",
    x = "Number of Neighbors (within 20 m)",
    y = "Number of Focal Trees"
  )

neigh_count_20


# Check neighbors by focal species 
neigh_spp_20 <- ggplot(neighbor_summary_20, aes(x = focal_species, y = n_neighbors, fill = focal_species)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values=all_hosts, 
                    name="Host Species",
                    breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                    labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(
    title = "",
    x = "",
    y = "Number of Neighbors"
  ) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

neigh_spp_20 

# Test for significant differences between species 
aov_neigh_20 <- aov(n_neighbors ~ focal_species, data = neighbor_summary_20)
summary(aov_neigh_20)
# p = 0.177

tuk_neigh_20 <- TukeyHSD(aov_neigh_20)
tuk_neigh_20

# No significant differences between species


# Check for relationship between tree size and number of neighbors 
tree_size_20 <- ggplot(neighbor_summary_20, aes(x = focal_DBH, y = n_neighbors)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    title = "Relationship Between Focal Tree DBH and Number of Neighbors",
    x = "Focal Tree DBH (cm)",
    y = "Number of Neighbors"
  )

tree_size_20 

size_mod_20 <- lm(focal_DBH ~ n_neighbors, data = neighbor_summary_20)

summary(size_mod_20)

# No relationship here, p = 0.883 


# Check for relationship between focal tree size and neighbor size 
neigh_size_20 <- ggplot(neighbor_summary_20, aes(x = focal_DBH, y = mean_neighbor_DBH)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    title = "Relationship Between Focal Tree Size and Average Neighbor Size",
    x = "Focal Tree DBH (cm)",
    y = "Mean Neighbor DBH (cm)")

neigh_size_20

size_mod2_20 <- lm(focal_DBH ~ mean_neighbor_DBH, data = neighbor_summary_20)

summary(size_mod2_20)

# No relationship here, p = 0.749 


### Overall consistent results for both 9m and 20m neighborhoods, 
# with no significant differences or relationships

#################################################################################

##################################### -- 
# (3) CALCULATE NEIGHBORHOOD METRICS
##################################### -- 

# Need to calculate distances between the focal tree to each of its neighbors

# The trees_geo df has the predicted and corrected UTMs for each of the trees, so I can 
# subset this, then merge it with the neighbors data 

### Explore 9m data ### -- 

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
    crowding = (neigh_DBH^2) / (distance^2)
  )


# Summarize crowding 
crowding_summary <- neighbors %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
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

# Significant relationship, p = 0.018, Multiple R-squared:  0.09544,	Adjusted R-squared:  0.07928 
# bigger trees are more crowded (though this study was focused on 
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
  theme(axis.text.x = element_text(hjust = 1))

spp_crowd 


# Test for significant differences between species 
aov_crowd <- aov(crowding_index ~ focal_species, data = crowding_summary)
summary(aov_crowd)
# p = 0.667

tuk_crowd <- TukeyHSD(aov_crowd)
tuk_crowd

# No significant differences between species


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
    title = "Spatial Distribution of 9m Neighborhood Crowding")

space_crowd 

# Doesn't seem like a consistent trend with crowding across the plot, except for maybe a bit 
# higher on the eastern side 


### Explore 20m data ### -- 

neighbors_20 <- merge(neighbors_20, geo_sub, by = "neigh_stem_tag")


# Each neighbor tree now has UTM coordinates, so can calculate it's distance to the focal tree


# Make sure coordinates are numeric
neighbors_20 <- neighbors_20 %>%
  mutate(UTM_X_pred = as.numeric(UTM_X_pred),
         UTM_Y_pred = as.numeric(UTM_Y_pred))

# Make a few things factors 
neighbors_20$focal_cell <- as.factor(neighbors_20$focal_cell)
neighbors_20$focal_species <- as.factor(neighbors_20$focal_species)
neighbors_20$neigh_species <- as.factor(neighbors_20$neigh_species)


# Calculate distance and then calculate the crowding index per neighbor
neighbors_20 <- neighbors_20 %>%
  mutate(
    # Step 1: Euclidean distance between focal and neighbor trees
    distance = sqrt((UTM_X_pred - UTM_X)^2 + (UTM_Y_pred - UTM_Y)^2),
    
    # Step 2: Neighborhood crowding index contribution per neighbor
    crowding = (neigh_DBH^2) / (distance^2)
  )


# Summarize crowding 
crowding_summary_20 <- neighbors_20 %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")

# A high crowding index means the focal tree is surrounded by many large, nearby neighbors 
# A low crowding index means few, small, or distant neighbors 

### Visualize ### 

# Relationship between focal tree size and crowding 
size_crowd_20 <- ggplot(crowding_summary_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding")

size_crowd_20 


crowd_mod1_20 <- lm(focal_DBH ~ crowding_index, data = crowding_summary_20)

summary(crowd_mod1_20)

# Significant relationship, p = 0.018, Multiple R-squared:  0.09608,	Adjusted R-squared:  0.07994 
# bigger trees are more crowded (though this study was focused on 
# trees 10 - 20 cm DBH so this is only a small snapshot)

# Crowding across different species 
spp_crowd_20 <- ggplot(crowding_summary_20, aes(x = focal_species, y = crowding_index, fill = focal_species)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    x = "Focal Species",
    y = "Neighborhood Crowding Index",
    title = "Variation in Crowding by Focal Species"
  ) +
  theme(axis.text.x = element_text(hjust = 1))

spp_crowd_20 


# Test for significant differences between species 
aov_crowd_20 <- aov(crowding_index ~ focal_species, data = crowding_summary_20)
summary(aov_crowd_20)
# p = 0.669

tuk_crowd_20 <- TukeyHSD(aov_crowd_20)
tuk_crowd_20

# No significant differences between species


# Merge back in UTM info to be able to plot spatially 
geo_sub_focal_20 <- dplyr::select(neighbors_20, focal_stem_tag, UTM_X, UTM_Y) %>% unique()

crowding_summary_20 <- merge(crowding_summary_20, geo_sub_focal_20, by = "focal_stem_tag")


# Spatial map of crowding across the plot 
space_crowd_20 <- ggplot(crowding_summary_20, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of 20m Neighborhood Crowding")

space_crowd_20 

# Doesn't seem like a consistent trend with crowding across the plot 


### Overall consistent results for both 9m and 20m neighborhoods, 
# even the spatial map of crowding across the plot looks identical, even though 
# the crowding values for the neighborhoods are certainly different.

#################################################################################

####################################################### -- 
# (4) CALCULATE CON- AND HETEROSPECIFIC NEIGHBORHOODS  
####################################################### --

# Want to separate out the neighborhoods for all conspecific trees to the focal 
# tree, and then all heterospecific trees to the focal tree to see how these 
# differ. 

### Explore 9m data ### -- 

# 'neighbors' df has 'focal_species' and 'neigh_species' columns 

neighborhoods <- neighbors %>%
  mutate(type = if_else(as.character(neigh_species) == as.character(focal_species),
                        "conspecific", "heterospecific"))


# Calculate distance and then calculate the crowding index per neighbor
neighborhoods <- neighborhoods %>%
  mutate(
    # Step 1: Euclidean distance between focal and neighbor trees
    distance = sqrt((UTM_X_pred - UTM_X)^2 + (UTM_Y_pred - UTM_Y)^2),
    
    # Step 2: Neighborhood crowding index contribution per neighbor
    crowding = (neigh_DBH^2) / (distance^2)
  )

# Summarize values for the conspecific and heterospecific neighborhoods 
# Split into separate dataframes 
neighborhoods$type <- as.factor(neighborhoods$type)

conspecific <- neighborhoods %>% 
  filter(type == "conspecific")

heterospecific <- neighborhoods %>% 
  filter(type == "heterospecific")


# Get some summary stats for these neighborhoods 

summary(conspecific$neigh_DBH)

# Calculate Mean
mean_val <- mean(conspecific$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(conspecific$neigh_DBH)

# Function to calculate Standard Error
se <- function(x) {
  sd(x) / sqrt(length(x))
}

# Calculate Standard Error
se_val <- se(conspecific$neigh_DBH)

# conspecific: count = 300 trees, mean = 14.571, SE = 1.11


summary(heterospecific$neigh_DBH)

# Calculate Mean
mean_val <- mean(heterospecific$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(heterospecific$neigh_DBH)

# Calculate Standard Error
se_val <- se(heterospecific$neigh_DBH)

# heterospecific: count = 1,789 trees, mean = 9.75, SE = 0.44 


# These are different sizes because they have the same number of focal trees, and 
# the same number of focal cells, but the number of neighbors is very different as 
# there are many more heterospecific neighbors than there are conspecific 

con_summary <- conspecific %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")

het_summary <- heterospecific %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")


## Visualize 

# Relationship between focal tree size and crowding 
size_crowd_con <- ggplot(con_summary, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Conspecific Neighborhoods")

size_crowd_con 


crowd_mod1_con <- lm(focal_DBH ~ crowding_index, data = con_summary)

summary(crowd_mod1_con)

# p = 0.044 

# Spatial map of crowding across the plot 
space_crowd_con <- ggplot(con_summary, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Conspecific Neighborhood Crowding")

space_crowd_con 


size_crowd_het <- ggplot(het_summary, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Heterospecific Neighborhoods")

size_crowd_het 


crowd_mod1_het <- lm(focal_DBH ~ crowding_index, data = het_summary)

summary(crowd_mod1_het) # Not significant 

# This has one high outlier of a THPL in H30 that is surrounded by four large diameter TSHE 

# Get rid of this one and see how things change 

# Remove this focal tree and test how it changes the relationship  
het_summary_no_outlier <- het_summary %>% filter(focal_cell != "H30") %>% droplevels()

# Replot and retest 


size_crowd_het_no_outlier <- ggplot(het_summary_no_outlier, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Heterospecific Neighborhoods")

size_crowd_het_no_outlier 


crowd_mod1_het_no_outlier <- lm(focal_DBH ~ crowding_index, data = het_summary_no_outlier)

summary(crowd_mod1_het_no_outlier) 

# When the one high outlier is removed, there is still no significant relationship between the 
# focal tree size and the neighborhood crowding index 


# Spatial map of crowding across the plot 
space_crowd_het <- ggplot(het_summary, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Heterospecific Neighborhood Crowding")

space_crowd_het 


### Explore 20m data ### -- 

# 'neighbors_20' df has 'focal_species' and 'neigh_species' columns 

neighborhoods_20 <- neighbors_20 %>%
  mutate(type = if_else(as.character(neigh_species) == as.character(focal_species),
                        "conspecific", "heterospecific"))


# Calculate distance and then calculate the crowding index per neighbor
neighborhoods_20 <- neighborhoods_20 %>%
  mutate(
    # Step 1: Euclidean distance between focal and neighbor trees
    distance = sqrt((UTM_X_pred - UTM_X)^2 + (UTM_Y_pred - UTM_Y)^2),
    
    # Step 2: Neighborhood crowding index contribution per neighbor
    crowding = (neigh_DBH^2) / (distance^2)
  )

# Summarize values for the conspecific and heterospecific neighborhoods 
# Split into separate dataframes 
neighborhoods_20$type <- as.factor(neighborhoods_20$type)

conspecific_20 <- neighborhoods_20 %>% 
  filter(type == "conspecific")

heterospecific_20 <- neighborhoods_20 %>% 
  filter(type == "heterospecific")


# Get some summary stats for these neighborhoods 

summary(conspecific_20$neigh_DBH)

# Calculate Mean
mean_val <- mean(conspecific_20$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(conspecific_20$neigh_DBH)

# Calculate Standard Error
se_val <- se(conspecific_20$neigh_DBH)

# conspecific: count = 895 trees, mean = 17.89, SE = 0.77


summary(heterospecific_20$neigh_DBH)

# Calculate Mean
mean_val <- mean(heterospecific_20$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(heterospecific_20$neigh_DBH)

# Calculate Standard Error
se_val <- se(heterospecific_20$neigh_DBH)

# heterospecific: count = 8,240 trees, mean = 10.78, SE = 0.22



# These are different sizes because they have the same number of focal trees, and 
# the same number of focal cells, but the number of neighbors is very different as 
# there are many more heterospecific neighbors than there are conspecific 

con_summary_20 <- conspecific_20 %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")

het_summary_20 <- heterospecific_20 %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")


## Visualize 

# Relationship between focal tree size and crowding 
size_crowd_con_20 <- ggplot(con_summary_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in 20m Conspecific Neighborhoods")

size_crowd_con_20 


crowd_mod1_con_20 <- lm(focal_DBH ~ crowding_index, data = con_summary_20)

summary(crowd_mod1_con_20)

# p = 0.043 

# Spatial map of crowding across the plot 
space_crowd_con_20 <- ggplot(con_summary_20, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Conspecific 20m Neighborhood Crowding")

space_crowd_con_20 


size_crowd_het_20 <- ggplot(het_summary_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Heterospecific 20m Neighborhoods")

size_crowd_het_20 


crowd_mod1_het_20 <- lm(focal_DBH ~ crowding_index, data = het_summary_20)

summary(crowd_mod1_het_20) # Not significant 

# This has one high outlier of a THPL in H30 that is surrounded by four large diameter TSHE, 
# same issue as in the 9m neighborhoods 

# Get rid of this one and see how things change 

# Remove this focal tree and test how it changes the relationship  
het_summary_no_outlier_20 <- het_summary_20 %>% filter(focal_cell != "H30") %>% droplevels()

# Replot and retest 


size_crowd_het_no_outlier_20 <- ggplot(het_summary_no_outlier_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Heterospecific 20m Neighborhoods")

size_crowd_het_no_outlier_20 


crowd_mod1_het_no_outlier_20 <- lm(focal_DBH ~ crowding_index, data = het_summary_no_outlier_20)

summary(crowd_mod1_het_no_outlier_20) 

# When the one high outlier is removed, there is still no significant relationship between the 
# focal tree size and the neighborhood crowding index 


# Spatial map of crowding across the plot 
space_crowd_het_20 <- ggplot(het_summary_20, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Heterospecific 20m Neighborhood Crowding")

space_crowd_het_20 


### Overall consistent results for both 9m and 20m neighborhoods, 
# there are many more heterospecific neighbors than there are conspecific, and that 
# pattern holds for both neighborhood sizes. There is one high THPL outlier in the 
# heterospecific results, but it doesn't change the significance of the results when 
# it is removed. 

#################################################################################

########################################################################### -- 
# (5) CALCULATE SAME- AND DIFFERENT MYCORRHIZAL ASSOCIATION NEIGHBORHOODS  
########################################################################### --

# Want to separate out the neighborhoods for all trees that have the same 
# mycorrhizal association as the focal tree, and then all trees that have a 
# different mycorrhizal association to the focal tree to see how these 
# differ. I'm saying tree here, but there are some neighbors that are shrubs, and 
# I am leaving those in the analyses for now. The possible mycorrhizal associations 
# are AM, DUAL (AM and EM), EM, and ERM (ericoid mycorrhizal)

### Explore 9m data ### -- 

# Add to neighborhoods df up above a column that determines in the focal tree 
# and neighbor tree have the same or different mycorrhizal association 

# For dual hosts, I am going to mark that they have the same if the neighbor is either 
# AM, EM, or DUAL
neighborhoods <- neighbors %>%
  mutate(myco_match = case_when(focal_myco == neigh_myco ~ "same",
      focal_myco == "DUAL" & neigh_myco %in% c("AM", "EM") ~ "same",
      neigh_myco == "DUAL" & focal_myco %in% c("AM", "EM") ~ "same",
      TRUE ~ "different"))

# Check assignments 
neighborhoods %>%
  count(focal_myco, neigh_myco, myco_match) %>%
  arrange(focal_myco, neigh_myco)

# This looks good!


# Summarize crowding values that were already calculated for the same and different 
# mycorrhizal neighborhoods 
# Split into separate dataframes 
neighborhoods$myco_match <- as.factor(neighborhoods$myco_match)

same_myco_assoc <- neighborhoods %>% 
  filter(myco_match == "same")

diff_myco_assoc <- neighborhoods %>% 
  filter(myco_match == "different")


# Get some summary stats for these neighborhoods 

summary(same_myco_assoc$neigh_DBH)

# Calculate Mean
mean_val <- mean(same_myco_assoc$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(same_myco_assoc$neigh_DBH)

# Calculate Standard Error
se_val <- se(same_myco_assoc$neigh_DBH)

# same myco assoc: count = 1,385 trees, mean = 11.00, SE = 0.47


summary(diff_myco_assoc$neigh_DBH)

# Calculate Mean
mean_val <- mean(diff_myco_assoc$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(diff_myco_assoc$neigh_DBH)

# Calculate Standard Error
se_val <- se(diff_myco_assoc$neigh_DBH)

# diff myco assoc: count = 704 trees, mean = 9.34, SE = 0.80


# These are different sizes because they have the same number of focal trees, and 
# the same number of focal cells, but the number of neighbors is very different as 
# there are many more neighbors with the same myco association as the focal tree, than 
# there are with different associations. The different dataset will also contain all of
# the ericoid mycorrhizal shrubs, because these would be considered different than 
# all other mycorrhizal types. 

same_myco_summary <- same_myco_assoc %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")

diff_myco_summary <- diff_myco_assoc %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")


## Visualize 

# Relationship between focal tree size and crowding 
size_crowd_same_myco <- ggplot(same_myco_summary, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Same Mycorrizal Neighborhoods")

size_crowd_same_myco 


crowd_mod1_same_myco <- lm(focal_DBH ~ crowding_index, data = same_myco_summary)

summary(crowd_mod1_same_myco)

# p = 0.043 

# Spatial map of crowding across the plot 
space_crowd_same_myco <- ggplot(same_myco_summary, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Same Mycorrhizal Neighborhood Crowding")

space_crowd_same_myco 


size_crowd_diff_myco <- ggplot(diff_myco_summary, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Different Mycorrhizal Neighborhoods")

size_crowd_diff_myco 


crowd_mod1_diff_myco <- lm(focal_DBH ~ crowding_index, data = diff_myco_summary)

summary(crowd_mod1_diff_myco) # Not significant 

# This has the same high outlier of a THPL in H30 that is surrounded by four large diameter TSHE 

# Get rid of this one and see how things change 

# Remove this focal tree and test how it changes the relationship  
diff_myco_summary_no_outlier <- diff_myco_summary %>% filter(focal_cell != "H30") %>% droplevels()

# Replot and retest 
size_crowd_diff_myco_no_outlier <- ggplot(diff_myco_summary_no_outlier, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Different Mycorrhizal Neighborhoods")

size_crowd_diff_myco_no_outlier


crowd_mod1_diff_myco_no_outlier <- lm(focal_DBH ~ crowding_index, data = diff_myco_summary_no_outlier)

summary(crowd_mod1_diff_myco_no_outlier) 

# When the one high outlier is removed, this becomes a significant relationship 
# Multiple R-squared:  0.2408,	Adjusted R-squared:  0.2243, F-statistic: 14.59 on 1 and 46 DF,  p-value: 0.0003989
# Smaller focal trees have more crowding in different mycorrhizal neighborhoods 


# Spatial map of crowding across the plot 
space_crowd_diff_myco <- ggplot(diff_myco_summary, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Different Mycorrhizal Neighborhood Crowding")

space_crowd_diff_myco


# Map this with the high outlier removed too 
space_crowd_diff_myco_no_outlier <- ggplot(diff_myco_summary_no_outlier, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Different Mycorrhizal Neighborhood Crowding, with outlier removed")

space_crowd_diff_myco_no_outlier

# This shows much more variation. Here the crowding index is generally very low, because a lot of the 
# different mycorrhizal neighbors to the focal trees are the small ericoid shrubs. In general we see 
# more crowding in the wetter part of the plot. 


### Explore 20m data ### -- 

# Add to neighborhoods_20 df up above a column that determines in the focal tree 
# and neighbor tree have the same or different mycorrhizal association 

# For dual hosts, I am going to mark that they have the same if the neighbor is either 
# AM, EM, or DUAL
neighborhoods_20 <- neighbors_20 %>%
  mutate(myco_match = case_when(focal_myco == neigh_myco ~ "same",
                                focal_myco == "DUAL" & neigh_myco %in% c("AM", "EM") ~ "same",
                                neigh_myco == "DUAL" & focal_myco %in% c("AM", "EM") ~ "same",
                                TRUE ~ "different"))

# Check assignments 
neighborhoods_20 %>%
  count(focal_myco, neigh_myco, myco_match) %>%
  arrange(focal_myco, neigh_myco)

# This looks good!


# Summarize crowding values that were already calculated for the same and different 
# mycorrhizal neighborhoods 
# Split into separate dataframes 
neighborhoods_20$myco_match <- as.factor(neighborhoods_20$myco_match)

same_myco_assoc_20 <- neighborhoods_20 %>% 
  filter(myco_match == "same")

diff_myco_assoc_20 <- neighborhoods_20 %>% 
  filter(myco_match == "different")



# Get some summary stats for these neighborhoods 

summary(same_myco_assoc_20$neigh_DBH)

# Calculate Mean
mean_val <- mean(same_myco_assoc_20$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(same_myco_assoc_20$neigh_DBH)

# Calculate Standard Error
se_val <- se(same_myco_assoc_20$neigh_DBH)

# same myco assoc: count = 5,723 trees, mean = 12.95, SE = 0.28


summary(diff_myco_assoc_20$neigh_DBH)

# Calculate Mean
mean_val <- mean(diff_myco_assoc_20$neigh_DBH)

# Calculate Standard Deviation
sd_val <- sd(diff_myco_assoc_20$neigh_DBH)

# Calculate Standard Error
se_val <- se(diff_myco_assoc_20$neigh_DBH)

# diff myco assoc: count = 3,412 trees, mean = 9.01, SE = 0.34


# These are different sizes because they have the same number of focal trees, and 
# the same number of focal cells, but the number of neighbors is very different as 
# there are many more neighbors with the same myco association as the focal tree, than 
# there are with different associations. The different dataset will also contain all of
# the ericoid mycorrhizal shrubs, because these would be considered different than 
# all other mycorrhizal types. 

same_myco_summary_20 <- same_myco_assoc_20 %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")

diff_myco_summary_20 <- diff_myco_assoc_20 %>%
  group_by(focal_stem_tag, focal_species, focal_DBH, focal_cell, UTM_X, UTM_Y) %>%
  summarise(
    mean_neighbor_DBH = mean(neigh_DBH, na.rm = TRUE),
    num_neighbors = n(),
    crowding_index = sum(crowding, na.rm = TRUE),
    .groups = "drop")


## Visualize 

# Relationship between focal tree size and crowding 
size_crowd_same_myco_20 <- ggplot(same_myco_summary_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Same Mycorrizal 20m Neighborhoods")

size_crowd_same_myco_20 


crowd_mod1_same_myco_20 <- lm(focal_DBH ~ crowding_index, data = same_myco_summary_20)

summary(crowd_mod1_same_myco_20)

# p = 0.038

# Spatial map of crowding across the plot 
space_crowd_same_myco_20 <- ggplot(same_myco_summary_20, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Same Mycorrhizal 20m Neighborhood Crowding")

space_crowd_same_myco_20 


size_crowd_diff_myco_20 <- ggplot(diff_myco_summary_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Different Mycorrhizal 20m Neighborhoods")

size_crowd_diff_myco_20 


crowd_mod1_diff_myco_20 <- lm(focal_DBH ~ crowding_index, data = diff_myco_summary_20)

summary(crowd_mod1_diff_myco_20) # Not significant 

# This has the same high outlier of a THPL in H30 that is surrounded by four large diameter TSHE 

# Get rid of this one and see how things change 

# Remove this focal tree and test how it changes the relationship  
diff_myco_summary_no_outlier_20 <- diff_myco_summary_20 %>% filter(focal_cell != "H30") %>% droplevels()

# Replot and retest 
size_crowd_diff_myco_no_outlier_20 <- ggplot(diff_myco_summary_no_outlier_20, aes(x = crowding_index, y = focal_DBH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(
    x = "Neighborhood Crowding Index (Σ(DBH² / distance²))",
    y = "Focal Tree DBH (cm)",
    title = "Relationship between Focal Tree Size and Local Crowding in Different Mycorrhizal 20m Neighborhoods")

size_crowd_diff_myco_no_outlier_20


crowd_mod1_diff_myco_no_outlier_20 <- lm(focal_DBH ~ crowding_index, data = diff_myco_summary_no_outlier_20)

summary(crowd_mod1_diff_myco_no_outlier_20) 

# When the one high outlier is removed, this becomes a significant relationship 
# Multiple R-squared:   0.26,	Adjusted R-squared:  0.2454, F-statistic: 17.91 on 1 and 51 DF,  p-value: 9.633e-05
# Smaller focal trees have more crowding in different mycorrhizal neighborhoods 


# Spatial map of crowding across the plot 
space_crowd_diff_myco_20 <- ggplot(diff_myco_summary_20, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Different Mycorrhizal 20m Neighborhood Crowding")

space_crowd_diff_myco_20


# Map this with the high outlier removed too 
space_crowd_diff_myco_no_outlier_20 <- ggplot(diff_myco_summary_no_outlier_20, aes(x = UTM_X, y = UTM_Y, color = crowding_index)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    x = "UTM Easting",
    y = "UTM Northing",
    color = "Crowding Index",
    title = "Spatial Distribution of Different Mycorrhizal 20m Neighborhood Crowding, with outlier removed")

space_crowd_diff_myco_no_outlier_20

# This shows much more variation. Here the crowding index is generally very low, because a lot of the 
# different mycorrhizal neighbors to the focal trees are the small ericoid shrubs. In general we see 
# more crowding in the wetter part of the plot. 


### Overall consistent results for both 9m and 20m neighborhoods, 
# there are many more neighbors with the same mycorrhizal associations than there are different associations, 
# and that pattern holds for both neighborhood sizes. There is one high THPL outlier in the 
# different myco results, same as in the heterospecific neighbors results, and it changes the significance 
# of the results when it is removed, which is different than the con-het neighborhoods comparison 

# This outlier pops up again because the THPL is being crowded by TSHEs that are EM while THPL is AM. 
# It will be interesting to see if the THPL has lower growth than any other THPLs that aren't being 
# crowded as much. 

#################################################################################

################################################################ -- 
# (6) SAVE NEIGHBORHOOD SUMMARY FILES FOR DEMOGRAPHIC ANALYSES  
################################################################ --

# There are a lot of summary files because these were made for both the 9 and 20m 
# radius neighborhoods 

# save crowding_summary file to explore relationships with focal tree growth 
# 9m 
write.csv(crowding_summary, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/combined_crowd_summary_9m.csv")

# 20m 
write.csv(crowding_summary_20, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/combined_crowd_summary_20m.csv")


# save con_summary for conspecific neighborhoods 
# 9m 
write.csv(con_summary, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/conspecific_crowd_summary_9m.csv")

# 20m 
write.csv(con_summary_20, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/conspecific_crowd_summary_20m.csv")


# save het_summary for heterospecific neighborhoods 
# 9m 
write.csv(het_summary, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/heterospecific_crowd_summary_9m.csv")

# 20m 
write.csv(het_summary_20, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/heterospecific_crowd_summary_20m.csv")


# save same_myco_summary for same mycorrhizal neighborhoods 
# 9m 
write.csv(same_myco_summary, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/same_myco_crowd_summary_9m.csv")

# 20m 
write.csv(same_myco_summary_20, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/same_myco_crowd_summary_20m.csv")


# save diff_myco_summary for different mycorrhizal neighborhoods 
# 9m 
write.csv(diff_myco_summary, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/diff_myco_crowd_summary_9m.csv")

# 20m 
write.csv(diff_myco_summary_20, 
          "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Neighborhood_files/diff_myco_crowd_summary_20m.csv")



## -- END -- ## 

