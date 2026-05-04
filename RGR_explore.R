# -----------------------------------------------------------------------------#
# Explore Tree RGR in WFDP
# Original Author: L. McKinley Nevins 
# March 19, 2026
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     sf v 1.0.19
#                     raster v 3.6.32
#                     cowplot v 1.2.0
#                     emmeans v 2.0.1
#                     lme4 v 1.1.38
#                     gstat v 2.1.4
#                     car v 3.1.3
#                     FNN v 1.1.4.1
#                     spdep v 1.4.1
#                     ade4 v 1.7.23
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(sf); packageVersion("sf")
library(raster); packageVersion("raster")
library(cowplot); packageVersion("cowplot")
library(emmeans); packageVersion("emmeans")
library(lme4); packageVersion("lme4")
library(gstat); packageVersion("gstat")
library(car); packageVersion("car")
library(FNN); packageVersion("FNN")
library(spdep); packageVersion("spdep")
library(ade4); packageVersion("ade4")

#################################################################################
#                               Main workflow                                   #
#  Explore the growth data for all of the trait trees sampled at WFDP. Assess   #
#  for differences between species, or any relationships to environmental       #
#  variables so far.                                                            # 
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)


## Tree Demography
# Read in tree demographic data 
growth <- read.csv("stems_WFDP_20250206_trimmed.csv")

# each tree is identified by its WFDP stem_tag. There are three census years - 2011, 
# 2016, and 2021.


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
    .groups = "drop")


diams$WFDP_Code <- diams$Stem_Tag


## Environmental/ Topographic Data 

# This can be updated once the krieged environmental data is in hand 

env_data <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# Trim to just relevant columns 
env_data <- dplyr::select(env_data, Cell, EM_Sample_Name, slope, aspect, elevation_m, Association, WFDP_Code)

# Make association a factor 
env_data$Association <- as.factor(env_data$Association)


# Merge enviro data to the existing files according to the WFDP_Code
growth_env <- merge(diams, env_data, by = "WFDP_Code")


# make dataframe of full species names 

full <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
          "T. plicata", "T. heterophylla")

Species <- c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")


taxa <- data.frame(full, Species)

# Merge to demo data by Species 
growth_env <- merge(growth_env, taxa, by = "Species")



#########################################################################################

## PREP SPATIAL DATA

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


###############################

##SET PLOTTING SPECS 

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# Define shapes for species 

# ABAM, ABGR, ALRU, CONU, TABR, THPL, TSHE  
species_shapes <- c(15, 16, 17, 18, 7, 8, 9)



spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
growth_env$full <- factor(growth_env$full, levels = spp_order)

#########################################################################################################################

##################################################### -- 
# (2) RGR VARIATION BETWEEN SPECIES AND ACROSS PLOT 
##################################################### -- 


# Visualize mean RGR between the focal tree species 

spp_rgr <- ggplot(growth_env, aes(y = mean_RGR, x = full, fill = full)) +
  geom_boxplot() +
  geom_point(alpha = 0.6) +
  theme_bw() +
  scale_fill_manual(values=all_hosts, 
                    name="Focal Species",
                    breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                             "T. plicata", "T. heterophylla"),
                    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                             "T. plicata", "T. heterophylla")) +
  labs(x = "", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 14, colour="black", face = "italic"),
    axis.text.y = element_text(size = 16, colour="black"),
    axis.title.y = element_text(size = 16, colour="black"),
    axis.title.x = element_text(size = 16, colour="black")) +
  theme(legend.text = element_text(size = 14, colour="black"), 
        legend.title = element_text(size = 14, face = "bold", colour="black")) +
  theme(legend.position = "none")

spp_rgr


# Save Figure 2 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/RGR_by_species.png", 
       plot = spp_rgr, width = 10, height = 8, units = "in", dpi = 300)


#summaries 

rgr_spp_aov <- aov(mean_RGR ~ Species, data = growth_env)

summary(rgr_spp_aov)

# RGR significantly differs between species 
# 
#           Df   Sum Sq   Mean Sq F value         Pr(>F)    
#     Species      6 0.002929 0.0004882   5.553 0.00016 ***
#   Residuals   53 0.004660 0.0000879                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


tuk_rgr_spp <- TukeyHSD(rgr_spp_aov)
tuk_rgr_spp

# 
# diff           lwr           upr     p adj
# ABGR-ABAM -6.952969e-03 -0.0217905642  0.0078846261 0.7799461
# ALRU-ABAM  1.022609e-02 -0.0055115590  0.0259637335 0.4322334
# CONU-ABAM -4.838083e-03 -0.0176878172  0.0080116516 0.9079049
# TABR-ABAM -1.493780e-02 -0.0281396491 -0.0017359550 0.0170683
# THPL-ABAM  4.200038e-03 -0.0086496967  0.0170497720 0.9514257
# TSHE-ABAM -4.903199e-03 -0.0177529337  0.0079465350 0.9024389
# ALRU-ABGR  1.717906e-02 -0.0002195662  0.0345776789 0.0550804
# CONU-ABGR  2.114886e-03 -0.0127227089  0.0169524814 0.9994224
# TABR-ABGR -7.984833e-03 -0.0231283902  0.0071587242 0.6729593
# THPL-ABGR  1.115301e-02 -0.0036845885  0.0259906019 0.2616727
# TSHE-ABGR  2.049770e-03 -0.0127878255  0.0168873649 0.9995173
# CONU-ALRU -1.506417e-02 -0.0308018163  0.0006734762 0.0690773
# TABR-ALRU -2.516389e-02 -0.0411903238 -0.0091374548 0.0002476
# THPL-ALRU -6.026050e-03 -0.0217636959  0.0097115966 0.9009760
# TSHE-ALRU -1.512929e-02 -0.0308669329  0.0006083596 0.0669937
# TABR-CONU -1.009972e-02 -0.0233015663  0.0031021278 0.2429691
# THPL-CONU  9.038120e-03 -0.0038116139  0.0218878548 0.3367882
# TSHE-CONU -6.511657e-05 -0.0129148509  0.0127846178 1.0000000
# THPL-TABR  1.913784e-02  0.0059359926  0.0323396868 0.0008547
# TSHE-TABR  1.003460e-02 -0.0031672444  0.0232364498 0.2497722
# TSHE-THPL -9.103237e-03 -0.0219529714  0.0037464973 0.3283815


# TABR is different from ABAM, ALRU, THPL


# Visualize in a way to be able to zoom in on which trees were just growing particularly faster 
# or slower, and what their conditions were 


# merge diams and env_data by WFDP_Code to get the subplot IDs as Cell 

growth_env <- growth_env %>% arrange(desc(mean_RGR)) %>% 
  mutate(WFDP_Code = factor(WFDP_Code, levels = WFDP_Code))


spp_rgr_plot <- ggplot(growth_env, aes(x = WFDP_Code, y = mean_RGR, colour = full, shape = full)) + 
  geom_point(alpha = 1, size = 3) +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
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
  labs(x = "Tree Stemtag", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black", hjust = 1, angle = 45),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 14, colour="black")) +
  theme(legend.text = element_text(size = 14, colour="black", face = "italic"), 
        legend.title = element_text(size = 14, face = "bold", colour="black")) +
  theme(legend.position = "bottom")

spp_rgr_plot



# Plot variation in mean RGR spatially as a representation across the plot 

# merge the rotated focal_trees df with the plot_df dataframe to get the data for plotting

trees_spatial_df <- merge(focal_trees_rotated, growth_env, by = "Stem_Tag")


spp_rgr_spatial <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = trees_spatial_df, aes(color = mean_RGR, shape = full), size = 4) +
  scale_color_gradient2(low = "red3", mid = "tan", high = "blue4", 
                        midpoint = 0.0085, limits = c(-0.012, 0.043), 
                        name = expression("Mean Relative\nGrowth Rate ("*yr^{-1}*")")) +
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(colour="black", size=14, face="bold"),
    legend.text  = element_text(colour="black", size=14, face = "italic"),
    axis.text    = element_text(size=1, colour="white"), # hide the coordinates 
    axis.title   = element_text(size=14, colour="black"))

# removing legends for plotting too 

spp_rgr_spatial


## It's too hard to see the colors with the different shapes, so going to do it just where the growth rate 
# is reflected in the size of the point? 


spp_rgr_spatial2 <- ggplot() +
  geom_sf(data = wfdp_poly_rotated, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = trees_spatial_df, aes(colour = full, size = mean_RGR))+
  scale_size_continuous(limits = c(-0.012, 0.043), 
                        name = expression("Mean Relative\nGrowth Rate ("*yr^{-1}*")"), 
  ) +
  scale_color_manual(values=all_hosts, 
                     name="Focal Species",
                     breaks=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla"),
                     labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
                              "T. plicata", "T. heterophylla")) +
  coord_sf(expand = FALSE) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=14, face="bold"),
    legend.text  = element_text(colour="black", size=14, face = "italic"),
    axis.text    = element_text(size=1, colour="white"), # hide the coordinates 
    axis.title   = element_text(size=14, colour="black"))

# removing legends for plotting too 

spp_rgr_spatial2


# Save as an additional panel for Figure 2? 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/RGR_by_species_spatial.png", 
       plot = spp_rgr_spatial2, width = 10, height = 6, units = "in", dpi = 300)



#### Test for spatial structuring of growth rate across the plot #### -- 


## PREP

# Get focal tree coords 
tree_coords <- st_coordinates(trees_spatial_df) %>% as.data.frame()

# Add Stem_tag column 
tree_coords$Stem_Tag <- trees_spatial_df$Stem_Tag # only doing it this way because I know they are in the same order 


# Merge back to trees_spatial_df by Stem_Tag
trees_spatial_df <- merge(trees_spatial_df, tree_coords, by = "Stem_Tag")

# get spatial dist for all focal trees 
tree.dists <- dist(cbind(trees_spatial_df$X, trees_spatial_df$Y))


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


## TESTS 

# Perform nearest-neighbor comparison of the data 
nn <- get.knn(tree.dists, k = 2)

trees_spatial_df$nn_dist <- nn$nn.dist[,1]
trees_spatial_df$nn_diff <- abs(
  trees_spatial_df$mean_RGR -
    trees_spatial_df$mean_RGR[nn$nn.index[,1]])


nn_plot <- ggplot2::ggplot(trees_spatial_df, aes(nn_dist, nn_diff)) +
  geom_point() +
  theme_minimal() +
  labs(
    x = "Distance to nearest tree (m)",
    y = "Difference in mean RGR")

nn_plot 


# get euclidian distance for mean_RGR for trees
env_dist_mean_RGR <- dist(trees_spatial_df$mean_RGR)


# Run mantel correlogram 
mc_mean_RGR <- ade4::mantel.randtest(env_dist_mean_RGR, tree.dists, nrepet = 9999)

mc_mean_RGR


# Monte-Carlo test
# Call: ade4::mantel.randtest(m1 = env_dist_mean_RGR, m2 = tree.dists, 
#                             nrepet = 9999)
# 
# Observation: 0.02121855 
# 
# Based on 9999 replicates
# Simulated p-value: 0.2841 
# Alternative hypothesis: greater 
# 
# Std.Obs  Expectation     Variance 
# 0.5103218465 0.0002924255 0.0016814701 


# No spatial autocorrelation in tree RGR


mantel_mean_RGR <- ggplot_mantel(mc_mean_RGR)

mantel_mean_RGR

# Throws some errors, but this plots the same as if you just did plot() so it's okay 

# Vertical line reflects observed value, and the frequency is the distribution of the 
# simulated values. If the distribution is skewed to one side or the other of the vertical 
# line, then this reflects a significant trend in spatial autocorrelation, either positively 
# or negatively 

#########################################################################################################################

###################################### -- 
# (3) RGR RELATIONSHIP TO TOPOGRAPHY
###################################### -- 

# Exploring separate slopes per species here, because we would expect that growth will vary between species, 
# and may also vary individually with topography. Also testing the linear regressions with an interaction between 
# the topographic variable and species, to get the individual slopes for each species for the relationship 
# between growth and the topographic variable. 


# If desired, should update these to reflect the full species names instead 


### Explore RGR relationships 

## Slope
slope2 <- ggplot(growth_env, aes(x = slope, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Slope", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

slope2

# test relationships 
lm_slope2 <- lm(mean_RGR ~ slope * Species, data = growth_env)

emtrends(lm_slope2, ~ Species, var = "slope") %>% test(adjust = "fdr")


# Species slope.trend      SE df t.ratio p.value
# ABAM       0.000417 0.00275 46   0.152  0.9599
# ABGR       0.000810 0.00894 46   0.091  0.9599
# ALRU      -0.013660 0.00571 46  -2.391  0.1467
# CONU      -0.000265 0.00524 46  -0.051  0.9599
# TABR       0.000967 0.00256 46   0.378  0.9599
# THPL       0.002027 0.00195 46   1.039  0.9599
# TSHE       0.001395 0.00366 46   0.381  0.9599


anova(lm_slope2)

# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# slope          1 0.0001564 0.00015639  1.7828 0.1883771    
# Species        6 0.0027952 0.00046587  5.3106 0.0003132 ***
#   slope:Species  6 0.0006022 0.00010036  1.1441 0.3525130    
# Residuals     46 0.0040353 0.00008772                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



## Aspect
aspect2 <- ggplot(growth_env, aes(x = aspect, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Aspect", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

aspect2


# test relationships 
lm_aspect2 <- lm(mean_RGR ~ aspect * Species, data = growth_env)

emtrends(lm_aspect2, ~ Species, var = "aspect") %>% test(adjust = "fdr")


# Species aspect.trend       SE df t.ratio p.value
# ABAM       -6.95e-05 9.96e-05 46  -0.697  0.9126
# ABGR       -4.72e-05 1.00e-04 46  -0.471  0.9126
# ALRU        2.47e-05 1.82e-04 46   0.135  0.9126
# CONU        1.37e-05 5.91e-05 46   0.233  0.9126
# TABR        7.05e-05 1.50e-04 46   0.469  0.9126
# THPL        1.77e-05 3.29e-05 46   0.537  0.9126
# TSHE        1.39e-05 1.26e-04 46   0.110  0.9126


anova(lm_aspect2)

# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# aspect          1 0.0000018 0.00000181  0.0184 0.8926854    
# Species         6 0.0029389 0.00048981  4.9722 0.0005358 ***
#   aspect:Species  6 0.0001169 0.00001948  0.1978 0.9758017    
# Residuals      46 0.0045315 0.00009851                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## Elevation
elev2 <- ggplot(growth_env, aes(x = elevation_m, y = mean_RGR, colour = Species)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Species, scales = "free_x") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Elevation (m)", y = expression("Mean Relative Growth Rate ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "none")

elev2


# test relationships 
lm_elev2 <- lm(mean_RGR ~ elevation_m * Species, data = growth_env)

emtrends(lm_elev2, ~ Species, var = "elevation_m") %>% test(adjust = "fdr")

# Species elevation_m.trend       SE df t.ratio p.value
# ABAM            -6.68e-04 0.000356 46  -1.874  0.4714
# ABGR             1.54e-04 0.000452 46   0.340  0.9988
# ALRU            -3.75e-04 0.000532 46  -0.704  0.9988
# CONU            -8.62e-07 0.000550 46  -0.002  0.9988
# TABR             8.30e-05 0.000427 46   0.195  0.9988
# THPL             4.52e-04 0.000406 46   1.113  0.9507
# TSHE            -5.75e-05 0.000333 46  -0.173  0.9988

anova(lm_elev2)

# Response: mean_RGR
# Df    Sum Sq    Mean Sq F value    Pr(>F)    
# elevation_m          1 0.0001623 0.00016235  1.7918 0.1872884    
# Species              6 0.0027895 0.00046492  5.1312 0.0004158 ***
#   elevation_m:Species  6 0.0004694 0.00007823  0.8634 0.5288464    
# Residuals           46 0.0041679 0.00009061                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Outcome: Same general trends between diam_diff and RGR, species differ in their growth rates, 
# but these are not being driven by any of the topographic variables. 



## -- END -- ##


