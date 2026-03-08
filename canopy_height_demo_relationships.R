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

#################################################################################

###################################################
# (2) VISUALIZE CANOPY VARIATION ACROSS THE PLOT 
###################################################

# Need to remove all of the coordinates from the plots 

# 1. Mean Canopy Height - 9m radius 
plot_mean_9m <- ggplot() +
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = all_canopy_growth_sf, aes(color = X9m_mean), size = 4) +
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
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = all_canopy_growth_sf, aes(color = X20m_mean), size = 6) +
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
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = all_canopy_growth_sf, aes(color = pr_2_20), size = 6) +
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
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = all_canopy_growth_sf, aes(color = pr_5_20), size = 6) +
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
  geom_sf(data = wfdp_poly, fill = NA, color = "black", linewidth = 0.6) +
  geom_sf(data = all_canopy_growth_sf, aes(color = p_10_20), size = 6) +
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



# Are there spatial relationships with canopy height? 

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")

# with easting in 9m neighborhoods 
east_height_9m <- ggplot(all_canopy_growth, aes(x = UTM_X, y = X9m_mean, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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
north_height_9m <- ggplot(all_canopy_growth, aes(x = UTM_Y, y = X9m_mean, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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

#################################################################################

#########################################################################
# (3) MODEL RELATIONSHIPS OF NEIGHBORHOOD CANOPY HEIGHT AND TREE GROWTH 
#########################################################################

# 1. Mean canopy height 

# Relationship between focal tree RGR and mean canopy height in 9m radius 
RGR_height_9m <- ggplot(all_canopy_growth, aes(x = X9m_mean, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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
RGR_height_20m <- ggplot(all_canopy_growth, aes(x = X20m_mean, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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
RGR_max_height_9m <- ggplot(all_canopy_growth, aes(x = X9m_max, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Maximum Canopy Height (m)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

RGR_max_height_9m


# The regression model
RGR_mod2_09 <- lm(mean_RGR ~ X9m_max, data = all_canopy_growth)

summary(RGR_mod2_09)

# No Significant relationship
# Multiple R-squared:  0.02009,	Adjusted R-squared:  0.003196 
# F-statistic: 1.189 on 1 and 58 DF,  p-value: 0.28



# Relationship between focal tree RGR and max canopy height in 20m radius 
RGR_max_height_20m <- ggplot(all_canopy_growth, aes(x = X20m_max, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(x = "Maximum Canopy Height (m)", y = expression("Mean RGR ("*yr^{-1}*")")) +
  theme(
    axis.text.x = element_text(size = 11, colour="black"),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 11, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

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
RGR_2open_20m <- ggplot(all_canopy_growth, aes(x = pr_2_20, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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
RGR_5open_20m <- ggplot(all_canopy_growth, aes(x = pr_5_20, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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
RGR_10open_20m <- ggplot(all_canopy_growth, aes(x = p_10_20, y = mean_RGR, colour = Species.x)) +
  geom_point(alpha = 1, cex = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Focal Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
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


############################################## -- 
# (4) GATHER SIGNIFICANT RESULTS OF INTEREST  
############################################## -- 

# Gather up the most interesting plots and organize the plots to show the canopy height variation 
# across the plot, and then to summarize the relationships to focal tree growth


# Spatial plots 
spatial_plots <- plot_grid(plot_2open_20m, plot_mean_9m, plot_5open_20m, plot_mean_20m, plot_10open_20m,
                           ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)'))

spatial_plots


# Organize and number plots for the growth relationships with mean canopy height 
canopy_plots <- plot_grid(RGR_2open_20m, RGR_height_9m, RGR_5open_20m, RGR_height_20m, RGR_10open_20m,
                          ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)', '(e)'))

canopy_plots


## -- END -- ## 
