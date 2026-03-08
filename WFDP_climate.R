# -----------------------------------------------------------------------------#
# Process NOAA 30 year Normals (1990-2020) data for WFDP
# Original Author: L. McKinley Nevins 
# January 5, 2026
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     cowplot v 1.2.0
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")

#################################################################################
#                               Main workflow                                   #
#  Use data of monthly 30 year normalized temperature and precipitation for     # 
#  WFDP, collected from the NOAA station at the Carson Fish Hatchery. Calculate #
#  summary values and plot variation over the year to reflect summer drought    #
#  conditions.                                                                  #
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/"
setwd(wd)

# Load in monthly NOAA data 

clim <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/NOAA_monthly_30yr_normals.csv")

# This doesn't contain data for the full 30 year period, there is only 24 years of temperature data and 15 years of precip data 


# STATION	NAME                            	LATITUDE	   LONGITUDE	   ELEVATION	
# USC00451160	CARSON FISH HATCHERY, WA US	  45.8678	     -121.9733	   345.6	 


# ANN-PRCP-NORMAL  	years_ANN-PRCP-NORMAL	ANN-TAVG-NORMAL		years_ANN-TAVG-NORMAL	  DJF-PRCP-NORMAL	 	years_DJF-PRCP-NORMAL	
#       2258.1	             	15	              9.1	                    24	               1056.4	    	        15	

# JJA-PRCP-NORMAL	years_JJA-PRCP-NORMAL	  MAM-PRCP-NORMAL	years_MAM-PRCP-NORMAL	SON-PRCP-NORMAL	   years_SON-PRCP-NORMAL						
#       86.6	           	25	                    489.2	         19	                    625.9	                19						



# Columns of interest are MLY.TMAX.NORMAL, MLY.TMIN.NORMAL, MLY.PRCP.NORMAL

# There is no STDDEV info for precip, so I'm just going to ignore those columns for temp. 

sub_clim <- dplyr::select(clim, month, MLY.TMAX.NORMAL, MLY.TAVG.NORMAL, MLY.TMIN.NORMAL, MLY.PRCP.NORMAL)


# Add a column of the month abbreviation for clearer plotting 
sub_clim$month_code <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

sub_clim$month_code <- as.factor(sub_clim$month_code)


# Convert temp from F to C, and precip from in to mm

sub_clim <- sub_clim %>%
  dplyr::mutate(
    temp_max_C   = (MLY.TMAX.NORMAL - 32) * 5/9,
    temp_min_C   = (MLY.TMIN.NORMAL - 32) * 5/9,
    temp_mean_C   = (MLY.TAVG.NORMAL - 32) * 5/9,
    precip_mm = MLY.PRCP.NORMAL * 25.4)

mean(sub_clim$temp_mean_C)

sum(sub_clim$precip_mm)

# Mean annual precip from 1991-2020 was 2,258.06 mm 


# Calculate Walter–Lieth drought threshold 
sub_clim <- sub_clim %>%
  dplyr::mutate(
    drought_threshold = 2 * temp_mean_C)



# Reshape data for plotting
clim_long <- sub_clim %>%
  dplyr::select(month, temp_max_C, temp_min_C, temp_mean_C, precip_mm) %>%
  tidyr::pivot_longer(
    cols = -month,
    names_to = "variable",
    values_to = "value"
  )

#################################################################################

############################# --
# (2) PLOT MONTHLY NORMALS 
############################# --

# add shading for growing season 
grow_start <- 5
grow_end   <- 10

plot2 <- ggplot(sub_clim, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.5, xmax = grow_end + 0.5,
    ymin = -Inf, ymax = Inf, fill = "darkgreen", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted", linewidth = 1) +
  geom_vline(
    xintercept = c(grow_start - 0.5, grow_end + 0.5), linetype = "dashed",
    color = "darkgreen", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
    labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "blue4", "Mean Temperature" = "darkorange3",
      "Max Temperature" = "darkred", "Min Temperature" = "orange2")) +
  labs(x = "Month", color = "") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

plot2


# -- END -- # 
