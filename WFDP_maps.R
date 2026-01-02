# -----------------------------------------------------------------------------#
# Maps for WFDP Location 
# Original Author: L. McKinley Nevins 
# December 29, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     ggplot2 v 3.5.1
#                     maps v 3.4.3
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(ggplot2); packageVersion("ggplot2")
library(maps); packageVersion("maps")

#################################################################################
#                               Main workflow                                   #
#  Generate several maps at different scales to show the location of the WFDP   #
#  in North America and in the Northwestern US. Use as insets for figure 1.     #
#                                                                               #
#################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/"
setwd(wd)


################ --
# (2) MAPPING
################ --

# Broad view
maps::map(database = 'world', boundary = TRUE, xlim=c(-140,-50.0), ylim=c(15.0,70.0), 
          col = "gray", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red3', cex = 1.5)
map.axes()

# Regional view 
maps::map(database = 'state', region = c('washington', 'oregon', 'idaho', 'montana'), 
          boundary = TRUE, xlim=c(-126,-115), ylim=c(45.0,49.0), col = "gray", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red3', cex = 1.5)
map.axes()

# Counties - probably not relevant 
maps::map('county', 'washington', boundary = TRUE, xlim=c(-124,-120), ylim=c(45.3,47.0), 
          col = "gray", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red3', cex = 2)
map.axes()

