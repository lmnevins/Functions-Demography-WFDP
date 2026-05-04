# -----------------------------------------------------------------------------#
# Exploring WFDP tree root and leaf trait data 
# Original Author: L. McKinley Nevins 
# January 9, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 4.0.1
#                     rstatix v 0.7.3
#                     vegan v 2.7.2
#                     rcompanion v 2.5.1
#                     ggfortify v 0.4.19
#                     gginnards v 0.2.0.2
#                     ggrepel v 0.9.6
#                     corrplot v 0.95
#                     car v 3.1.3
#                     multcomp v 1.4.29
#                     multcompView v 0.1.10
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(vegan); packageVersion("vegan")
library(rcompanion); packageVersion("rcompanion")
library(ggfortify); packageVersion("ggfortify")
library(gginnards); packageVersion("gginnards")
library(ggrepel); packageVersion("ggrepel")
library(corrplot); packageVersion("corrplot")
library(car); packageVersion("car")
library(multcomp); packageVersion("multcomp")
library(multcompView); packageVersion("multcompView")

#################################################################################
#                               Main workflow                                   #
#  Explore the leaf and root trait data. Determine if there are differences     #
#  between species, and major PCA axes for the separate and combined traits.    #
#  Explore relationships between traits typically considered to be equivalent   #
#  above and belowground.                                                       #
#                                                                               #
#################################################################################

############### -- 
# (1) DATA PREP
############### -- 

wd <- "~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/"
setwd(wd)

# Load in trait datasets and clean a bit 

##Leaves 
leaf <- read.csv("WFDP_leaf_traits.csv")

# make species a factor 
leaf$Host_ID <- as.factor(leaf$Host_ID)

# filter out PSME as a host since there were so few sampled 
leaf <- leaf %>% filter(Host_ID != 'PSME')

# filter out T-TABR-03 as a host since it had no fungal community data 
leaf <- leaf %>% filter(code != 'T-TABR-03')

## Roots
root <- read.csv("WFDP_root_traits.csv")

root$Host_ID <- as.factor(root$Host_ID)

root <- root %>% filter(Host_ID != 'PSME')
root <- root %>% filter(code != 'T-TABR-03')


# combine into one dataset for the trees 
traits <- merge(leaf, root, by = 'code')

# subset just variables of interest, excluding petiole data since it is absent for
# needle-leaf species 

traits <- dplyr::select(traits, code, WFDP_Code = WFDP_Code.x, sub_plot = sub_plot.x, Host_ID = Host_ID.x, SLA_leaf, LDMC_leaf, LMA_leaf, 
                 leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N, specific_root_length, specific_root_area,
                 root_dry_matter_cont, root_CN, root_15N, avg_root_dia, root_pct_N, root_pct_C)

## UPDATE: removed C isotope traits 


# make dataframe of full species names 

sci_name <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
          "T. plicata", "T. heterophylla")

Host_ID <- c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")


taxa <- data.frame(sci_name, Host_ID)

# Merge to traits data by Species 
traits <- merge(traits, taxa, by = "Host_ID")

leaf <- merge(leaf, taxa, by = "Host_ID")

root <- merge(root, taxa, by = "Host_ID")

############################################## -- 
## create species x site matrix 

# load in all WFDP environmental data 

# This can be updated when krieged data is in hand 

env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# grab just sub_plot (referring to the subplot in WFDP - the site level) and Host_ID

env$sub_plot <- env$Cell

env <- dplyr::select(env, Host_ID, sub_plot)


# create a presence-absence matrix from these 
sites_tree <- env %>%
  mutate(presence = 1) %>%
  pivot_wider(names_from = Host_ID, values_from = presence, values_fill = list(presence = 0))

print(sites_tree)

sites_tree$sub_plot <- as.factor(sites_tree$sub_plot)

#move subplot into the rownames 
sites_tree <- sites_tree %>%
  column_to_rownames(var = "sub_plot")

############################################## -- 
## create species x trait matrix 

# need species as the rownames and columns to only be the traits 

# remove code column 
traits_tree <- subset(traits, select = -c(code, WFDP_Code, sub_plot, Host_ID, sci_name))



#################################################################################

########################################## -- 
# (2) ASSESS VARIATION IN RAW TRAIT DATA
########################################## -- 

# set colors for hosts 
# ABAM      ABGR      ALRU        CONU     TABR        THPL       TSHE        
all_hosts <- c("#FFD373", "#FD8021", "#E05400", "#0073CC","#003488", "#001D59", "#001524")


# Define shapes for species 

# ABAM, ABGR, ALRU, CONU, TABR, THPL, TSHE  
species_shapes <- c(15, 16, 17, 18, 7, 8, 9)


# pause and look at variation in raw traits between species 


# Grab environmental data again for some association labels 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

env <- dplyr::select(env, WFDP_Code, Association, Host_ID)


#subset a few columns from traits 

traits_sub <- dplyr::select(traits, -c("code", "WFDP_Code", "sub_plot", "sci_name"))


# Convert to long format for ggplot
long_traits <- pivot_longer(traits_sub, cols = -(Host_ID), names_to = "trait", values_to = "value")

# Boxplot
trait_plot_OG <- ggplot(long_traits, aes(x = Host_ID, y = value, fill = Host_ID)) +
  geom_boxplot() +
  facet_wrap(~trait, scales = "free_y") +  # Separate plots for each trait
  theme_bw() +
  scale_colour_manual(values=all_hosts, 
                      name="Host Species",
                      breaks=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE"),
                      labels=c("ABAM", "ABGR", "ALRU", "CONU", "TABR", "THPL", "TSHE")) +
  labs(title = "Trait Variation Across Species", y = "Trait Value") +
 theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 11)) + 
  theme(
    axis.text.x = element_text(size = 11, colour="black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 11, colour="black"))

trait_plot_OG


## Get summary table of species-level means and standard error 

# Standard error helper function
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

trait_summary <- traits_sub %>%
  group_by(Host_ID) %>%
  summarise(across(
    where(is.numeric), 
    list(mean = ~mean(.x, na.rm = TRUE), se = ~se(.x)),
    .names = "{.col}_{.fn}"))





# SRA across species 

# add association 
root$Association <- env$Association


SRA_plot <- ggplot(root, aes(x = sci_name, y = specific_root_area, fill = Association)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\n Association") + 
  theme_bw() +
  labs(title = "", y = "Specific Root Area", x = "") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black", face = "italic"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "bottom")

SRA_plot


# Test for differences between species 
aov_SRA <- aov(specific_root_area ~ Host_ID, data = root)
summary(aov_SRA)

tuk_SRA <- TukeyHSD(aov_SRA)
tuk_SRA

# Df  Sum Sq   Mean Sq F value Pr(>F)  
# Host_ID      6 0.01267 0.0021124   2.921 0.0155 *
#   Residuals   53 0.03833 0.0007232      


# Only between ALRU-ABGR  0.051200000  0.001301082 0.101098918 0.0407861


# Test for differences between association
aov_SRA2 <- aov(specific_root_area ~ Association, data = root)
summary(aov_SRA2)

tuk_SRA2 <- TukeyHSD(aov_SRA2)
tuk_SRA2


# Marginally significant 
# Df  Sum Sq   Mean Sq F value Pr(>F)  
# Association  2 0.00507 0.0025357   3.147 0.0505 .
# Residuals   57 0.04593 0.0008058

# and sig diff between EM and AM 
# diff         lwr           upr     p adj
# DUAL-AM -0.006466165 -0.03052652  0.0175941906 0.7949653
# EM-AM   -0.020598441 -0.04105375 -0.0001431274 0.0480612 * 
# EM-DUAL -0.014132275 -0.03662965  0.0083651040 0.2931193


# SRL across species 
SRL_plot <- ggplot(root, aes(x = sci_name, y = specific_root_length, fill = Association)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\n Association") + 
  theme_bw() +
  labs(title = "", y = "Specific Root Length", x = "") +
  theme(legend.position = "right")  +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black", face = "italic"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

SRL_plot


# Test for differences between species 
aov_SRL <- aov(specific_root_length ~ Host_ID, data = root)
summary(aov_SRL)

tuk_SRL <- TukeyHSD(aov_SRL)
tuk_SRL


# $Host_ID
# diff         lwr         upr     p adj
# ABGR-ABAM -0.10673333 -0.38867912  0.17521246 0.9057196
# ALRU-ABAM  0.26380000 -0.03524867  0.56284867 0.1174988
# CONU-ABAM -0.00490000 -0.24907222  0.23927222 1.0000000
# TABR-ABAM -0.08628889 -0.33715199  0.16457422 0.9384496
# THPL-ABAM  0.04590000 -0.19827222  0.29007222 0.9972452
# TSHE-ABAM -0.04980000 -0.29397222  0.19437222 0.9956850
# ALRU-ABGR  0.37053333  0.03992259  0.70114408 0.0187135
# CONU-ABGR  0.10183333 -0.18011246  0.38377912 0.9232717
# TABR-ABGR  0.02044444 -0.26731527  0.30820416 0.9999902
# THPL-ABGR  0.15263333 -0.12931246  0.43457912 0.6455301
# TSHE-ABGR  0.05693333 -0.22501246  0.33887912 0.9959140
# CONU-ALRU -0.26870000 -0.56774867  0.03034867 0.1050654
# TABR-ALRU -0.35008889 -0.65462515 -0.04555263 0.0145950 **
# THPL-ALRU -0.21790000 -0.51694867  0.08114867 0.2960614
# TSHE-ALRU -0.31360000 -0.61264867 -0.01455133 0.0340552 *
# TABR-CONU -0.08138889 -0.33225199  0.16947422 0.9530983
# THPL-CONU  0.05080000 -0.19337222  0.29497222 0.9951910
# TSHE-CONU -0.04490000 -0.28907222  0.19927222 0.9975622
# THPL-TABR  0.13218889 -0.11867422  0.38305199 0.6736149
# TSHE-TABR  0.03648889 -0.21437422  0.28735199 0.9993513
# TSHE-THPL -0.09570000 -0.33987222  0.14847222 0.8906700


## STOP: Carbon isotopes have been removed now, so this is outdated. Can add back in at the top 
# if needed to rerun this 


# Take a closer look at just the isotopic results 
isotopes <- dplyr::select(traits, sci_name, leaf_15N, root_15N)

# Merge with env to get associations 
isotopes$Association <- env$Association

# Convert to long format for ggplot
isotopes <- pivot_longer(isotopes, cols = -c(sci_name, Association), names_to = "trait", values_to = "value")


# Create a dataframe to label the isotopes better for faceting 

trait <- c("leaf_15N", "root_15N")

label <- c("leaf_15N", "root_15N")

isotope_labels <- data.frame(trait, label)

isotopes <- merge(isotopes, isotope_labels, by = "trait")


# Put species in order for plotting 
spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
isotopes$sci_name <- factor(isotopes$sci_name, levels = spp_order)


# Boxplot
isotopes_plot <- ggplot(isotopes, aes(x = sci_name, y = value, fill = Association)) +
  geom_boxplot() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\n Association") + 
  facet_wrap(~label, scales = "free_y") +  # Separate plots for each trait
  theme_bw() +
  labs(title = "", y = "Isotope per mil", x = "") +
theme(
    axis.text.x = element_text(size = 12, colour="black", face = "italic", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 14, colour="black")) +
  theme(legend.text = element_text(size = 14, colour="black"), 
        legend.title = element_text(size = 14, face = "bold", colour="black")) +
  theme(legend.position = "none")

isotopes_plot

## Can do some stats here if desired to test for differences between the species for 
# these values 

isotopes$trait <- as.factor(isotopes$trait)

leaf_13C <- isotopes %>%
  filter(trait == "leaf_13C") %>%
  droplevels()

# Test for significant differences between species 
aov_13C_leaf <- aov(value ~ Host_ID, data = leaf_13C)
summary(aov_13C_leaf)

tuk_13C_leaf <- TukeyHSD(aov_13C_leaf)
tuk_13C_leaf

# Df Sum Sq Mean Sq F value   Pr(>F)    
# Host_ID      6  60.83   10.14   17.48 4.92e-11 ***
#   Residuals   53  30.74    0.58                     
# ---


# diff        lwr        upr     p adj
# ABGR-ABAM  0.2486667 -0.9564932  1.4538265 0.9954029
# ALRU-ABAM  0.6330000 -0.6452650  1.9112650 0.7332349
# CONU-ABAM  0.5390000 -0.5046990  1.5826990 0.6937297
# TABR-ABAM  2.6325556  1.5602567  3.7048544 0.0000000 ***
# THPL-ABAM  2.4380000  1.3943010  3.4816990 0.0000001 ***
# TSHE-ABAM  0.7680000 -0.2756990  1.8116990 0.2850233
# ALRU-ABGR  0.3843333 -1.0288418  1.7975085 0.9802568
# CONU-ABGR  0.2903333 -0.9148265  1.4954932 0.9894530
# TABR-ABGR  2.3838889  1.1538778  3.6139000 0.0000046 ***
# THPL-ABGR  2.1893333  0.9841735  3.3944932 0.0000175 ***
# TSHE-ABGR  0.5193333 -0.6858265  1.7244932 0.8393164
# CONU-ALRU -0.0940000 -1.3722650  1.1842650 0.9999880
# TABR-ALRU  1.9995556  0.6978342  3.3012769 0.0003530 ***
# THPL-ALRU  1.8050000  0.5267350  3.0832650 0.0012452 **
# TSHE-ALRU  0.1350000 -1.1432650  1.4132650 0.9998984
# TABR-CONU  2.0935556  1.0212567  3.1658544 0.0000039 ***
# THPL-CONU  1.8990000  0.8553010  2.9426990 0.0000170 ***
# TSHE-CONU  0.2290000 -0.8146990  1.2726990 0.9935873
# THPL-TABR -0.1945556 -1.2668544  0.8777433 0.9977375
# TSHE-TABR -1.8645556 -2.9368544 -0.7922567 0.0000410 ***
# TSHE-THPL -1.6700000 -2.7136990 -0.6263010 0.0001809 ***


# Test for significant differences between association
aov_13C_leaf2 <- aov(value ~ Association, data = leaf_13C)
summary(aov_13C_leaf2)

tuk_13C_leaf2 <- TukeyHSD(aov_13C_leaf2)
tuk_13C_leaf2


# Df Sum Sq Mean Sq F value   Pr(>F)    
# Association  2  24.25  12.124   10.27 0.000156 ***
#   Residuals   57  67.32   1.181    


# diff        lwr        upr     p adj
# DUAL-AM  0.4651128 -0.4560297  1.3862553 0.4492664
# EM-AM   -1.0337232 -1.8168479 -0.2505985 0.0067019
# EM-DUAL -1.4988360 -2.3601404 -0.6375315 0.0002877


# EM distinct from both AM and dual. 



root_13C <- isotopes %>%
  filter(trait == "root_13C") %>%
  droplevels()

# Test for significant differences between species 
aov_13C_root <- aov(value ~ Host_ID, data = root_13C)
summary(aov_13C_root)

tuk_13C_root <- TukeyHSD(aov_13C_root)
tuk_13C_root


# Df Sum Sq Mean Sq F value Pr(>F)  
# Host_ID      6   26.2   4.367   1.986 0.0841 .
# Residuals   53  116.5   2.199 

# No differences between species 



# Test for significant differences between association
aov_13C_root2 <- aov(value ~ Association, data = root_13C)
summary(aov_13C_root2)

tuk_13C_root2 <- TukeyHSD(aov_13C_root2)
tuk_13C_root2


# Df Sum Sq Mean Sq F value Pr(>F)  
# Association  2  18.43   9.214   4.225 0.0195 *
#   Residuals   57 124.31   2.181  

# diff       lwr        upr     p adj
# DUAL-AM -0.2520677 -1.503780  0.9996447 0.8788740
# EM-AM   -1.2053216 -2.269486 -0.1411575 0.0228109
# EM-DUAL -0.9532540 -2.123654  0.2171463 0.1315723

# AM and EM are different 



leaf_15N <- isotopes %>%
  filter(trait == "leaf_15N") %>%
  droplevels()

# Test for significant differences between species 
aov_15N_leaf <- aov(value ~ Host_ID, data = leaf_15N)
summary(aov_15N_leaf)

tuk_15N_leaf <- TukeyHSD(aov_15N_leaf)
tuk_15N_leaf

# 
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Host_ID      6 169.86  28.310   20.45 3.19e-12 ***
#   Residuals   53  73.39   1.385   



# diff        lwr        upr     p adj
# ABGR-ABAM -0.9246667 -2.7867196  0.9373863 0.7307170
# ALRU-ABAM  2.4210000  0.4459946  4.3960054 0.0073917 *
# CONU-ABAM  2.7250000  1.1124149  4.3375851 0.0000696 ***
# TABR-ABAM  3.1692222  1.5124485  4.8259959 0.0000061 ***
# THPL-ABAM  0.5600000 -1.0525851  2.1725851 0.9357069
# TSHE-ABAM -1.3230000 -2.9355851  0.2895851 0.1751707
# ALRU-ABGR  3.3456667  1.1622161  5.5291173 0.0003672 ***
# CONU-ABGR  3.6496667  1.7876137  5.5117196 0.0000036 ***
# TABR-ABGR  4.0938889  2.1934391  5.9943387 0.0000004 ***
# THPL-ABGR  1.4846667 -0.3773863  3.3467196 0.2015198
# TSHE-ABGR -0.3983333 -2.2603863  1.4637196 0.9944070
# CONU-ALRU  0.3040000 -1.6710054  2.2790054 0.9991041
# TABR-ALRU  0.7482222 -1.2630248  2.7594692 0.9125904
# THPL-ALRU -1.8610000 -3.8360054  0.1140054 0.0770615
# TSHE-ALRU -3.7440000 -5.7190054 -1.7689946 0.0000073 ***
# TABR-CONU  0.4442222 -1.2125515  2.1009959 0.9816339
# THPL-CONU -2.1650000 -3.7775851 -0.5524149 0.0024635 *
# TSHE-CONU -4.0480000 -5.6605851 -2.4354149 0.0000000 ***
# THPL-TABR -2.6092222 -4.2659959 -0.9524485 0.0002357 ***
# TSHE-TABR -4.4922222 -6.1489959 -2.8354485 0.0000000 ***
# TSHE-THPL -1.8830000 -3.4955851 -0.2704149 0.0124536


# Test for significant differences between association
aov_15N_leaf2 <- aov(value ~ Association, data = leaf_15N)
summary(aov_15N_leaf2)

tuk_15N_leaf2 <- TukeyHSD(aov_15N_leaf2)
tuk_15N_leaf2

# 
# Df Sum Sq Mean Sq F value   Pr(>F)    
# Association  2  132.8   66.41   34.28 1.68e-10 ***
#   Residuals   57  110.4    1.94      


# $Association
# diff         lwr       upr     p adj
# DUAL-AM  1.229211  0.04946539  2.408956 0.0393263
# EM-AM   -2.328752 -3.33173245 -1.325772 0.0000020
# EM-DUAL -3.557963 -4.66107106 -2.454855 0.0000000


# All groups distinct from eachother 




root_15N <- isotopes %>%
  filter(trait == "root_15N") %>%
  droplevels()

# Test for significant differences between species 
aov_15N_root <- aov(value ~ Host_ID, data = root_15N)
summary(aov_15N_root)

tuk_15N_root <- TukeyHSD(aov_15N_root)
tuk_15N_root


# Df Sum Sq Mean Sq F value  Pr(>F)   
# Host_ID      6  37.82   6.303   4.299 0.00133 **
#   Residuals   53  77.71   1.466   


# diff        lwr        upr     p adj
# ABGR-ABAM -0.46666667 -2.3827600  1.4494266 0.9888308
# ALRU-ABAM  0.61100000 -1.4213239  2.6433239 0.9674694
# CONU-ABAM -1.12900000 -2.7883855  0.5303855 0.3763640
# TABR-ABAM -0.80277778 -2.5076343  0.9020787 0.7760961
# THPL-ABAM  1.13900000 -0.5203855  2.7983855 0.3657760
# TSHE-ABAM -0.83600000 -2.4953855  0.8233855 0.7174414
# ALRU-ABGR  1.07766667 -1.1691519  3.3244852 0.7609892
# CONU-ABGR -0.66233333 -2.5784266  1.2537600 0.9370463
# TABR-ABGR -0.33611111 -2.2917157  1.6194934 0.9983289
# THPL-ABGR  1.60566667 -0.3104266  3.5217600 0.1569138
# TSHE-ABGR -0.36933333 -2.2854266  1.5467600 0.9968356
# CONU-ALRU -1.74000000 -3.7723239  0.2923239 0.1395988
# TABR-ALRU -1.41377778 -3.4833951  0.6558395 0.3715344
# THPL-ALRU  0.52800000 -1.5043239  2.5603239 0.9843740
# TSHE-ALRU -1.44700000 -3.4793239  0.5853239 0.3225546
# TABR-CONU  0.32622222 -1.3786343  2.0310787 0.9969609
# THPL-CONU  2.26800000  0.6086145  3.9273855 0.0019468 *
# TSHE-CONU  0.29300000 -1.3663855  1.9523855 0.9980561
# THPL-TABR  1.94177778  0.2369213  3.6466343 0.0160028 **
# TSHE-TABR -0.03322222 -1.7380787  1.6716343 1.0000000
# TSHE-THPL -1.97500000 -3.6343855 -0.3156145 0.0101980 **


# Test for significant differences between association
aov_15N_root2 <- aov(value ~ Association, data = root_15N)
summary(aov_15N_root2)

tuk_15N_root2 <- TukeyHSD(aov_15N_root2)
tuk_15N_root2


#              Df Sum Sq Mean Sq F value Pr(>F)
# Association  2   0.74   0.370   0.184  0.833
# Residuals   57 114.78   2.014   
# 
# 
# $Association
# diff       lwr       upr     p adj
# DUAL-AM -0.20496241 -1.407753 0.9978284 0.9116268
# EM-AM   -0.25136452 -1.273937 0.7712081 0.8252054
# EM-DUAL -0.04640212 -1.171059 1.0782546 0.9945805

# No differences between association 




## Get summary table of myco group-level means and standard error 
leaf_15N <- dplyr::select(traits, sci_name, leaf_15N)

root_15N <- dplyr::select(traits, sci_name, root_15N)

# Merge with env to get associations 
leaf_15N$Association <- env$Association

# Merge with env to get associations 
root_15N$Association <- env$Association


leaf_15N_summary <- leaf_15N %>%
  group_by(Association) %>%
  summarise(across(
    where(is.numeric), 
    list(mean = ~mean(.x, na.rm = TRUE), se = ~se(.x)),
    .names = "{.col}_{.fn}"))


root_15N_summary <- root_15N %>%
  group_by(Association) %>%
  summarise(across(
    where(is.numeric), 
    list(mean = ~mean(.x, na.rm = TRUE), se = ~se(.x)),
    .names = "{.col}_{.fn}"))










# Save figure 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Demography/Figures/isotope_plots.png", 
       plot = isotopes_plot, width = 9, height = 9, units = "in", dpi = 300)


# Interpretation: For 13C, leaf values are very distinct between species and myco associations, which 
# shows some differences that likely aren't just from canopy position. Though we expect that understory 
# species would have less negative values because straight soil respiration is -12 to -15. Root values are 
# not different between species, but they are between AM and EM associations, and I wonder if this has anything 
# to do with the carbon being allocated from the host to the fungi? Seems like the root values are about 3 permil 
# less negative than the leaf values, though not always. 

# For 15N, leaf values have many differences between species, and all associations are different from eachother. 
# The values also have a greater spread, compared to the root values that have less differentiation between species, 
# and no differences between associations. So the important pathway here seems to be that the myco types are bringing 
# in different types of nitrogen and then as it is processed and translocated up to the leaves of the tree, the 
# 15N signal changes due to translocation. 



################################### -- 
# (2) EXPLORE TRAIT RELATIONSHIPS 
################################### -- 

# Just doing visualizations right now but could test for significant trends if desired 

# can use just the 'traits' data 
traits <- merge(traits, env, by = "WFDP_Code")

# lazy fix for naming 
traits$Host_ID <- traits$Host_ID.x


# Define shapes for species 

            # ABAM, ABGR, ALRU, CONU, TABR, THPL, TSHE  
species_shapes <- c(15, 16, 17, 18, 7, 8, 9)

# leaf C and N
leaf_CN <- ggplot(traits, aes(x = leaf_pct_N, y = leaf_pct_C, color = Association)) +
  geom_point(size = 3, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

leaf_CN

# shows clustering of ALRU and CONU away from other species, makes sense since they 
# are the only broadleaf species 


# root C and N
root_CN <- ggplot(traits, aes(x = root_pct_N, y = root_pct_C, color = Association)) +
  geom_point(size = 3, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  labs(title = "", x = "Root N (%)", y = "Root C (%)") + 
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

root_CN

# no visible trends 


# whole plant C 

plant_C <- ggplot(traits, aes(x = leaf_pct_C, y = root_pct_C, color = Association)) +
  geom_point(size = 3, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  labs(title = "", x = "Leaf C (%)", y = "Root C (%)") + 
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

plant_C


# whole plant N 

plant_N <- ggplot(traits, aes(x = leaf_pct_N, y = root_pct_N, color = Association)) +
  geom_point(size = 3, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  labs(title = "", x = "Leaf N (%)", y = "Root N (%)") + 
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

plant_N


# SLA x SRL 

SLA_SRL <- ggplot(traits, aes(x = SLA_leaf, y = specific_root_length, color = Association)) +
  geom_point(size = 3, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  labs(title = "", x = "SLA (mm2/mg)", y = "SRL (cm/mg)") + 
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

SLA_SRL


# LDMC x RDMC

LDMC_RDMC <- ggplot(traits, aes(x = LDMC_leaf, y = root_dry_matter_cont, color = Association)) +
  geom_point(size = 3, aes(shape = sci_name)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "EM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "", x = "Leaf N (%)", y = "Leaf C (%)") +
  labs(title = "", x = "LDMC (mg/g)", y = "RDMC (mg/g)") + 
  theme(legend.title = element_text(colour="black", size=12, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(
    axis.text.x = element_text(size = 12, colour="black"),
    axis.text.y = element_text(size = 12, colour="black"),
    axis.title.y = element_text(size = 12, colour="black"),
    axis.title.x = element_text(size = 12, colour="black"), 
    strip.text = element_text(size = 12, colour="black")) +
  theme(legend.text = element_text(size = 12, colour="black"), 
        legend.title = element_text(size = 12, face = "bold", colour="black")) +
  theme(legend.position = "right")

LDMC_RDMC



## Can do more here, and could switch these to be elipses instead to show more species variation, 
# rather than doing it by mycorrhizal association. 

# LDMC x RDMC, for example, has a lot of clear clustering between species 



###### OPTIONAL #########

#check the spread of the raw trait data 
summary(traits_tree)

#need to get average traits for each species 
traits_tree <- traits_tree %>% group_by(Host_ID) %>%
  summarise(across(everything(), mean))


# traits do not have the same sample sizes, and should be standardized so they 
# are comparable 

#move species into the rownames 
traits_tree <- traits_tree %>%
  column_to_rownames(var = "Host_ID")


#scaling
traits_tree_sc <- scale(traits_tree)
summary(traits_tree_sc)

traits_tree_sc <- as.data.frame(traits_tree_sc)


# #could also scale all traits between zero and 1 
# 
# min_values <- as.numeric(lapply(as.data.frame(traits_tree), min))
# max_values <- as.numeric(lapply(as.data.frame(traits_tree), max))
# 
# traits_tree_minmax <- apply(traits_tree, 1, function(x) {
#   (x - min_values)/(max_values - min_values)
# })
# traits_tree_minmax <- t(traits_tree_minmax)
# summary(traits_tree_minmax)


# make a scaled dataset that doesn't have rownames 

scale_4_plots <- traits_tree_sc %>%
  rownames_to_column(var = "Host_ID")


# backtrack formatting a bit 
# move rownames back into a column 
long_traits_tree_sc <- traits_tree_sc %>%
  rownames_to_column(var = "Host_ID")

# Convert to long format for ggplot
long_traits_tree_sc <- pivot_longer(long_traits_tree_sc, cols = -Host_ID, names_to = "trait", values_to = "value")


#################################################################################

##################################### -- 
# (3) PRINCIPLE COMPONENT ANALYSIS
##################################### -- 

# perform a PCA to assess the axes of trait variation for the trees 
    # using unscaled data because the PCA does the scaling internally, and not using means 

# "traits" data has info columns in the first four positions, so need to exclude those 

#All Traits and individual tree positions 
alltraits.pca = prcomp(traits[5:19], center = T, scale = T)

sd.alltraits = alltraits.pca$sdev
loadings.alltraits = alltraits.pca$rotation
trait.names.alltraits = colnames(traits[5:19])
scores.alltraits = as.data.frame(alltraits.pca$x)
scores.alltraits$WFDP_Code = traits$WFDP_Code
summary(alltraits.pca)

# Save loadings for all traits
write.csv(loadings.alltraits, "./PCA/PCA_loadings_alltraits_no13C_tree.csv", row.names = TRUE)

#Save species scores
write.csv(scores.alltraits, "./PCA/PCA_scores_alltraits_no13C_tree.csv")


# PCA scores are 'scores.alltraits' with column for Host_ID

loadings.alltraits <- as.data.frame(loadings.alltraits)

# get proportion of variance explained to add to each axis label 
pca_var <- alltraits.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

#Merge in mycorrhizal association for plotting 
scores.alltraits <- merge(scores.alltraits, env, by = "WFDP_Code")

# Merge in scientific name for plotting 
scores.alltraits <- merge(scores.alltraits, taxa, by = "Host_ID")


# Change loadings names to something cleaner 
new_loadings <- c("SLA", "LDMC", "LMA", "Leaf_PctN", "Leaf_PctC", "Leaf_CN", "Leaf_d15N", 
                  "SRL", "SRA", "RDMC", "Root_CN", 
                  "Root_d15N", "RD", "Root_PctN", "Root_PctC")
rownames(loadings.alltraits) <- new_loadings

PCA_plot <- ggplot(scores.alltraits, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3.5, aes(shape = sci_name)) +
  geom_segment(data = loadings.alltraits, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.alltraits, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.alltraits)),
                  color = "black", size = 5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  guides(color = guide_legend(nrow = 3, ncol = 1), shape = guide_legend(nrow = 3, ncol = 3)) + 
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)"), 
       color = "Focal Species") +
  theme(legend.position = "bottom")  +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 14, face = "italic")) + 
  theme(
    axis.text.x = element_text(size = 14, colour="black"),
    axis.text.y = element_text(size = 14, colour="black"),
    axis.title.y = element_text(size = 14, colour="black"),
    axis.title.x = element_text(size = 14, colour="black"))

PCA_plot

## beautiful 


##Broken-Stick test for the significance of the loadings
print(alltraits.pca)

plot(alltraits.pca, type = "l")

ev = alltraits.pca$sdev^2

evplot = function(ev) {
  # Broken stick model (MacArthur 1957)
  n = length(ev)
  bsm = data.frame(j=seq(1:n), p=0)
  bsm$p[1] = 1/n
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p = 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

evplot(ev)


summary(alltraits.pca)

#can safely retain first two - PC1, PC2, and maybe the third one too 

# figure out which traits were most strongly retained on the first two PC's

# Get top 4 traits for PC1
top_PC1 <- loadings.alltraits[order(abs(loadings.alltraits$PC1), decreasing = TRUE), ][1:4, ]

# Get top 4 traits for PC2
top_PC2 <- loadings.alltraits[order(abs(loadings.alltraits$PC2), decreasing = TRUE), ][1:4, ]


# PC1 is showing a spread of Leaf_PctN, Leaf_CN, SLA, and LMA. This reflects conservative vs 
# acquisitive strategies, with more conservative needle-leaf species with more negative values, 
# and more acquisitive broadleaf species with more positive values. 

# PC2 is showing a spread of SRA, SRL, Root_CN, and Root_PctN. This is pretty much just the equivalent 
# of PC1 reflecting leaf strategies. This also shows a "Do-it-yourself" vs "Outsourcing" strategies trade-off. 

###

# Could get the eigenvectors for each host tree's position along these axes, which could 
# summarize its position in this functional space 

# Then could regress these with environmental properties in the plot. 
# Could also look at fungal community composition dispersion to see if this is 
# related to the traits of the host tree 

# Grab PC1 and PC2 for each tree 
alltrees_PCs <- dplyr::select(scores.alltraits, PC1, PC2, WFDP_Code, Host_ID, sci_name)

# Save this file for later 
write.csv(alltrees_PCs, "~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/tree_PC_no13C__scores.csv")


################################### -- 

# Perform separate PCAs for leaf and root traits 

# Grab environmental data again for some association labels 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

env <- dplyr::select(env, WFDP_Code, Association)

## Use the original leaf and root trait datasets 

## Leaves first: 

# Subset to traits excluding petiole 
leaf_sub <- dplyr::select(leaf, code, WFDP_Code, sub_plot, Host_ID, SLA_leaf, LDMC_leaf, LMA_leaf, 
                   leaf_pct_N, leaf_pct_C, leaf_CN, leaf_15N)

# Then roots: 
root_sub <- dplyr::select(root, code, WFDP_Code, sub_plot, Host_ID, specific_root_length, 
                   specific_root_area, root_dry_matter_cont, root_CN, root_15N, 
                   avg_root_dia, root_pct_N, root_pct_C)


################################### -- 

# PCAs

## Leaves
leaf.pca = prcomp(leaf_sub[5:11], center = T, scale = T)

sd.leaf = leaf.pca$sdev
loadings.leaf = leaf.pca$rotation
trait.names.leaf = colnames(leaf_sub[5:11])
scores.leaf = as.data.frame(leaf.pca$x)
scores.leaf$WFDP_Code = leaf_sub$WFDP_Code
scores.leaf$Host_ID = leaf_sub$Host_ID
summary(leaf.pca)

# Save loadings for leaf traits
write.csv(loadings.leaf, "./PCA/PCA_loadings_leaf_traits_no13C.csv", row.names = TRUE)

#Save species scores
write.csv(scores.leaf, "./PCA/PCA_scores_leaf_traits_no13C.csv")


# PCA scores are 'scores.leaf' with column for Host_ID

loadings.leaf <- as.data.frame(loadings.leaf)

# get proportion of variance explained to add to each axis label 
pca_var <- leaf.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage

#Merge in mycorrhizal association for plotting 
scores.leaf <- merge(scores.leaf, env, by = "WFDP_Code")

# Change loadings names to something cleaner 
new_loadings <- c("SLA", "LDMC", "LMA", "PctN", "PctC", "C:N", "d15N")
rownames(loadings.leaf) <- new_loadings


# Merge in scientific name for plotting 
scores.leaf <- merge(scores.leaf, taxa, by = "Host_ID")


## Removing legend from all plots to save for formatting 

# Visualize
PCA_plot_leaf <- ggplot(scores.leaf, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3.5, aes(shape = sci_name)) +
  geom_segment(data = loadings.leaf, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.leaf, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.leaf)),
                  color = "black", size = 6.5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  guides(shape = guide_legend(nrow = 4, ncol = 2)) + 
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)")) +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16, face = "italic")) + 
  theme(
  axis.text.x = element_text(size = 18, colour="black"),
  axis.text.y = element_text(size = 18, colour="black"),
  axis.title.y = element_text(size = 18, colour="black"),
  axis.title.x = element_text(size = 18, colour="black"))

PCA_plot_leaf

  
# Get top 3 traits for PC1
top_PC1_leaf <- loadings.leaf[order(abs(loadings.leaf$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_leaf <- loadings.leaf[order(abs(loadings.leaf$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in SLA, C:N ratio, and LMA, which matches what 
# would be expected for broad leaf vs needle-leaf species 

# PC2 is being driven by the leaf C content, d15N, and N content.  

## Together the axes explain 87.4% of the variation 

### Calculate average PC1 score for each species and plot, compare statistically 

PC1_leaf <- dplyr::select(scores.leaf, PC1, WFDP_Code, Host_ID, sci_name)

PC1_leaf_summary <- PC1_leaf %>%
  group_by(sci_name) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1   = sd(PC1, na.rm = TRUE),
    n        = n(),
    se_PC1   = sd_PC1 / sqrt(n)
  )


# Test for significant differences between species 
aov_PC1_leaf <- aov(PC1 ~ Host_ID, data = PC1_leaf)
summary(aov_PC1_leaf)

tuk_PC1_leaf <- TukeyHSD(aov_PC1_leaf)
tuk_PC1_leaf

# Significant differences between species

# 
#         Df Sum Sq Mean Sq F value Pr(>F)    
#    Host_ID      6 278.05   46.34   73.16 <2e-16 ***
#   Residuals   53  33.57    0.63                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# diff         lwr        upr     p adj
# ABGR-ABAM  0.1662578 -1.09318295  1.4256985 0.9996285
# ALRU-ABAM  5.7050870  4.36924837  7.0409256 0.0000000
# CONU-ABAM  5.6015496  4.51084195  6.6922573 0.0000000
# TABR-ABAM  2.1807431  1.06014752  3.3013387 0.0000042
# THPL-ABAM  0.7583286 -0.33237910  1.8490362 0.3504788
# TSHE-ABAM  1.3339967  0.24328905  2.4247044 0.0075804
# ALRU-ABGR  5.5388292  4.06200406  7.0156544 0.0000000
# CONU-ABGR  5.4352918  4.17585113  6.6947326 0.0000000
# TABR-ABGR  2.0144853  0.72907405  3.2998966 0.0002554
# THPL-ABGR  0.5920708 -0.66736993  1.8515115 0.7774055
# TSHE-ABGR  1.1677389 -0.09170177  2.4271797 0.0859092
# CONU-ALRU -0.1035374 -1.43937597  1.2323012 0.9999836
# TABR-ALRU -3.5243439 -4.88469531 -2.1639924 0.0000000
# THPL-ALRU -4.9467584 -6.28259702 -3.6109198 0.0000000
# TSHE-ALRU -4.3710903 -5.70692887 -3.0352517 0.0000000
# TABR-CONU -3.4208065 -4.54140209 -2.3002109 0.0000000
# THPL-CONU -4.8432211 -5.93392871 -3.7525134 0.0000000
# TSHE-CONU -4.2675529 -5.35826056 -3.1768453 0.0000000
# THPL-TABR -1.4224146 -2.54301014 -0.3018190 0.0049446
# TSHE-TABR -0.8467464 -1.96734199  0.2738492 0.2560444
# TSHE-THPL  0.5756682 -0.51503950  1.6663758 0.6719600


## Removing the 13C did not change any of the significant differences between species. 



# Put species in order for plotting 
spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
PC1_leaf_summary$sci_name <- factor(PC1_leaf_summary$sci_name, levels = spp_order)


# Add in mycorrhizal associations 
PC1_leaf_summary$Association <- c("ECM", "ECM", "DUAL", "AM", "DUAL", "ECM", "AM")


# Visualize 
PC1_leaf_plot <- ggplot(PC1_leaf_summary, aes(x = sci_name, y = mean_PC1, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  labs(title = "", x = "", y = "PCA Axis 1 - LES") +
  theme(legend.position = "none")  +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 18, colour="black", face = "italic", hjust = 1, angle = 45),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"))

PC1_leaf_plot


### Calculate average PC2 score for each species and plot, compare statistically 

PC2_leaf <- dplyr::select(scores.leaf, PC2, WFDP_Code, sci_name)

PC2_leaf_summary <- PC2_leaf %>%
  group_by(sci_name) %>%
  summarise(
    mean_PC2 = mean(PC2, na.rm = TRUE),
    sd_PC2   = sd(PC2, na.rm = TRUE),
    n        = n(),
    se_PC2   = sd_PC2 / sqrt(n)
  )



# Test for significant differences between species 
aov_PC2_leaf <- aov(PC2 ~ sci_name, data = PC2_leaf)
summary(aov_PC2_leaf)

tuk_PC2_leaf <- TukeyHSD(aov_PC2_leaf)
tuk_PC2_leaf

# Significant differences between species
# 

#                Df Sum Sq Mean Sq F value   Pr(>F)    
#     sci_name     6  37.82   6.303   29.39 3.27e-15 ***
#   Residuals   53  11.37   0.214                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# diff         lwr         upr     p adj
# A. grandis-A. amabilis         0.4700106 -0.26278680  1.20280803 0.4481112
# A. rubra-A. amabilis          -0.6265695 -1.40381858  0.15067949 0.1911490
# C. nuttallii-A. amabilis       1.0178844  0.38326321  1.65250557 0.0001737 ***
# T. brevifolia-A. amabilis     -1.4398115 -2.09182281 -0.78780025 0.0000002 *** 
# T. heterophylla-A. amabilis    0.7673514  0.13273027  1.40197262 0.0086005 ***
# T. plicata-A. amabilis         0.1615578 -0.47306341  0.79617895 0.9859313
# A. rubra-A. grandis           -1.0965802 -1.95586130 -0.23729903 0.0046387 ***
# C. nuttallii-A. grandis        0.5478738 -0.18492364  1.28067119 0.2674907
# T. brevifolia-A. grandis      -1.9098221 -2.65773038 -1.16191392 0.0000000 ***
# T. heterophylla-A. grandis     0.2973408 -0.43545659  1.03013824 0.8737177
# T. plicata-A. grandis         -0.3084528 -1.04125026  0.42434457 0.8535222
# C. nuttallii-A. rubra          1.6444539  0.86720490  2.42170297 0.0000006 ***
# T. brevifolia-A. rubra        -0.8132420 -1.60475366 -0.02173032 0.0403433 *
# T. heterophylla-A. rubra       1.3939210  0.61667196  2.17117002 0.0000226 ***
# T. plicata-A. rubra            0.7881273  0.01087828  1.56537634 0.0448494 *
# T. brevifolia-C. nuttallii    -2.4576959 -3.10970720 -1.80568464 0.0000000 ***
# T. heterophylla-C. nuttallii  -0.2505329 -0.88515412  0.38408823 0.8873040 
# T. plicata-C. nuttallii       -0.8563266 -1.49094780 -0.22170544 0.0023067 **
# T. heterophylla-T. brevifolia  2.2071630  1.55515170  2.85917426 0.0000000 ***
# T. plicata-T. brevifolia       1.6013693  0.94935802  2.25338058 0.0000000 ***
# T. plicata-T. heterophylla    -0.6057937 -1.24041485  0.02882750 0.0704307


## UPDATE: Lots of new differences between species 


# Put species in order for plotting 
spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
PC2_leaf_summary$sci_name <- factor(PC2_leaf_summary$sci_name, levels = spp_order)


# Add in mycorrhizal associations 
PC2_leaf_summary$Association <- c("ECM", "ECM", "DUAL", "AM", "DUAL", "ECM", "AM")


# Visualize 
PC2_leaf_plot <- ggplot(PC2_leaf_summary, aes(x = sci_name, y = mean_PC2, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  labs(title = "", x = "", y = "PCA Axis 2") +
  theme(legend.position = "none")  +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 18, colour="black", face = "italic", hjust = 1, angle = 45),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"), 
    strip.text = element_text(size = 18, colour="black"))

PC2_leaf_plot

############### -- 

## Roots
root.pca = prcomp(root_sub[5:12], center = T, scale = T)

sd.root = root.pca$sdev
loadings.root = root.pca$rotation
trait.names.root = colnames(root_sub[5:12])
scores.root = as.data.frame(root.pca$x)
scores.root$WFDP_Code = root_sub$WFDP_Code
scores.root$Host_ID = root_sub$Host_ID
summary(root.pca)

# Save loadings for root traits
write.csv(loadings.root, "./PCA/PCA_loadings_root_traits_no13C.csv", row.names = TRUE)

#Save species scores
write.csv(scores.root, "./PCA/PCA_scores_root_traits_no13C.csv")


# PCA scores are 'scores.root' with column for Host_ID

# The PC1 axis is the equivalent of the leaf conservative-acquisitive axis, but right now the direction is 
# flipped so it's less intuitive to interpret them together. Going to multiply the PC1 values by -1 to 
# flip the orientation, and this doesn't do anything to the actual PCA, or the PC2 axis. 

# Flip PC1 scores
root.pca$x[, "PC1"] <- -1 * root.pca$x[, "PC1"]

# Flip PC1 loadings
root.pca$rotation[, "PC1"] <- -1 * root.pca$rotation[, "PC1"]

# Load back into items for plotting 
scores.root   <- as.data.frame(root.pca$x)
scores.root$WFDP_Code = root_sub$WFDP_Code
scores.root$Host_ID = root_sub$Host_ID

loadings.root <- as.data.frame(root.pca$rotation)

# Change loadings names to something cleaner 
new_loadings_root <- c("SRL", "SRA", "RDMC", "C:N", "d15N", "RD", "PctN", "PctC")
rownames(loadings.root) <- new_loadings_root

#Merge in mycorrhizal association for plotting 
scores.root <- merge(scores.root, env, by = "WFDP_Code")


# Merge in scientific name for plotting 
scores.root <- merge(scores.root, taxa, by = "Host_ID")


# get proportion of variance explained to add to each axis label 
pca_var <- root.pca$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


# Visualize
PCA_plot_root <- ggplot(scores.root, aes(x = PC1, y = PC2, color = Association)) +
  geom_point(size = 3.5, aes(shape = sci_name)) +
  geom_segment(data = loadings.root, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.root, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.root)),
                  color = "black", size = 6.5, max.overlaps = 10) +
  theme_minimal() +
  scale_color_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal\nAssociation") + 
  scale_shape_manual(
    values = species_shapes, 
    breaks = c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla"), 
    name = "Focal Species",  
    labels=c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", "T. plicata", "T. heterophylla")) +
  labs(title = "",
       x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(
    axis.text.x = element_text(size = 18, colour="black"),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"), 
    strip.text = element_text(size = 18, colour="black"))

PCA_plot_root


# Get top 3 traits for PC1
top_PC1_root <- loadings.root[order(abs(loadings.root$PC1), decreasing = TRUE), ][1:3, ]

# Get top 3 traits for PC2
top_PC2_root <- loadings.root[order(abs(loadings.root$PC2), decreasing = TRUE), ][1:3, ]


# PC1 is being driven by variation in Specific root area, specific root length, and root N content. 

# PC2 is being driven by root C content, d15N, and root diameter. 

## Together the axes explain 70.9% of the variation 


### Calculate average PC1 score for each species and plot, compare statistically 

PC1_root <- dplyr::select(scores.root, PC1, WFDP_Code, Host_ID, sci_name)

PC1_root_summary <- PC1_root %>%
  group_by(sci_name) %>%
  summarise(
    mean_PC1 = mean(PC1, na.rm = TRUE),
    sd_PC1   = sd(PC1, na.rm = TRUE),
    n        = n(),
    se_PC1   = sd_PC1 / sqrt(n)
  )


# Test for significant differences between species 
aov_PC1_root <- aov(PC1 ~ Host_ID, data = PC1_root)
summary(aov_PC1_root)

tuk_PC1_root <- TukeyHSD(aov_PC1_root)
tuk_PC1_root


# Df Sum Sq Mean Sq F value  Pr(>F)   
# Host_ID      6   80.3  13.384   4.002 0.00224 **
#   Residuals   53  177.2   3.344                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# $Host_ID
# diff        lwr        upr     p adj
# ABGR-ABAM -0.68740481 -3.5813023  2.2064927 0.9902084
# ALRU-ABAM  2.93084009 -0.1386018  6.0002819 0.0702911
# CONU-ABAM  0.11120033 -2.3949884  2.6173891 0.9999994
# TABR-ABAM -1.02174948 -3.5966137  1.5531147 0.8848469
# THPL-ABAM  1.25183901 -1.2543498  3.7580278 0.7253808
# TSHE-ABAM -0.98788271 -3.4940715  1.5183061 0.8880183
# ALRU-ABGR  3.61824490  0.2248493  7.0116405 0.0295168 *
# CONU-ABGR  0.79860514 -2.0952924  3.6925027 0.9787438
# TABR-ABGR -0.33434467 -3.2879165  2.6192271 0.9998477
# THPL-ABGR  1.93924382 -0.9546537  4.8331414 0.3947193
# TSHE-ABGR -0.30047790 -3.1943754  2.5934196 0.9999081
# CONU-ALRU -2.81963976 -5.8890816  0.2498021 0.0913045
# TABR-ALRU -3.95258957 -7.0783562 -0.8268230 0.0051721 **
# THPL-ALRU -1.67900108 -4.7484429  1.3904408 0.6344123
# TSHE-ALRU -3.91872280 -6.9881647 -0.8492809 0.0046157 **
# TABR-CONU -1.13294981 -3.7078140  1.4419144 0.8258306
# THPL-CONU  1.14063868 -1.3655501  3.6468275 0.8021808
# TSHE-CONU -1.09908304 -3.6052718  1.4071057 0.8280304
# THPL-TABR  2.27358849 -0.3012757  4.8484527 0.1168157
# TSHE-TABR  0.03386677 -2.5409974  2.6087310 1.0000000
# TSHE-THPL -2.23972172 -4.7459105  0.2664671 0.1086208


## UPDATE: Removing 13C did not change any of the significant differences between species 


# Put species in order for plotting 
spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
PC1_root_summary$sci_name <- factor(PC1_root_summary$sci_name, levels = spp_order)


# Add in mycorrhizal associations 
PC1_root_summary$Association <- c("ECM", "ECM", "DUAL", "AM", "DUAL", "ECM", "AM")



# Visualize 
PC1_root_plot <- ggplot(PC1_root_summary, aes(x = sci_name, y = mean_PC1, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC1 - se_PC1, ymax = mean_PC1 + se_PC1), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 1 - RES") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(
    axis.text.x = element_text(size = 18, colour="black", face = "italic", hjust = 1, angle = 45),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"), 
    strip.text = element_text(size = 18, colour="black"))

PC1_root_plot



### Calculate average PC2 score for each species and plot, compare statistically 

PC2_root <- dplyr::select(scores.root, PC2, WFDP_Code, Host_ID, sci_name)

PC2_root_summary <- PC2_root %>%
  group_by(sci_name) %>%
  summarise(
    mean_PC2 = mean(PC2, na.rm = TRUE),
    sd_PC2   = sd(PC2, na.rm = TRUE),
    n        = n(),
    se_PC2   = sd_PC2 / sqrt(n)
  )


# Test for significant differences between species 
aov_PC2_root <- aov(PC2 ~ Host_ID, data = PC2_root)
summary(aov_PC2_root)

tuk_PC2_root <- TukeyHSD(aov_PC2_root)
tuk_PC2_root

# SIGNIFICANT differences between species

# Df Sum Sq Mean Sq F value   Pr(>F)    
# Host_ID      6  29.39   4.898   5.486 0.000178 ***
#   Residuals   53  47.32   0.893                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# $Host_ID
# diff        lwr         upr     p adj
# ABGR-ABAM -1.17977226 -2.6749467  0.31540218 0.2116233
# ALRU-ABAM -0.51600821 -2.1018802  1.06986376 0.9524521
# CONU-ABAM  0.61083601 -0.6840230  1.90569505 0.7746276
# TABR-ABAM  0.55119178 -0.7791494  1.88153298 0.8625169
# THPL-ABAM -0.87735666 -2.1722157  0.41750238 0.3813392
# TSHE-ABAM  0.74446677 -0.5503923  2.03932581 0.5788524
# ALRU-ABGR  0.66376405 -1.0894834  2.41701148 0.9056858
# CONU-ABGR  1.79060826  0.2954338  3.28578270 0.0095445 **
# TABR-ABGR  1.73096403  0.2049580  3.25697005 0.0166594 *
# THPL-ABGR  0.30241560 -1.1927588  1.79759003 0.9958773
# TSHE-ABGR  1.92423902  0.4290646  3.41941346 0.0041904 **
# CONU-ALRU  1.12684422 -0.4590278  2.71271619 0.3249093
# TABR-ALRU  1.06719999 -0.5477730  2.68217295 0.4117011
# THPL-ALRU -0.36134845 -1.9472204  1.22452352 0.9921491
# TSHE-ALRU  1.26047498 -0.3253970  2.84634695 0.2045507
# TABR-CONU -0.05964423 -1.3899854  1.27069697 0.9999994
# THPL-CONU -1.48819266 -2.7830517 -0.19333362 0.0146296 *
# TSHE-CONU  0.13363076 -1.1612283  1.42848980 0.9999113
# THPL-TABR -1.42854843 -2.7588896 -0.09820723 0.0277404 *
# TSHE-TABR  0.19327499 -1.1370662  1.52361619 0.9993557
# TSHE-THPL  1.62182342  0.3269644  2.91668247 0.0057828 **

# Now a lot of interesting differences between species 


# Put species in order for plotting 
spp_order <- c("A. amabilis", "A. grandis", "A. rubra", "C. nuttallii", "T. brevifolia", 
               "T. plicata", "T. heterophylla")


# Convert the column to a factor with the specified levels
PC2_root_summary$sci_name <- factor(PC2_root_summary$sci_name, levels = spp_order)


# Add in mycorrhizal associations 
PC2_root_summary$Association <- c("ECM", "ECM", "DUAL", "AM", "DUAL", "ECM", "AM")


# Visualize 
PC2_root_plot <- ggplot(PC2_root_summary, aes(x = sci_name, y = mean_PC2, fill = Association)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_PC2 - se_PC2, ymax = mean_PC2 + se_PC2), width = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("AM" = "#dc267f", "DUAL" = "#648fff", "ECM" = "#ffb000"), name = "Mycorrhizal Association") + 
  labs(title = "", x = "", y = "PCA Axis 2") +
  theme(legend.position = "none")  +
  theme(legend.title = element_text(colour="black", size=16, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) + 
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(
    axis.text.x = element_text(size = 18, colour="black", face = "italic", hjust = 1, angle = 45),
    axis.text.y = element_text(size = 18, colour="black"),
    axis.title.y = element_text(size = 18, colour="black"),
    axis.title.x = element_text(size = 18, colour="black"), 
    strip.text = element_text(size = 18, colour="black"))

PC2_root_plot


#################################################################################

############################################ -- 
# (4) ENVIRONMENTAL RELATIONSHIPS WITH PCAs
############################################ -- 

# Read back in environmental data for the trees 

# Right now this is just slope, aspect, and elevation, but will be able to pull in 
    # krieged data when it exists 
env <- read.csv("~/Dropbox/WSU/WFDP_Chapter_3_Project/Enviro_Data/WFDP_enviro_data_all.csv")

# subset to data of interest 

env <- dplyr::select(env, Cell, WFDP_Code, slope, aspect, elevation_m)

tree_PCs_env <- merge(env, alltrees_PCs, by = "WFDP_Code")

#### Assess relationships between PC1 and PC2 values and slope, aspect, and elevation 

### PC1

## Slope
slope <- ggplot(tree_PCs_env, aes(x = slope, y = PC1)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope

# test relationships 
lm_slope <- lm(PC1 ~ slope, data = tree_PCs_env)
summary(lm_slope) # NOT SIGNIFICANT


## Aspect
aspect <- ggplot(tree_PCs_env, aes(x = aspect, y = PC1)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect

# test relationships 
lm_aspect <- lm(PC1 ~ aspect, data = tree_PCs_env)
summary(lm_aspect) # NOT SIGNIFICANT


## Elevation
elev <- ggplot(tree_PCs_env, aes(x = elevation_m, y = PC1)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
elev

# test relationships 
lm_elev <- lm(PC1 ~ elevation_m, data = tree_PCs_env)
summary(lm_elev) # NOT SIGNIFICANT


### PC2

## Slope
slope2 <- ggplot(tree_PCs_env, aes(x = slope, y = PC2)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
slope2

# test relationships 
lm_slope2 <- lm(PC2 ~ slope, data = tree_PCs_env)
summary(lm_slope2) # SIGNIFICANT

# Adjusted R-squared:  0.086, p-value: 0.013


## Aspect
aspect2 <- ggplot(tree_PCs_env, aes(x = aspect, y = PC2)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
aspect2

# test relationships 
lm_aspect2 <- lm(PC1 ~ aspect, data = tree_PCs_env)
summary(lm_aspect2) # NOT SIGNIFICANT


## Elevation
elev2 <- ggplot(tree_PCs_env, aes(x = elevation_m, y = PC2)) +
  geom_point(aes(color = Host_ID)) +
  geom_smooth(method = "lm") +
  theme_minimal()
elev2

# test relationships 
lm_elev2 <- lm(PC2 ~ elevation_m, data = tree_PCs_env)
summary(lm_elev2) # NOT SIGNIFICANT


## RESULT: Only significant relationship was between slope and PC2 values 

# This is just for the overall PCA, could also do it with the separate leaf and root 
# PCAs if desired 

#################################################################################

## Score csvs are saved for leaf and root traits already. 

# Save final plots

# All traits plot 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/all_traits_biplot.png", 
       plot = PCA_plot, width = 10, height = 8, units = "in", dpi = 300)

# leaf PCA
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/leaf_traits_biplot.png", 
       plot = PCA_plot_leaf, width = 8, height = 6.5, units = "in", dpi = 300)

# PC1 leaf bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC1_leaf_barplot.png", 
       plot = PC1_leaf_plot, width = 7, height = 6.5, units = "in", dpi = 300)

# PC2 leaf bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC2_leaf_barplot.png", 
       plot = PC2_leaf_plot, width = 7, height = 6.5, units = "in", dpi = 300)

# root PCA
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/root_traits_biplot.png", 
       plot = PCA_plot_root, width = 8, height = 6.5, units = "in", dpi = 300)

# PC1 root bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC1_root_barplot.png", 
       plot = PC1_root_plot, width = 7, height = 6.5, units = "in", dpi = 300)

# PC2 root bar chart 
ggsave("~/Dropbox/WSU/WFDP_Chapter_3_Project/Trait_Data/PCA/PC2__root_barplot.png", 
       plot = PC2_root_plot, width = 7, height = 6.5, units = "in", dpi = 300)


## -- END -- ## 

