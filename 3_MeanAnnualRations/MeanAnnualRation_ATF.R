# This script file includes code necessary to calculate mean annual rations for Arrowtooth Flounder (g per g per yr). Bioenergetics model parameters used to estimate maximum daily consumption (Cmax; g per g per d), the temperature scaling functions, and mean relative foraging rates (RFR) were obtained from Holsman and Aydin (2015). Mean region- and size-specific number of foraging days per year were also specified in Holsman and Aydin (2015). 

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Beaudreau, A. H., and T. E. Essington. 2009. Development of a new field-based approach for estimating consumption rates of fishes and comparison with a bioenergetics model for lingcod (Ophiodon elongatus). Canadian Journal of Fisheries and Aquatic Sciences 66:565−578.
# Holsman, K. K., and K. Aydin. 2015. Comparative methods for evaluating climate change impacts on the foraging ecology of Alaskan groundfish. Marine Ecology Progress Series 521:217–235.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##########################################################
### INITIAL DATA PREPARATION ###
##########################################################
# Read in and format food habits data for Arrowtooth Flounder:
FHdata = read.csv("1_Data/GOA_RawPP_ATFdiets.csv") 

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to survey data below):
FHdata$Haul_Join = paste(FHdata$Vessel, FHdata$Cruise, FHdata$Haul, sep="")

# Exclude data prior to 1990, when survey methods were standardized:
FHdata = subset(FHdata, Year >= 1990)

# Select predator sizes of interest (ATF >= 19 cm):
FHdata = subset(FHdata, Pred_len >= 19)

##########################################################
### CONSTRUCT UNIFORM GRID FOR SPATIALLY-EXPLICIT ANALYSES 
##########################################################
# Create a uniform grid spanning the study area (to aggregate diet data):
require(tidyr)
require(raster)
require(rgeos)
require(rgbif)
require(viridis)
require(gridExtra)
require(rasterVis)
require(purrr)
require(mapproj)
require(devtools)
require(stringr)
require(maptools)
require(rgdal)
require(ggplot2)
require(ggmap)
require(dplyr)

# Read in and prepare shapefile for INPFC statistical areas (610 to 650). 
INPFC = readOGR(".", "GOA_Shapes")
INPFC_Pr = spTransform(INPFC, CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs")) 
  # UTM projection maintains cosistent grid cell area regardless of geographic location.

INPFCdata = data.frame()
INPFCdata = rbind(INPFCdata, INPFC_Pr@data)
INPFCdata$OBJECTID = as.character(INPFCdata$OBJECTID)
INPFCdata$REP_AREA = as.character(INPFCdata$REP_AREA)
INPFCdata$Region = NA; INPFCdata$Region = "GOA"

INPFC_Pr@data$OBJECTID = as.character(INPFC_Pr@data$OBJECTID)
INPFC_Pr@data = full_join(INPFC_Pr@data, INPFCdata, by = "REP_AREA")
row.names(INPFC_Pr) = row.names(INPFC_Pr@data)
INPFC_Pr = spChFIDs(INPFC_Pr, row.names(INPFC_Pr))
INPFC_Pr = gUnaryUnion(INPFC_Pr, id = INPFC_Pr@data$Region)
row.names(INPFC_Pr) = as.character(1:length(INPFC_Pr))

INPFCdata = unique(INPFCdata$Region)
INPFCdata = as.data.frame(INPFCdata)
colnames(INPFCdata) = "Region"  
INPFC_Pr = SpatialPolygonsDataFrame(INPFC_Pr, INPFCdata)

# Set size for square grid cells (50 km x 50 km): 
my.interval=50 

# Specify grid boundaries (UTM projection):
lon_min = -875
lon_max = 1975
lat_min = 5075
lat_max = 7975

# Compile series of points for uniform grid:
mygrd = expand.grid(
  LON = seq(lon_min, lon_max, by=my.interval),
  LAT = seq(lat_min, lat_max, by=my.interval)) %>% 
  mutate(my.z=1:n()) %>% 
  data.frame

# Convert mygrd to a spatial dataframe:
coordinates(mygrd) = ~ LON + LAT

# Convert SpatialPoints to a SpatialPixelsDataFrame:
mygrd = as(SpatialPixelsDataFrame(mygrd, mygrd@data, tolerance=0.00086), "SpatialPolygonsDataFrame")

# Project, clip, and reproject mygrid to INPFC statistical areas:
proj4string(mygrd) = proj4string(INPFC_Pr)
INPFC_Pr = gBuffer(INPFC_Pr, byid=T, width=0)
mygrd = gBuffer(mygrd, byid=T, width=0)
clip_INPFC = gIntersection(INPFC_Pr, mygrd, byid = T, drop_lower_td = T)
proj4string(clip_INPFC) = proj4string(INPFC_Pr)
plot(clip_INPFC, col="grey") # check

# Convert grid to a data frame (in decimal degrees):
clip_INPFC_dd = clip_INPFC
clip_INPFC_dd = spTransform(clip_INPFC_dd, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
 
# Find where food habits data overlap with clipped grid and convert back to a data frame:
coordinates(FHdata) = ~ Rlong + Rlat 
proj4string(FHdata) = proj4string(clip_INPFC_dd)

tempdat = data.frame(myrows=names(over(FHdata, clip_INPFC_dd)), mygrid=over(FHdata, clip_INPFC_dd))
FHdata$id2 = over(FHdata, clip_INPFC_dd)
FHdata = as.data.frame(FHdata)

# Create a data frame with grid cell identifiers (id2) and coordinates (vertices and centers):
mycenter_INPFC = as.data.frame(gCentroid(clip_INPFC, byid=TRUE)) %>% 
  mutate(id2=1:n(),
         EEZgrid=rownames(.))
colnames(mycenter_INPFC)[3] = "id"
xy=mycenter_INPFC[,c(1,2)]
mycenter_INPFC_sp = SpatialPointsDataFrame(coords=xy, data=mycenter_INPFC, proj4string=CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs"))
mycenter_INPFC_Pr = spTransform(mycenter_INPFC_sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
mycenter_INPFC_df = as.data.frame(mycenter_INPFC_Pr) 
names(mycenter_INPFC_df) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")

# Join food habits and geographic data:
FHdata_grid = mycenter_INPFC_df %>% right_join(FHdata)

##########################################################
### BIOENGERTICS CALCULATIONS ### 
##########################################################
FHdata_grid = subset(FHdata_grid, Gear_temp > 0) 

# Parameter estimates from Holsman and Aydin (2015):
Ca = 0.125
Cb = -0.199
Cq = 2.497
Tco = 20.512
Tcm = 26

Y = (log(Cq)) * (Tcm - Tco + 2)
Z = (log(Cq)) * (Tcm - Tco)
X = (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
FHdata_grid$V = (Tcm - FHdata_grid$Gear_temp)/(Tcm - Tco)

# Compute temperature scaling function:
FHdata_grid$fT = with(FHdata_grid, V^X * (exp(X*(1-V))))

# Calculate maximum daily consumption (Cmax; g/g/d):
FHdata_grid$Cmax_ggd = Ca * (FHdata_grid$Pred_wt^Cb) * FHdata_grid$fT

# Assume fish are feeding at their theoretical maximum (RFR = 1): 
    # RFR, juv(< 40 cm): 0.79 # Holsman and Aydin 2015
    # RFR, adult(>= 40 cm): 1.07 # Holsman and Aydin 2015

# Multiply Cmax by effective number of foraging days (EFD; Holsman and Aydin 2015) to scale to year:
  # EFD ATF <  40 cm = 346 d
  # EFD ATF >= 40 cm = 306 d
FHdata_grid$EFD = ifelse(FHdata_grid$Pred_len < 40, 346, 306)
FHdata_grid$Cmax_ggy = FHdata_grid$Cmax_ggd * FHdata_grid$EFD
FHdata_grid = subset(FHdata_grid, Cmax_ggy !="Inf")

# Remove INPFC statistical areas 649 and 659 (inside waters of Southeast Alaska):
FHdata_grid = subset(FHdata_grid, INPFC_Area != "649" & INPFC_Area != "659")
saveRDS(FHdata_grid, "3_MeanAnnualRations/Rations_ATF.rds")

##########################################################
### MEAN ANNUAL RATIONS ###
##########################################################
# Calculate mean and sd Cmax (g/g/y) for each survey year:
Cmax_ggy_ATF_mean = FHdata_grid %>%
  group_by(Year) %>%
  summarise(meanCmax_ggy = mean(Cmax_ggy)) 
Cmax_ggy_ATF_sd = FHdata_grid %>%
  group_by(Year) %>%
  summarise(sdCmax_ggy = sd(Cmax_ggy)) 
Cmax_ggy_ATF = merge(Cmax_ggy_ATF_mean, Cmax_ggy_ATF_sd)

# Calculate 95% confidence intervals:
Cmax_ggy_ATF$CI_Hi = Cmax_ggy_ATF$meanCmax_ggy + (1.96*Cmax_ggy_ATF$sdCmax_ggy)
Cmax_ggy_ATF$CI_Lo = Cmax_ggy_ATF$meanCmax_ggy - (1.96*Cmax_ggy_ATF$sdCmax_ggy)

# Plot results (part of Fig. 2):
MeanAnnualRation_ATF = ggplot(Cmax_ggy_ATF, aes(x = Year, y = meanCmax_ggy)) +
  geom_line(col="lemonchiffon", lwd=1) +
  geom_point(size = 3.5, col="yellow2", stroke=2, shape=3, alpha = 1) +
  # geom_errorbar(aes(ymin=meanC_ggy-CIlo, ymax=meanC_ggy+CIhi), color="goldenrod1", width=0.5) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "vertical", legend.position = "right", legend.box.margin = margin(0,0,0,5), legend.background = element_rect(fill="transparent")) +
  theme(legend.position = c(0.0445, 0.865)) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(1.0, "mm")) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(legend.title = element_text(family="Arial", color="black", size=11)) +
  theme(legend.text = element_text(family="Arial", color="black", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=12)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=4, size=12)) +
  xlab("") +
  ylab("Mean Annual Ration (g per g per y)") +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=c(0,2,4,6,8)) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

MeanAnnualRation_ATF