# This script file includes code necessary to estimate relative densities of Sablefish in the Gulf of Alaska. Methods are described in Barnes et al. (2018) and were modified from Hunsicker et al. (2013) and Shelton et al. (2017). Relative predator densities (kg per ha) were multiplied by total predator biomass, mean annual rations, and age-specific proportions of pollock consumed to estimate year-, area-, and predator-specific predation.

# Standardized survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# Below, we use generalized additive mixed models (GAMMs) to quantify probabilities of occurrence and predator densities across a uniform grid spanning the study area. We included survey year, tow location (latitude and longitude), depth (m), and bottom temperature (degrees C) as model covariates and allowed a spatial autocorrelation term to vary by year. We multiplied year- and grid cell-specific probabilities of occurrence and predicted densities from best-fit models to estimate predator biomass. These values were then normalized to estimate relative predator density in each unique combination of survey year and grid cell.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Barnes, C. L., A. H. Beaudreau, M. E. Hunsicker, and L. Ciannelli (2018). Assessing the potential for competition between Pacific Halibut (Hippoglossus stenolepis) and Arrowtooth Flounder (Atheresthes stomias) in the Gulf of Alaska. PLoS ONE 13(12):e0209402.
# Hanselman, D. H., C. J. Rodgveller, C. R. Lunsford, and K. H. Fenske. 2017. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report 327–502.
# Hunsicker, M. E., L. Ciannelli, K. M. Bailey, S. Zador, and L. Stige. 2013. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLoS ONE 8(6):e66025.
# Shelton, A. O., M. E. Hunsicker, E. J. Ward, B. E. Feist, R. Blake, C. L. Ward, et al. 2017. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES Journal of Marince Science doi:10.1093/icesjms/fsx079.
# Sigler, M. F., and H. H. Zenger, Jr. 1989. Assessment of Gulf of Alaska Sablefish and other groundfish based on the domestic longline survey, 1987. NOAA Technical Memorandum NMFS-AFSC Report 169. 
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
#########################################################
### INITIAL DATA PREPARATION ###
#########################################################
# Prepare and format AFSC longline survey data:
# These data include all sets conducted between 1979 and 2017.
LL = read.csv(unz("1_Data/SBL_catch_summary_view_with_nulls.csv.zip", "SBL_catch_summary_view_with_nulls.csv"), skip=6, header=T, na.strings=c("","NA")) 

# Eliminate extraneous survey years (for comparability with AFSC bottom trawl survey data):
LLsub = subset(LL, Year >= 1990 & Year < 2018)

require(dplyr)
# Eliminate areas outside of the Gulf of Alaska:
LLsub = subset(LLsub, subset = NMFS.Mgmt.Area %in% c("610", "620", "630", "640", "650"))

# Remove remaining Japanese surveys:
LLsub = subset(LLsub, Survey.Country != "Japan")

# Create unique identifier by concatenating year, station, and depth stratum: 
LLsub$UniqueID = paste(LLsub$Year, LLsub$Station.Number, sep="_")

# Calculate mean depth per set (m):
LLsub$Starting.Depth..m. = as.numeric(as.character(LLsub$Starting.Depth..m.))
LLsub$Ending.Depth..m. = as.numeric(as.character(LLsub$Ending.Depth..m.))
LLsub$Depth_m = (LLsub$Starting.Depth..m. + LLsub$Ending.Depth..m.) / 2

# Set NAs to 0:
LLsub$SumFreq = LLsub$SumFreq
LLsub$SumFreq[is.na(LLsub$SumFreq)] = 0

# Remove unnecessary columns and rename:
LLsub = (LLsub[,c("Year", "Station.Number", "NMFS.Mgmt.Area", "Start.Latitude..DD.", "Start.Longitude..DD.", "Depth.Stratum", "SumFreq", "UniqueID", "Depth_m")])
names(LLsub)[2] = "Station"
names(LLsub)[3] = "StatArea"

# Sum numbers of fish and calculate mean depth (m) by Unique ID:
LL_enviro = LLsub %>%
  group_by(UniqueID) %>%
  mutate(sumFish = sum(SumFreq)) %>%
  mutate(Lat = mean(Start.Latitude..DD.)) %>%
  mutate(Lon = mean(Start.Longitude..DD.)) %>%
  mutate(Depth = mean(Depth_m))

LLenviro = as.data.frame(distinct(LL_enviro[,c("Year", "Station", "StatArea", "UniqueID", "sumFish", "Lat", "Lon", "Depth")]))

#########################################################
### INCREASE COMPARABILITY WITH TOTAL BIOMASS ESTIMATES 
#########################################################
# Adjust set-specific CPUE estimates (kg per station) to exclude fish encompassed within assessment-based estimates of total biomass (Sablefish < 2 yr or < 45 cm). 

# Read in and format length data from AFSC longline survey:
lengths = read.csv(unz("1_Data/SBL_length_summary_view.csv.zip", "SBL_length_summary_view.csv"), header=T, skip=6) 

# Eliminate extraneous survey years and locations:
lengths = subset(lengths, Year >= 1990 & Year < 2018)
lengths = subset(lengths, subset = NMFS.Area.Code %in% c("610", "620", "630", "640", "650"))

# Remove remaining Japanese surveys:
lengths = subset(lengths, Country != "Japan")

# Select and rename columns:
lengths = lengths[,c("Year", "Station.Number", "NMFS.Area.Code", "Depth.Stratum", "Sex", "Length", "Frequencey")]
names(lengths) = c("Year", "Station", "StatArea", "Stratum", "Sex", "Length.cm", "No_Fish")

# Create unique identifier by concatenating Year, Station, and Depth Stratum:
lengths$UniqueID = paste(lengths$Year, lengths$Station, sep="_")

# Create fork length bins (cm):
lengths$Length.cm = as.numeric(as.character(lengths$Length.cm))
lengths$FLBin = cut(lengths$Length.cm, breaks = c(0,14,24,34,44,54,64,74,84,94,120))
levels(lengths$FLBin) = c("<=14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-94", ">=95")

#########################################################
# Calculate proportions of fish (by weight) sampled in each length bin and set:

# Estimate weight (kg) from length (cm) (parameter estimates from Hanselman et al. 2007 - Table 5; 1 = male, 2 = female):
lengths$WT_est = with(lengths, 
ifelse(Sex == "1", (0.0000124 * Length.cm^2.960),
ifelse(Sex == "2", (0.0000101 * Length.cm^3.015), NA)))

# Calculate mean weight at length for fish with no sex identified:
lengths = lengths %>%
  group_by(Length.cm) %>%
  mutate(WT_calc = mean(na.omit(WT_est)))

lengths$WT = lengths$WT_est
lengths$WT = with(lengths,
  ifelse(is.na(WT), WT_calc, WT_est))
lengths = na.omit(lengths)

# Calculate total weight (kg; all size classes) per year, station, and depth strata:
lengths$No_Fish = as.numeric(as.character(lengths$No_Fish))
LL_all = lengths %>%
  group_by(UniqueID) %>%
  mutate(WT_bin = WT * No_Fish) %>%
  summarise(Total_WT = sum(WT_bin))

# Limit data to sizes of interest (SBL >= 45 cm):
lengths_red = subset(lengths, Length.cm >= 45)

LL_45 = lengths_red %>%
  group_by(UniqueID) %>%
  mutate(WT_bin = WT * No_Fish) %>%
  summarise(Recruit_WT = sum(WT_bin))

# Calculate proportion of set measuring >= 45 cm:
SBLprop = as.data.frame(LL_all %>% full_join(LL_45)); SBLprop[is.na(SBLprop)] = 0
SBLprop$SBLprop = SBLprop$Recruit_WT / SBLprop$Total_WT

#########################################################
# Join setline survey and proportional length data:
LL_prop = LLenviro %>% left_join(SBLprop)

# Assign weights with NAs to zero when no fish were caught:
LL_prop$Total_WT[is.na(LL_prop$Total_WT) & LL_prop$sumFish==0] = 0
LL_prop$Recruit_WT[is.na(LL_prop$Recruit_WT) & LL_prop$sumFish==0] = 0
LL_prop$SBLprop[is.na(LL_prop$SBLprop) & LL_prop$sumFish==0] = 0

# Remove sets with weights, but no fish:
LL_prop = subset(LL_prop, subset =!(sumFish == 0 & Total_WT > 0))

# Remove stations without lengthed fish:
LL_prop = na.omit(LL_prop)

# Adjust CPUE (kg per station) based on proportion of catch >= 45 cm:
LL_prop$adjCPUE = LL_prop$Total_WT * LL_prop$SBLprop

# Remove extreme outlier (n = 1, depth > 1000 m):
LL_prop = subset(LL_prop, Depth < 1000)

#########################################################
### MODEL FITTING AND PLOTTING ###
#########################################################
require(lme4)
require(MuMIn)
  options(na.action = "na.fail") 
require(sp)
require(maps)
require(mapdata)
require(visreg)
require(ggplot2)
require(PBSmapping)
require(mgcv)
  
#########################################################
### Model Presence-Absence ### 
#########################################################
# Label each station as present or absent for Sablefish:
LL_prop$SBLpa = as.numeric(LL_prop$adjCPUE > 0)
length(LL_prop$SBLpa) # total number of hauls conducted
sum(LL_prop$SBLpa) # total number of hauls that captured SBL
  # There was only one station that did not sample Sablefish. Do not run separate presence-absence model.
LL_prop$Year = as.factor(LL_prop$Year)

#########################################################
### Model CPUE (where present) ### 
#########################################################
# Subset stations to include only those that sampled Sablefish:
SBL = subset(LL_prop, adjCPUE > 0)
SBL$logSBL = log(SBL$adjCPUE) # log-transform CPUE data

# Run the full model, without spatial autocorrelation:
SBL.cpue.gam_full = gam(logSBL ~ Year + s(Lon, Lat) + s(Depth), data = SBL, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
SBL.cpue.gam_select = dredge(SBL.cpue.gam_full, beta=FALSE, evaluate=T, rank="AIC", trace=FALSE)
  print(SBL.cpue.gam_select, abbrev.names=FALSE, warnings=T) # Table S2 
  summary(SBL.cpue.gam_full); sum(SBL.cpue.gam_full$edf) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
SBL.cpue.gamm_full = gamm(logSBL ~ Year + s(Lon, Lat) + s(Depth), correlation = corGaus(form = ~ Lon + Lat | Year), data = SBL, family = gaussian(link=identity))
  summary(SBL.cpue.gamm_full)
  summary(SBL.cpue.gamm_full$gam)
  summary(SBL.cpue.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
SBL.cpue.gamm_full_noSA = gamm(logSBL ~ Year + s(Lon, Lat) + s(Depth), data = SBL, family = gaussian(link=identity))
  AIC(SBL.cpue.gamm_full, SBL.cpue.gamm_full_noSA)
    # GAMM without spatial autocorrelation = best-fit model
      summary(SBL.cpue.gam_full$gam); sum(SBL.cpue.gam_full$gam$edf) # Table S4
      save(SBL.cpue.gam_full, file = "SBL.cpue.gam_full.rda")
      plot(SBL.cpue.gam_full$lme, xlab="", ylab="", col="red", pch=1) # Fig. S2

#########################################################
# Plot partial effects of each model covariate on CPUE of Sablefish, where present (Fig. S3):
lon_min = -172
lon_max = -130
lat_min = 52
lat_max = 62

data(worldHiresMapEnv)
SBL.cpue.gamm_full$gam$data = SBL

# Latitude x Longitude:
vis.gam(SBL.cpue.gamm_full$gam, c("Lon", "Lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Sablefish, Partial Effect on log-CPUE (kg per hectare)", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")

# Survey Year:
visreg(SBL.cpue.gamm_full$gam, xvar="Year", line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), band=T, partial=FALSE, rug=FALSE, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Sablefish") +
  theme(plot.title = element_text(size=18, hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Survey Year", y="Partial Effect on log-CPUE (kg per hectare)") +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(limits=c(0,10), breaks=c(0,2, 4, 6, 8, 10))

# Depth (m):
visreg(SBL.cpue.gamm_full$gam, xvar="Depth",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Sablefish") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on log-CPUE (kg per hectare)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000), breaks=c(0,150,300,450,600,750,900)) +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_y_continuous(limits=c(0,10), breaks=c(0,2, 4, 6, 8, 10)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter") 

#########################################################
### MODEL PREDICTIONS ###
#########################################################
# Create a uniform grid spanning the study area for model predictions:
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

# Read in and prepare shapefile for INPFC statistical areas (610 to 650): 
INPFC = readOGR(".", "GOA_Shapes")
INPFC_Pr = spTransform(INPFC, CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs")) 
    # UTM projection maintains consistent area within grid cells regardless of geographic location.

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

# Select range of coordinates for grid boundaries (UTM projection).
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

# Convert grid to a data frame (in decimal degrees) for plotting:
clip_INPFC_dd = clip_INPFC
clip_INPFC_dd = spTransform(clip_INPFC_dd, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
clip_INPFC_dd = fortify(clip_INPFC_dd)

#########################################################
# Prepare survey data for model predictions:

# Convert bottom trawl survey data to spatial data frame:
clipLL = clip_INPFC
clipLL = spTransform(clipLL, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
coordinates(LL_enviro) = c("Lon", "Lat")
proj4string(LL_enviro) = proj4string(clipLL)

# Find where survey data overlap with the clipped grid and convert back to a data frame:
tempdat = data.frame(myrows=names(over(LL_enviro, clipLL)), mygrid=over(LL_enviro, clipLL))
LL_enviro$id2 = over(LL_enviro, clipLL)
LL_enviro = as.data.frame(LL_enviro)

# Create data frame with id values and grid cell coordinates:
mycenter_INPFC = as.data.frame(gCentroid(clip_INPFC, byid=T)) %>% 
  mutate(id2=1:n(),
EEZgrid=rownames(.))

colnames(mycenter_INPFC)[3] = "id"
xy=mycenter_INPFC[,c("x", "y")]
mycenter_allp = SpatialPointsDataFrame(coords=xy, data=mycenter_INPFC, proj4string=CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs"))
mycenter_INPFC_Pr = spTransform(mycenter_allp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
mycenter_INPFC_df = as.data.frame(mycenter_INPFC_Pr) 
names(mycenter_INPFC_df) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")

# Join survey and grid cell data:
mycenter_all = mycenter_INPFC_df %>% left_join(LL_enviro)
mycenter_all = unique(mycenter_all)

#########################################################
# Compile input data for models:
geographData = na.omit(unique(mycenter_all[,c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")]))

# Remove few grid cells in inside waters of Southeast Alaska:
geographData = subset(geographData, id2 != 330 & id2 != 398 & id2 != 431 & id2 != 432)

# Add survey years:
geographData = geographData[rep(seq_len(nrow(geographData)), each = 13),]
geographData$YEAR = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017)
geographData$YEAR = as.factor(geographData$YEAR)

# Calculate year- and grid cell-specific mean temperatures:
HaulData = unique(mycenter_all[,c("id2","Year","Depth")])
LL_depth = HaulData %>%
  group_by(id2) %>%
  mutate(meanDepth = mean(Depth, na.omit=T))
LL_depth = as.data.frame(na.omit(unique(LL_depth[,c("id2", "meanDepth")])))

# Merge geographic and depth data:
EnviroData = unique(merge(LL_depth, geographData))

# Rename columns to match initial input data:
names(EnviroData) = c("id2", "Depth", "x.UTM", "y.UTM", "EEZgrid", "Lon", "Lat", "Year")

##########################################################
### Predict CPUE for Sablefish ###
  # Probability of occurrence is equal to one (not separately modeled).

### Predict CPUE for Sablefish ###
SBLcpue_predict = predict.gam(SBL.cpue.gamm_full$gam, newdata=EnviroData, type="response", se.fit=T)
SBLcpue_pred = cbind(EnviroData, SBLcpue_predict)

# Calculate 95% confidence intervals:
SBLcpue_predCI = within(SBLcpue_pred, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(SBLcpue_predCI) = c("id2", "Depth", "x.UTM", "y.UTM", "EEZgrid", "Lon", "Lat", "Year", "SBLcpue_fit", "SBLcpue_se.fit", "SBLcpue_upperCI", "SBLcpue_lowerCI")

# Normalize so that all grid cells in a given survey year sum to one:
normAbun_SBL = SBLcpue_predCI %>%
  group_by(Year) %>%
  mutate(relAbun = SBLcpue_fit/sum(SBLcpue_fit))

# Save relative predator densities:
saveRDS(normAbun_SBL, "Predation/TrophicStability/normSBLabun.rds")

#########################################################
### PLOT PLOT PREDATOR DENSITIES, Sablefish ###
#########################################################
# Load map data and set coordinate boundaries (DD):
data(nepacLLhigh)
lon_min = -172.5
lon_max = -129.5
lat_min = 49
lat_max = 62

# Clip maps to coordinate boundaries:
# World:
world = fortify(nepacLLhigh)
world2 = clipPolys(world, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))

# Canada:
Canada = raster::getData("GADM", country = "CAN", level = 0)
Canada = fortify(Canada)
names(Canada) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Canada$PID = as.numeric(Canada$PID)
Canada = clipPolys(Canada, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))

# INPFC Statistical Areas:
INPFC_shape = readOGR(".", "GOA_Shapes")
INPFC_shape = spTransform(INPFC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
INPFC_plot = fortify(INPFC_shape)
names(INPFC_plot) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
INPFC_plot$PID = as.numeric(INPFC_plot$PID)
INPFC_plot = clipPolys(INPFC_plot, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))

# Join summary information and spatial data: 
goa.df = fortify(clip_INPFC_dd, region='id')
colnames(goa.df)[colnames(goa.df) == "id"] = "EEZgrid"
plot_SBLpredictions = goa.df %>% right_join(SBLpredictions)
plot_SBLnormAbun = goa.df %>% right_join(normAbun_SBL)

# Plot average relative density in each grid cell (all years combined; Fig. S4):
textAbun = data.frame(text = c("Relative Density"))
rDensity_SBL = ggplot() +
  geom_polygon(data=plot_SBLnormAbun, aes(x=long, y=lat, group=group, fill=relAbun), col="black", lwd = 0.25) + 
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,0.0157), breaks=c(0.00,0.005,0.010,0.015)) +
  geom_polygon(data=INPFC_plot, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
  geom_text(data=textAbun, aes(label=text, x=-132.37, y=61.65, size=12), show.legend=FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.962, 0.837)) +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(family="Arial", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(8.0, "mm")) +
  theme(legend.spacing.x = unit(1.5, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=12)) +
  theme(axis.text.x = element_text(family="Arial", size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(hjust=0.46, size=12)) +
  labs(x="Longitude", y="Latitude") +
  scale_x_continuous(limits = c(lon_min, lon_max), expand = c(0,0)) +
  scale_y_continuous(limits = c(lat_min, lat_max), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

rDensity_SBL