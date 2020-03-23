# This script file includes code necessary to estimate relative densities of Pacific Halibut in the Gulf of Alaska. Methods are described in Barnes et al. (2018) and were modified from Hunsicker et al. (2013) and Shelton et al. (2017). Relative predator densities (kg per ha) were multiplied by total predator biomass, mean annual rations, and age-specific proportions of pollock consumed to estimate year-, area-, and predator-specific predation.

# Standardized survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# Below, we use generalized additive mixed models (GAMMs) to quantify probabilities of occurrence and predator densities across a uniform grid spanning the study area. We included survey year, tow location (latitude and longitude), depth (m), and bottom temperature (degrees C) as model covariates and allowed a spatial autocorrelation term to vary by year. We multiplied year- and grid cell-specific probabilities of occurrence and predicted densities from best-fit models to estimate predator biomass. These values were then normalized to estimate relative predator density in each unique combination of survey year and grid cell.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Barnes, C. L., A. H. Beaudreau, M. E. Hunsicker, and L. Ciannelli (2018). Assessing the potential for competition between Pacific Halibut (Hippoglossus stenolepis) and Arrowtooth Flounder (Atheresthes stomias) in the Gulf of Alaska. PLoS ONE 13(12):e0209402.
# Clark, W. G., and S. R. Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. IPHC Scientific Report 83. 
# Hunsicker, M. E., L. Ciannelli, K. M. Bailey, S. Zador, and L. Stige. 2013. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLoS ONE 8(6):e66025.
# Shelton, A. O., M. E. Hunsicker, E. J. Ward, B. E. Feist, R. Blake, C. L. Ward, et al. 2017. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES Journal of Marince Science doi:10.1093/icesjms/fsx079.
# Stewart, I., and A. Hicks. 2017. Assessment of the Pacific halibut (Hippoglossus stenolepis) stock at the end of 2017. International Pacific Halibut Commission IPHC-2018-AM094-10.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Technical Memorandum NMFS-AFSC-325. 

#########################################################
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##########################################################
### INITIAL DATA PREPARATION ###
##########################################################
# Prepare and format IPHC setline survey data:
# These data include all sets conducted between 1998 and 2017.
setline = read.csv("1_Data/IPHC_SetLineData.csv", header=T) 
setline$StatArea = as.numeric(as.character(setline$Statarea))

# Group current statistical areas (i.e., those defined in 2003) by IPHC regulatory area:
setline$RegArea = with(setline, 
ifelse(StatArea >= 006 & StatArea <= 050, "2A",
ifelse(StatArea >= 060 & StatArea <= 135, "2B",
ifelse(StatArea >= 140 & StatArea <= 184, "2C", 
ifelse(StatArea >= 185 & StatArea <= 281, "3A",
ifelse(StatArea >= 290 & StatArea <= 340, "3B",
ifelse(StatArea >= 350 & StatArea <= 395, "4A",
ifelse(StatArea >= 400 & StatArea <= 510, "4B", "Other"))))))))
setline$RegArea[is.na(setline$RegArea)] = "EBS"

# Select only the regulatory areas encompassed by the Gulf of Alaska:
setline = subset(setline, RegArea != "EBS" & RegArea != "2A" & RegArea != "2B" & RegArea != "4B" & RegArea != "Other")

# Remove irregular year designations:
setline = subset(setline, Year < 3000)

# Treat survey year as a factor:
setline$Year = as.factor(setline$Year)

# Convert from lb to kg and remove commas from data frame:
setline$O32_weight = as.numeric(gsub(",","", setline$O32_weight)) 
setline$O32_kg = setline$O32_weight * 0.453592

# Calculate mean depths (fm) for each set and convert to km:
setline$Depth_fm = (setline$Max_depth_fm + setline$Min_depth_fm) / 2
setline$Depth_m = setline$Depth_fm * 1.8288

# Ensure that the same tows are included in each model by removing those with incomplete environmental data (i.e., hauls with missing depths):
setline = subset(setline, !is.na(Depth_m))

##########################################################
### MODEL FITTING AND PLOTTING ###
##########################################################
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
# Label each set as being present or absent for Pacific Halibut:
setline$PHpa = as.numeric(setline$O32_kg > 0)
length(setline$PHpa) # total number of skates fished
sum(setline$PHpa) # total number of skates that captured halibut

# Convert from latitude and longitude from decimal minutes to decimal degrees:
require(measurements)
setline$Lat = as.numeric(as.character(measurements::conv_unit(setline$Lat_fished, from = 'deg_dec_min', to = 'dec_deg')))
setline$Lon = - as.numeric(as.character(measurements::conv_unit(setline$Lon_fished, from = 'deg_dec_min', to = 'dec_deg')))

require(dplyr)
# Remove duplicate locations:
setline$LocCheck = setline$Lat * setline$Lon
SLsub = setline %>%
  group_by(LocCheck) %>%
  mutate(meanO32 = mean(O32_kg))
SLsub$meanO32 = round(SLsub$meanO32, digits=1)
SLsub = SLsub[,c("Year", "Lat", "Lon", "LocCheck", "Depth_m", "PHpa", "meanO32")]
SLsub = na.omit(unique(SLsub))

# Run the full model, without spatial autocorrelation:
PH.pa.gam_full = gam(PHpa ~ Year + s(Lon, Lat) + s(Depth_m), data = SLsub, family = binomial(link=logit), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PH.pa.gam_select = dredge(PH.pa.gam_full, beta=FALSE, evaluate=T, rank="AIC", trace=FALSE)
  print(PH.pa.gam_select, abbrev.names=FALSE, warnings=T) # Table S2 
  summary(PH.pa.gam_full) # full model = best-fit model
    
# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
  # GAMM for Pacific Halibut does not converge. Use GAM for probabilities of occurrence.
    save(PH.pa.gam_full, file = "2_RelativePredatorDensities/PH.pa.gam_full.rda")
    summary(PH.pa.gam_full); sum(PH.pa.gam_full$edf) # Table S4
    plot.lme(PH.pa.gam_full, xlab="", ylab="", col="red", pch=1) # Fig. S2

##########################################################
# Plot partial effects of each model covariate on the probability of occurrence for Pacific Halibut (Fig. S3):
lon_min = -172
lon_max = -130
lat_min = 52
lat_max = 62

data(worldHiresMapEnv)
PH.pa.gam_full$data = SLsub

# Latitude x Longitude:
vis.gam(PH.pa.gam_full, view=c("Lon", "Lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Pacific Halibut, Partial Effect on Presence or Absence", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")

# Survey Year:
visreg(fit=PH.pa.gam_full, xvar="Year", band=T, partial=FALSE, rug=FALSE, line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), trans=binomial()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(size=18, hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Survey Year", y="Partial effect on presence-absence") +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25, 0.50, 0.75, 1)) +
  theme(legend.background = element_rect(fill="transparent"))

# Depth (m):
visreg(PH.pa.gam_full, xvar="Depth_m",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on Presence (1) or Absence (0)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000), breaks=c(0,150,300,450,600,750,900)) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25, 0.50, 0.75, 1)) +
  theme(legend.background = element_rect(fill="transparent")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")

# Limited bottom temperature data (temperature covariate not included in models).

#########################################################
### Model CPUE (where present) ### 
#########################################################
# Subset skates to include only those that caught Pacific Halibut:
PH = subset(SLsub, meanO32 > 0)
PH$logPH = log(PH$meanO32) # log-transform CPUE data

# Run the full model, without spatial autocorrelation:
PH.cpue.gam_full = gam(logPH ~ Year + s(Lon, Lat) + s(Depth_m), data = PH, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PH.cpue.gam_select = dredge(PH.cpue.gam_full, beta=FALSE, evaluate=T, rank="AIC", trace=FALSE)
  print(PH.cpue.gam_select, abbrev.names=FALSE, warnings=T)
  summary(PH.cpue.gam_full); sum(PH.cpue.gam_full$edf) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
PH.cpue.gamm_full = gamm(logPH ~ Year + s(Lon, Lat) + s(Depth_m), correlation = corGaus(form = ~ Lon + Lat | Year), data = PH, family = gaussian(link=identity))
  summary(PH.cpue.gamm_full)
  summary(PH.cpue.gamm_full$gam)
  summary(PH.cpue.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
PH.cpue.gamm_full_noSA = gamm(logPH ~ Year + s(Lon, Lat) + s(Depth_m), data = PH, family = gaussian(link=identity))
  AIC(PH.cpue.gamm_full, PH.cpue.gamm_full_noSA)
    # GAMM with spatial autocorrelation = best-fit model
      summary(PH.cpue.gamm_full$gam); sum(PH.cpue.gamm_full_noSA$gam$edf) # Table S4
      save(PH.cpue.gamm_full, file = "PH.cpue.gamm_full.rda")
      plot(PH.cpue.gamm_full$lme, xlab="", ylab="", col="red", pch=1) # Fig. S2
      
##########################################################
# Plot partial effects of each model covariate on CPUE of Pacific Halibut, where present (Fig. S3):
PH.cpue.gamm_full$gam$data = PH

# Latitude x Longitude:
vis.gam(PH.cpue.gamm_full$gam, c("Lon", "Lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Pacific Halibut, Partial Effect on log-CPUE (kg per hectare)", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")

# Survey Year:
visreg(PH.cpue.gamm_full$gam, xvar="Year", line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), band=T, partial=FALSE, rug=FALSE, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
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
visreg(PH.cpue.gamm_full$gam, xvar="Depth_m",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
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

# Limited bottom temperature data (temperature covariate not included in models).

##########################################################
### MODEL PREDICTIONS ###
##########################################################
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

# Select range of coordinates for grid boundaries (UTM projection):
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

########################################################### 
# Prepare survey data for model predictions:

# Convert bottom trawl survey data to spatial data frame:
SL_enviro = SLsub
clipSL = clip_INPFC
clipSL = spTransform(clipSL, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
coordinates(SL_enviro) = c("Lon", "Lat")
proj4string(SL_enviro) = proj4string(clipSL)

# Find where survey data overlap with the clipped grid and convert back to a data frame:
tempdat = data.frame(myrows=names(over(SL_enviro, clipSL)), mygrid=over(SL_enviro, clipSL))
SL_enviro$id2 = over(SL_enviro, clipSL)
SL_enviro = as.data.frame(SL_enviro)

# Create data frame with id values and grid cell coordinates:
mycenter_INPFC = as.data.frame(gCentroid(clip_INPFC, byid=T)) %>% 
  mutate(id2=1:n(),
EEZgrid=rownames(.))

colnames(mycenter_INPFC)[3] = "id"
xy=mycenter_INPFC[,c("x", "y")]
mycenter_INPFC_sp = SpatialPointsDataFrame(coords=xy, data=mycenter_INPFC, proj4string=CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs"))
mycenter_INPFC_Pr = spTransform(mycenter_INPFC_sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
mycenter_INPFC_df = as.data.frame(mycenter_INPFC_Pr) 
names(mycenter_INPFC_df) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")

# Join survey and grid cell data:
mycenter_all = mycenter_INPFC_df %>% right_join(SL_enviro)
mycenter_all = unique(mycenter_all)

#########################################################
# Compile input data for models:
geographData = na.omit(unique(mycenter_all[,c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")]))

# Remove few grid cells in inside waters of Southeast Alaska:
geographData = subset(geographData, id2 != 330 & id2 != 398 & id2 != 431 & id2 != 432)

# Add survey years:
geographData = geographData[rep(seq_len(nrow(geographData)), each = 10),]
geographData$YEAR = c(1999,2001,2003,2005,2007,2009,2011,2013,2015,2017)
geographData$YEAR = as.factor(geographData$YEAR)

# Calculate year- and grid cell-specific mean temperatures:
HaulData = unique(mycenter_all[,c("id2","Year","Depth_m")])
SL_depth = HaulData %>%
  group_by(id2) %>%
  summarize(meanDepth = mean(Depth_m, na.rm=T))
SL_depth = as.data.frame(na.omit(unique(SL_depth[,c("id2", "meanDepth")])))

# Merge geographic and depth data:
EnviroData = unique(merge(SL_depth, geographData))

# Rename columns to match initial input data:
names(EnviroData) = c("id2", "Depth_m", "x.UTM", "y.UTM", "EEZgrid", "Lon", "Lat", "Year")

##########################################################
### Predict presence-absence for Pacific Halibut ###
PHpa_predict = predict.gam(PH.pa.gam_full, newdata=EnviroData, type="response", se.fit=T)
PHpa_pred = cbind(EnviroData, PHpa_predict)

# Calculate 95% confidence intervals:
PHpa_predCI = within(PHpa_pred, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(PHpa_predCI) = c("id2", "Depth_m", "x.UTM", "y.UTM", "EEZgrid", "Lon", "Lat", "Year", "PHpa_fit", "PHpa_se.fit", "PHpa_upperCI", "PHpa_lowerCI")

### Predict CPUE for Pacific Halibut ###
PHcpue_predict = predict.gam(PH.cpue.gamm_full$gam, newdata=EnviroData, type="response", se.fit=T)
PHcpue_pred = cbind(EnviroData, PHcpue_predict)

# Calculate 95% confidence intervals:
PHcpue_predCI = within(PHcpue_pred, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(PHcpue_predCI) = c("id2", "Depth_m", "x.UTM", "y.UTM", "EEZgrid", "Lon", "Lat", "Year", "PHcpue_fit", "PHcpue_se.fit", "PHcpue_upperCI", "PHcpue_lowerCI")

# Merge predictions:
PHpredictions = PHpa_predCI %>% full_join(PHcpue_predCI)

# Calculate predator density, accounting for presence-absence:
PHpredictions$PHpredAbun = PHpredictions$PHpa_fit * PHpredictions$PHcpue_fit

# Normalize so that all grid cells in a given survey year sum to one:
normAbun_PH = PHpredictions %>%
  group_by(Year) %>%
  mutate(relAbun = PHpredAbun/sum(PHpredAbun))

# Save relative predator densities:
saveRDS(normAbun_PH, "2_RelativePredatorDensities/normPHabun.rds")

##########################################################
### PLOT PREDATOR DENSITIES, Pacific Halibut ###
##########################################################
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
plot_PHpredictions = goa.df %>% right_join(PHpredictions)
plot_PHnormAbun = goa.df %>% right_join(normAbun_PH)

# Plot average relative density in each grid cell (all years combined; Fig. S4):
textAbun = data.frame(text = c("Relative Density"))
rDensity_PH = ggplot() +
  geom_polygon(data=plot_PHnormAbun, aes(x=long, y=lat, group=group, fill=relAbun), col="black", lwd = 0.25) + 
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,0.007), breaks=c(0.00,0.002,0.004,0.006)) +
  geom_polygon(data=INPFC_plot, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
  geom_text(data=textAbun, aes(label=text, x=-132.36, y=61.65, size=12), show.legend=FALSE) +
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

rDensity_PH