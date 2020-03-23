# This script file includes code necessary to estimate relative densities of Pacific Cod in the Gulf of Alaska. Methods are described in Barnes et al. (2018) and were modified from Hunsicker et al. (2013) and Shelton et al. (2017). Relative predator densities (kg per ha) were multiplied by total predator biomass, mean annual rations, and age-specific proportions of pollock consumed to estimate year-, area-, and predator-specific predation.

# Standardized survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# Below, we use generalized additive mixed models (GAMMs) to quantify probabilities of occurrence and predator densities across a uniform grid spanning the study area. We included survey year, tow location (latitude and longitude), depth (m), and bottom temperature (degrees C) as model covariates and allowed a spatial autocorrelation term to vary by year. We multiplied year- and grid cell-specific probabilities of occurrence and predicted densities from best-fit models to estimate predator biomass. These values were then normalized to estimate relative predator density in each unique combination of survey year and grid cell.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Barnes, C. L., A. H. Beaudreau, M. E. Hunsicker, and L. Ciannelli (2018). Assessing the potential for competition between Pacific Halibut (Hippoglossus stenolepis) and Arrowtooth Flounder (Atheresthes stomias) in the Gulf of Alaska. PLoS ONE 13(12):e0209402.
# Barbeaux, S., K. Aydin, B. Fissel, K. Holsman, and W. Palsson. 2017. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 189–332.
# Hunsicker, M. E., L. Ciannelli, K. M. Bailey, S. Zador, and L. Stige. 2013. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLoS ONE 8(6):e66025.
# Shelton, A. O., M. E. Hunsicker, E. J. Ward, B. E. Feist, R. Blake, C. L. Ward, et al. 2017. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES Journal of Marince Science doi:10.1093/icesjms/fsx079.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
#########################################################
### INITIAL DATA PREPARATION ###
#########################################################
# Prepare and format AFSC bottom trawl survey data:
# These data include all survey tows conducted between 1984 and 2017.
trawl = read.csv("1_Data/AFSC_TrawlData_1984_2017.csv") 

# Manually assign statistical areas defined by the International North Pacific Fisheries Commission based on survey strata. Note: The second number from the right corresponds with individual statistical areas (e.g., STRATUM entry '120' = INPFC StatArea '620'; STRATUM entry '251' = StatArea '650').
strata = c(unique(trawl$STRATUM))
stat.area = c("640", "640", "640", "610", "650", "610", "650", "610", "610", "610", "610", "610", "610", "610", "610", "620", "620",
           "630", "630", "630", "640", "630", "630", "620", "620", "620", "630", "630", "620", "630", "630", "640", "640", "640",
           "650", "650", "650", "640", "630", "630", "620", "620", "630", "650", "650", "620", "640", "650", "650", "640", "640",
           "620", "630", "640", "630", "620", "630", "630", "610") 
trawl$StatArea = stat.area[match(trawl$STRATUM, strata)]

# Exclude data from 1984 and 1987 (survey methods were standardized in 1990):
trawl = subset(trawl, YEAR >= 1990)

# Treat survey year as a factor:
trawl$YEAR = as.factor(trawl$YEAR)

# Create a unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to proportional length data below):
trawl$Haul_Join = paste(trawl$VESSEL, trawl$CRUISE, trawl$HAUL, sep="")

# Relabel species codes: 10110 = Arrowtooth Flounder (ATF), 10120 = Pacific Halibut (PH), 20510 = Sablefish (SBL), 21720 = Pacific Cod (PC), 21740 = Walleye Pollock (WEP):
trawl$SPECIES_CODE = as.factor(trawl$SPECIES_CODE)
trawl$Species = trawl$SPECIES_CODE
levels(trawl$Species) = list(ATF="10110", PH="10120", SBL="20510", PC="21720", WEP="21740")

# Select the species of interest for this analysis (i.e., Pacific Cod):
trawl_PC = subset(trawl, Species=="PC")

# Reshape data frame (covert long to wide format):
require(reshape2)
trawl_wide = dcast(trawl_PC, YEAR + HAULJOIN + VESSEL + CRUISE + Haul_Join + STRATUM + DISTANCE_FISHED + NET_WIDTH + STATIONID + START_LATITUDE + START_LONGITUDE + GEAR_DEPTH + GEAR_TEMPERATURE + StatArea ~ Species, value.var="WGTCPUE")

#########################################################
### INCREASE COMPARABILITY WITH TOTAL BIOMASS ESTIMATES 
#########################################################
# Adjust haul-specific CPUE (kg per ha) to exclude fish not encompassed within assessment-based estimates of total biomass (PC < 0 yr or < 0 cm). Note: 100 to 200 fish were subsampled for length measurements per species and haul. 

# Read in and format length data from AFSC bottom trawl survey:
lengths = read.csv("1_Data/race_length_by_haul_PC.csv", header=T, skip=7)

# Remove and rename columns:
lengths = lengths[,c(1:10, 13:14, 25:27, 28:31)]
colnames(lengths) = c("Survey", "Year", "Cruise.Join", "Haul.Join", "Catch.Join", "Cruise.Number", "Vessel.Number", "Haul.Number", "START_LAT", "START_LON", "Stratum", "INPFCArea", "GEAR_DEPTH", "BOT_DEPTH", "SpeciesCode", "Common_Name", "Sci_Name", "Length..mm.", "Frequency")

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to survey data, prepared above):
lengths$Haul_Join = paste(lengths$Vessel.Number, lengths$Cruise.Number, lengths$Haul.Number, sep="")

# Exclude data prior to 1990, when survey methods were standardized:
lengths$Year = as.numeric(as.character(lengths$Year))
lengths = na.omit(subset(lengths, Year >= 1990))

# Convert fork length measurements (mm to cm):
lengths$Length..mm. = as.numeric(as.character(lengths$Length..mm.))
lengths$PredLength = lengths$Length..mm. / 10

# Create fork length bins (cm):
lengths$FLBin = cut(lengths$PredLength, breaks = c(0,10,20,30,40,50,60,70,80,90,100,110))
levels(lengths$FLBin) = c("<=10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", "101-110")

#######################################################
# Calculate proportions of fish (by weight) sampled by length bin and haul:

# Predict weight (kg) from length (cm) (parameter estimates from Barbeaux et al. 2017):
lengths$PredWT = 0.000005631 * (lengths$PredLength ^ 3.1306)

# Calculate total weight (kg) per haul:
lengths$Total_WT = as.numeric(as.character(lengths$Frequency)) * lengths$PredWT

require(dplyr)
# Calculate total weight (all size classes) of PC, by haul:
PC_all = lengths %>%
  group_by(Haul_Join) %>%
  summarise(Total_WT = sum(PredWT))

# Limit data to sizes of interest (PC >= 0 cm):
lengths_red = subset(lengths, PredLength >= 0)

# Calculate total weight (>= 0 cm) of PC, by haul:
PC_0 = lengths_red %>%
  group_by(Haul_Join) %>%
  summarise(Recruit_WT = sum(PredWT))

# Calculate proportion of haul measuring >= 0 cm:
PCprop = merge(PC_all, PC_0, all=T); PCprop[is.na(PCprop)] = 0
PCprop$PCprop = PCprop$Recruit_WT / PCprop$Total_WT

#######################################################
# Join bottom trawl survey and proportional length data:
trawl_prop = trawl_wide %>% left_join(PCprop)

# Assign mean proportion of >= 0 cm fish to hauls without length data:
trawl_prop$PCprop[is.na(trawl_prop$PCprop)] = mean(trawl_prop$PCprop, na.rm=T)

# Adjust CPUE (kg per ha) based on proportion of catch >= 0 cm:
trawl_prop$adjCPUE = trawl_prop$PC * trawl_prop$PCprop

# Ensure that the same tows are included in each model by removing those with incomplete environmental data (i.e., hauls with missing depths or bottom temperatures):
trawl_comp = subset(trawl_prop, !is.na(GEAR_DEPTH))
trawl_comp = subset(trawl_comp, !is.na(GEAR_TEMPERATURE))

################################################################
### MODEL FITTING AND PLOTTING ###
################################################################
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
### Model, Presence-Absence ### 
#########################################################
# Label each haul as present or absent for Pacific Cod:
trawl_comp$PCpa = as.numeric(trawl_comp$adjCPUE > 0)
length(trawl_comp$PCpa) # total number of hauls conducted
sum(trawl_comp$PCpa) # total number of hauls that captured cod

# Run the full model, without spatial autocorrelation:
PC.pa.gam_full = gam(PCpa ~ YEAR + s(START_LONGITUDE, START_LATITUDE) + s(GEAR_DEPTH) + s(GEAR_TEMPERATURE, k=4), data = trawl_comp, family = binomial(link=logit), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PC.pa.gam_select = dredge(PC.pa.gam_full, beta=FALSE, evaluate=T, rank="AIC", trace=FALSE)
  print(PC.pa.gam_select, abbrev.names=FALSE, warnings=T) # Table S2
    summary(PC.pa.gam_full); sum(PC.pa.gam_full$edf) # full model = best-fit model 

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
PC.pa.gamm_full = gamm(PCpa ~ YEAR + s(START_LONGITUDE, START_LATITUDE) + s(GEAR_DEPTH) + s(GEAR_TEMPERATURE, k=4), correlation = corGaus(form = ~ START_LONGITUDE + START_LATITUDE | YEAR), data = trawl_comp, family = binomial(link=logit))
  summary(PC.pa.gamm_full)
  summary(PC.pa.gamm_full$gam)
  summary(PC.pa.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
PC.pa.gamm_full_noSA = gamm(PCpa ~ YEAR + s(START_LONGITUDE, START_LATITUDE) + s(GEAR_DEPTH) + s(GEAR_TEMPERATURE, k=4), data = trawl_comp, family = binomial(link=logit))
  AIC(PC.pa.gamm_full, PC.pa.gamm_full_noSA) # Table S3
    # GAMM with spatial autocorrelation = best-fit model
      summary(PC.pa.gamm_full$gam); sum(PC.pa.gamm_full$gam$edf) # Table S4
      save(PC.pa.gamm_full, file = "PC.pa.gamm_full.rda")
      plot(PC.pa.gamm_full$lme, xlab="", ylab="", col="red", pch=1) # Fig. S2
    
#################################################################
# Plot partial effects of each model covariate on the probability of occurrence for Pacific Cod (Fig. S3):
lon_min = -172
lon_max = -130
lat_min = 52
lat_max = 62

data(worldHiresMapEnv)
PC.pa.gamm_full$gam$data = trawl_comp

# Latitude x Longitude:
vis.gam(PC.pa.gamm_full$gam, view=c("START_LONGITUDE", "START_LATITUDE"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Pacific Cod, Partial Effect on Presence or Absence", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")

# Survey Year:
visreg(fit=PC.pa.gamm_full$gam, xvar="YEAR", band=T, partial=FALSE, rug=FALSE, line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), trans=binomial()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Cod") +
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
visreg(PC.pa.gamm_full$gam, xvar="GEAR_DEPTH",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Cod") +
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

# Bottom Temperature (degrees C):
visreg(PC.pa.gamm_full$gam, xvar="GEAR_TEMPERATURE",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Cod") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x=expression(paste("Bottom Temperature (",degree,"C)")), y="Partial Effect on Presence (1) or Absence (0)") +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25, 0.50, 0.75, 1)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")

#########################################################
### Model CPUE (where present) ### 
#########################################################
# Subset hauls to include only those that sampled Pacific Cod:
PC = subset(trawl_comp, adjCPUE > 0)
PC$logPC = log(PC$adjCPUE)  # log-transform CPUE data

# Run the full model, without spatial autocorrelation:
PC.cpue.gam_full = gam(logPC ~ YEAR + s(START_LONGITUDE, START_LATITUDE) + s(GEAR_DEPTH) + s(GEAR_TEMPERATURE, k=4), data = PC, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PC.cpue.gam_select = dredge(PC.cpue.gam_full, beta=FALSE, evaluate=T, rank="AIC", trace=FALSE)
  print(PC.cpue.gam_select, abbrev.names=FALSE, warnings=T) # Table S2 
  summary(PC.cpue.gam_full); sum(PC.cpue.gam_full$edf) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
PC.cpue.gamm_full = gamm(logPC ~ YEAR + s(START_LONGITUDE, START_LATITUDE) + s(GEAR_DEPTH) + s(GEAR_TEMPERATURE, k=4), correlation = corGaus(form = ~ START_LONGITUDE + START_LATITUDE | YEAR), data = PC, family = gaussian(link=identity))
  summary(PC.cpue.gamm_full)
  summary(PC.cpue.gamm_full$gam)
  summary(PC.cpue.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
PC.cpue.gamm_full_noSA = gamm(logPC ~ YEAR + s(START_LONGITUDE, START_LATITUDE) + s(GEAR_DEPTH) + s(GEAR_TEMPERATURE, k=4), data = PC, family = gaussian(link=identity))
  AIC(PC.cpue.gamm_full, PC.cpue.gamm_full_noSA)   
    # GAMM with spatial autocorrelation = best-fit model
    summary(PC.cpue.gamm_full$gam); sum(PC.cpue.gamm_full$gam$edf) # Table S4
    save(PC.cpue.gamm_full, file = "PC.cpue.gamm_full.rda")
    plot(PC.cpue.gamm_full$lme, xlab="", ylab="", col="red", pch=1) # Fig. S2
  
#################################################################
# Plot partial effects of each model covariate on CPUE of Pacific Cod, where present (Fig. S3):
PC.cpue.gamm_full$gam$data = PC

# Latitude x Longitude:
vis.gam(PC.cpue.gamm_full$gam, c("START_LONGITUDE", "START_LATITUDE"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Pacific Cod, Partial Effect on log-CPUE (kg per hectare)", too.far=0.025, n.grid=500, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max))
maps::map('worldHires', fill=T, xlim=c(lon_min, lon_max), ylim=c(lat_min, lat_max), add=T, col="lightgrey")

# Survey Year:
visreg(PC.cpue.gamm_full$gam, xvar="YEAR", line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), band=T, partial=FALSE, rug=FALSE, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Cod") +
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
visreg(PC.cpue.gamm_full$gam, xvar="GEAR_DEPTH",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Cod") +
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

# Bottom Temperature (degrees C):
visreg(PC.cpue.gamm_full$gam, "GEAR_TEMPERATURE",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=T, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=T) +
  theme_bw() +
  ggtitle("Pacific Cod") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x=expression(paste("Bottom Temperature (",degree,"C)")), y="Partial Effect on log-CPUE (kg per hectare)") +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +
  scale_y_continuous(limits=c(0,10), breaks=c(0,2, 4, 6, 8, 10)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")

################################################################
### MODEL PREDICTIONS ###
################################################################
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

########################################################## 
# Prepare survey data for model predictions:

# Convert bottom trawl survey data to a spatial data frame:
trawl_enviro = trawl
clipTrawl = clip_INPFC
clipTrawl = spTransform(clipTrawl, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
coordinates(trawl_enviro) = c("START_LONGITUDE", "START_LATITUDE")
proj4string(trawl_enviro) = proj4string(clipTrawl)

# Find where survey data overlap with the clipped grid and convert back to a data frame:
tempdat = data.frame(myrows=names(over(trawl_enviro, clipTrawl)), mygrid = over(trawl_enviro, clipTrawl))
trawl_enviro$id2 = over(trawl_enviro, clipTrawl)
trawl_enviro = as.data.frame(trawl_enviro)

# Create a data frame with grid cell id values and coordinates:
mycenter_INPFC = as.data.frame(gCentroid(clip_INPFC, byid = T)) %>%   mutate(id2 = 1:n(),
EEZgrid = rownames(.))

colnames(mycenter_INPFC)[3] = "id"
xy = mycenter_INPFC[,c("x", "y")]
mycenter_INPFC_sp = SpatialPointsDataFrame(coords = xy, data=mycenter_INPFC, proj4string = CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs"))
mycenter_INPFC_Pr = spTransform(mycenter_INPFC_sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
mycenter_INPFC_df = as.data.frame(mycenter_INPFC_Pr) 
names(mycenter_INPFC_df) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")

# Join survey and grid cell data:
mycenter_all = mycenter_INPFC_df %>% right_join(trawl_enviro)
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
HaulData = unique(mycenter_all[,c("id2","YEAR","GEAR_DEPTH","GEAR_TEMPERATURE", "StatArea", "Haul_Join")])
trawl_BT = HaulData %>%
  group_by(YEAR, id2) %>%
  mutate(meanTemp = mean(GEAR_TEMPERATURE, na.rm=T)) 

# Join geographic and temperature data:
BTdata = geographData %>% full_join(trawl_BT)

# Calculate grid cell-specific (all years combined) mean depths:
trawl_depth = HaulData  %>%
  group_by(id2) %>%
  mutate(meanDepth = mean(GEAR_DEPTH, na.rm=T))
trawl_depth = as.data.frame(na.omit(unique(trawl_depth[,c("id2", "meanDepth")])))

# Merge geographic and depth data:
depthData = geographData %>% right_join(trawl_depth)
depthData = na.omit(depthData)

# Merge temperature and depth data:
EnviroDataPre = BTdata %>% right_join(depthData)

# Manually assign INPFC statistical areas:
EnviroDataPre$StatArea = with(EnviroDataPre, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), 610, 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), 620, 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), 630, 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), 640, 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), 650, 0))))))
EnviroDataPre = subset(EnviroDataPre, StatArea != 0)

# Calculate year- and grid-specific mean bottom temperatures:
BT = HaulData %>%
  group_by(id2) %>% 
  summarise(meanBT = mean(GEAR_TEMPERATURE, na.rm=T))

# Replace missing temperature data with mean values:
EnviroData = BT %>% right_join(EnviroDataPre)
EnviroData$inputBT = EnviroData$meanTemp
EnviroData$inputBT[is.na(EnviroData$inputBT)] = EnviroData$meanBT[is.na(EnviroData$inputBT)]
EnviroData = unique(EnviroData)

# Select and rename columns to match initial input data:
EnviroData = EnviroData[,c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat", "YEAR", "meanDepth", "inputBT")]
colnames(EnviroData) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "START_LONGITUDE", "START_LATITUDE", "YEAR", "GEAR_DEPTH", "GEAR_TEMPERATURE")

###############################################################
### Predict presence-absence for Pacific Cod ###
PCpa_predict = predict.gam(PC.pa.gamm_full$gam, newdata=EnviroData, type="response", se.fit=T)
PCpa_pred = cbind(EnviroData, PCpa_predict)

# Calculate 95% confidence intervals:
PCpa_predCI = within(PCpa_pred, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(PCpa_predCI) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "START_LONGITUDE", "START_LATITUDE", "YEAR", "GEAR_DEPTH", "GEAR_TEMPERATURE", "PCpa_fit", "PCpa_se.fit", "PCpa_upperCI", "PCpa_lowerCI")

### Predict CPUE for Pacific Cod ###
PCcpue_predict = predict.gam(PC.cpue.gamm_full$gam, newdata=EnviroData, type="response", se.fit=T)
PCcpue_pred = cbind(EnviroData, PCcpue_predict)

# Calculate 95% confidence intervals:
PCcpue_predCI = within(PCcpue_pred, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(PCcpue_predCI) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "START_LONGITUDE", "START_LATITUDE", "YEAR", "GEAR_DEPTH", "GEAR_TEMPERATURE", "PCcpue_fit", "PCcpue_se.fit", "PCcpue_upperCI", "PCcpue_lowerCI")

# Merge predictions:
PCpredictions = merge(PCpa_predCI, PCcpue_predCI)

# Calculate predator density, accounting for presence-absence:
PCpredictions$PCpredAbun = PCpredictions$PCpa_fit * PCpredictions$PCcpue_fit
PCpredictions = unique(PCpredictions)

# Reassign negative values as zero:
PCpredictions$PCpredAbun = ifelse(PCpredictions$PCpredAbun < 0, 0, PCpredictions$PCpredAbun)

# Normalize so that all grid cells in a given survey year sum to one:
normAbun_PC = PCpredictions %>%
  group_by(YEAR) %>%
  mutate(PCnormAbun = PCpredAbun/sum(PCpredAbun))

# Save relative predator densities:
saveRDS(normAbun_PC, "Predation/TrophicStability/normPCabun.rds")

#########################################################
### PLOT PREDATOR DENSITIES, Pacific Cod ###
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
plot_PCpredictions = goa.df %>% right_join(PCpredictions)
plot_PCnormAbun = goa.df %>% right_join(normAbun_PC)

# Plot average relative density in each grid cell (all years combined; Fig. S4):
textAbun = data.frame(text = c("Relative Density"))
rDensity_PC = ggplot() +
  geom_polygon(data=plot_PCnormAbun, aes(x=long, y=lat, group=group, fill=PCnormAbun), col="black", lwd = 0.25) + 
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,0.01), breaks=c(0.00,0.0025,0.005,0.0075,0.010)) +
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

rDensity_PC