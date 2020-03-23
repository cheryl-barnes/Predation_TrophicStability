# This script file includes code necessary to calculate fork length- and biomass-weighted proportions of pollock (Chipps and Garvey 2007) consumed by Arrowtooth Flounder. The analyses below rely on standardized bottom trawl surveys (Gulf of Alaska, 1990 to 2015) conducted by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (AFSC, National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]). Food habits data (Gulf of Alaska, 1990 to 2015) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php. See von Szalay and Raring (2016) and Livingston et al. (2017) for data collection methods.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Chipps, S. R., and J. E. Garvey. 2007. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. C. S. Guy and M. L. Brown, editors. Bethesda, MD. American Fisheries Society 473–514.
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environmental Biology of Fishes. 100(4):443–470.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##########################################################
### INITIAL DATA PREPARATION ###
##########################################################
# Read in and format food habits data from subsampled stomachs:
preyData = read.csv("1_Data/GOA_RawPP_ATFdiets.csv") 

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to AFSC bottom trawl survey data below):
preyData$Haul_Join = paste(preyData$Vessel, preyData$Cruise, preyData$Haul, sep="")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
preyData = subset(preyData, Year >= 1990)

# Remove statistical areas other than 610, 620, 630, 640, and 650:
preyData = subset(preyData, INPFC_Area != 649 & INPFC_Area != 659)

require(reshape2)
require(dplyr)
# Remove all empty stomachs:
preyWT = dcast(preyData, Year + INPFC_Area + Haul_Join + Rlat + Rlong + Gear_depth + Gear_temp + Pred_name + Pred_len ~ Prey_Name, value.var = "Prey_twt", sum) # convert data, long to wide form

preyWT = subset(preyWT,select = -c(Empty)) # exclude 'Empty' prey ID
preyWT = preyWT %>% 
    mutate(sumWT = rowSums(select(., Agonidae:Zoarcoidae), na.rm = TRUE))
preyWTcontents = subset(preyWT, sumWT > 0) # remove hidden empties

# Create 10-cm predator length bins:
preyWTcontents$FLBin = cut(preyWTcontents$Pred_len, breaks = c(0,18,28,38,48,58,68,78,88))
levels(preyWTcontents$FLBin) = c("<=18", "19-28", "29-38", "39-48", "49-58", "59-68", "69-78", "79-88")
preyWTcontents = rename(preyWTcontents, PredLength = Pred_len)

##########################################################
# Read in and format predator length data from AFSC bottom trawls:
lengths = read.csv("1_Data/race_length_by_haul_ATF.csv", header=T, skip=7)

# Select and rename columns to match survey data:
lengths = subset(lengths, select = c(Survey:Starting.Longitude..dd., Stratum:Stratum.INPFC.Area,Gear.Depth:Species.Code, Common.Name:Frequency))
colnames(lengths) = c("Survey", "Year", "Cruise.Join", "Haul.Join", "Catch.Join", "Cruise.Number", "Vessel.Number", "Haul.Number", "START_LAT", "START_LON", "Stratum", "INPFCArea", "GEAR_DEPTH", "BOT_DEPTH", "SpeciesCode", "Common_Name", "Sci_Name", "Length..mm.", "Frequency")

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to food habits data prepared above):
lengths$Haul_Join = paste(lengths$Vessel.Number, lengths$Cruise.Number, lengths$Haul.Number, sep="")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
lengths$Year = as.numeric(as.character(lengths$Year))
lengths = na.omit(subset(lengths, Year >= 1990))

# Convert predator lengths from mm to cm:
lengths$PredLength = as.numeric(as.character(lengths$Length..mm.))/10

# Create 10-cm predator length bins:
lengths$FLBin = cut(lengths$PredLength, breaks = c(0,18,28,38,48,58,68,78,88,100))
levels(lengths$FLBin) = c("<=18", "19-28", "29-38", "39-48", "49-58", "59-68", "69-78", "79-88", ">88")

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

# Read in and prepare shape file for INPFC statistical areas (610 to 650):
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/1_Data/")
INPFC = readOGR(".", "GOA_Shapes")
INPFC_Pr = spTransform(INPFC, CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs")) 
  # UTM projection maintains consistent grid cell area regardless of geographic location.

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
lonmin = -875
lonmax = 1975
latmin = 5075
latmax = 7975

# Compile series of points for uniform grid:
mygrd = expand.grid(
  LON = seq(lonmin, lonmax, by=my.interval),
  LAT = seq(latmin, latmax, by=my.interval)) %>% 
  mutate(my.z=1:n()) %>% 
  data.frame

# Convert mygrd to a spatial dataframe:
coordinates(mygrd) = ~ LON + LAT

# Convert SpatialPoints to a SpatialPixelsDataFrame:
mygrd = (as(SpatialPixelsDataFrame(mygrd, mygrd@data, tolerance=.00086), "SpatialPolygonsDataFrame"))

# Project, clip, and reproject mygrid to INPFC statistical areas:
proj4string(mygrd) = proj4string(INPFC_Pr)
INPFC_Pr = gBuffer(INPFC_Pr, byid=TRUE, width=0)
mygrd = gBuffer(mygrd, byid=TRUE, width=0)
clip_INPFC = gIntersection(INPFC_Pr, mygrd, byid = TRUE, drop_lower_td = TRUE)
proj4string(clip_INPFC) = proj4string(INPFC_Pr)

# Convert grid to a data frame (in decimal degrees):
clip_INPFC_dd = clip_INPFC
clip_INPFC_dd = spTransform(clip_INPFC_dd, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Find where food habits data overlap with clipped grid and convert back to a data frame:
coordinates(preyWTcontents) = ~ Rlong + Rlat 
proj4string(preyWTcontents) = proj4string(clip_INPFC_dd)

tempdat = data.frame(myrows=names(over(preyWTcontents, clip_INPFC_dd)), mygrid=over(preyWTcontents, clip_INPFC_dd))
preyWTcontents$id2 = over(preyWTcontents, clip_INPFC_dd)
preyWTcontents = as.data.frame(preyWTcontents)

# Create data frame with grid cell identifiers (id2) and coordinates (vertices and centers):
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
preyWTgrid = mycenter_INPFC_df %>% right_join(preyWTcontents)

# Find where predator length data overlap with clipped grid and convert back to a data frame:
lengths$START_LON = as.numeric(as.character(lengths$START_LON))
lengths$START_LAT = as.numeric(as.character(lengths$START_LAT))
lengths = na.omit(lengths)

coordinates(lengths) = ~ START_LON + START_LAT
proj4string(lengths) = proj4string(clip_INPFC_dd)

tempdat = data.frame(myrows=names(over(lengths, clip_INPFC_dd)), mygrid=over(lengths, clip_INPFC_dd))
lengths$id2 = over(lengths, clip_INPFC_dd)
lengths = as.data.frame(lengths)

# Join predator length and geographic data:
lengthsGrid = mycenter_INPFC_df %>% right_join(lengths)

##########################################################
### SUMMARIZE FOOD HABITS DATA ###
##########################################################
# Reshape data from wide to long and select columns of interest:
colnames(preyWTgrid[16]); colnames(preyWTgrid[87]); ncol(preyWTgrid)
preyWTlong = melt(preyWTgrid, id.vars = c("id2", "EEZgrid", "Year", "INPFC_Area", "Haul_Join", "Rlong", "Rlat", "Pred_name", "PredLength", "FLBin"), measure.vars = 16:87, variable.name = "Prey_Name", value.name = "preyWT_g")
preyWTlong = na.omit(preyWTlong)
preyWTlong$Pred_name = "Arrowtooth Flounder"
preyWTlong$INPFC_Area = as.factor(preyWTlong$INPFC_Area)
levels(preyWTlong$INPFC_Area) = c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern")

# Recategorize prey taxa:
levels(preyWTlong$Prey_Name) = list(Gadus.chalcogrammus="Walleye.pollock", Misc.Crabs="Tanner.Crab", Misc.Shrimps="Pandalidae..shrimp.", Misc.Crabs="Paguridae", Teleostei="Misc.Teleost", Trachiniformes="Ammodytidae", Scorpaeniformes="Atka.Mackeral", Euphausiacea="Euphausiacea", Annelida="Polychaeta", Offal="Offal", Misc.Shrimps="Crangonidae..shrimp.", Misc.Crabs="Misc.Anomura", Osmeriformes="Osmerid", Pleuronectiformes="Misc.Flatfish", Misc.Crabs="Misc.Majidae", Octopoda="Octopoda", Scorpaeniformes="Cottid", Misc.Shrimps="Misc.Shrimp", Pleuronectiformes="Arrowtooth.flounder", Perciformes="Stichaeidae", Misc.Crabs="Cancridea", Misc.Shrimps="Hippolytidae..shrimp.", Crustacea="Gammaridea", Scorpaeniformes="Sebastes", Gadiformes="Misc.Gadidae", Pleuronectiformes="Flathead.sole", Annelida="Misc.Worm", Misc.Crabs="Misc.Crab", Clupeiformes="Clupeoidei", Misc.Crabs="Chionoecetes.spp.", Scorpaeniformes="Cyclopteridae", Teuthida="Teuthida", Trachiniformes="Pacific.sandfish", Misc.Crabs="Misc.Lithodidae", Mollusca="Bivalvia", Mollusca="Gastropod", Other="Misc.Org", Other="Misc.Invert", Salmoniformes="Salmonidae", Scorpaeniformes="Agonidae", Gadus.macrocephalus="Pacific.Cod", Chondrichthyes="Rajadae", Cnidaria="Cnidaria", Crustacea="Cumacea", Mollusca="Misc.Cephalopoda", Crustacea="Misc.Crustacea", Crustacea="Isopoda", Pleuronectiformes="Lepidopsetta.sp", Scorpaeniformes="Unid.Rockfish", Echinodermata="Brittle.Star", Mollusca="Misc.Mollusca", Echinodermata="Sea.Cucumber", Misc.Crabs="Misc.Brachyura", Scorpaeniformes="Sebastelobus", Echinodermata="Misc.Echinoderm", Teleostei="Fish.Eggs", Crustacea="Copepoda", Misc.Shrimps="Mysidacea", Ctenophora="Ctenophora", Mollusca="Pteropoda", Misc.Crabs="Opilio.Crab", Echinodermata="Sea.Urchin", Perciformes="Pholidae", Misc.Shrimps="Hyperiidea", Tunicata="Larvacea", Other="Misc.Bird", Other="Unid.Eggs", Misc.Shrimps="Capreillidea", Pleuronectiformes="Kamchatka.flounder", Misc.Shrimps="Misc.Decapoda", Tunicata="Tunicate", Pleuronectiformes="Northern.rock.sole", Pleuronectiformes="Southern.rock.sole", Chondrichthyes="Misc.Non.teleost.fish", Myctophiformes="Myctophidae", Other="Misc", Crustacea="Misc.Amphipoda", Pleuronectiformes="Pacific.halibut", Other="Chaetognatha", Gadiformes="Macrouridae", Echinodermata="Sand.Dollar", Scorpaeniformes="Misc.Hexagrammidae")

# Order prey taxa by phylogeny:
preyWTlong$Prey_Name = ordered(preyWTlong$Prey_Name, levels = c("Ctenophora", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Annelida", "Crustacea", "Euphausiacea", "Misc.Shrimps", "Misc.Crabs", "Echinodermata", "Tunicata", "Chondrichthyes", "Teleostei", "Clupeiformes", "Salmoniformes", "Osmeriformes", "Myctophiformes", "Gadiformes", "Gadus.chalcogrammus", "Gadus.macrocephalus", "Perciformes", "Trachiniformes", "Scorpaeniformes", "Pleuronectiformes", "Offal", "Other"))

##########################################################
### CALCULATE WEIGHTED PROPORTIONS ###
##########################################################
# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (food habits):
sumSizeDiets = preyWTgrid %>%
  group_by(Year, id2, FLBin) %>% 
  summarise(sumFreq_Diets = length(PredLength))

propSizeDiets = sumSizeDiets %>%
  group_by(Year, id2) %>%
  mutate(propFreq_Diets = sumFreq_Diets/sum(sumFreq_Diets)) 
propSizeDiets = na.omit(propSizeDiets)

# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (bottom trawl survey):
sumSizeLengths = lengths %>%
  group_by(Year, id2, FLBin) %>% 
  summarise(sumFreq_Lengths = sum(as.numeric(as.character(Frequency)) ))

propSizeLengths = sumSizeLengths %>%
  group_by(Year, id2) %>% 
  mutate(propFreq_Lengths = sumFreq_Lengths/sum(sumFreq_Lengths)) 
propSizeLengths = na.omit(propSizeLengths)

# Calculate sample weights from dividing length-based proportions of fish caught by those subsampled for gut content analysis in each survey year and grid cell:
LengthProp = merge(propSizeDiets, propSizeLengths)
LengthProp[is.na(LengthProp)] = 0

LengthProp$WT = LengthProp$propFreq_Lengths/LengthProp$propFreq_Diets
LengthProp = subset(LengthProp, select = c(Year:FLBin,propFreq_Diets,propFreq_Lengths:WT))

# Account for length-structured stomach sampling - Multiply prey mass (g) by sample weights according to predator length frequencies:
preyWTweighted = preyWTlong %>% left_join(LengthProp)
preyWTweighted = na.omit(preyWTweighted)
preyWTweighted = subset(preyWTweighted, preyWT_g > 0)
preyWTweighted$preyFLwt = preyWTweighted$preyWT_g * preyWTweighted$WT

# Select only sizes of fish encompassed in total predator biomass estimates:
preyWTweighted19 = subset(preyWTweighted, PredLength >= 19)

# Account for unequal sampling across the Gulf of Alaska - Multiply prey mass (g) by sample weights according to grid-specific biomass:
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
ATFabunPred = readRDS("2_RelativePredatorDensities/normATFabun.rds")
ATFabunPred = rename(ATFabunPred, Year = YEAR) # match other index components
ATFabunPred$Year = as.numeric(as.character(ATFabunPred$Year))

BioWT = ATFabunPred %>%
  group_by(Year) %>% 
  mutate(BioWT = ATFpredAbun / mean(ATFpredAbun))
BioWT = as.data.frame(BioWT)

BioWT = unique(BioWT[,c("id2", "Year", "BioWT")])

WeightedPreyData = preyWTweighted%>% left_join(BioWT)
WeightedPreyData$WTdiet = WeightedPreyData$preyFLwt * WeightedPreyData$BioWT

# Select only sizes of fish encompassed in total predator biomass estimates:
WeightedPreyData_19 = subset(WeightedPreyData, PredLength >= 19)
WeightedPreyData_19 = na.omit(WeightedPreyData_19)

# Save results:
saveRDS(WeightedPreyData_19, "4_ProportionsPollockConsumed/propWEP_ATF_19.rds")

# Plot proportions of pollock consumed, by survey year:
preyWT_YR = ggplot(WeightedPreyData_19, aes(x = Year, y = WTdiet, fill = Prey_Name, order = -as.numeric(Prey_Name))) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "horizontal", legend.position = "right", legend.box.margin = margin(-20,2.5,-5,-5), legend.background = element_rect(fill="transparent")) +
  theme(legend.text.align=0) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(axis.text = element_text(family="Arial", size=11)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=0, size=12)) +
  theme(axis.title.y = element_text(vjust=2.5, size=12)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=12)) +
  scale_fill_manual(values = c("Ctenophora"="lavenderblush2", "Cnidaria"="lightpink1", "Mollusca"="brown4", "Teuthida"="firebrick3", "Octopoda"="red", "Annelida"="chocolate4", "Crustacea"="chocolate", "Euphausiacea"="chocolate1", "Misc.Shrimps"="orange", "Misc.Crabs"="gold", "Echinodermata"="yellow", "Tunicata"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupeiformes"="forestgreen", "Salmoniformes"="lightcyan1", "Osmeriformes"="cadetblue1", "Argentiniformes"="deepskyblue3", "Myctophiformes"="deepskyblue4", "Gadiformes"="blue", "Gadus.chalcogrammus"="mediumblue", "Gadus.macrocephalus"="blue4", "Perciformes"="mediumpurple1", "Trachiniformes"="mediumorchid3", "Scorpaeniformes"="darkviolet", "Pleuronectiformes"="purple4", "Offal"="gray91", "Other"="azure4"), name="", labels = c("Ctenophora", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Annelida", "Crustacea", "Euphausiacea", "Misc. Shrimps", "Misc. Crabs", "Echinodermata", "Tunicata", "Chondrichthyes", "Teleostei", "Clupeiformes", "Salmoniformes", "Osmeriformes", "Argentiniformes", "Myctophiformes", "Gadiformes", expression(paste(italic("   Gadus chalcogrammus"))), expression(paste(italic("   Gadus macrocephalus"))), "Perciformes", "Trachiniformes", "Scorpaeniformes", "Pleuronectiformes", "Offal", "Other"), breaks = c("Ctenophora", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Annelida", "Crustacea", "Euphausiacea", "Misc.Shrimps", "Misc.Crabs", "Echinodermata", "Tunicata", "Chondrichthyes", "Teleostei", "Clupeiformes", "Salmoniformes", "Osmeriformes", "Argentiniformes", "Myctophiformes", "Gadiformes", "Gadus.chalcogrammus", "Gadus.macrocephalus", "Perciformes", "Trachiniformes", "Scorpaeniformes", "Pleuronectiformes", "Offal", "Other"), drop=FALSE) +
  xlab("") +
  ylab("Proportion of Prey by Weight") +
  scale_x_continuous(breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01))

preyWT_YR