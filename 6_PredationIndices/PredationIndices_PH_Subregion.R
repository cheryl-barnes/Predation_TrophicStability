# This script file includes code necessary to calculate an index of predation for Walleye Pollock in the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predator: Pacific Halibut (Hippoglossus stenolepis). 

# When using any portion of the code herein, please cite: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

# Data Sources: All the data necessary to complete the following analyses can be found in the 'Data' folder on github. 
  # Total biomass estimates were obtained from the most recent stock assessment for Pacific Halibut (Stewart and Hicks 2017). 
  # Bottom trawl survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm (for methods, see von Szalay et al. 2016). Setline survey data (1998 to 2017) were collected by the International Pacific Halibut Commission (IPHC) and are publicly available at: https://iphc.int/data/fiss-data-query (for methods, see Clark and Hare 2006). 
  # Food habits data (a1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 
  
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Clark, W. G., and S. R. Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. International Pacific Halibut Commission Scientific Report 83. 
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##################################################################
# Input total biomass estimates from most recent stock assessment (PH = age 8+, coast-wide [EBS, GOA, AI, BC, West Coast], Mlb - converted to MT):
# Coast-wide biomass estimates were adjusted by proportions of over 32 inch biomass that were sampled in the GOA as part of IPHC's setline survey (IPHC Regulatory Areas 4A, 3B, 3A, and 2C). No setline survey data were available for 1990, 1993, or 1996. Thus, biomass estimates for those years were calculated (from 1996) based on trends from the bottom trawl survey - which were highly correlated with setline survey biomass.
PHpredationWEP = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), "TotBio_t" = c(
  461359,  
  653478, 
  625352, 
  615092,
  465851,
  428400,
  368729,
  340579,
  256951,
  237279,
  252087,
  223213))

# Convert MT to g (comparable with units for annual ration and proportions of pollock consumed):
PHpredationWEP$TotBio_g = PHpredationWEP$TotBio_t * 1000000

#########################################################
# Input relative predator densities (unitless) predicted from GAMMs by survey year:
RelPredDensity = readRDS("2_RelativePredatorDensities/normPHabun.rds") # rename to match other data frames
RelPredDensity$Year = as.numeric(as.character(RelPredDensity$Year))

# Assign GOA subregions based on grid cell id:
RelPredDensity$SubRegion = with(RelPredDensity, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Western", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Central", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Central", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Eastern", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Eastern", "Other"))))))

# Eliminate locations outside of the study area (western, central, and eastern GOA):
RelPredDensity$SubRegion = as.factor(RelPredDensity$SubRegion)
RelPredDensity = as.data.frame(RelPredDensity)
RelPredDensity = subset(RelPredDensity, SubRegion != "Other")
RelPredDensity$SubRegion = ordered(RelPredDensity$SubRegion, levels=c("Western", "Central", "Eastern"))

# Sum year-specific relative predator densities: 
require(dplyr)
RPD_sum = RelPredDensity %>%
  group_by(Year, SubRegion) %>%
  summarise(NormAbun = sum(relAbun))
RPD_sum = as.data.frame(RPD_sum)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Assign mean relative densities (entire study area) to missing survey years (1990, 1993, and 1996):

# Western GOA:
RPDw = subset(RPD_sum, SubRegion == "Western")
meanRPDw = mean(RPDw$NormAbun)

RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Western", meanRPDw)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Western", meanRPDw)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Western", meanRPDw)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Central GOA:
RPDc = subset(RPD_sum, SubRegion == "Central")
meanRPDc = mean(RPDc$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Central", meanRPDc)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Central", meanRPDc)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Central", meanRPDc)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Eastern GOA:
RPDe = subset(RPD_sum, SubRegion == "Eastern")
meanRPDe = mean(RPDe$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Eastern", meanRPDe)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Eastern", meanRPDe)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Eastern", meanRPDe)

# Remove survey years not sampled by the AFSC bottom trawl survey:
RPD_sum = RPD_sum[RPD_sum$Year %in% c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015), ]
RPD_sum$Year = as.numeric(RPD_sum$Year)

# Join to predator biomass data:
PHpredationWEP = PHpredationWEP %>% left_join(RPD_sum)

#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models:
Ration_PH = readRDS("3_MeanAnnualRations/Rations_PH.rds")

# Assign GOA subregions based on grid cell id:
Ration_PH$SubRegion = with(Ration_PH, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Western", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Central", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Central", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Eastern", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Eastern", "Other"))))))

# Eliminate locations outside of the study area (western, central, and eastern GOA):
Ration_PH$SubRegion = as.factor(Ration_PH$SubRegion)
Ration_PH = subset(Ration_PH, SubRegion != "Other")
Ration_PH$SubRegion = ordered(as.factor(Ration_PH$SubRegion), levels=c("Western", "Central", "Eastern"))

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
Ration_PH_Cmax = Ration_PH %>%
  group_by(Year, SubRegion) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm=TRUE))

# Assume that the mean annual ration in 1990 in the western GOA (missing temperature data) was the same as in 1993:
West93 = subset(Ration_PH_Cmax, Year == 1993 & SubRegion == "Western")
Ration_PH_Cmax = as.data.frame(Ration_PH_Cmax)
Ration_PH_Cmax = add_row(Ration_PH_Cmax)
Ration_PH_Cmax[nrow(Ration_PH_Cmax), ] = c(1990, "Western", West93[1,3])

# Join to total predator biomass and relative predator density estimates:
PHpredationWEP = PHpredationWEP %>% left_join(Ration_PH_Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless):
propWEP_PH = readRDS("4_ProportionsPollockConsumed/propWEP_PH_82.rds")
propWEP_PH$Year = as.numeric(as.character(propWEP_PH$Year))

# Assign GOA subregions based on grid cell id:
propWEP_PH$SubRegion = with(propWEP_PH, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Western", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Central", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Central", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Eastern", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Eastern", "Other"))))))

# Eliminate locations outside of the study area (western, central, and eastern GOA):
propWEP_PH$SubRegion = as.factor(propWEP_PH$SubRegion)
propWEP_PH = as.data.frame(propWEP_PH)
propWEP_PH = subset(propWEP_PH, SubRegion != "Other")
propWEP_PH$SubRegion = ordered(as.factor(propWEP_PH$SubRegion), levels=c("Western", "Central", "Eastern"))

# Calculate year-specific proportions of prey consumed:
propWEP_PH_mean = propWEP_PH %>%
  group_by(Year, SubRegion, Prey_Name) %>% 
  summarise(sumWT = sum(WTdiet)) %>%
  mutate(propWT = sumWT / sum(sumWT))

# Select only proportions of Walleye Pollock:
propWEP_PHsub = subset(propWEP_PH_mean, Prey_Name == "Gadus.chalcogrammus")

# Add zeros for when stomachs were sampled (n >= 5) but no pollock were consumed:
propWEP_PHsub = as.data.frame(propWEP_PHsub)
propWEP_PHsub = add_row(as.data.frame(propWEP_PHsub))
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(2003, "Eastern", "Gadus.chalcogrammus", 0, 0)
propWEP_PHsub$Year = as.numeric(propWEP_PHsub$Year)

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
PHpredWEP = PHpredationWEP %>% right_join(propWEP_PHsub)

# Remove subregions and survey years without stomach samples (Eastern, 1996 to 2001):
PHpredWEP = subset(PHpredWEP, SubRegion != "Eastern" | Year < 1996 | Year > 2001)

#########################################################
# Calculate year-specific predation by Arrowtooth Flounder (g; all pollock age classes combined):
PHpredWEP$NormAbun = as.numeric(PHpredWEP$NormAbun)
PHpredWEP$propWT = as.numeric(PHpredWEP$propWT)
PHpredWEP$predWEP_g_max = with(PHpredWEP, TotBio_g * NormAbun * Annual_Cmax * propWT)

# Convert g of pollock consumed to metric tons:
PHpredWEP$predWEP_t_max = round((PHpredWEP$predWEP_g_max / 1000000), digits=0)

#########################################################
# Calculate age-specific predation on pollock by Arrowtooth Flounder (g):

# Input proportions of pollock consumed by age class:
propWEP_age_PH = readRDS("5_AgeCompositionsPrey/propWEP_age_PH.rds")

require(reshape2)
propWEP_Age_PH = dcast(propWEP_age_PH, Year ~ AgeClass, value.var = "propWT", fun.aggregrate = sum) # long to wide format
colnames(propWEP_Age_PH) = c("Year", "prop_Age0", "prop_Age1", "prop_Age2", "prop_Age3plus")
propWEP_Age_PH$Year = as.numeric(as.character(propWEP_Age_PH$Year))

# Join to total predator biomasses, relative predator densities, mean annual rations, and proportions of pollock consumed:
PHpredWEP_age = PHpredWEP %>% left_join(propWEP_Age_PH)

# Calculate age-specific consumption:
PHpredWEP_age$predWEP_t_max_Age0 = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age0, digits=0)
PHpredWEP_age$predWEP_t_max_Age1 = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age1, digits=0)
PHpredWEP_age$predWEP_t_max_Age2 = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age2, digits=0)
PHpredWEP_age$predWEP_t_max_Age3plus = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age3plus, digits=0)

#########################################################
# Convert dataframe back to long format for plotting:
PHpredWEP_age_t_max = melt(PHpredWEP_age, id.vars = c("Year", "SubRegion"), measure.vars = c("predWEP_t_max", "predWEP_t_max_Age0", "predWEP_t_max_Age1", "predWEP_t_max_Age2", "predWEP_t_max_Age3plus"), variable.name = "AgeClass", value.name = "MT")
levels(PHpredWEP_age_t_max$AgeClass) = c(" All Ages", " 0", " 1", " 2", " 3+")
PHpredWEP_age_t_max$AgeClass = ordered(PHpredWEP_age_t_max$AgeClass, levels = c(" All Ages", " 3+", " 2", " 1", " 0"))
PHpredWEP_age_t_max$millMT = round((PHpredWEP_age_t_max$MT / 1000000), digits=2)