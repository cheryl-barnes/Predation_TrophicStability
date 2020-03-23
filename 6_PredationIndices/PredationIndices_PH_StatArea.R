# This script file includes code necessary to calculate statistical area indices of predation for Walleye Pollock in the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predator: Pacific Halibut (Hippoglossus stenolepis). 

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
RelPredDensity = readRDS("2_RelativePredatorDensities/normPHabun.rds")
RelPredDensity$Year = as.numeric(as.character(RelPredDensity$Year))

# Assign GOA statistical areas based on grid cell id:
RelPredDensity$StatArea = with(RelPredDensity, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the study area: 
RelPredDensity$StatArea = as.factor(RelPredDensity$StatArea)
RelPredDensity = as.data.frame(RelPredDensity)
RelPredDensity = subset(RelPredDensity, StatArea != "Other")
RelPredDensity$StatArea = ordered(RelPredDensity$StatArea, levels=c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Sum year-specific relative predator densities: 
require(dplyr)
RPD_sum = RelPredDensity %>%
  group_by(Year, StatArea) %>%
  summarise(NormAbun = sum(relAbun))
RPD_sum = as.data.frame(RPD_sum)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Assign mean relative densities (entire study area) to missing survey years (1990, 1993, and 1996):

# Shumagin:
RPDshum = subset(RPD_sum, StatArea == "Shumagin")
meanRPDshum = mean(RPDshum$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Shumagin", meanRPDshum)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Shumagin", meanRPDshum)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Shumagin", meanRPDshum)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Chirikof:
RPDchir = subset(RPD_sum, StatArea == "Chirikof")
meanRPDchir = mean(RPDchir$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Chirikof", meanRPDchir)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Chirikof", meanRPDchir)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Chirikof", meanRPDchir)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Kodiak:
RPDkod = subset(RPD_sum, StatArea == "Kodiak")
meanRPDkod = mean(RPDkod$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Kodiak", meanRPDkod)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Kodiak", meanRPDkod)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Kodiak", meanRPDkod)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Yakutat:
RPDyak = subset(RPD_sum, StatArea == "Yakutat")
meanRPDyak = mean(RPDyak$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Yakutat", meanRPDyak)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Yakutat", meanRPDyak)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Yakutat", meanRPDyak)
RPD_sum$NormAbun = as.numeric(RPD_sum$NormAbun)

# Southeastern:
RPDse = subset(RPD_sum, StatArea == "Southeastern")
meanRPDse = mean(RPDse$NormAbun)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1990, "Southeastern", meanRPDse)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1993, "Southeastern", meanRPDse)
RPD_sum = add_row(RPD_sum)
RPD_sum[nrow(RPD_sum), 1:3] = c(1996, "Southeastern", meanRPDse)

# Remove survey years not sampled by the AFSC bottom trawl survey:
RPD_sum = RPD_sum[RPD_sum$Year %in% c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015), ]
RPD_sum$Year = as.numeric(RPD_sum$Year)

# Join to predator biomass data:
PHpredationWEP = PHpredationWEP %>% left_join(RPD_sum)

#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models:
Ration_PH = readRDS("3_MeanAnnualRations/Rations_PH.rds")

# Assign GOA statistical areas based on grid cell id:
Ration_PH$StatArea = with(Ration_PH, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the study area:  
Ration_PH$StatArea = as.factor(Ration_PH$StatArea)
Ration_PH = subset(Ration_PH, StatArea != "Other")
Ration_PH$StatArea = ordered(Ration_PH$StatArea, levels=c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
Ration_PH_Cmax = Ration_PH %>%
  group_by(Year, StatArea) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm=TRUE))

# Assume that mean annual rations in Yakutat in 2003 and 2005 were the same as in 2007:
Ration_PH_Cmax = as.data.frame(Ration_PH_Cmax)
Yak07 = subset(Ration_PH_Cmax, Year == 2007 & StatArea == "Yakutat")
Ration_PH_Cmax = add_row(Ration_PH_Cmax)
Ration_PH_Cmax[nrow(Ration_PH_Cmax), ] = c(2003, "Yakutat", Yak07[1,3])
Ration_PH_Cmax = add_row(Ration_PH_Cmax)
Ration_PH_Cmax[nrow(Ration_PH_Cmax), ] = c(2005, "Yakutat", Yak07[1,3])

# Assume that the mean annual ration in 1990 in the Shumagin statistical area (missing temperature data) was the same as in 1993:
Shum93 = subset(Ration_PH_Cmax, Year == 1993 & StatArea == "Shumagin")
Ration_PH_Cmax = as.data.frame(Ration_PH_Cmax)
Ration_PH_Cmax = add_row(Ration_PH_Cmax)
Ration_PH_Cmax[nrow(Ration_PH_Cmax), ] = c(1990, "Shumagin", Shum93[1,3])
Ration_PH_Cmax$Year = as.numeric(Ration_PH_Cmax$Year)

# Join to total predator biomass and relative predator density estimates:
PHpredationWEP = PHpredationWEP %>% left_join(Ration_PH_Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless):
propWEP_PH = readRDS("4_ProportionsPollockConsumed/propWEP_PH_82.rds")
propWEP_PH$Year = as.numeric(as.character(propWEP_PH$Year))

# Assign GOA statistical areas based on grid cell id:
propWEP_PH$StatArea = with(propWEP_PH, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the study area: 
propWEP_PH$StatArea = as.factor(propWEP_PH$StatArea)
propWEP_PH = as.data.frame(propWEP_PH)
propWEP_PH = subset(propWEP_PH, StatArea != "Other")
propWEP_PH$StatArea = ordered(propWEP_PH$StatArea, levels=c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Calculate year-specific proportions of prey consumed:
propWEP_PH_mean = propWEP_PH %>%
  group_by(Year, StatArea, Prey_Name) %>% 
  summarise(sumWT = sum(WTdiet)) %>%
  mutate(propWT = sumWT / sum(sumWT)) 

# Select only proportions of Walleye Pollock:
propWEP_PHsub = subset(propWEP_PH_mean, Prey_Name == "Gadus.chalcogrammus")

propWEP_PHsub$Year = as.numeric(as.character(propWEP_PHsub$Year))
propWEP_PHsub$propWT = as.numeric(as.character(propWEP_PHsub$propWT))

# Add zeros for when stomachs were sampled (n >= 5) but no pollock were consumed:
propWEP_PHsub = as.data.frame(propWEP_PHsub)
propWEP_PHsub = add_row(propWEP_PHsub)
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(2005, "Kodiak", "Gadus.chalcogrammus", 0, 0)
propWEP_PHsub = add_row(propWEP_PHsub)
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(2009, "Yakutat", "Gadus.chalcogrammus", 0, 0)
propWEP_PHsub = add_row(propWEP_PHsub)
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(2003, "Southeastern", "Gadus.chalcogrammus", 0, 0)

# Assume proportions of pollock consumed in Chirikof in 1999 and in Yakutat in 2003 and 2005 were equal to area-specific means:
propWEP_PHsub = as.data.frame(propWEP_PHsub)
propWEP_PHsub$propWT = as.numeric(propWEP_PHsub$propWT)
ChirProp = subset(propWEP_PHsub, StatArea == "Chirikof")
propWEP_PHsub = add_row(propWEP_PHsub)
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(1999, "Chirikof", "Gadus.chalcogrammus", 0, mean(ChirProp$propWT))

propWEP_PHsub = as.data.frame(propWEP_PHsub)
propWEP_PHsub$propWT = as.numeric(propWEP_PHsub$propWT)
YakProp = subset(propWEP_PHsub, StatArea == "Yakutat")
propWEP_PHsub = add_row(propWEP_PHsub)
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(2003, "Yakutat", "Gadus.chalcogrammus", 0, mean(YakProp$propWT))
propWEP_PHsub = add_row(propWEP_PHsub)
propWEP_PHsub[nrow(propWEP_PHsub), ] = c(2005, "Yakutat", "Gadus.chalcogrammus", 0, mean(YakProp$propWT))
propWEP_PHsub$Year = as.numeric(propWEP_PHsub$Year)

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
PHpredWEP = PHpredationWEP %>% right_join(propWEP_PHsub)

# Remove subregions and survey years without stomach samples (Eastern, 1996 to 2001):
PHpredWEP = subset(PHpredWEP, StatArea != "Yakutat" | Year < 1996 | Year > 2001)

# Remove areas and years with < 5 stomach samples:
PHpredWEP = subset(PHpredWEP, StatArea != "Southeastern" | Year > 2003)

#########################################################
# Calculate year-specific predation by Pacific Halibut (g; all pollock age classes combined):
PHpredWEP$NormAbun = as.numeric(PHpredWEP$NormAbun)
PHpredWEP$Annual_Cmax = as.numeric(PHpredWEP$Annual_Cmax)
PHpredWEP$propWT = as.numeric(PHpredWEP$propWT)
PHpredWEP$predWEP_g_max = with(PHpredWEP, TotBio_g * NormAbun * Annual_Cmax * propWT)

# Convert g of pollock consumed to metric tons:
PHpredWEP$predWEP_t_max = round((PHpredWEP$predWEP_g_max / 1000000), digits=0)

#########################################################
# Calculate age-specific predation on pollock by Pacific Halibut (g):

# Input proportions of pollock consumed by age class:
propWEP_age_PH = readRDS("5_AgeCompositionsPrey/propWEP_age_PH.rds")

require(reshape2)
propWEP_Age_PH = dcast(propWEP_age_PH, Year ~ AgeClass, value.var = "propWT", fun.aggregrate = sum) # long to wide format
colnames(propWEP_Age_PH) = c("Year", "prop_Age0", "prop_Age1", "prop_Age2", "prop_Age3plus")
propWEP_Age_PH$Year = as.numeric(as.character(propWEP_Age_PH$Year))

# Join to total predator biomasses, relative predator densities, mean annual rations, and proportions of pollock consumed:
PHpredWEP_age = PHpredWEP %>% left_join(propWEP_Age_PH)

PHpredWEP_age$StatArea = as.factor(PHpredWEP_age$StatArea)
PHpredWEP_age$StatArea = ordered(PHpredWEP_age$StatArea, levels = c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Calculate age-specific consumption:
PHpredWEP_age$predWEP_t_max_Age0 = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age0, digits=0)
PHpredWEP_age$predWEP_t_max_Age1 = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age1, digits=0)
PHpredWEP_age$predWEP_t_max_Age2 = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age2, digits=0)
PHpredWEP_age$predWEP_t_max_Age3plus = round(PHpredWEP_age$predWEP_t_max * PHpredWEP_age$prop_Age3plus, digits=0)

#########################################################
# Convert dataframe back to long format for plotting:
PHpredWEP_age_t_max = melt(PHpredWEP_age, id.vars = c("Year", "StatArea"), measure.vars = c("predWEP_t_max", "predWEP_t_max_Age0", "predWEP_t_max_Age1", "predWEP_t_max_Age2", "predWEP_t_max_Age3plus"), variable.name = "AgeClass", value.name = "MT")
levels(PHpredWEP_age_t_max$AgeClass) = c(" All Ages", " 0", " 1", " 2", " 3+")
PHpredWEP_age_t_max$AgeClass = ordered(PHpredWEP_age_t_max$AgeClass, levels = c(" All Ages", " 3+", " 2", " 1", " 0"))
PHpredWEP_age_t_max$millMT = round((PHpredWEP_age_t_max$MT / 1000000), digits=2)