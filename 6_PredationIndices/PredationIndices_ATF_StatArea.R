# This script file includes code necessary to calculate statistical area indices of predation for Walleye Pollock in the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predator: Arrowtooth Flounder (Atheresthes stomias). 

# When using any portion of the code herein, please cite: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

# Data Sources: All the data necessary to complete the following analyses can be found in the 'Data' folder on github. 
  # Total biomass estimates were obtained from the most recent stock assessment for Arrowtooth Flounder (Spies et al. 2017). 
  # Bottom trawl survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm (for methods, see von Szalay et al. 2016). 
  # Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 
  
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.
# Spies, I., K. Aydin, J. N. Ianelli, and W. Palsson. 2017. Assessment of the Arrowtooth Flounder stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 749–846.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##################################################################
# Input total biomass estimates from most recent stock assessment (ATF = age 1+, Gulf-wide, MT):
ATFpredationWEP = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), "TotBio_t" = c(
  1660800, 
  1773450,
  1770270,
  1835310,
  1957130,
  2035310,
  2069910,
  2054040,
  1962540,
  1826620,
  1701770,
  1571460))
str(ATFpredationWEP)

# Convert MT to g (comparable with units for annual ration and proportions of pollock consumed):
ATFpredationWEP$TotBio_g = ATFpredationWEP$TotBio_t * 1000000

#########################################################
# Input relative predator densities (unitless) predicted from GAMMs by survey year:
RelPredDensity = readRDS("2_RelativePredatorDensities/normATFabun.rds")
colnames(RelPredDensity)[7] = c("Year") # rename to match other data frames
RelPredDensity$Year = as.numeric(as.character(RelPredDensity$Year))

# Assign GOA statistical areas based on grid cell id:
RelPredDensity$StatArea = with(RelPredDensity, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the study area:  
RelPredDensity = as.data.frame(RelPredDensity)
RelPredDensity$StatArea = as.factor(RelPredDensity$StatArea)
RelPredDensity = subset(RelPredDensity, StatArea != "Other")
RelPredDensity$StatArea = ordered(RelPredDensity$StatArea, levels=c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Sum year-specific relative predator densities: 
require(dplyr)
RPD_sum = RelPredDensity %>%
  group_by(Year, StatArea) %>%
  summarise(NormAbun = sum(ATFnormAbun))

# Join to predator biomass data:
ATFpredationWEP = ATFpredationWEP %>% left_join(RPD_sum)

#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models:
Ration_ATF = readRDS("3_MeanAnnualRations/Rations_ATF.rds")

# Assign GOA statistical areas based on grid cell id:
Ration_ATF$StatArea = with(Ration_ATF, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the study area:  
Ration_ATF$StatArea = as.factor(Ration_ATF$StatArea)
Ration_ATF = subset(Ration_ATF, StatArea != "Other")
Ration_ATF$StatArea = ordered(Ration_ATF$StatArea, levels=c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
Ration_ATF_Cmax = Ration_ATF %>%
  group_by(Year, StatArea) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm=TRUE))

# Assume that the mean annual ration in 1990 in the Shumagin statistical area (missing temperature data) was the same as in 1993:
Shum93 = subset(Ration_ATF_Cmax, Year == 1993 & StatArea == "Shumagin")
Ration_ATF_Cmax = as.data.frame(Ration_ATF_Cmax)
Ration_ATF_Cmax = add_row(Ration_ATF_Cmax)
Ration_ATF_Cmax[nrow(Ration_ATF_Cmax), ] = c(1990, "Shumagin", Shum93[1,3])

# Join to total predator biomass and relative predator density estimates:
ATFpredationWEP = ATFpredationWEP %>% left_join(Ration_ATF_Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless):
propWEP_ATF = readRDS("4_ProportionsPollockConsumed/propWEP_ATF_19.rds")
propWEP_ATF$Year = as.numeric(as.character(propWEP_ATF$Year))

# Assign GOA statistical areas based on grid cell id:
propWEP_ATF$StatArea = with(propWEP_ATF, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the study area:  
propWEP_ATF$StatArea = as.factor(propWEP_ATF$StatArea)
propWEP_ATF = subset(propWEP_ATF, StatArea != "Other")
propWEP_ATF$StatArea = ordered(propWEP_ATF$StatArea, levels=c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Calculate year-specific proportions of prey consumed:
propWEP_ATF_mean = propWEP_ATF %>%
  group_by(Year, StatArea, Prey_Name) %>% 
  summarise(sumWT = sum(WTdiet)) %>%
  mutate(propWT = sumWT / sum(sumWT)) 

# Select only proportions of Walleye Pollock:
propWEP_ATFsub = subset(propWEP_ATF_mean, Prey_Name == "Gadus.chalcogrammus")

# Add zeros for when stomachs were sampled (n >= 5) but no pollock were consumed:
propWEP_ATFsub = as.data.frame(propWEP_ATFsub)
propWEP_ATFsub = add_row(propWEP_ATFsub)
propWEP_ATFsub[nrow(propWEP_ATFsub), ] = c(2005, "Yakutat", "Gadus.chalcogrammus", 0, 0)
propWEP_ATFsub$Year = as.numeric(as.character(propWEP_ATFsub$Year))

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
ATFpredWEP = ATFpredationWEP %>% right_join(propWEP_ATFsub)

# Remove subregions and survey years without stomach samples (Eastern, 1996 to 2001):
ATFpredWEP = subset(ATFpredWEP, StatArea != "Yakutat" | Year < 1996 | Year > 2001)
ATFpredWEP = subset(ATFpredWEP, StatArea != "Southeastern" | Year > 2003)

#########################################################
# Calculate year-specific predation by Arrowtooth Flounder (g; all pollock age classes combined):
ATFpredWEP$propWT = as.numeric(as.character(ATFpredWEP$propWT))
ATFpredWEP$predWEP_g_max = with(ATFpredWEP, TotBio_g * NormAbun * Annual_Cmax * propWT)

# Convert g of pollock consumed to metric tons:
ATFpredWEP$predWEP_t_max = round((ATFpredWEP$predWEP_g_max / 1000000), digits=0)

#########################################################
# Calculate age-specific predation on pollock by Arrowtooth Flounder (g):

# Input proportions of pollock consumed by age class:
propWEP_age_ATF = readRDS("5_AgeCompositionsPrey/propWEP_age_ATF.rds")

require(reshape2)
propWEP_Age_ATF = dcast(propWEP_age_ATF, Year ~ AgeClass, value.var = "propWT", fun.aggregrate = sum) # long to wide format
colnames(propWEP_Age_ATF) = c("Year", "prop_Age0", "prop_Age1", "prop_Age2", "prop_Age3plus")
propWEP_Age_ATF$Year = as.numeric(as.character(propWEP_Age_ATF$Year))

# Join to total predator biomasses, relative predator densities, mean annual rations, and proportions of pollock consumed:
ATFpredWEP_age = ATFpredWEP %>% left_join(propWEP_Age_ATF)

# Calculate age-specific consumption:
ATFpredWEP_age$predWEP_t_max_Age0 = round(ATFpredWEP_age$predWEP_t_max * ATFpredWEP_age$prop_Age0, digits=0)
ATFpredWEP_age$predWEP_t_max_Age1 = round(ATFpredWEP_age$predWEP_t_max * ATFpredWEP_age$prop_Age1, digits=0)
ATFpredWEP_age$predWEP_t_max_Age2 = round(ATFpredWEP_age$predWEP_t_max * ATFpredWEP_age$prop_Age2, digits=0)
ATFpredWEP_age$predWEP_t_max_Age3plus = round(ATFpredWEP_age$predWEP_t_max * ATFpredWEP_age$prop_Age3plus, digits=0)

#########################################################
# Convert dataframe back to long format for plotting:
ATFpredWEP_age_t_max = melt(ATFpredWEP_age, id.vars = c("Year", "StatArea"), measure.vars = c("predWEP_t_max", "predWEP_t_max_Age0", "predWEP_t_max_Age1", "predWEP_t_max_Age2", "predWEP_t_max_Age3plus"), variable.name = "AgeClass", value.name = "MT")
levels(ATFpredWEP_age_t_max$AgeClass) = c(" All Ages", " 0", " 1", " 2", " 3+")
ATFpredWEP_age_t_max$AgeClass = ordered(ATFpredWEP_age_t_max$AgeClass, levels = c(" All Ages", " 3+", " 2", " 1", " 0"))
ATFpredWEP_age_t_max$millMT = round((ATFpredWEP_age_t_max$MT / 1000000), digits=2)