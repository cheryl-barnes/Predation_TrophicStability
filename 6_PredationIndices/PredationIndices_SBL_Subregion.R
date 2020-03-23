# This script file includes code necessary to calculate subregional indices of predation for Walleye Pollock in the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predator: Sablefish (Anoplopoma fimbria). 

# When using any portion of the code herein, please cite: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

# Data Sources: All the data necessary to complete the following analyses can be found in the 'Data' folder on github. 
  # Total biomass estimates were obtained from the most recent stock assessment for Sablefish ( Hanselman et al. 2017). 
  # Longline survey data (1990 to 2017) were collected by personnel at the Alaska Fisheries Science Center's Auke Bay Laboratories and can be found at https://www.afsc.noaa.gov/maps/longline/Map.php (for methods, see Sigler and Zenger 1989). 
  # Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 
  
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Hanselman, D. H., C. J. Rodgveller, C. R. Lunsford, and K. H. Fenske. 2017. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report 327–502.
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##################################################################
# Input total biomass estimates from most recent stock assessment (SBL = age 2+, Region-specific, MT):
  # provided by S. Barbeaux via email [Model 17.09.35_3.24U]
SBLpredationWEP = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), "TotBio_t" = c(
  251000,
  261000,
  200000,
  183000,
  182000,
  202000,
  197000,
  183000,
  164000,
  181000,
  157000,
  140000)) 

# Convert MT to g (comparable with units for annual ration and proportions of pollock consumed):
SBLpredationWEP$TotBio_g = SBLpredationWEP$TotBio_t * 1000000

#########################################################
# Input relative predator densities (unitless) predicted from GAMMs by survey year:
RelPredDensity = readRDS("2_RelativePredatorDensities/normSBLabun.rds")
colnames(RelPredDensity)[8] = c("Year") # rename to match other data frames
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
RelPredDensity = subset(RelPredDensity, SubRegion != "Other")
RelPredDensity$SubRegion = ordered(as.factor(RelPredDensity$SubRegion), levels=c("Western", "Central", "Eastern"))

# Sum year-specific relative predator densities: 
require(dplyr)
RPD_sum = RelPredDensity %>%
  group_by(Year, SubRegion) %>%
  summarise(NormAbun = sum(relAbun))

# Join to predator biomass data:
SBLpredationWEP = SBLpredationWEP %>% left_join(RPD_sum)

#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models:
Ration_SBL = readRDS("3_MeanAnnualRations/Rations_SBL.rds")

# Assign GOA subregions based on grid cell id:
Ration_SBL$SubRegion = with(Ration_SBL, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Western", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Central", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Central", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Eastern", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Eastern", "Other"))))))

# Eliminate locations outside of the study area (western, central, and eastern GOA):
Ration_SBL$SubRegion = as.factor(Ration_SBL$SubRegion)
Ration_SBL = subset(Ration_SBL, SubRegion != "Other")
Ration_SBL$SubRegion = ordered(as.factor(Ration_SBL$SubRegion), levels=c("Western", "Central", "Eastern"))

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
Ration_SBL_Cmax = Ration_SBL %>%
  group_by(Year, SubRegion) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy))
range(Ration_SBL_Cmax$Annual_Cmax)

# Assume that the mean annual ration in 1990 in the western GOA (missing temperature data) was the same as in 1993:
West93 = subset(Ration_SBL_Cmax, Year == 1993 & SubRegion == "Western")
Ration_SBL_Cmax = as.data.frame(Ration_SBL_Cmax)
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(1990, "Western", West93[1,3])

# Assume that the mean annual rations in 2013 and 2015 (missing Sablefish data) were the same as in 2011:

# Western GOA:
Ration_SBL_Cmax = as.data.frame(Ration_SBL_Cmax)
West11 = subset(Ration_SBL_Cmax, Year == 2011 & SubRegion == "Western")
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(2013, "Western", West11[1,3])
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(2015, "Western", West11[1,3])

# Central GOA:
Ration_SBL_Cmax = as.data.frame(Ration_SBL_Cmax)
Cent11 = subset(Ration_SBL_Cmax, Year == 2011 & SubRegion == "Central")
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(2013, "Central", Cent11[1,3])
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(2015, "Central", Cent11[1,3])

# Eastern GOA:
Ration_SBL_Cmax = as.data.frame(Ration_SBL_Cmax)
East11 = subset(Ration_SBL_Cmax, Year == 2011 & SubRegion == "Eastern")
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(2013, "Eastern", East11[1,3])
Ration_SBL_Cmax = add_row(Ration_SBL_Cmax)
Ration_SBL_Cmax[nrow(Ration_SBL_Cmax), ] = c(2015, "Eastern", East11[1,3])
Ration_SBL_Cmax$Year = as.numeric(Ration_SBL_Cmax$Year)

# Join to total predator biomass and relative predator density estimates:
SBLpredationWEP = SBLpredationWEP %>% left_join(Ration_SBL_Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless):
propWEP_SBL = readRDS("4_ProportionsPollockConsumed/propWEP_SBL_45.rds")
propWEP_SBL$Year = as.numeric(as.character(propWEP_SBL$Year))

# Assign GOA subregions based on grid cell id:
propWEP_SBL$SubRegion = with(propWEP_SBL, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Western", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Central", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Central", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Eastern", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Eastern", "Other"))))))

# Eliminate locations outside of the study area (western, central, and eastern GOA):
propWEP_SBL$SubRegion = as.factor(propWEP_SBL$SubRegion)
propWEP_SBL = subset(propWEP_SBL, SubRegion != "Other")
propWEP_SBL$SubRegion = ordered(as.factor(propWEP_SBL$SubRegion), levels=c("Western", "Central", "Eastern"))

# Calculate year-specific proportions of prey consumed:
propWEP_SBL_mean = propWEP_SBL %>%
  group_by(Year, SubRegion, Prey_Name) %>% 
  summarise(sumWT = sum(WTdiet)) %>%
  mutate(propWT = sumWT / sum(sumWT))

# Select only proportions of Walleye Pollock:
propWEP_SBLsub = subset(propWEP_SBL_mean, Prey_Name == "Gadus.chalcogrammus")

# Add zeros for when stomachs were sampled (n >= 5) but no pollock were consumed:
propWEP_SBLsub = as.data.frame(propWEP_SBLsub)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(1990, "Western", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(1999, "Western", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2003, "Western", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2007, "Western", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2009, "Western", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2011, "Western", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(1996, "Central", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(1993, "Eastern", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2005, "Eastern", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2007, "Eastern", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2009, "Eastern", "Gadus.chalcogrammus", 0, 0)
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2011, "Eastern", "Gadus.chalcogrammus", 0, 0)

# Assume proportions in 2003 (Eastern only), 2013, and 2015 were equal to the subregional mean:
propWEP_SBLsub = as.data.frame(propWEP_SBLsub)
propWEP_SBLsub$propWT = as.numeric(propWEP_SBLsub$propWT)

# Calculate mean proportions of pollock consumed by subregion:
meanSubRegionProp = propWEP_SBLsub %>%
  group_by(SubRegion) %>%
  summarise(meanProp = mean(propWT))

# Add zeros for when stomachs were sampled (n >= 5) but no pollock were consumed:
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2013, "Western", "Gadus.chalcogrammus", NA, meanSubRegionProp[1,2])
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2015, "Western", "Gadus.chalcogrammus", NA, meanSubRegionProp[1,2])
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2013, "Central", "Gadus.chalcogrammus", NA, meanSubRegionProp[2,2])
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2015, "Central", "Gadus.chalcogrammus", NA, meanSubRegionProp[2,2])
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2003, "Eastern", "Gadus.chalcogrammus", NA, meanSubRegionProp[3,2])
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2013, "Eastern", "Gadus.chalcogrammus", NA, meanSubRegionProp[3,2])
propWEP_SBLsub = add_row(propWEP_SBLsub)
propWEP_SBLsub[nrow(propWEP_SBLsub), ] = c(2015, "Eastern", "Gadus.chalcogrammus", NA, meanSubRegionProp[3,2])
propWEP_SBLsub = as.data.frame(propWEP_SBLsub)
propWEP_SBLsub$propWT = as.numeric(propWEP_SBLsub$propWT)
propWEP_SBLsub$Year = as.numeric(propWEP_SBLsub$Year)

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
SBLpredWEP = SBLpredationWEP %>% right_join(propWEP_SBLsub)

# Remove subregions and survey years without stomach samples (Eastern, 1996 to 2001):
SBLpredWEP = subset(SBLpredWEP, SubRegion != "Eastern" | Year < 1996 | Year > 2001)

#########################################################
# Calculate year-specific predation by Arrowtooth Flounder (g; all pollock age classes combined):
SBLpredWEP$Annual_Cmax = as.numeric(SBLpredWEP$Annual_Cmax)
SBLpredWEP$predWEP_g_max = with(SBLpredWEP, TotBio_g * NormAbun * Annual_Cmax * propWT)

# Convert g of pollock consumed to metric tons:
SBLpredWEP$predWEP_t_max = round((SBLpredWEP$predWEP_g_max / 1000000), digits=0)

#########################################################
# Calculate age-specific predation on pollock by Arrowtooth Flounder (g):

# Input proportions of pollock consumed by age class:
propWEP_age_SBL = readRDS("5_AgeCompositionsPrey/propWEP_age_SBL.rds")

require(reshape2)
propWEP_Age_SBL = dcast(propWEP_age_SBL, Year ~ AgeClass, value.var = "propWT", fun.aggregrate = sum) # long to wide format
colnames(propWEP_Age_SBL) = c("Year", "prop_Age0", "prop_Age1", "prop_Age2", "prop_Age3plus")
propWEP_Age_SBL$Year = as.numeric(as.character(propWEP_Age_SBL$Year))

# Join to total predator biomasses, relative predator densities, mean annual rations, and proportions of pollock consumed:
SBLpredWEP_age = SBLpredWEP %>% left_join(propWEP_Age_SBL)

# Calculate age-specific consumption:
SBLpredWEP_age$predWEP_t_max_Age0 = round(SBLpredWEP_age$predWEP_t_max * SBLpredWEP_age$prop_Age0, digits=0)
SBLpredWEP_age$predWEP_t_max_Age1 = round(SBLpredWEP_age$predWEP_t_max * SBLpredWEP_age$prop_Age1, digits=0)
SBLpredWEP_age$predWEP_t_max_Age2 = round(SBLpredWEP_age$predWEP_t_max * SBLpredWEP_age$prop_Age2, digits=0)
SBLpredWEP_age$predWEP_t_max_Age3plus = round(SBLpredWEP_age$predWEP_t_max * SBLpredWEP_age$prop_Age3plus, digits=0)

#########################################################
# Convert dataframe back to long format for plotting:
SBLpredWEP_age_t_max = melt(SBLpredWEP_age, id.vars = c("Year", "SubRegion"), measure.vars = c("predWEP_t_max","predWEP_t_max_Age0","predWEP_t_max_Age1", "predWEP_t_max_Age2", "predWEP_t_max_Age3plus"), variable.name = "AgeClass", value.name = "MT")
levels(SBLpredWEP_age_t_max$AgeClass) = c(" All Ages", " 0", " 1", " 2", " 3+")
SBLpredWEP_age_t_max$AgeClass = ordered(SBLpredWEP_age_t_max$AgeClass, levels = c(" All Ages", " 3+", " 2", " 1", " 0"))
SBLpredWEP_age_t_max$millMT = round((SBLpredWEP_age_t_max$MT / 1000000), digits=2)