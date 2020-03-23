# This script file includes code necessary to calculate an index of predation for Walleye Pollock within the assessment are of the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predator: Walleye Pollock (Gadus chalcogrammus) conspecifics. 

# When using any portion of the code herein, please cite: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

# Data Sources: All the data necessary to complete the following analyses can be found in the 'Data' folder on github. 
  # Total biomass estimates were obtained from the most recent stock assessment for Walleye Pollock (Dorn et al. 2017). 
  # Bottom trawl survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm (for methods, see von Szalay et al. 2016). 
  # Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 
  
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Dorn, M., K. Aydin, B. Fissel, D. Jones, A. McCarthy, W. Palsson, and K. Spalinger. 2017. Assessment of the Walleye Pollock stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 47–182.
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##################################################################
# Input total biomass estimates from most recent stock assessment (WEP = age 3+, GOA [-140 LON westward] & SE separate, MT):
WEPpredationWEP = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), "TotBio_t" = c(
  1479000,
  1748000,
  1013000,
  737000,
  625000,
  1021000,
  713000,
  580000,
  1170000,
  1330000,
  1277000,
  1771000))

# Convert MT to g (comparable with units for annual ration and proportions of pollock consumed):
WEPpredationWEP$TotBio_g = WEPpredationWEP$TotBio_t * 1000000

#########################################################
# Input relative predator densities (unitless) predicted from GAMMs by survey year:
RelPredDensity = readRDS("2_RelativePredatorDensities/normWEPabun.rds")
colnames(RelPredDensity)[7] = c("Year") # rename to match other data frames 
RelPredDensity$Year = as.numeric(as.character(RelPredDensity$Year))
RelPredDensity = as.data.frame(RelPredDensity)

# Assign INPFC statistical areas based on grid cell id:
RelPredDensity$StatArea = with(RelPredDensity, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the main (i.e., 'GOA') Walleye Pollock stock assessment (areas east of -140):
RelPredDensity = as.data.frame(RelPredDensity)
RPD_poll = subset(RelPredDensity, StatArea != "Southeastern" & StatArea != "Other")

# Sum year-specific relative predator densities: 
require(dplyr)
RPD_sum = RPD_poll %>%
  group_by(Year) %>%
  summarise(NormAbun = sum(WEPnormAbun))

# Join to predator biomass data:
WEPpredationWEP = WEPpredationWEP %>% left_join(RPD_sum)

#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models:
Ration_WEP = readRDS("3_MeanAnnualRations/Rations_WEP.rds")

# Assign INPFC statistical areas based on grid cell id:
Ration_WEP$StatArea = with(Ration_WEP, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the main (i.e., 'GOA') Walleye Pollock stock assessment (areas east of -140):
Ration_WEPpoll = subset(Ration_WEP, StatArea != "Southeastern" & StatArea != "Other")

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
Ration_WEP_Cmax = Ration_WEPpoll %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy))

# Join to total predator biomass and relative predator density estimates:
WEPpredationWEP = WEPpredationWEP %>% left_join(Ration_WEP_Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless):
propWEP_WEP = readRDS("4_ProportionsPollockConsumed/propWEP_WEP_37.rds")
propWEP_WEP$Year = as.numeric(as.character(propWEP_WEP$Year))

# Assign INPFC statistical areas based on grid cell id:
propWEP_WEP$StatArea = with(propWEP_WEP, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the main (i.e., 'GOA') Walleye Pollock stock assessment (areas east of -140):
propWEP_WEPpoll = subset(propWEP_WEP, StatArea != "Southeastern" & StatArea != "Other")

# Calculate year-specific proportions of prey consumed:
propWEP_WEP_mean = propWEP_WEPpoll %>%
  group_by(Year, Prey_Name) %>% 
  summarise(sumWT = sum(WTdiet)) %>%
  mutate(propWT = sumWT / sum(sumWT))

# Select only proportions of Walleye Pollock:
propWEP_WEPsub = subset(propWEP_WEP_mean, Prey_Name == "Gadus.chalcogrammus")

# Assign zeros to survey years with no pollock consumed:
propWEP_WEPsub = add_row(as.data.frame(propWEP_WEPsub))
propWEP_WEPsub[nrow(propWEP_WEPsub), ] = c(2005, "Gadus.chalcogrammus", 0, 0)
propWEP_WEPsub = add_row(as.data.frame(propWEP_WEPsub))
propWEP_WEPsub[nrow(propWEP_WEPsub), ] = c(2011, "Gadus.chalcogrammus", 0, 0)
propWEP_WEPsub$Year = as.numeric(as.character(propWEP_WEPsub$Year))
propWEP_WEPsub$propWT = as.numeric(propWEP_WEPsub$propWT)

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
WEPpredWEP = WEPpredationWEP %>% left_join(propWEP_WEPsub)

#########################################################
# Calculate year-specific predation by Walleye Pollock conspecifics (g; all pollock age classes combined):
WEPpredWEP$predWEP_g_max = with(WEPpredWEP, TotBio_g * NormAbun * Annual_Cmax * propWT)

# Convert g of pollock consumed to metric tons:
WEPpredWEP$predWEP_t_max = round((WEPpredWEP$predWEP_g_max / 1000000), digits=0)

#########################################################
# Calculate age-specific predation on pollock by Walleye Pollock conspecifics (g):

# Input proportions of pollock consumed by age class:
propWEP_age_WEP = readRDS("5_AgeCompositionsPrey/propWEP_age_WEP.rds")

require(reshape2)
propWEP_Age_WEP = dcast(propWEP_age_WEP, Year ~ AgeClass, value.var = "propWT", fun.aggregrate = sum) # convert long to wide format
colnames(propWEP_Age_WEP) = c("Year", "prop_Age0", "prop_Age1")
propWEP_Age_WEP$Year = as.numeric(as.character(propWEP_Age_WEP$Year))

# Join to total predator biomasses, relative predator densities, mean annual rations, and proportions of pollock consumed:
WEPpredWEP_age = WEPpredWEP %>% left_join(propWEP_Age_WEP)

# Calculate age-specific consumption (no ages 2 or 3+):
WEPpredWEP_age$predWEP_t_max_Age0 = round(WEPpredWEP_age$predWEP_t_max * WEPpredWEP_age$prop_Age0, digits=0)
WEPpredWEP_age$predWEP_t_max_Age1 = round(WEPpredWEP_age$predWEP_t_max * WEPpredWEP_age$prop_Age1, digits=0)

#########################################################
# Convert dataframe back to long format for plotting:
WEPpredWEP_age_t_max = melt(WEPpredWEP_age, id.vars = c("Year"), measure.vars = c("predWEP_t_max","predWEP_t_max_Age0","predWEP_t_max_Age1"), variable.name = "AgeClass", value.name = "MT")
WEPpredWEP_age_t_max$millMT = round((WEPpredWEP_age_t_max$MT / 1000000), digits=2)
levels(WEPpredWEP_age_t_max$AgeClass) = c(" All Ages", " 0", " 1", " 2", " 3+")
WEPpredWEP_age_t_max$AgeClass = ordered(WEPpredWEP_age_t_max$AgeClass, levels = c(" All Ages", " 3+", " 2", " 1", " 0"))

require(ggplot2)
# Plot age-specific predation on Walleye Pollock by Walleye Pollock conspecifics:
WEPpredWEP_age_t_max_stack = subset(WEPpredWEP_age_t_max, AgeClass != " All Ages")
WEPpredWEP_age_t_max_stacked = ggplot(WEPpredWEP_age_t_max_stack, aes(x = Year, y = millMT, group = AgeClass, fill = AgeClass)) +
  geom_area(position="stack", show.legend=T, col="gray20", lwd=0.1) +
  scale_fill_manual(values = c("lightskyblue", "white")) +
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
  theme(panel.spacing.y = unit(1.0, "lines")) +
  theme(legend.position = c(0.95,0.90)) +
  theme(legend.text = element_text(colour="black", family="Arial", size=10.5), legend.text.align = 0) +
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(7.5, "mm")) +
  theme(legend.key.height = unit(3.5, "mm")) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=11)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=11)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=10)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12.5)) +
  xlab("") +
  ylab("Walleye Pollock Predation on Walleye Pollock (mill MT)") +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=c(0,2,4,6,8)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

WEPpredWEP_age_t_max_stacked