# This script file includes code necessary to an calculate index of predation for Walleye Pollock within the assessment are of the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predator: Pacific Cod (Gadus macrocephalus). 

# When using any portion of the code herein, please cite: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

# Data Sources: All the data necessary to complete the following analyses can be found in the 'Data' folder on github. 
  # Total biomass estimates were obtained from the most recent stock assessment for Pacific Cod (Barbeaux et al. 2017). 
  # Bottom trawl survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm (for methods, see von Szalay et al. 2016). 
  # Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 
  
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Barbeaux, S., K. Aydin, B. Fissel, K. Holsman, and W. Palsson. 2017. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 189–332.
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
##################################################################
# Input total biomass estimates from most recent stock assessment (PC = age 0+, Gulf-wide, MT):
PCpredationWEP = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), "TotBio_t" = c(
  583841,
  516782,
  429292,
  320235,
  286165,
  292752,
  247481,
  246629,
  307285,
  345269,
  316926,
  312414)) # total biomass estimates provided by email [Barbeaux Model 17.09.35_3.24U] 

# Convert MT to g (comparable with units for annual ration and proportions of pollock consumed):
PCpredationWEP$TotBio_g = PCpredationWEP$TotBio_t * 1000000

#########################################################
# Input relative predator densities (unitless) predicted from GAMMs by survey year:
RelPredDensity = readRDS("2_RelativePredatorDensities/normPCabun.rds")
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
  summarise(NormAbun = sum(PCnormAbun))

# Join to predator biomass data:
PCpredationWEP = PCpredationWEP %>% left_join(RPD_sum)
#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models:
Ration_PC = readRDS("3_MeanAnnualRations/Rations_PC.rds")

# Assign INPFC statistical areas based on grid cell id:
Ration_PC$StatArea = with(Ration_PC, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the main (i.e., 'GOA') Walleye Pollock stock assessment (areas east of -140):
Ration_PCpoll = subset(Ration_PC, StatArea != "Southeastern" & StatArea != "Other")

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
Ration_PC_Cmax = Ration_PCpoll %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy))

# Join to total predator biomass and relative predator density estimates:
PCpredationWEP = PCpredationWEP %>% left_join(Ration_PC_Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless):
propWEP_PC = readRDS("4_ProportionsPollockConsumed/propWEP_PC_0.rds")
propWEP_PC$Year = as.numeric(as.character(propWEP_PC$Year))

# Assign INPFC statistical areas based on grid cell id:
propWEP_PC$StatArea = with(propWEP_PC, 
ifelse(id2 %in% c(123,150,151,152,153,154,155,156,178,179,180,181,182,183,184,185,186,187,188,189,190,206,20,208,209,210,211,212,213,214,215,236,237,238,239,240,241,242,268,269,270,271,272,299), "Shumagin", 
ifelse(id2 %in% c(216,217,218,243,244,245,246,247,273,274,275,276,277,278,279,300,301,302,303,304,305,306,333,334,335,336,337,338,369,370,371,372,403,404,405,406,437,438,439,470), "Chirikof", 
ifelse(id2 %in% c(280,307,308,339,340,341,373,374,375,376,407,408,409,410,411,440,441,442,443,444,445,446,471,472,473,474,475,476,477,499,500,501,502,503,504,505,525,526,527,528,529,530,531,532,533,551,552,553,554,555,556,557,558,578,579,580,582), "Kodiak", 
ifelse(id2 %in% c(534,540,541,559,560,561,562,563,564,565,566,583,584,585,586,587,588,589,590,604,605,608,609), "Yakutat", 
ifelse(id2 %in% c(327,328,329,361,362,363,364,395,396,397,429,430,461,462,491,492,493,515,516,517,518,519,542,543,544,545,567,568,569,591,611), "Southeastern", "Other"))))))

# Eliminate locations outside of the main (i.e., 'GOA') Walleye Pollock stock assessment (areas east of -140):
propWEP_PCpoll = subset(propWEP_PC, StatArea != "Southeastern" & StatArea != "Other")

# Calculate year-specific proportions of prey consumed:
propWEP_PC_mean = propWEP_PCpoll %>%
  group_by(Year, Prey_Name) %>% 
  summarise(sumWT = sum(WTdiet)) %>%
  mutate(propWT = sumWT / sum(sumWT))

# Select only proportions of Walleye Pollock:
propWEP_PCsub = subset(propWEP_PC_mean, Prey_Name == "Gadus.chalcogrammus")

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
PCpredWEP = PCpredationWEP %>% left_join(propWEP_PCsub)

#########################################################
# Calculate year-specific predation by Pacific Cod (g; all pollock age classes combined):
PCpredWEP$predWEP_g_max = with(PCpredWEP, TotBio_g * NormAbun * Annual_Cmax * propWT)

# Convert g of pollock consumed to metric tons:
PCpredWEP$predWEP_t_max = round((PCpredWEP$predWEP_g_max / 1000000), digits=0)

#########################################################
# Calculate age-specific predation on pollock by Pacific Cod (g):

# Input proportions of pollock consumed by age class:
propWEP_age_PC = readRDS("5_AgeCompositionsPrey/propWEP_age_PC.rds")

require(reshape2)
propWEP_Age_PC = dcast(propWEP_age_PC, Year ~ AgeClass, value.var = "propWT", fun.aggregrate = sum) # convert long to wide format
colnames(propWEP_Age_PC) = c("Year", "prop_Age0", "prop_Age1", "prop_Age2", "prop_Age3plus")
propWEP_Age_PC$Year = as.numeric(as.character(propWEP_Age_PC$Year))

# Join to total predator biomasses, relative predator densities, mean annual rations, and proportions of pollock consumed:
PCpredWEP_age = PCpredWEP %>% left_join(propWEP_Age_PC)

# Calculate age-specific consumption:
PCpredWEP_age$predWEP_t_max_Age0 = round(PCpredWEP_age$predWEP_t_max * PCpredWEP_age$prop_Age0, digits=0)
PCpredWEP_age$predWEP_t_max_Age1 = round(PCpredWEP_age$predWEP_t_max * PCpredWEP_age$prop_Age1, digits=0)
PCpredWEP_age$predWEP_t_max_Age2 = round(PCpredWEP_age$predWEP_t_max * PCpredWEP_age$prop_Age2, digits=0)
PCpredWEP_age$predWEP_t_max_Age3plus = round(PCpredWEP_age$predWEP_t_max * PCpredWEP_age$prop_Age3plus, digits=0)

#########################################################
# Convert dataframe back to long format for plotting:
PCpredWEP_age_t_max = melt(PCpredWEP_age, id.vars = c("Year"), measure.vars = c("predWEP_t_max","predWEP_t_max_Age0","predWEP_t_max_Age1","predWEP_t_max_Age2","predWEP_t_max_Age3plus"), variable.name = "AgeClass", value.name = "MT")
PCpredWEP_age_t_max$millMT = round((PCpredWEP_age_t_max$MT / 1000000), digits=2)
levels(PCpredWEP_age_t_max$AgeClass) = c(" All Ages", " 0", " 1", " 2", " 3+")
PCpredWEP_age_t_max$AgeClass = ordered(PCpredWEP_age_t_max$AgeClass, levels = c(" All Ages", " 3+", " 2", " 1", " 0"))

require(ggplot2)
# Plot age-specific predation on Walleye Pollock by Pacific Cod:
PCpredWEP_age_t_max_stack = subset(PCpredWEP_age_t_max, AgeClass != " All Ages")
PCpredWEP_age_t_max_stacked = ggplot(PCpredWEP_age_t_max_stack, aes(x = Year, y = millMT, group = AgeClass, fill = AgeClass)) +
  geom_area(position="stack", show.legend=T, col="gray20", lwd=0.1) +
  scale_fill_manual(values = c("chartreuse4", "chartreuse2", "olivedrab1", "white")) +
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
  ylab("Pacific Cod Predation on Walleye Pollock (mill MT)") +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=c(0,2,4,6,8)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

PCpredWEP_age_t_max_stacked