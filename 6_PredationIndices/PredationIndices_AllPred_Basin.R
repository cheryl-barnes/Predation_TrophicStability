# This script file includes code necessary to calculate indices of predation for Walleye Pollock in the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predators included Arrowtooth Flounder (Atheresthes stomias), Pacific Cod (Gadus macrocephalus), Pacific Halibut (Hippoglossus stenolepis), Sablefish (Anoplopoma fimbria), and Walleye Pollock (Gadus chalcogrammus) conspecifics. This script also includes estimates of synchrony (by way of variance ratio calculations) and portfolio effects to enable inferences about trophic stability among groundfishes within the study area.

# When using any portion of the code herein, please cite: Barnes, C. L., A. H. Beaudreau, M. W. Dorn, K. K. Holsman, and F. J. Mueter. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecological Applications. 

# Data Sources: All the data necessary to complete the following analyses can be found in the 'Data' folder on github. 
  # Total biomass estimates were obtained from the most recent stock assessment for each groundfish predator (Barbeaux et al. 2017, Dorn et al. 2017, Hanselman et al. 2017, Spies et al. 2017, Stewart and Hicks 2017). 
  # Bottom trawl survey data (all groundfish predators; 1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm (for methods, see von Szalay et al. 2016). Setline survey data (Pacific Halibut; 1998 to 2017) were collected by the International Pacific Halibut Commission and are publicly available at: https://iphc.int/data/fiss-data-query (for methods, see Clark and Hare 2006). Longline survey data (Sablefish; 1990 to 2017) were collected by personnel at the Alaska Fisheries Science Center's Auke Bay Laboratories and can be found at https://www.afsc.noaa.gov/maps/longline/Map.php (for methods, see Sigler and Zenger 1989). 
  # Food habits data (all groundfish predators; 1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 
  
# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Barbeaux, S., K. Aydin, B. Fissel, K. Holsman, and W. Palsson. 2017. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 189–332.
# Clark, W. G., and S. R. Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. International Pacific Halibut Commission Scientific Report 83. 
# Dorn, M., K. Aydin, B. Fissel, D. Jones, A. McCarthy, W. Palsson, and K. Spalinger. 2017. Assessment of the Walleye Pollock stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 47–182.
# Hanselman, D. H., C. J. Rodgveller, C. R. Lunsford, and K. H. Fenske. 2017. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report 327–502.
# Livingston, P. A., K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.
# Spies, I., K. Aydin, J. N. Ianelli, and W. Palsson. 2017. Assessment of the Arrowtooth Flounder stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report 749–846.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. NOAA Technical Memorandum NMFS-AFSC-325. 

##################################################################
# Read in predator-specific estimates of pollock consumption:
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_ATF_Basin.R")
  ATFpredData = ATFpredWEP_age_t_max
  ATFpredData$Species = "ATF"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_PC_Basin.R")
  PCpredData = PCpredWEP_age_t_max
  PCpredData$Species = "PC"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_PH_Basin.R")
  PHpredData = PHpredWEP_age_t_max
  PHpredData$Species = "PH"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_SBL_Basin.R")
  SBLpredData = SBLpredWEP_age_t_max
  SBLpredData$Species = "SBL"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_WEP_Basin.R")
  WEPpredData = WEPpredWEP_age_t_max
  WEPpredData$Species = "WEP"

PredData = rbind(ATFpredData, PCpredData, PHpredData, SBLpredData, WEPpredData)
PredData_all = subset(PredData, AgeClass == " All Ages")

##################################################################
# Estimate relative contributions of each predator to pollock predation mortality in the Gulf of Alaska (1990 to 2015; Table 3):
PredData_all = PredData_all %>%
  group_by(Year) %>%
  mutate(SumMT = sum(MT))

# Results in Table 3:
PredData_all %>%
  group_by(Species) %>%
  summarise(mean(MT / SumMT))
PredData_all %>%
  group_by(Species) %>%
  summarise(sd(MT / SumMT))

##################################################################
# Plot pollock consumption by predator (all pollock age classes combined):
PredWEP_t = PredData_all %>%
  group_by(Year, Species) %>%
  summarise(sumMillMT = sum(millMT))
PredWEP_t = as.data.frame(PredWEP_t)

PredWEP_t$Species = ordered(as.factor(PredWEP_t$Species), levels = c("ATF", "PC", "PH", "SBL", "WEP"))
levels(PredWEP_t$Species) = c(" Arrowtooth Flounder", " Pacific Cod", " Pacific Halibut", " Sablefish", " Walleye Pollock")

# Figure 3 (top):
range(PredWEP_t$sumMillMT)
PredWEP_t_stacked = ggplot(PredWEP_t, aes(x = Year, y = sumMillMT, group = Species, fill = Species)) +
  geom_area(position="stack", show.legend=T, col="gray50", lwd=0.1) +
  scale_fill_manual(values = c("deepskyblue4", "lightseagreen", "olivedrab3", "orange", "yellow")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(panel.spacing.y = unit(1.0, "lines")) +
  theme(legend.position = c(0.85,0.87)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(legend.text = element_text(colour="black", family="Arial", size=11), legend.text.align = 0) +
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(7.5, "mm")) +
  theme(legend.key.height = unit(3.5, "mm")) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks = element_line(color="black")) +
  theme(axis.text = element_text(family="Arial", color="black", size=11)) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Predator Species") +
  scale_x_continuous(limits = c(1990,2015), breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=c(0,2,4,6,8)) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

PredWEP_t_stacked

##################################################################
# Plot total pollock consumption (all predators combined) by pollock age class:
SumPredWEP_t = subset(PredData, AgeClass != " All Ages")
SumPredWEP_t = SumPredWEP_t %>%
  group_by(Year, AgeClass) %>%
  summarise(sumMillMT = sum(millMT))
SumPredWEP_t = as.data.frame(SumPredWEP_t)

# Figure 3 (bottom):
SumPredWEP_t_plot_stacked = ggplot(SumPredWEP_t, aes(x = Year, y = sumMillMT, group = AgeClass, fill = AgeClass)) +
  geom_area(position="stack", show.legend=T, col="gray20", lwd=0.1) +
  scale_fill_manual(values = c("mediumblue", "dodgerblue1", "lightskyblue1", "azure")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(panel.spacing.y = unit(1.0, "lines")) +
  theme(legend.position = c(0.92,0.87)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(legend.text = element_text(colour="black", family="Arial", size=11), legend.text.align = 0) +
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(7.5, "mm")) +
  theme(legend.key.height = unit(3.5, "mm")) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text = element_text(family="Arial", color="black", size=11)) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Age Class (yr)") +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=c(0,2,4,6,8)) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

SumPredWEP_t_plot_stacked

##################################################################
# Calculate correlations in consumption among predators:
PredData_wide = dcast(PredData_all, Year ~ Species, value.var = "MT", fun.aggregrate = sum) # convert long to wide format

require(Hmisc); require(corrplot)
CorrMat = rcorr(as.matrix(PredData_wide[2:6]), type="pearson")

# Results in Table 4:
corrplot(CorrMat$r, type="upper", order="hclust", diag=F, tl.pos="td", tl.cex=1, tl.col="black", p.mat = CorrMat$P, sig.level = 0.1, insig = "blank", method="number")

##################################################################
# Calculate variance ratios (sum of preadator-specific variances in consumption / variance in total pollock consumption):
require(zoo)

# Calculate variances in pollock consumption for each predator using a 5-year moving window:
VarSpp = PredData_all %>% 
  group_by(Species) %>%
  mutate(VarEst = rollapplyr(MT, width=5, FUN=var, align="left", partial=T))

# Sum predator-specific variance estimates for each window:
VarSpp_summ = VarSpp %>%
  group_by(Year) %>%
  summarise(VarSpp = sum(VarEst))

# Calculate variances in total pollock consumption (sum of consumption by all predators) using the same 5-year moving window:
SumPred_MT = PredData_all %>%
  group_by(Year) %>%
  summarise(sumMT = sum(MT))

VarTot = SumPred_MT %>%
  mutate(VarTot = rollapplyr(sumMT, width=5, FUN=var, align="left", partial=T))

# Combine data frames and calculate rolling variance ratios (VR; i.e., synchrony) and portfolio effects (PE; 1 - VR):
VarRatio = VarSpp_summ %>% left_join(VarTot)
  VarRatio = subset(VarRatio, Year <= 2007) # remove entries with fewer than five survey years in estimate

VarRatio = VarRatio %>%
  mutate(Sync = VarTot / VarSpp) %>%
  mutate(PE = 1 - Sync)

# Plot variance ratios (> 1: synchronous [neg. PE], 1: indep., < 1: asynchronous [pos. PE]); Figure 7 (top):
VarRatio$Area = as.factor("Basin") # dummy variable for plotting
SyncPlot = ggplot(VarRatio, aes(x=Year, y=Sync, group=Area, col=Area, fill=Area)) +
  geom_point(shape=16, size=2) +
  geom_line() +
  scale_color_manual(values = "black") +
  geom_hline(yintercept=1, linetype="dashed") +
  theme(panel.background = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(family="Arial", size=11), legend.text.align = 0) +
  theme(legend.position = c(0.742,0.95)) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.spacing.x = unit(0.05, 'cm')) +
  theme(legend.key = element_blank()) +
  theme(legend.key.height = unit(0.5, 'cm')) +
  theme(axis.text = element_text(family="Arial", color="black", size=11)) +
  theme(axis.title = element_text(family="Arial", color="black", hjust=0.46, vjust=1.5, size=12)) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  labs(y="Variance Ratio", x="Start of 5-yr Moving Window") +
  scale_x_continuous(limits = c(1990,2015), breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015)) +
  scale_y_continuous(limits = c(0.5,1.763), breaks = c(0.5, 1.0, 1.5)) 

SyncPlot