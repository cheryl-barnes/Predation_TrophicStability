# This script file includes code necessary to calculate indices of predation for Walleye Pollock in assessment area of the Gulf of Alaska (MT per year; 1990 to 2015). Pollock predators included Arrowtooth Flounder (Atheresthes stomias), Pacific Cod (Gadus macrocephalus), Pacific Halibut (Hippoglossus stenolepis), Sablefish (Anoplopoma fimbria), and Walleye Pollock (Gadus chalcogrammus) conspecifics. This script also includes estimates of synchrony (by way of variance ratio calculations) and portfolio effects to enable inferences about trophic stability among groundfishes in the study area.

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
source("6_PredationIndices/PredationIndices_ATF_PollAssess.R")
  ATFpredData = ATFpredWEP_age_t_max
  ATFpredData$Species = "ATF"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_PC_PollAssess.R")
  PCpredData = PCpredWEP_age_t_max
  PCpredData$Species = "PC"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_PH_PollAssess.R")
  PHpredData = PHpredWEP_age_t_max
  PHpredData$Species = "PH"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_SBL_PollAssess.R")
  SBLpredData = SBLpredWEP_age_t_max
  SBLpredData$Species = "SBL"
setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
source("6_PredationIndices/PredationIndices_WEP_PollAssess.R")
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
VarRatio$Area = as.factor("Assessment Area") # dummy variable for plotting
SyncPlot = ggplot(VarRatio, aes(x=Year, y=Sync, group=Area, col=Area, fill=Area)) +
  geom_point(shape=1, size=2) +
  geom_line() +
  scale_color_manual(values = "blue") +
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

##################################################################
# Plot age 3+ consumption of pollock against age-3+ pollock biomass estimates from the GOA stock assessment (Dorn et al. 2017):

# Sum consumption of age-3+ pollock (all predators):
PredData_3plus = subset(PredData, AgeClass == " 3+")
PredWEP_3t = PredData_3plus %>%
  group_by(Year) %>%
  summarise(cons_millMT = sum(millMT))
PredWEP_3t = as.data.frame(PredWEP_3t)

# Add column with age-3+ pollock biomass (from assessment):
PredWEP_3t$WEPbio = c(
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
  1771000)

# Convert to pollock biomass to million MT:
PredWEP_3t$WEP_millMT = PredWEP_3t$WEPbio / 1000000

# Calculate ratio of pollock consumption to pollock biomass:
PredWEP_3t$PredBioRatio = round(PredWEP_3t$cons_millMT / PredWEP_3t$WEP_millMT, digits=2)

# Test for correlation between pollock consumption and pollock biomass:
  # Full time series (1990 to 2015):
cor.test(PredWEP_3t$cons_millMT, PredWEP_3t$WEP_millMT, method ="pearson")
  # Early portion of the time series (1990 to 2003):
corr1990 = subset(PredWEP_3t, Year < 2005)
cor.test(corr1990$cons_millMT, corr1990$WEP_millMT, method ="pearson")
  # Late portion of the time series (2005 to 2015):
corr2005 = subset(PredWEP_3t, Year > 2003)
cor.test(corr2005$cons_millMT, corr2005$WEP_millMT, method ="pearson")

# Plot results (Figure l6):
WEP_3plusRatio = ggplot(PredWEP_3t, aes(x = Year)) +
  geom_line(aes(y=cons_millMT), linetype="solid", lwd=1.25, col="mediumblue", show.legend = T) +
  geom_line(aes(y=WEP_millMT), col="black", linetype="12345678", lwd=0.5, show.legend = T) +
  geom_point(aes(y=WEP_millMT), shape=5, size=1.5, stroke=0.75, col="black", fill="black", show.legend = T) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(plot.background = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.position = c(0.87,0.90)) +
  theme(legend.text = element_text(colour="black", family="Arial", size=10.5), legend.text.align = 0) +
  theme(legend.key.width = unit(12.75, "mm")) +
  theme(legend.key.height = unit(4.50, "mm")) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text = element_text(family="Arial", color="black", size=11)) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  xlab("") +
  ylab("Predation on age-3+ pollock (mill MT)") +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3), breaks=c(0,1,2,3)) +
  theme(legend.background = element_rect(fill="transparent")) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

WEP_3plusRatio