# This script file includes code necessary to identify age compositions of pollock consumed by Arrowtooth Flounder, Pacific Cod, Pacific Halibut, Sablefish, and Walleye Pollock conspecifics in the Gulf of Alaska (1990 to 2015). First, we used survey data to estimate von Bertalanffy (1938) growth parameters to predict pollock age from length. We then quantified bias-corrected length-weight relationships (Brodziak 2012) for all pollock measured from stomach contents of our focal predators. We used multinomial logistic regression to estimates year-specific age compositions of pollock consumed. Small sample sizes precluded spatially-explicit age compositions.

# Bottom trawl survey data were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# Food habits data were obtained from the Resource Ecology and Ecosystem Modeling (REEM) program of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See Livingston et al. (2017) for methods details. 

# References:
# Brodziak, J. 2012. Fitting length-weight relationships with linear regression using the log-transformed allometric model with bias-correction. NOAA Technical Memorandum PIFSC-H-12-03. 
# Livingston, P. A. , K. Aydin, T. W. Buckley, G. M. Lang, M-S. Yang, and B. S. Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environmental Biology of Fishes 100(4):443–470.
# von Bertalanffy, L. 1938. A quantitative theory of organic growth. Human Biology 10:181–213.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Technical Memorandum NMFS-AFSC-325. 

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

setwd("~/Documents/UAF/Dissertation/GitHub/Predation_TrophicStability/")
#########################################################
### QUANTIFY LENGTH-AGE RELATIONSHIPS, Walleye Pollock
#########################################################
# Input pollock length, weight, and age data from the AFSC bottom trawl survey:
WEPdata = read.csv("1_Data/WEPagesSPC.csv")

# Manually remove obvious outliers (i.e., age 0 fish > 400 mm):
WEPdata = WEPdata[-c(1679, 1759, 1861, 1871),]

# Set starting values for VBGF (based on raw data):
theta = c(650, 0.2, 0) 

# Estimate von Bertalanffy growth parameters:
SSQ = function(theta, x) {
  Linf = theta[1]
  K = theta[2]
  t0 = theta[3]
  epsilon = rep(0, length(WEPdata$AGE))
  lpred = rep(0, length(WEPdata$AGE))
  for (i in 1:length(WEPdata$AGE)) {
    lpred[i] = Linf * (1 - exp(-K * (WEPdata$AGE[i] - t0)))
    epsilon[i] = (WEPdata$LENGTH[i] - lpred[i])^2
  }
  ssq = sum(epsilon)
  return(ssq)
}

# Solve using least squares:
out = optim(theta, fn = SSQ, method = "BFGS", x = WEPdata$AGE, hessian = TRUE)
out$V = solve(out$hessian)  # solve the hessian
out$SE = sqrt(diag(out$V))  # estimate SE
out$R = out$V/(out$SE %o% out$SE)  # Correlation
out$par; out$SE

# Estimate pollock length from age:
WEPlengths = read.csv("1_Data/WEPSizeComps.csv")
WEPlengths = subset(WEPlengths, Year >= 1990)

Linf = out$par[1]
K = out$par[2]
t0 = out$par[3]
WEPlengths$predAge = ((-(log(-(WEPlengths$Length..mm./Linf)+1)))/K)+t0

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
WEPlengths$AgeClass = with(WEPlengths, 
    ifelse(predAge < 1, "0", 
    ifelse(predAge < 2, "1", 
    ifelse(predAge < 3, "2", "3+")))) 

#########################################################
### QUANTIFY LENGTH-WEIGHT RELATIONSHIPS, Walleye Pollock
#########################################################
# Allometric model with bias-correction (Brodziak 2012):
WEPdata$logW = log(WEPdata$WEIGHT)
WEPdata$logL = log(WEPdata$LENGTH)
WEPdata = subset(WEPdata, logW !="NA")

L_W = lm(logW ~ logL, data=WEPdata)
summary(L_W)
a = L_W$coefficients[1]
b = L_W$coefficients[2]

# Calculate correction factor for predicting weight on original scale:
syx = summary(L_W)$sigma
cf = exp((syx^2)/2) 
WEPdata$predWlog = predict(L_W, data.frame(logL=log(WEPdata$LENGTH)), interval="c") 
WEPdata$predW_unbi = cf *(exp(WEPdata$predWlog))

# Estimate pollock weight from length:
WEPlengths$logL = log(WEPlengths$Length..mm.)
WEPlengths$predWlog = predict(L_W, data.frame(logL=log(WEPlengths$Length..mm.)), interval="c") 
WEPlengths$predW_unbi = cf *(exp(WEPlengths$predWlog)) 

#########################################################
### CALCULATE PROPORTIONS OF POLLOCK CONSUMED, BY AGE CLASS 
#########################################################
# Predator: Arrowtooth Flounder #
#########################################################
# Input food habits data for Arrowtooth Flounder:
WEPpreyL_ATF = read.csv("1_Data/GOA_RawPL_ATFdiets.csv", check.names=F)

# Remove fish not included in assessment-based estimates of total predator biomass (ATF < 19 cm):
WEPpreyL_ATF = subset(WEPpreyL_ATF, Pred_len >= 19)

# Select pollock as prey:
WEPpreyL_ATF = subset(WEPpreyL_ATF, Prey_Name=="Walleye.pollock.Gadus.chalcogrammus")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
WEPpreyL_ATF$Year = as.factor(WEPpreyL_ATF$Year)
WEPpreyL_ATF = subset(WEPpreyL_ATF, Year != "1981" & Year != "1984" & Year != "1987")

# Predict pollock age from length:
WEPpreyL_ATF$predAge = (log(-(WEPpreyL_ATF$Prey_sz1/Linf)+1))/-K

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
WEPpreyL_ATF$AgeClass = with(WEPpreyL_ATF, 
  ifelse(predAge < 1, "0", 
     ifelse(predAge < 2, "1", 
         ifelse(predAge < 3, "2", "3+")))) 
WEPpreyL_ATF$AgeClass = as.factor(WEPpreyL_ATF$AgeClass)

# Predict pollock weight from length:
WEPpreyL_ATF$predWlog = predict(L_W, data.frame(logL=log(WEPpreyL_ATF$Prey_sz1)), interval="c") 
WEPpreyL_ATF$predWlog = as.vector(WEPpreyL_ATF$predWlog[,1])
WEPpreyL_ATF$predW_unbi = cf *(exp(WEPpreyL_ATF$predWlog))

require(reshape2)
require(VGAM) 
require(ggplot2)

# Estimate age compositions of pollock consumed, by survey year:
WEP_age_ATFdata = dcast(WEPpreyL_ATF, Year ~ AgeClass, value.var = "predW_unbi", sum)
AgeMatrix_ATF = as.matrix(round(WEP_age_ATFdata[,2:5], digits=0)) # convert long to wide format
rownames(AgeMatrix_ATF) = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015) 
WEP_age_ATFdata$Year = as.numeric(as.character(WEP_age_ATFdata$Year))

ATFfit_gam = vgam(AgeMatrix_ATF ~ s(Year), multinomial, data=WEP_age_ATFdata)
  summary(ATFfit_gam) 

predProp_ATF_gam = as.data.frame(predict(ATFfit_gam, type="response"))
predProp_ATF_gam$Year = c("1990", "1993", "1996", "1999", "2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015")

ATF_multi_gam = melt(predProp_ATF_gam, id.vars = c("Year"), measure.vars = 1:4, variable.name = "AgeClass", value.name = "propWT")
ATF_multi_gam$Predator = "ATF"

# Save results:
saveRDS(ATF_multi_gam, "5_AgeCompositionsPrey/propWEP_age_ATF.rds")

# Plot (portion of Fig. 2):
ATF_multi_gam$AgeClass = ordered(ATF_multi_gam$AgeClass, levels = c("3+", "2", "1", "0"))
WEP_Age_ATF = ggplot(ATF_multi_gam, aes(x = as.numeric(as.character(Year)), y = propWT, fill = AgeClass, order = -as.numeric(AgeClass))) + 
  geom_bar(position = "fill", stat = "identity",color="black") +
  scale_fill_manual(values = c("orange", "yellow", "lemonchiffon", "white"), name=" Age", drop=FALSE) +
  theme_bw() +
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
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "vertical", legend.position = "right", legend.margin = margin(0,0,0,5), legend.background = element_rect(fill="transparent")) +
  theme(legend.position = c(0.0385, 0.865)) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(1.0, "mm")) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(legend.title = element_text(family="Arial", color="black", size=11)) +
  theme(legend.text = element_text(family="Arial", color="black", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=12)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  xlab("") +
  ylab("Proportion of Pollock Consumed") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) +
  scale_x_continuous(breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01))

WEP_Age_ATF

#########################################################
# Predator: Pacific Cod #
#########################################################
# Input food habits data for Pacific Cod:
WEPpreyL_PC = read.csv("1_Data/GOA_RawPL_PCdiets.csv", check.names=F)

# Remove fish not included in assessment-based estimates of total predator biomass (PC < 0 cm):
WEPpreyL_PC = subset(WEPpreyL_PC, Pred_len >= 0)

# Select pollock as prey:
WEPpreyL_PC = subset(WEPpreyL_PC, Prey_Name=="Walleye pollock Gadus chalcogrammus")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
WEPpreyL_PC$Year = as.factor(WEPpreyL_PC$Year)
WEPpreyL_PC = subset(WEPpreyL_PC, Year != "1981" & Year != "1984" & Year != "1987")

# Predict age from length:
WEPpreyL_PC$predAge = (log(-(WEPpreyL_PC$Prey_sz1/Linf)+1))/-K

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
WEPpreyL_PC$AgeClass = with(WEPpreyL_PC, 
    ifelse(predAge < 1, "0", 
      ifelse(predAge < 2, "1", 
        ifelse(predAge < 3, "2", "3+"))))

# Predict pollock weight from length:
WEPpreyL_PC$predWlog = predict(L_W, data.frame(logL=log(WEPpreyL_PC$Prey_sz1)), interval="c") 
WEPpreyL_PC$predWlog = as.vector(WEPpreyL_PC$predWlog[,1])
WEPpreyL_PC$predW_unbi = cf *(exp(WEPpreyL_PC$predWlog))

# Estimate age compositions of pollock consumed, by survey year:
WEP_age_PCdata = dcast(WEPpreyL_PC, Year ~ AgeClass, value.var = "predW_unbi", sum)
AgeMatrix_PC = as.matrix(round(WEP_age_PCdata[,2:5], digits=0))
rownames(AgeMatrix_PC) = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015) 
WEP_age_PCdata$Year = as.numeric(as.character(WEP_age_PCdata$Year))

PCfit_gam = vgam(AgeMatrix_PC ~ s(Year), multinomial, data=WEP_age_PCdata)
  summary(PCfit_gam) 

predProp_PC_gam = as.data.frame(predict(PCfit_gam, type="response"))
predProp_PC_gam$Year = c("1990", "1993", "1996", "1999", "2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015")

PC_multi_gam = melt(predProp_PC_gam, id.vars = c("Year"), measure.vars = 1:4, variable.name = "AgeClass", value.name = "propWT")
PC_multi_gam$Predator = "PC"

# Save results:
saveRDS(PC_multi_gam, "5_AgeCompositionsPrey/propWEP_age_PC.rds")

# Plot (portion of Fig. 2):
PC_multi_gam$AgeClass = ordered(PC_multi_gam$AgeClass, levels = c("3+", "2", "1", "0"))
WEP_Age_PC = ggplot(PC_multi_gam, aes(x = as.numeric(as.character(Year)), y = propWT, fill = AgeClass, order = -as.numeric(AgeClass))) + 
  geom_bar(position = "fill", stat = "identity",color="black") +
  scale_fill_manual(values = c("chartreuse4", "chartreuse2", "olivedrab1", "white"), name=" Age", drop=FALSE) +
  theme_bw() +
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
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "vertical", legend.position = "right", legend.margin = margin(0,0,0,5), legend.background = element_rect(fill="transparent")) +
  theme(legend.position = c(0.0385, 0.865)) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(1.0, "mm")) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(legend.title = element_text(family="Arial", color="black", size=11)) +
  theme(legend.text = element_text(family="Arial", color="black", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=12)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  xlab("") +
  ylab("Proportion of Pollock Consumed") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) +
  scale_x_continuous(breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01))

WEP_Age_PC

#########################################################
# Predator: Pacific Halibut #
#########################################################
# Input food habits data for Pacific Halibut:
WEPpreyL_PH = read.csv("1_Data/GOA_RawPL_PHdiets.csv", check.names=F)

# Remove fish not included in assessment-based estimates of total predator biomass (PH < 82 cm):
WEPpreyL_PH = subset(WEPpreyL_PH, Pred_len >= 82)

# Select pollock as prey:
WEPpreyL_PH = subset(WEPpreyL_PH, .Prey_Name=="Walleye.pollock.Gadus.chalcogrammus")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
WEPpreyL_PH$Year = as.factor(WEPpreyL_PH$Year)
WEPpreyL_PH = subset(WEPpreyL_PH, Year != "1981" & Year != "1984" & Year != "1987")

# Predict pollock age from length:
WEPpreyL_PH$predAge = (log(-(WEPpreyL_PH$Prey_sz1/Linf)+1))/-K

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
WEPpreyL_PH$AgeClass = with(WEPpreyL_PH, 
  ifelse(predAge < 1, "0", 
    ifelse(predAge < 2, "1", 
      ifelse(predAge < 3, "2", "3+"))))

# Predict pollock weight from length:
WEPpreyL_PH$predWlog = predict(L_W, data.frame(logL=log(WEPpreyL_PH$Prey_sz1)), interval="c") 
WEPpreyL_PH$predWlog = as.vector(WEPpreyL_PH$predWlog[,1])
WEPpreyL_PH$predW_unbi = cf *(exp(WEPpreyL_PH$predWlog))

# Estimate age compositions of pollock consumed, by survey year:
WEP_age_PHdata = dcast(WEPpreyL_PH, Year ~ AgeClass, value.var = "predW_unbi", sum)
AgeMatrix_PH = as.matrix(round(WEP_age_PHdata[,2:5], digits=0))
rownames(AgeMatrix_PH) = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015) 
WEP_age_PHdata$Year = as.numeric(as.character(WEP_age_PHdata$Year))

PHfit_gam = vgam(AgeMatrix_PH ~ s(Year), multinomial, data=WEP_age_PHdata)
  summary(PHfit_gam)   
  
predProp_PH_gam = as.data.frame(predict(PHfit_gam, type="response"))
predProp_PH_gam$Year = c("1990", "1993", "1996", "1999", "2001", "2003", "2005", "2007", "2009", "2011", "2013", "2015")

PH_multi_gam = melt(predProp_PH_gam, id.vars = c("Year"), measure.vars = 1:4, variable.name = "AgeClass", value.name = "propWT")
PH_multi_gam$Predator = "PH"

# Save results:
saveRDS(PH_multi_gam, "5_AgeCompositionsPrey/propWEP_age_PH.rds")

# Plot (portion of Fig. 2):
PH_multi_gam$AgeClass = ordered(PH_multi_gam$AgeClass, levels = c("3+", "2", "1", "0"))
WEP_Age_PH = ggplot(PH_multi_gam, aes(x = as.numeric(as.character(Year)), y = propWT, fill = AgeClass, order = -as.numeric(AgeClass))) + 
  geom_bar(position = "fill", stat = "identity",color="black") +
  scale_fill_manual(values = c("indianred2", "lightpink2", "thistle1", "white"), name=" Age", drop=FALSE) +
  theme_bw() +
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
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "vertical", legend.position = "right", legend.margin = margin(0,0,0,5), legend.background = element_rect(fill="transparent")) +
  theme(legend.position = c(0.0385, 0.865)) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(1.0, "mm")) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(legend.title = element_text(family="Arial", color="black", size=11)) +
  theme(legend.text = element_text(family="Arial", color="black", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=12)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  xlab("") +
  ylab("Proportion of Pollock Consumed") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) +
  scale_x_continuous(breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01))

WEP_Age_PH

#########################################################
# Predator: Sablefish #
#########################################################
# Input food habits data for Sablefish:
WEPpreyL_SBL = read.csv("1_Data/GOA_RawPL_SBLdiets.csv", check.names=F)

# Remove fish not included in assessment-based estimates of total predator biomass (SBL < 45 cm):
WEPpreyL_SBL = subset(WEPpreyL_SBL, Pred_len >= 45)

# Select pollock as prey:
WEPpreyL_SBL = subset(WEPpreyL_SBL, Prey_Name=="Walleye pollock Gadus chalcogrammus")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
WEPpreyL_SBL$Year = as.factor(WEPpreyL_SBL$Year)
WEPpreyL_SBL = subset(WEPpreyL_SBL, Year != "1981" & Year != "1984" & Year != "1987")

# Predict pollock age from length:
WEPpreyL_SBL$predAge = (log(-(WEPpreyL_SBL$Prey_sz1/Linf)+1))/-K

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
WEPpreyL_SBL$AgeClass = with(WEPpreyL_SBL, 
    ifelse(predAge < 1, "0", 
        ifelse(predAge < 2, "1", 
            ifelse(predAge < 3, "2", "3+"))))

# Predict pollock weight from length:
WEPpreyL_SBL$predWlog = predict(L_W, data.frame(logL=log(WEPpreyL_SBL$Prey_sz1)), interval="c") 
WEPpreyL_SBL$predWlog = as.vector(WEPpreyL_SBL$predWlog[,1])
WEPpreyL_SBL$predW_unbi = cf *(exp(WEPpreyL_SBL$predWlog))

# Estimate age compositions of pollock consumed, by survey year:
WEP_age_SBLdata = dcast(WEPpreyL_SBL, Year ~ AgeClass, value.var = "predW_unbi", sum)
AgeMatrix_SBL = as.matrix(round(WEP_age_SBLdata[,2:5], digits=0))
rownames(AgeMatrix_SBL) = c(1990,1993,1996,1999,2001,2003,2007,2009,2011) 
WEP_age_SBLdata$Year = as.numeric(as.character(WEP_age_SBLdata$Year))

SBLfit_gam = vgam(AgeMatrix_SBL ~ s(Year), multinomial, data=WEP_age_SBLdata)
  summary(SBLfit_gam)   

predProp_SBL_gam = as.data.frame(predict(SBLfit_gam, type="response"))
  
# Assign overall mean age compositions when data were missing:
predProp_SBL_gam[10,] = c(mean(predProp_SBL_gam$`0`), mean(predProp_SBL_gam$`1`), mean(predProp_SBL_gam$`2`), mean(predProp_SBL_gam$`3+`))   
predProp_SBL_gam[11,] = c(mean(predProp_SBL_gam$`0`), mean(predProp_SBL_gam$`1`), mean(predProp_SBL_gam$`2`), mean(predProp_SBL_gam$`3+`))
predProp_SBL_gam[12,] = c(mean(predProp_SBL_gam$`0`), mean(predProp_SBL_gam$`1`), mean(predProp_SBL_gam$`2`), mean(predProp_SBL_gam$`3+`))

predProp_SBL_gam$Year = c("1990", "1993", "1996", "1999", "2001", "2003", "2007", "2009", "2011", "2005", "2013", "2015")

SBL_multi_gam = melt(predProp_SBL_gam, id.vars = c("Year"), measure.vars = 1:4, variable.name = "AgeClass", value.name = "propWT")
SBL_multi_gam$Predator = "SBL"

# Save results:
saveRDS(SBL_multi_gam, "5_AgeCompositionsPrey/propWEP_age_SBL.rds")

# Plot (portion of Fig. 2):
SBL_multi_gam$AgeClass = ordered(SBL_multi_gam$AgeClass, levels = c("3+", "2", "1", "0"))
WEP_Age_SBL = ggplot(SBL_multi_gam, aes(x = as.numeric(as.character(Year)), y = propWT, fill = AgeClass, order = -as.numeric(AgeClass))) + 
  geom_bar(position = "fill", stat = "identity",color="black") +
  scale_fill_manual(values = c("chocolate1", "orange1", "gold", "white"), name=" Age", drop=FALSE) +
  theme_bw() +
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
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "vertical", legend.position = "right", legend.margin = margin(0,0,0,5), legend.background = element_rect(fill="transparent")) +
  theme(legend.position = c(0.0385, 0.865)) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(1.0, "mm")) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(legend.title = element_text(family="Arial", color="black", size=11)) +
  theme(legend.text = element_text(family="Arial", color="black", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=12)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  xlab("") +
  ylab("Proportion of Pollock Consumed") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) +
  scale_x_continuous(breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01))

WEP_Age_SBL

#########################################################
# Predator: Walleye Pollock #
#########################################################
# Input food habits data for Walleye Pollock:
WEPpreyL_WEP = read.csv("1_Data/GOA_RawPL_WEPdiets.csv", check.names=F)

# Remove fish not included in assessment-based estimates of total predator biomass (WEP < 37 cm):
WEPpreyL_WEP = subset(WEPpreyL_WEP, Pred_len >= 37)

# Select pollock as prey:
WEPpreyL_WEP = subset(WEPpreyL_WEP, Prey_Name=="Walleye pollock Gadus chalcogrammus")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
WEPpreyL_WEP$Year = as.factor(WEPpreyL_WEP$Year)
WEPpreyL_WEP = subset(WEPpreyL_WEP, Year != "1981" & Year != "1984" & Year != "1987")

# Predict pollock age from length:
WEPpreyL_WEP$predAge = (log(-(WEPpreyL_WEP$Prey_sz1/Linf)+1))/-K

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
WEPpreyL_WEP$AgeClass = with(WEPpreyL_WEP, 
    ifelse(predAge < 1, "0", 
        ifelse(predAge < 2, "1", 
            ifelse(predAge < 3, "2", "3+"))))

# Predict pollock weight from length:
WEPpreyL_WEP$predWlog = predict(L_W, data.frame(logL=log(WEPpreyL_WEP$Prey_sz1)), interval="c") 
WEPpreyL_WEP$predWlog = as.vector(WEPpreyL_WEP$predWlog[,1])
WEPpreyL_WEP$predW_unbi = cf *(exp(WEPpreyL_WEP$predWlog))

# Estimate age compositions of pollock consumed, by survey year:
WEP_age_WEPdata = dcast(WEPpreyL_WEP, Year ~ AgeClass, value.var = "predW_unbi", sum)

# Remove year with (essentially) zero pollock of any age class:
WEP_age_WEPdata = subset(WEP_age_WEPdata, Year != "2015")

AgeMatrix_WEP = as.matrix(round(WEP_age_WEPdata[,2:3], digits=0))
rownames(AgeMatrix_WEP) = c(1990,1993,1996,1999,2001,2003,2007,2009,2013) 
WEP_age_WEPdata$Year = as.numeric(as.character(WEP_age_WEPdata$Year))

WEPfit_gam = vgam(AgeMatrix_WEP ~ s(Year), multinomial, data=WEP_age_WEPdata)
  summary(WEPfit_gam) 

predProp_WEP_gam = as.data.frame(predict(WEPfit_gam, type="response"))

# Assign overall mean when missing data:
predProp_WEP_gam[10,] = c(mean(predProp_WEP_gam$`0`), mean(predProp_WEP_gam$`1`))   
predProp_WEP_gam[11,] = c(mean(predProp_WEP_gam$`0`), mean(predProp_WEP_gam$`1`))   
predProp_WEP_gam[12,] = c(mean(predProp_WEP_gam$`0`), mean(predProp_WEP_gam$`1`))   

predProp_WEP_gam$Year = c("1990", "1993", "1996", "1999", "2001", "2003", "2007", "2009", "2013", "2005", "2011", "2015")

WEP_multi_gam = melt(predProp_WEP_gam, id.vars = c("Year"), measure.vars = 1:2, variable.name = "AgeClass", value.name = "propWT")
WEP_multi_gam$Predator = "WEP"

# Save results:
saveRDS(WEP_multi_gam, "5_AgeCompositionsPrey/propWEP_age_WEP.rds")

# Plot (portion of Fig. 2):
WEP_multi_gam$AgeClass = ordered(WEP_multi_gam$AgeClass, levels = c("3+", "2", "1", "0"))
WEP_Age_WEP = ggplot(WEP_multi_gam, aes(x = as.numeric(as.character(Year)), y = propWT, fill = AgeClass, order = -as.numeric(AgeClass))) + 
  geom_bar(position = "fill", stat = "identity",color="black") +
  scale_fill_manual(values = c("deepskyblue4", "deepskyblue3", "lightskyblue", "white"), name=" Age", drop=FALSE) +
  theme_bw() +
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
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.direction = "vertical", legend.position = "right", legend.margin = margin(0,0,0,5), legend.background = element_rect(fill="transparent")) +
  theme(legend.position = c(0.0385, 0.865)) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(1.0, "mm")) +
  guides(fill=guide_legend(ncol=c(1))) +
  theme(legend.title = element_text(family="Arial", color="black", size=11)) +
  theme(legend.text = element_text(family="Arial", color="black", size=11), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(axis.line = element_line(color="black", linetype="solid")) +
  theme(axis.ticks=element_line(color="black")) +
  theme(axis.text.y = element_text(family="Arial", color="black", size=12)) +
  theme(axis.text.x = element_text(family="Arial", color="black",  size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(family="Arial", color="black", hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(family="Arial", color="black", vjust=2, size=12)) +
  xlab("") +
  ylab("Proportion of Pollock Consumed") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) +
  scale_x_continuous(breaks = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015), expand = c(0.01, 0.01))

WEP_Age_WEP