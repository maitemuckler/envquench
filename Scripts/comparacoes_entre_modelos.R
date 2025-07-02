setwd("~/Work/Research/Astronomy/Projects/envquenching/")

# Setting:
TType_lim     <- 0
critval       <- 1.96 
Ma            <- 12.3

LIMass <- read.csv("inputdata_zmax0.03_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min9.17.csv")
IHMass <- read.csv("inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")

# Satélites

# LI-Mass
LIMass$SF_char <- factor(LIMass$SF)
LIMass$SF      <- ifelse(LIMass$SF_char == "Star-forming", 1, 0)

LIMass$LT <- ifelse(LIMass$TType >= TType_lim, 1, 0)
LIMass$LT <- factor(LIMass$LT)

LIMass.s <- subset(LIMass, LIMass$type == "Satellite")
LIMass.c <- subset(LIMass, LIMass$type == "Central")

# IH-Mass
IHMass$SF_char <- factor(IHMass$SF)
IHMass$SF      <- ifelse(IHMass$SF_char == "Star-forming", 1, 0)

IHMass$LT <- ifelse(IHMass$TType >= TType_lim, 1, 0)
IHMass$LT <- factor(IHMass$LT)

IHMass.s <- subset(IHMass, IHMass$type == "Satellite")
IHMass.c <- subset(IHMass, IHMass$type == "Central")

length(which(LIMass.s$logMgroup < 14))/nrow(LIMass.s)
length(which(LIMass.s$logMgroup >= 14))/nrow(LIMass.s)

length(which(IHMass.s$logMgroup < 14))/nrow(IHMass.s)
length(which(IHMass.s$logMgroup >= 14))/nrow(IHMass.s)

# Modelos com Mgroup
fSFG_LIMass <- glm(SF ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = LIMass.s)
fSFG_IHMass <- glm(SF ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = IHMass.s)

fLTG_LIMass <- glm(LT ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = LIMass.s)
fLTG_IHMass <- glm(LT ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = IHMass.s)

summary(fSFG_LIMass)
summary(fSFG_IHMass)

summary(fLTG_LIMass)
summary(fLTG_IHMass)

# Modelos com Mstar
fSFG_LIMass <- glm(SF ~ logMstar + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = LIMass.s)
fSFG_IHMass <- glm(SF ~ logMstar + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = IHMass.s)

fLTG_LIMass <- glm(LT ~ logMstar + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = LIMass.s)
fLTG_IHMass <- glm(LT ~ logMstar + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = IHMass.s)

summary(fSFG_LIMass)
summary(fSFG_IHMass)

summary(fLTG_LIMass)
summary(fLTG_IHMass)

# Modelos com sigma e Mstar
fSFG_LIMass <- glm(SF ~ logMstar + logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = LIMass.s)
fSFG_IHMass <- glm(SF ~ logMstar + logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = IHMass.s)

fLTG_LIMass <- glm(LT ~ logMstar + logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = LIMass.s)
fLTG_IHMass <- glm(LT ~ logMstar + logvelDisp_e + logRproj_rvir + logMgroup, family = binomial(link = "logit"), data = IHMass.s)

summary(fSFG_LIMass)
summary(fSFG_IHMass)

summary(fLTG_LIMass)
summary(fLTG_IHMass)
