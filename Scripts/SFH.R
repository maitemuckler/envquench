## Diretórios ----
wdcode <- "Scripts/"
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"
wdfigs <- "~/Work/Research/Astronomy/Projects/envquench/Figures/"

# Bibliotecas ----
library(data.table)
library(binom)
library(ggplot2)
library(InformationValue)
library(caret)
library(dplyr)
library(scales)
library(viridis)
library(ggthemes)
library(pROC)
library(ggtext)
library(tidyr)
library(patchwork)

# Códigos extras ----
source("Scripts/Themes/ggplot_theme_Publication-2.R")
source("Scripts/Themes/my_theme.R")

# Opções ----
options(scipen = 999)

critval   <- 1.96 
TType_lim <- 0

width_figs  <- 11
height_figs <- 7

# Dados ----

input_data <- paste0("SFH_inputdata_zmax0.03_Rlim2.5_Ma12.3_logMstar_min9.17.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

df$SF_char <- df$SF
df$SF      <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType_label == "LTG", 1, 0)
df$LT <- as.factor(df$LT)

data.s_LIM <- subset(df, df$type == "Satellite")

input_data <- paste0("SFH_inputdata_zmax0.1_Rlim2.5_Ma12.3_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

df$SF_char <- df$SF
df$SF      <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType_label == "LTG", 1, 0)
df$LT <- as.factor(df$LT)

data.s_IHM <- subset(df, df$type == "Satellite")

rm(df)

# -------------------------------------------

# LI-M*

data.s_LIM <- subset(data.s_LIM, velDisp_e <= 200)

# Categorias para sigma 
breaks_sigma  <- c(30, 60, 80, 200)
levels_sigma  <- c("30 ≤ σₑ (km/s) < 60", "60 ≤ σₑ (km/s) < 80", "80 ≤ σₑ (km/s) ≤ 200")

data.s_LIM$sigma_char <- cut(10^data.s_LIM$logvelDisp_e, breaks = breaks_sigma, 
                             labels = levels_sigma, include.lowest = TRUE, right = FALSE)
data.s_LIM$sigma_char <- factor(data.s_LIM$sigma_char, levels = levels_sigma)

ggplot(data.s_LIM) + 
  geom_point(aes(x = logRproj_rvir, y = t_90, color = SF_char)) + 
  facet_grid(.~sigma_char)
  theme_Publication()
