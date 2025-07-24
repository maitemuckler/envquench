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
input_data <- paste0("inputdata_zmax0.03_Rlim2.5_Ma12.3_logMstar_min9.17.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

df$SF_char <- df$SF
df$SF      <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType_label == "LTG", 1, 0)
df$LT <- as.factor(df$LT)

data.s_LIM <- subset(df, df$type == "Satellite")

input_data <- paste0("inputdata_zmax0.1_Rlim2.5_Ma12.3_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

df$SF_char <- df$SF
df$SF      <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType_label == "LTG", 1, 0)
df$LT <- as.factor(df$LT)

data.s_IHM <- subset(df, df$type == "Satellite")

# Exploração!

# Combinação das categorias
data.s_LIM$tipo <- paste0(data.s_LIM$SF_char, "_", data.s_LIM$TType_label)
data.s_LIM$tipo <- as.factor(data.s_LIM$tipo)

table(data.s_LIM$tipo)/nrow(data.s_LIM) * 100

data.s_IHM$tipo <- paste0(data.s_IHM$SF_char, "_", data.s_IHM$TType_label)
data.s_IHM$tipo <- as.factor(data.s_IHM$tipo)

table(data.s_IHM$tipo)/nrow(data.s_IHM) * 100



# nr de grupos
length(unique(data.s$groupID))

# SF/Q
table(data.s$SF_char)

# LTG/ETG
table(data.s$TType_label)



paleta <- c("Quiescent_ETG" = "#ef233c",
            "Quiescent_LTG" = "#fca311",
            "Star-forming_ETG" = "#7E508D",
            "Star-forming_LTG" = "#00509d")





ggplot(data.s) + 
  geom_density(aes(x = logvelDisp_e, fill = SF_char, color = SF_char, linetype = TType_label), alpha = 0.1) + 
  theme_Publication()

ggplot(data.s) + 
  geom_density(aes(x = logMgroup, fill = SF_char, color = SF_char, linetype = TType_label), alpha = 0.1) + 
  facet_grid(.~SF_char) + 
  theme_Publication()


ggplot(data.s) + 
  geom_point(aes(x = logMgroup, y = logvelDisp_e, color = tipo)) + 
  scale_color_manual(values = paleta)  +
  theme_Publication()


