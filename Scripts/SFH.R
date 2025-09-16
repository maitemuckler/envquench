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
data.s_LIM <- subset(data.s_LIM, data.s_LIM$SF == 0)

# Categorias para sigma 
breaks_sigma  <- c(30, 60, 80, 200)
levels_sigma  <- c("30 ≤ σₑ (km/s) < 60", "60 ≤ σₑ (km/s) < 80", "80 ≤ σₑ (km/s) ≤ 200")

data.s_LIM$sigma_char <- cut(10^data.s_LIM$logvelDisp_e, breaks = breaks_sigma, 
                             labels = levels_sigma, include.lowest = TRUE, right = FALSE)
data.s_LIM$sigma_char <- factor(data.s_LIM$sigma_char, levels = levels_sigma)

# Categorias para logRproj
breaks_Rproj_rvir  <- c(-Inf, 0.5, 2, Inf)
levels_Rproj_rvir  <- c("0.5 < Rproj/rvir", "0.5 ≤ Rproj/rvir < 2", "Rproj/rvir ≥ 2")

data.s_LIM$Rproj_rvir_char <- cut(data.s_LIM$Rproj_rvir, breaks = breaks_Rproj_rvir, 
                                  labels = levels_Rproj_rvir, include.lowest = TRUE, right = FALSE)
data.s_LIM$Rproj_rvir_char <- factor(data.s_LIM$Rproj_rvir_char, levels = levels_Rproj_rvir)

# Calcular a mediana por painel
medianas_LIM <- data.s_LIM %>%
  group_by(Rproj_rvir_char, sigma_char) %>%
  summarise(mediana_t90 = median(t_90, na.rm = TRUE),
            mediana_t80 = median(t_80, na.rm = TRUE),
            .groups = "drop")

# Adicionar uma coordenada y razoável para o texto (ajuste se necessário)
medianas_LIM <- medianas_LIM %>%
  mutate(y_pos = 0.75,  # altura no gráfico para o texto
         label = round(mediana_t90, 2))  # arredondar o valor para exibição

# Criar o gráfico
ggplot(data.s_LIM) + 
  geom_density(aes(x = t_90)) + 
  facet_grid(Rproj_rvir_char ~ sigma_char) + 
  geom_vline(data = medianas_LIM, aes(xintercept = mediana_t90), linetype = "dashed", color = "red") +
  geom_text(data = medianas_LIM, aes(x = mediana_t90, y = y_pos, label = label), 
            color = "red", hjust = -0.1, size = 3.5) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) + 
  labs(x = "t_90", y = "density") + 
  ggtitle("LI-M*") + 
  theme_Publication()

ggsave(path = wdfigs,
       filename = paste0("LIM_t90.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs,
       units = "in", dpi = 600)

# -------------------------------------------

# IH-M*

data.s_IHM <- subset(data.s_IHM, data.s_IHM$SF == 0)

# Categorias para sigma 
breaks_sigma  <- c(30, 100, 150, Inf) # máximo de sigma round(10^(max(data.s$logvelDisp_e)))
levels_sigma  <- c("30 ≤ σₑ (km/s) < 100", "100 ≤ σₑ (km/s) < 150", "150 ≤ σₑ (km/s) < 317")

data.s_IHM$sigma_char <- cut(10^data.s_IHM$logvelDisp_e, breaks = breaks_sigma, 
                             labels = levels_sigma, include.lowest = TRUE, right = FALSE)
data.s_IHM$sigma_char <- factor(data.s_IHM$sigma_char, levels = levels_sigma)

# Categorias para logRproj
breaks_Rproj_rvir  <- c(-Inf, 0.5, 2, Inf)
levels_Rproj_rvir  <- c("0.5 < Rproj/rvir", "0.5 ≤ Rproj/rvir < 2", "Rproj/rvir ≥ 2")

data.s_IHM$Rproj_rvir_char <- cut(data.s_IHM$Rproj_rvir, breaks = breaks_Rproj_rvir, 
                                  labels = levels_Rproj_rvir, include.lowest = TRUE, right = FALSE)
data.s_IHM$Rproj_rvir_char <- factor(data.s_IHM$Rproj_rvir_char, levels = levels_Rproj_rvir)

# Calcular a mediana por painel
medianas_IHM <- data.s_IHM %>%
  group_by(Rproj_rvir_char, sigma_char) %>%
  summarise(mediana_t90 = median(t_90, na.rm = TRUE),
            mediana_t80 = median(t_80, na.rm = TRUE),
            .groups = "drop")

# Adicionar uma coordenada y razoável para o texto (ajuste se necessário)
medianas_IHM <- medianas_IHM %>%
  mutate(y_pos = 0.75,  # altura no gráfico para o texto
         label = round(mediana_t90, 2))  # arredondar o valor para exibição

# Criar o gráfico
ggplot(data.s_IHM) + 
  geom_density(aes(x = t_90)) + 
  facet_grid(Rproj_rvir_char ~ sigma_char) + 
  geom_vline(data = medianas_IHM, aes(xintercept = mediana_t90), linetype = "dashed", color = "red") +
  geom_text(data = medianas_IHM, aes(x = mediana_t90, y = y_pos, label = label), 
            color = "red", hjust = -0.1, size = 3.5) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) + 
  labs(x = "t_90", y = "density") + 
  ggtitle("IH-M*") + 
  theme_Publication()

ggsave(path = wdfigs,
       filename = paste0("IHM_t90.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs,
       units = "in", dpi = 600)

