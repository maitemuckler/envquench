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

data.s_IHM <- subset(df, df$type == "Satellite")

input_data <- paste0("SFH_inputdata_zmax0.1_Rlim2.5_Ma12.3_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

df$SF_char <- df$SF
df$SF      <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType_label == "LTG", 1, 0)
df$LT <- as.factor(df$LT)

data.s_IHM <- subset(df, df$type == "Satellite")

rm(df)

#-------------------------------------------------------------------------


# Função auxiliar para calcular proporções
calc_props <- function(df, vars, sample_name) {
  res <- lapply(vars, function(v) {
    prop_low  <- round(sum(df[[v]] <= 0.5) / nrow(df) * 100)
    prop_high <- round(sum(df[[v]] > 0.5) / nrow(df) * 100)
    c(prop_low, prop_high)
  })
  
  res_df <- as.data.frame(do.call(rbind, res))
  names(res_df) <- c("<=0.5", ">0.5")
  res_df$variavel <- vars
  res_df$amostra <- sample_name
  res_df
}

# Lista de variáveis
vars <- c("P_bulge", "P_disk", "P_bar_GZ2", "P_bar_Nair10", 
          "P_merg", "P_cigar", "P_S0")

# Calcular para as duas amostras
df_IHM <- calc_props(data.s_IHM, vars, "IHM")
df_LIM <- calc_props(data.s_LIM, vars, "LIM")

# Juntar em um só data.frame
df_final <- rbind(df_IHM, df_LIM)


# Filtrar apenas proporções > 0.5
df_high <- df_final %>%
  select(variavel, amostra, `>0.5`) %>%
  rename(proporcao = `>0.5`)

# Gráfico de barras comparando amostras
ggplot(df_high, aes(x = variavel, y = proporcao, fill = amostra)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("red", "blue")) + 
  scale_y_continuous(breaks = pretty_breaks(n = 11)) + 
  labs(x = "", y = "Proportion (%)", fill = "Sample",
       title = "Proportion of galaxies with a probability > 0.5 of exhibiting each feature") +
  theme_Publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "inside",
        legend.position.inside = c(0.15, 0.9))

ggsave(path = wdfigs,
       filename = paste0("prop_morfologicas.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs,
       units = "in", dpi = 600)

# ------------ fSF vs fLT em bins de Rvir

# LIM
# Categorias para logRproj
breaks_Rproj_rvir  <- c(-Inf, 0.5, 1, 2, Inf)
levels_Rproj_rvir  <- c("0.5 < Rproj/rvir", "0.5 ≤ Rproj/rvir < 1", "1 ≤ Rproj/rvir < 2", "Rproj/rvir ≥ 2")

data.s_LIM$Rproj_rvir_char <- cut(data.s_LIM$Rproj_rvir, breaks = breaks_Rproj_rvir, 
                                  labels = levels_Rproj_rvir, include.lowest = TRUE, right = FALSE)
data.s_LIM$Rproj_rvir_char <- factor(data.s_LIM$Rproj_rvir_char, levels = levels_Rproj_rvir)


data.s_LIM$Rproj_rvir_char
data.s_LIM$LT <- as.numeric(as.character(data.s_LIM$LT))

dentro <- data.s_LIM[data.s_LIM$Rproj_rvir_char == "0.5 < Rproj/rvir",]
sum(dentro$SF) / nrow(dentro)
sum(dentro$LT) / nrow(dentro)

meio1 <- data.s_LIM[data.s_LIM$Rproj_rvir_char == "0.5 ≤ Rproj/rvir < 1",]
sum(meio1$SF) / nrow(meio1)
sum(meio1$LT) / nrow(meio1)

meio2 <- data.s_LIM[data.s_LIM$Rproj_rvir_char == "1 ≤ Rproj/rvir < 2",]
sum(meio2$SF) / nrow(meio2)
sum(meio2$LT) / nrow(meio2)

fora <- data.s_LIM[data.s_LIM$Rproj_rvir_char == "Rproj/rvir ≥ 2",]
sum(fora$SF) / nrow(fora)
sum(fora$LT) / nrow(fora)


pontos <- data.frame(Rproj_bin = levels(data.s_LIM$Rproj_rvir_char),
           fSFG = c(sum(dentro$SF) / nrow(dentro), 
                    sum(meio1$SF) / nrow(meio1), 
                    sum(meio2$SF) / nrow(meio2), 
                    sum(fora$SF) / nrow(fora)),
           fLTG = c(sum(dentro$LT) / nrow(dentro), 
                    sum(meio1$LT) / nrow(meio1), 
                    sum(meio2$LT) / nrow(meio2), 
                    sum(fora$LT) / nrow(fora)))


plot(pontos$fSFG, pontos$fLTG)
