## Pacotes ----
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)

## Diretórios ----
wdcode <- "/home/muckler/Work/Research/Astronomy/Projects/envquenching/Scripts/"
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"
wdfigs <- "~/Work/Research/Astronomy/Projects/envquenching/Figures/"

## Outros códigos e funções ----
source(file.path(wdcode, "Themes/my_theme.R"))
source(file.path(wdcode, "MyFunctions/sigma_gal_correction.R"))
source(file.path(wdcode, "MyFunctions/SF_Q_class.R"))
source(file.path(wdcode, "MyFunctions/distance_to_line_calc.R"))

# Substituindo valores inválidos por NA
replace_na_values <- function(data, columns, threshold) {
  data %>% mutate(across(all_of(columns), ~ if_else(.x <= threshold, NA_real_, .x)))
}

## Configuração inicial ----
set.seed(123)
figs_width  <- 14
figs_height <- 8

## Definindo input e output files ----
input_file  <- "SDSS_DR18_Legacy_MGS_QSO_kcorr+GSWLC-X2_MGS+Simard11+DS18.csv"
output_file <- paste0(wddata, "clean_", input_file)

## Lendo os dados ----
df <- fread(paste0(wddata, input_file))

## Tratamento de dados ----

na_replacements <- list(
  list(columns = c("magPetro_u", "magPetro_g", "magPetro_r", "magPetro_i", "magPetro_z",
                   "magFiber_u", "magFiber_g", "magFiber_r", "magFiber_i", "magFiber_z",
                   "petroR90_r", "petroR50_r", "lgm_tot_p50", "lgm_tot_p16", "lgm_tot_p84", 
                   "sfr_tot_p50", "sfr_tot_p16", "sfr_tot_p84"), threshold = -9999),
  list(columns = c("logMstar", "logMstar_err", "logSFR_SED", "logSFR_SED_err"), threshold = -99)
)

for (replacement in na_replacements) {
  df <- replace_na_values(df, replacement$columns, replacement$threshold)
}

## Removendo NAs de colunas essenciais
df <- df %>% drop_na(absPetro_r, logMstar, logSFR_SED)

## Correção da dispersão de velocidades ----
df$velDisp_e    <- sigma_gal_correction(df$velDisp, df$Rhlr, df$Scale)
df$logvelDisp_e <- log10(df$velDisp_e)

# Removendo valores inválidos de sigma
df <- df %>% drop_na(velDisp_e) %>% filter(!is.infinite(logvelDisp_e))

# Corte em 30 km/s < velDisp_e < 320 km/s
limits <- c(30, 320)
df <- df %>% filter(between(velDisp_e, limits[1], limits[2]))

# Verificando se ficou algum NA
data.frame(NA_count = sapply(df, function(x) sum(is.na(x))))

# Criando variáveis --------------------

## Classificação SF/Q ----
# Tudo usando dados do SALIM!!!!
df$logsSFR  <- df$logSFR_SED - df$logMstar

df$SF <- SF_Q_class(df$logsSFR,  df$logMstar, slope = -0.45, intercept = -6.55)
df$SF <- factor(df$SF, levels = c("Star-forming", "Quiescent"))
table(df$SF)

## Criando variável dicotômica para T-Type ----
TType_lim      <- 0
df$TType_label <- factor(ifelse(df$TType <= TType_lim, "ETG", "LTG"))
table(df$TType_label)

# Gráfico de distribuição de T-Type
ggplot(data = df, aes(x = TType)) + 
  geom_histogram(aes(y = after_stat(density), 
                     fill = ifelse(TType <= TType_lim, "Early-type", "Late-type")), 
                 colour = "white", bins = 50) +
  scale_fill_manual(values = c('red', 'blue'), name = "Morphology") +
  labs(x = "T-Type", y = "Densidade") +
  theme_Publication()

## Adicionando labels para variáveis factor ----
df$bptClass <- factor(df$bptClass, 
                      labels = c("unclassifiable", "star-forming", "low S/N star-forming", "composite", "AGN (excluding liners)", "low S/N LINER"))
df$flag_sed <- factor(df$flag_sed, 
                      levels = c(0, 1, 2, 3, 5), 
                      labels = c("OK", "broad-line spectrum", "X²_r > 30", "Não encontrei label", "missing SDSS photometry"))
df$flag_uv  <- factor(df$flag_uv, 
                      labels = c("no UV", "FUV only", "NUV only", "both"))
df$flag_midir <- factor(df$flag_midir, 
                        labels = c("no mid-IR", "LIR based on 12 μm", "LIR based on 22 μm", "LIR corrected for mid-IR AGN emission", "Não encontrei label 6", "Não encontrei label 7"))

## Conversão de variáveis
lapply(df, class)
colunas_para_converter <- c("bptClass", "flag_sed", "flag_uv", "flag_midir", "SF", "TType_label")
df <- df %>% mutate(across(all_of(colunas_para_converter), factor))

## Salvando os dados ----
write.csv(df, output_file, row.names = FALSE)
