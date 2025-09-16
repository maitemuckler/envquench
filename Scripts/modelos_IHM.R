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
input_data <- paste0("inputdata_zmax0.1_Rlim2.5_Ma12.3_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

df$SF_char <- df$SF
df$SF      <- ifelse(df$SF_char == "Star-forming", 1, 0)

df$LT <- ifelse(df$TType_label == "LTG", 1, 0)
df$LT <- as.factor(df$LT)

data.s <- subset(df, df$type == "Satellite")
data.c <- subset(df, df$type == "Central")

# nr de grupos
length(unique(data.s$groupID))

# SF/Q
table(data.s$SF_char)
round(table(data.s$SF_char)/nrow(data.s)*100)

# LTG/ETG
table(data.s$TType_label)
round(table(data.s$TType_label)/nrow(data.s)*100)

# Ajusta modelos 
fit_fSFG <- glm(SF ~ logvelDisp_e + logMgroup + logRproj_rvir, family = binomial("logit"), data = data.s)
fit_fLTG <- glm(LT ~ logvelDisp_e + logMgroup + logRproj_rvir, family = binomial("logit"), data = data.s)

# Predição

data.s$pred_fSFG <- predict(fit_fSFG, data.s, type = "response")
optimal_fSFG     <- optimalCutoff(data.s$SF, data.s$pred_fSFG)[1]
misClass_fSFG    <- misClassError(data.s$SF, data.s$pred_fSFG, threshold = optimal_fSFG)
acuracia_fSFG    <- 1 - misClass_fSFG
acuracia_fSFG

data.s$pred_fLTG <- predict(fit_fLTG, data.s, type = "response")
optimal_fLTG     <- optimalCutoff(data.s$LT, data.s$pred_fLTG)[1]
misClass_fLTG    <- misClassError(data.s$LT, data.s$pred_fLTG, threshold = optimal_fLTG)
acuracia_fLTG    <- 1 - misClass_fLTG
acuracia_fLTG

### FIGURA --------------------------

# Categorias para sigma 
breaks_sigma  <- c(30, 100, 150, Inf) # máximo de sigma round(10^(max(data.s$logvelDisp_e)))
levels_sigma  <- c("30 ≤ σₑ (km/s) < 100", "100 ≤ σₑ (km/s) < 150", "150 ≤ σₑ (km/s) < 317")

data.s$sigma_char <- cut(10^data.s$logvelDisp_e, breaks = breaks_sigma, labels = levels_sigma, include.lowest = TRUE, right = FALSE)
data.s$sigma_char <- factor(data.s$sigma_char, levels = levels_sigma)

# Categorias para logMgroup
breaks_logMgroup  <- c(-Inf, 14, Inf)
levels_logMgroup  <- c("< 14", "≥ 14")
facet_labels_logMgroup <- c("< 14" = "log[10](M[h]/M['\\u2609']) < 14", "≥ 14" = "log[10](M[h]/M['\\u2609']) >= 14")

data.s$logMgroup_char <- cut(data.s$logMgroup, breaks = breaks_logMgroup, labels = levels_logMgroup, include.lowest = TRUE, right = FALSE)

# Inicializa data frames
modelos <- bins_total <- data.frame()

# Define domínio fixo de Rproj (global para todos os painéis)
range_logRproj     <- range(data.s$logRproj_rvir)
bins_logRproj_rvir <- seq(range_logRproj[1], range_logRproj[2], by = 0.01)

# Loop por sigma_char e logMgroup_char
for (sigma_label in levels(data.s$sigma_char)) {
  
  data_sigma <- subset(data.s, sigma_char == sigma_label)
  
  for (grupo in levels(data.s$logMgroup_char)) {
    
    data_sub <- subset(data_sigma, logMgroup_char == grupo)
    if (nrow(data_sub) == 0) next
    
    modelo <- data.frame(
      logvelDisp_e   = median(data_sub$logvelDisp_e),
      logRproj_rvir  = bins_logRproj_rvir,
      logMgroup      = median(data_sub$logMgroup)
    )
    
    pred_fSFG <- predict(fit_fSFG, type = "link", se.fit = TRUE, newdata = modelo)
    pred_fLTG <- predict(fit_fLTG, type = "link", se.fit = TRUE, newdata = modelo)
    
    modelo$pred_fSFG_prob <- fit_fSFG$family$linkinv(pred_fSFG$fit)
    modelo$prob_fSFG_upr  <- fit_fSFG$family$linkinv(pred_fSFG$fit + critval * pred_fSFG$se.fit)
    modelo$prob_fSFG_lwr  <- fit_fSFG$family$linkinv(pred_fSFG$fit - critval * pred_fSFG$se.fit)
    
    modelo$pred_fLTG_prob <- fit_fLTG$family$linkinv(pred_fLTG$fit)
    modelo$prob_fLTG_upr  <- fit_fLTG$family$linkinv(pred_fLTG$fit + critval * pred_fLTG$se.fit)
    modelo$prob_fLTG_lwr  <- fit_fLTG$family$linkinv(pred_fLTG$fit - critval * pred_fLTG$se.fit)
    
    modelo$sigma_char      <- sigma_label
    modelo$logMgroup_char  <- grupo
    
    modelos <- rbind(modelos, modelo)
    
    # Criação dos bins de logRproj_rvir e cálculo das frações
    aux <- quantile(data_sub$logRproj_rvir, probs = seq(0, 1, by = 0.2))
    for (i in 1:(length(aux)-1)) {
      bin <- data_sub$logRproj_rvir >= aux[i] & data_sub$logRproj_rvir < aux[i+1]
      bin_data <- data_sub[bin, ]
      if (nrow(bin_data) == 0) next
      
      meio_bin <- median(c(aux[i], aux[i+1]))
      ngal     <- nrow(bin_data)
      
      ngal_SF <- sum(bin_data$SF == 1)
      ngal_LT <- sum(bin_data$LT == 1)
      
      fSFG <- ngal_SF / ngal
      fLTG <- ngal_LT / ngal
      
      fSFG_up <- binom.confint(ngal_SF, ngal, conf.level = 0.95, method = "bayes")$upper
      fSFG_lw <- binom.confint(ngal_SF, ngal, conf.level = 0.95, method = "bayes")$lower
      fLTG_up <- binom.confint(ngal_LT, ngal, conf.level = 0.95, method = "bayes")$upper
      fLTG_lw <- binom.confint(ngal_LT, ngal, conf.level = 0.95, method = "bayes")$lower
      
      bins <- data.frame(
        sigma_char        = sigma_label,
        meio_bin          = meio_bin,
        fSFG              = fSFG, fSFG_up = fSFG_up, fSFG_lw = fSFG_lw,
        fLTG              = fLTG, fLTG_up = fLTG_up, fLTG_lw = fLTG_lw,
        logMgroup_char    = grupo,
        ngal              = ngal
      )
      
      bins_total <- rbind(bins_total, bins)
    }
  }
}

# Fatores ordenados
modelos$sigma_char        <- factor(modelos$sigma_char, levels = levels_sigma)
bins_total$sigma_char     <- factor(bins_total$sigma_char, levels = levels_sigma)
modelos$logMgroup_char    <- factor(modelos$logMgroup_char, levels = levels_logMgroup)
bins_total$logMgroup_char <- factor(bins_total$logMgroup_char, levels = levels_logMgroup)

# 1. Calcula os valores extremos diretamente de modelos
rvir_min <- 0.02

extremos_df <- modelos %>%
  group_by(sigma_char, logMgroup_char) %>%
  summarise(
    x_min = log10(rvir_min),
    x_max = max(logRproj_rvir),
    
    pred_fSFG_prob_min = pred_fSFG_prob[which.min(abs(logRproj_rvir - x_min))],
    pred_fSFG_prob_max = pred_fSFG_prob[which.max(logRproj_rvir)],
    
    pred_fLTG_prob_min = pred_fLTG_prob[which.min(abs(logRproj_rvir - x_min))],
    pred_fLTG_prob_max = pred_fLTG_prob[which.max(logRproj_rvir)],
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(pred_fSFG_prob_min, pred_fSFG_prob_max, pred_fLTG_prob_min, pred_fLTG_prob_max),
    names_to = "tipo",
    values_to = "y"
  ) %>%
  mutate(
    modelo = ifelse(grepl("SFG", tipo), "fSFG", "fLTG"),
    x = ifelse(grepl("min", tipo), x_min, x_max),
    tipo_ponto = ifelse(grepl("min", tipo), "min", "max")
  )


x_offset <- 0.4  # mais conservador
extremos_df <- extremos_df %>%
  mutate(
    painel_full_id = paste0(sigma_char, "___", logMgroup_char),
    
    vjust_label = case_when(
      painel_full_id == "30 ≤ σₑ (km/s) < 100___< 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "30 ≤ σₑ (km/s) < 100___< 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "30 ≤ σₑ (km/s) < 100___< 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.1,
      painel_full_id == "30 ≤ σₑ (km/s) < 100___< 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "100 ≤ σₑ (km/s) < 150___< 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "100 ≤ σₑ (km/s) < 150___< 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "100 ≤ σₑ (km/s) < 150___< 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "100 ≤ σₑ (km/s) < 150___< 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "150 ≤ σₑ (km/s) < 317___< 14"  & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "150 ≤ σₑ (km/s) < 317___< 14"  & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "150 ≤ σₑ (km/s) < 317___< 14"  & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "150 ≤ σₑ (km/s) < 317___< 14"  & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "30 ≤ σₑ (km/s) < 100___≥ 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.2,
      painel_full_id == "30 ≤ σₑ (km/s) < 100___≥ 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.4,
      painel_full_id == "30 ≤ σₑ (km/s) < 100___≥ 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.8,
      painel_full_id == "30 ≤ σₑ (km/s) < 100___≥ 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.5,
      
      painel_full_id == "100 ≤ σₑ (km/s) < 150___≥ 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "100 ≤ σₑ (km/s) < 150___≥ 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "100 ≤ σₑ (km/s) < 150___≥ 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "100 ≤ σₑ (km/s) < 150___≥ 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "150 ≤ σₑ (km/s) < 317___≥ 14"  & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "150 ≤ σₑ (km/s) < 317___≥ 14"  & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "150 ≤ σₑ (km/s) < 317___≥ 14"  & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "150 ≤ σₑ (km/s) < 317___≥ 14"  & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      TRUE ~ 0
    ),
    
    x_label = ifelse(tipo_ponto == "min", x - x_offset, x + x_offset)
  )

# 3. Cria padding_df para forçar expansão por painel com pontos invisíveis
limites_x <- modelos %>%
  group_by(sigma_char, logMgroup_char) %>%
  summarise(
    x_min = min(logRproj_rvir),
    x_max = max(logRproj_rvir),
    .groups = "drop"
  ) %>%
  mutate(
    x_min_exp = x_min - 0.5,
    x_max_exp = x_max + 0.5
  )

padding_df <- limites_x %>%
  mutate(
    y = NA,
    modelo = "fSFG"
  ) %>%
  pivot_longer(cols = c(x_min_exp, x_max_exp), names_to = NULL, values_to = "x") %>%
  select(sigma_char, logMgroup_char, x, y, modelo)

label_df <- data.frame(
  x = -1.4,
  y = 0.1,
  label = 'italic(IH-M["★"]~sample)',
  sigma_char = levels_sigma[1],
  logMgroup_char = levels_logMgroup[2]
  
)

label_df$sigma_char     <- factor(label_df$sigma_char, levels = levels_sigma)
label_df$logMgroup_char <- factor(label_df$logMgroup_char, levels = levels_logMgroup)

# Gráfico final 
ggplot() + 
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fSFG_prob, color = "fSFG", group = interaction(sigma_char, logMgroup_char)), linewidth = 1) + 
  geom_ribbon(data = modelos, aes(ymin = prob_fSFG_lwr, ymax = prob_fSFG_upr, x = logRproj_rvir, fill = "fSFG", group = interaction(sigma_char, logMgroup_char)), alpha = 0.5, show.legend = FALSE) + 
  
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fLTG_prob, color = "fLTG", group = interaction(sigma_char, logMgroup_char)), linewidth = 1) + 
  geom_ribbon(data = modelos, aes(ymin = prob_fLTG_lwr, ymax = prob_fLTG_upr, x = logRproj_rvir, fill = "fLTG", group = interaction(sigma_char, logMgroup_char)), alpha = 0.5, show.legend = FALSE) + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fSFG, color = "fSFG", group = interaction(sigma_char, logMgroup_char)), size = 2) +
  geom_errorbar(data = bins_total, aes(ymin = fSFG_lw, ymax = fSFG_up, x = meio_bin, color = "fSFG", group = interaction(sigma_char, logMgroup_char)), width = 0.1) + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fLTG, color = "fLTG", group = interaction(sigma_char, logMgroup_char)), size = 2) +
  geom_errorbar(data = bins_total, aes(ymin = fLTG_lw, ymax = fLTG_up, x = meio_bin, color = "fLTG", group = interaction(sigma_char, logMgroup_char)), width = 0.1) + 
  
  # Pontos nos extremos das curvas
  geom_point(data = extremos_df, 
             aes(x = x, y = y, color = modelo, group = sigma_char),
             size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Rótulos percentuais nos extremos
  geom_text(data = extremos_df,
            aes(x = x_label, y = y, label = percent(y, accuracy = 1), 
                color = modelo, group = sigma_char, 
                vjust = vjust_label),
            size = 4, show.legend = FALSE) +
  
  # Número de galáxias por bin
  geom_text(data = bins_total,
            aes(x = meio_bin, y = 1.2, label = ngal, group = interaction(sigma_char, logMgroup_char)),
            inherit.aes = FALSE,
            color = "black", size = 3,
            angle = -90) +
  
  # Adiciona nome da amostra
  geom_text(data = label_df,
            aes(x = x, y = y, label = label),
            parse = TRUE,
            inherit.aes = FALSE,
            hjust = 0,
            size = 5) +

  facet_grid(logMgroup_char ~ sigma_char,
             scales = "free_x",
             labeller = labeller(
               logMgroup_char = as_labeller(facet_labels_logMgroup, label_parsed))) + 
  
  # Padding invisível para forçar expansão do eixo X em cada facet
  geom_blank(data = padding_df, aes(x = x, group = interaction(sigma_char, logMgroup_char))) +
  
  ylab("Fraction") + 
  xlab(expression(log[10](R[proj]/r[vir]))) + 
  
  scale_color_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F"), 
                     name = "",
                     labels = c("fSFG" = expression(f[SFG]), "fLTG" = expression(f[LTG]))) + 
  scale_fill_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F")) + 
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.23)) +
  #scale_x_continuous(sec.axis = sec_axis(~ . , expression(paste(sigma[e], " (km/s)")), breaks = NULL, labels = NULL)) + 
  
  theme_Publication() + 
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.047, 0.075),
    legend.box.background = element_rect(color = "black", linewidth = 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
    strip.text.x = element_text(face = "plain"),
    strip.text.y = element_text(face = "plain")
  )

ggsave(path = wdfigs,
       filename = paste0("logistic_IHM.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs,
       units = "in", dpi = 600)

# Table com Delta_r
# f_esperado seria longe, 
# f_obs seria perto (o que é "observado de efeito ambiental")
# (f_esperado - f_obs) / f_esperado 

delta_r <- extremos_df[,c("modelo", "sigma_char", "logMgroup_char", "tipo_ponto", "y")]
delta_r$y_perc      <- round(delta_r$y * 100)

delta_r <- delta_r %>%
  group_by(modelo, sigma_char, logMgroup_char) %>%
  mutate(delta_y = (y[tipo_ponto == "max"] - y[tipo_ponto == "min"]) / y[tipo_ponto == "max"]) %>%
  mutate(delta_y_perc = (y_perc[tipo_ponto == "max"] - y_perc[tipo_ponto == "min"]) / y_perc[tipo_ponto == "max"])

delta_r$delta_y      <- round(delta_r$delta_y * 100)
delta_r$delta_y_perc <- round(delta_r$delta_y_perc * 100)

# tabela que vai no artigo
tab_delta_r <- 
  delta_r %>% 
  select(modelo, sigma_char, logMgroup_char, tipo_ponto, y_perc, delta_y_perc)


### Tabelas de coeficientes, OR, etc -------------

# Modelo fSFG
summary(fit_fSFG)
confint(fit_fSFG)
# Odds Ratio
OR_fSFG <- exp(coef(fit_fSFG))
OR_fSFG

# IC para Odds Ratio
conf_int_OR_fSFG <- exp(confint(fit_fSFG))
conf_int_OR_fSFG

# Tabela para OR
OR_table_fSFG <- data.frame(
  Variable = names(OR_fSFG),
  OR = OR_fSFG,
  Lower95 = conf_int_OR_fSFG[,1],
  Upper95 = conf_int_OR_fSFG[,2]
)

rownames(OR_table_fSFG) <- NULL

OR_table_fSFG

# Modelo fLTG
summary(fit_fLTG)
confint(fit_fLTG)
# Odds Ratio
OR_fLTG <- exp(coef(fit_fLTG))
OR_fLTG

# IC para Odds Ratio
conf_int_OR_fLTG <- exp(confint(fit_fLTG))
conf_int_OR_fLTG

# Tabela para OR
OR_table_fLTG <- data.frame(
  Variable = names(OR_fLTG),
  OR = OR_fLTG,
  Lower95 = conf_int_OR_fLTG[,1],
  Upper95 = conf_int_OR_fLTG[,2]
)

rownames(OR_table_fLTG) <- NULL

OR_table_fLTG

# Tabelas adicionais
# Valores medianos dos modelos

modelos %>% 
  distinct(sigma_char, logMgroup_char, round(10^logvelDisp_e, 2), round(logMgroup, 2))

# Densidades
label_df_density <- data.frame(
  label = 'italic(IH-M["★"]~sample)',
  sigma_char = levels_sigma[1],
  logMgroup_char = levels_logMgroup[1]
)

label_df_density$sigma_char     <- factor(label_df_density$sigma_char, levels = levels_sigma)
label_df_density$logMgroup_char <- factor(label_df_density$logMgroup_char, levels = levels_logMgroup)

# Gráfico 1: dispersão de velocidades
p1 <- ggplot(data.s, aes(x = logvelDisp_e)) + 
  geom_histogram(aes(y = after_stat(density), color = sigma_char), fill = 'white', bins = 8) + 
  geom_density(aes(fill = sigma_char, color = sigma_char), alpha = 0.5) +
  scale_color_brewer(palette = "Accent") + 
  scale_fill_brewer(palette = "Accent") + 
  facet_grid(. ~ sigma_char, scales = "free_x") + 
  xlab(label_logvelDisp_e) + 
  geom_text(data = label_df_density,
            aes(x = 1.5, y = 7.5, label = label),
            parse = TRUE, hjust = 0, size = 5) +
  theme_Publication() + 
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "plain"),
    strip.text.y = element_text(face = "plain")
  )

p1

# Gráfico 2: massa do halo
p2 <- ggplot(data.s, aes(x = logMgroup)) + 
  geom_histogram(aes(y = after_stat(density), color = logMgroup_char), fill = 'white', bins = 8) + 
  geom_density(aes(fill = logMgroup_char, color = logMgroup_char), alpha = 0.5) +
  scale_color_brewer(palette = "Accent") + 
  scale_fill_brewer(palette = "Accent") + 
  facet_grid(. ~ logMgroup_char, scales = "free_x") + 
  xlab(label_logMgroup) + 
  geom_text(data = label_df_density,
            aes(x = 12.3, y = 1.2, label = label),
            parse = TRUE, hjust = 0, size = 5) +
  theme_Publication() + 
  theme(
    legend.position = "none",
    strip.text.x = element_text(face = "plain"),
    strip.text.y = element_text(face = "plain")
  )

p2

p1 / p2  

# ggsave(path = wdfigs,
#        filename = paste0("densidades_IHM.pdf"),
#        device = cairo_pdf, width = width_figs, height = height_figs,
#        units = "in", dpi = 600)
