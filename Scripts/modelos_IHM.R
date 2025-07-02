## Diretórios ----
wdcode <- "Scripts/"
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"
wdfigs <- "~/Work/Research/Astronomy/Projects/envquenching/Figures/"

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

# Códigos extras ----
source("Scripts/Themes/my_theme.R")

# Opções ----
options(scipen = 999)

critval   <- 1.96 
TType_lim <- 0

width_figs  <- 11
height_figs <- 7

# Dados ----
input_data <- paste0("inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

data.s <- subset(df, df$type == "Satellite")
data.c <- subset(df, df$type == "Central")

# nr de grupos
length(unique(data.s$groupID))

data.s$SF_char <- data.s$SF
data.s$SF <- ifelse(data.s$SF == "Star-forming", 1, 0)

fit_fSFG  <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

data.s$LT <- ifelse(data.s$TType >= TType_lim, 1, 0)
data.s$LT <- as.factor(data.s$LT)

fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

# Categorização de logMgroup
data.s$logMgroup_char <- ifelse(data.s$logMgroup < 14, "< 14", "≥ 14")
data.s$logMgroup_char <- factor(data.s$logMgroup_char, levels = c("< 14", "≥ 14"))

# Ajuste dos modelos com logMgroup incluído
fit_fSFG <- glm(SF ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial("logit"), data = data.s)
fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir + logMgroup, family = binomial("logit"), data = data.s)

# Medidas do modelo:

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

# ------------- FIGURA:

# Define os limites dos bins de dispersão em km/s
bins_velDisp <- c(30, 100, 150, 10^(max(data.s$logvelDisp_e)))
# Converte para log10(velDisp) para comparar com logvelDisp_e
vetor_prob <- log10(bins_velDisp)

# Define domínio fixo
range_logRproj <- range(data.s$logRproj_rvir)
bins_logRproj_rvir <- seq(range_logRproj[1], range_logRproj[2], by = 0.01)

modelos <- bins_total <- data.frame()

# Loop sobre dispersão + logMgroup
for (a in 1:(length(vetor_prob)-1)) {
  lw <- vetor_prob[a]
  up <- vetor_prob[a+1]
  
  faixa <- paste0("[", round(10^lw), ", ", round(10^up), ifelse(up < log10(200), ")", "]"))
  data_painel <- subset(data.s, logvelDisp_e >= lw & logvelDisp_e < up)
  
  for (grupo in levels(data.s$logMgroup_char)) {
    data_sub <- subset(data_painel, logMgroup_char == grupo)
    if (nrow(data_sub) == 0) next
    
    modelo <- data.frame(
      logvelDisp_e = median(data_sub$logvelDisp_e),
      logRproj_rvir = bins_logRproj_rvir,
      logMgroup = median(data_sub$logMgroup)
    )
    
    pred_fSFG <- predict(fit_fSFG, type = "link", se.fit = TRUE, newdata = modelo)
    pred_fLTG <- predict(fit_fLTG, type = "link", se.fit = TRUE, newdata = modelo)
    
    modelo$pred_fSFG_prob <- fit_fSFG$family$linkinv(pred_fSFG$fit)
    modelo$prob_fSFG_upr  <- fit_fSFG$family$linkinv(pred_fSFG$fit + critval * pred_fSFG$se.fit)
    modelo$prob_fSFG_lwr  <- fit_fSFG$family$linkinv(pred_fSFG$fit - critval * pred_fSFG$se.fit)
    
    modelo$pred_fLTG_prob <- fit_fLTG$family$linkinv(pred_fLTG$fit)
    modelo$prob_fLTG_upr  <- fit_fLTG$family$linkinv(pred_fLTG$fit + critval * pred_fLTG$se.fit)
    modelo$prob_fLTG_lwr  <- fit_fLTG$family$linkinv(pred_fLTG$fit - critval * pred_fLTG$se.fit)
    
    modelo$painel_logvelDisp_e <- faixa
    modelo$logMgroup_char <- grupo
    
    modelos <- rbind(modelos, modelo)
    
    aux <- quantile(data_sub$logRproj_rvir, probs = seq(0, 1, by = 0.2))
    
    for (i in 1:(length(aux)-1)) {
      bin <- data_sub$logRproj_rvir >= aux[i] & data_sub$logRproj_rvir < aux[i+1]
      bin_data <- data_sub[bin, ]
      if (nrow(bin_data) == 0) next
      
      meio_bin <- median(c(aux[i], aux[i+1]))
      ngal <- nrow(bin_data)
      
      fSFG <- sum(bin_data$SF == 1) / ngal
      fLTG <- sum(bin_data$LT == 1) / ngal
      
      fSFG_up <- binom.confint(sum(bin_data$SF == 1), ngal, conf.level = 0.95, method = "bayes")$upper
      fSFG_lw <- binom.confint(sum(bin_data$SF == 1), ngal, conf.level = 0.95, method = "bayes")$lower
      fLTG_up <- binom.confint(sum(bin_data$LT == 1), ngal, conf.level = 0.95, method = "bayes")$upper
      fLTG_lw <- binom.confint(sum(bin_data$LT == 1), ngal, conf.level = 0.95, method = "bayes")$lower
      
      bins <- data.frame(
        painel_logvelDisp_e = faixa,
        meio_bin = meio_bin,
        fSFG = fSFG, fSFG_up = fSFG_up, fSFG_lw = fSFG_lw,
        fLTG = fLTG, fLTG_up = fLTG_up, fLTG_lw = fLTG_lw,
        logMgroup_char = grupo,
        ngal = ngal
      )
      
      bins_total <- rbind(bins_total, bins)
    }
  }
}

# Define a ordem desejada dos rótulos
ordem_bins <- c(
  "[30, 100)", 
  "[100, 150)", 
  paste0("[150, ", round(10^vetor_prob[length(vetor_prob)]), "]")
)

# Converte fatores
modelos$painel_logvelDisp_e     <- factor(modelos$painel_logvelDisp_e, levels = ordem_bins)
bins_total$painel_logvelDisp_e  <- factor(bins_total$painel_logvelDisp_e, levels = ordem_bins)
modelos$logMgroup_char <- factor(modelos$logMgroup_char, levels = c("< 14", "≥ 14"))
bins_total$logMgroup_char <- factor(bins_total$logMgroup_char, levels = c("< 14", "≥ 14"))

# 1. Calcula os valores extremos diretamente de `modelos`
extremos_df <- modelos %>%
  group_by(painel_logvelDisp_e, logMgroup_char) %>%
  summarise(
    x_min = min(logRproj_rvir),
    x_max = max(logRproj_rvir),
    
    pred_fSFG_prob_min = pred_fSFG_prob[which.min(logRproj_rvir)],
    pred_fSFG_prob_max = pred_fSFG_prob[which.max(logRproj_rvir)],
    
    pred_fLTG_prob_min = pred_fLTG_prob[which.min(logRproj_rvir)],
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
    painel_full_id = paste0(painel_logvelDisp_e, "___", logMgroup_char),
    
    vjust_label = case_when(
      painel_full_id == "[30, 60)___< 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "[30, 60)___< 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "[30, 60)___< 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.1,
      painel_full_id == "[30, 60)___< 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "[60, 80)___< 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "[60, 80)___< 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "[60, 80)___< 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "[60, 80)___< 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "[80, 200]___< 14"  & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "[80, 200]___< 14"  & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "[80, 200]___< 14"  & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "[80, 200]___< 14"  & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "[30, 60)___≥ 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.2,
      painel_full_id == "[30, 60)___≥ 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.4,
      painel_full_id == "[30, 60)___≥ 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.8,
      painel_full_id == "[30, 60)___≥ 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.5,
      
      painel_full_id == "[60, 80)___≥ 14"   & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "[60, 80)___≥ 14"   & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "[60, 80)___≥ 14"   & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "[60, 80)___≥ 14"   & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      painel_full_id == "[80, 200]___≥ 14"  & modelo == "fSFG" & tipo_ponto == "min" ~ -1.0,
      painel_full_id == "[80, 200]___≥ 14"  & modelo == "fLTG" & tipo_ponto == "min" ~  1.5,
      painel_full_id == "[80, 200]___≥ 14"  & modelo == "fSFG" & tipo_ponto == "max" ~ -0.7,
      painel_full_id == "[80, 200]___≥ 14"  & modelo == "fLTG" & tipo_ponto == "max" ~  1.6,
      
      TRUE ~ 0
    ),
    
    x_label = ifelse(tipo_ponto == "min", x - x_offset, x + x_offset)
  )

# 3. Cria padding_df para forçar expansão por painel com pontos invisíveis
limites_x <- modelos %>%
  group_by(painel_logvelDisp_e, logMgroup_char) %>%
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
  select(painel_logvelDisp_e, logMgroup_char, x, y, modelo)


label_df <- data.frame(
  x = -0.5,
  y = 0.22,
  label = 'italic(LI-M["★"]~sample)',
  painel_logvelDisp_e = modelos$painel_logvelDisp_e[1]
)

facet_labels_logMgroup <- c(
  "< 14" = "log[10](M[h]/h^-1*M['\\u2609']) < 14",
  "≥ 14" = "log[10](M[h]/h^-1*M['\\u2609']) >= 14"
)

ggplot() + 
  # Curvas para fSFG
  geom_line(data = modelos, 
            aes(x = logRproj_rvir, y = pred_fSFG_prob, 
                color = "fSFG", 
                group = interaction(painel_logvelDisp_e, logMgroup_char)), 
            linewidth = 1) + 
  geom_ribbon(data = modelos, 
              aes(x = logRproj_rvir, ymin = prob_fSFG_lwr, ymax = prob_fSFG_upr, 
                  fill = "fSFG", 
                  group = interaction(painel_logvelDisp_e, logMgroup_char)), 
              alpha = 0.4, show.legend = FALSE) + 
  
  # Curvas para fLTG
  geom_line(data = modelos, 
            aes(x = logRproj_rvir, y = pred_fLTG_prob, 
                color = "fLTG", 
                group = interaction(painel_logvelDisp_e, logMgroup_char)), 
            linewidth = 1) + 
  geom_ribbon(data = modelos, 
              aes(x = logRproj_rvir, ymin = prob_fLTG_lwr, ymax = prob_fLTG_upr, 
                  fill = "fLTG", 
                  group = interaction(painel_logvelDisp_e, logMgroup_char)), 
              alpha = 0.4, show.legend = FALSE) + 
  
  # Pontos observacionais para fSFG
  geom_point(data = bins_total, 
             aes(x = meio_bin, y = fSFG, 
                 color = "fSFG", 
                 group = interaction(painel_logvelDisp_e, logMgroup_char)), 
             size = 2) +
  geom_errorbar(data = bins_total, 
                aes(x = meio_bin, ymin = fSFG_lw, ymax = fSFG_up, 
                    color = "fSFG", 
                    group = interaction(painel_logvelDisp_e, logMgroup_char)), 
                width = 0.08) +
  
  # Pontos observacionais para fLTG
  geom_point(data = bins_total, 
             aes(x = meio_bin, y = fLTG, 
                 color = "fLTG", 
                 group = interaction(painel_logvelDisp_e, logMgroup_char)), 
             size = 2) +
  geom_errorbar(data = bins_total, 
                aes(x = meio_bin, ymin = fLTG_lw, ymax = fLTG_up, 
                    color = "fLTG", 
                    group = interaction(painel_logvelDisp_e, logMgroup_char)), 
                width = 0.08) +
  
  
  # Pontos nos extremos das curvas
  geom_point(data = extremos_df, 
             aes(x = x, y = y, color = modelo, group = painel_logvelDisp_e),
             size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  # Rótulos percentuais nos extremos
  geom_text(data = extremos_df,
            aes(x = x_label, y = y, label = percent(y, accuracy = 1), 
                color = modelo, group = painel_logvelDisp_e, 
                vjust = vjust_label),
            size = 4, show.legend = FALSE) +
  
  # Número de galáxias por bin
  geom_text(data = bins_total,
            aes(x = meio_bin, y = 1.2, label = ngal, group = interaction(painel_logvelDisp_e, logMgroup_char)),
            inherit.aes = FALSE,
            color = "black", size = 3,
            angle = -90) +
  
  # geom_text(data = bins_total,
  #           aes(x = meio_bin, y = 1.1, label = ngal, group = interaction(painel_logvelDisp_e, logMgroup_char)),
  #           inherit.aes = FALSE,
  #           color = "#7FD23F", size = 3.5) + 
  
  #facet_grid(logMgroup_char ~ painel_logvelDisp_e, scales = "free_x") + 
  facet_grid(logMgroup_char ~ painel_logvelDisp_e, 
             scales = "free_x",
             labeller = labeller(
               logMgroup_char = as_labeller(facet_labels_logMgroup, label_parsed)
             )) + 
  
  # Padding invisível para forçar expansão do eixo X em cada facet
  geom_blank(data = padding_df, aes(x = x, group = interaction(painel_logvelDisp_e, logMgroup_char))) +
  
  # Eixos e cores
  ylab("Fraction") + 
  xlab(expression(log[10](R[proj]/r[vir]))) + 
  
  scale_color_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F"),
                     name = "",
                     labels = c("fSFG" = expression(f[SFG]), "fLTG" = expression(f[LTG]))) +
  scale_fill_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F")) +
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.23)) +
  scale_x_continuous(
    sec.axis = sec_axis(~ ., name = expression(sigma[e] ~ "(km/s)"), breaks = NULL)
  ) +
  
  # Tema e legendas
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
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
  )

ggsave(path = wdfigs,
       filename = paste0("logistic_IHM.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs, 
       units = "in", dpi = 600)
