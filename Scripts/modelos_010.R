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

width_figs  <- 9
height_figs <- 5

# Dados ----
input_data <- paste0("inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")
df         <- fread(paste0(wddata, "inputModel/", input_data))

data.s <- subset(df, df$type == "Satellite")
data.c <- subset(df, df$type == "Central")

# nr de grupos
length(unique(data.s$groupID))

# Modelo para fSFG ----

data.s$SF_char <- data.s$SF
data.s$SF <- ifelse(data.s$SF == "Star-forming", 1, 0)

fit_fSFG  <- glm(SF ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

# Resumo do modelo
summary(fit_fSFG)

# Intervalo de confiança para os coeficientes
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

# Importância das covariáveis
importance_fSFG <- as.data.frame(varImp(fit_fSFG))
importance_fSFG <- data.frame(overall = importance_fSFG$Overall, names = rownames(importance_fSFG))
importance_fSFG[order(importance_fSFG$overall, decreasing = T),]

# Predição
data.s$pred_fSFG <- predict(fit_fSFG, data.s, type = "response")
optimal_fSFG     <- optimalCutoff(data.s$SF, data.s$pred_fSFG)[1]
misClass_fSFG    <- misClassError(data.s$SF, data.s$pred_fSFG, threshold = optimal_fSFG)
acuracia_fSFG    <- 1 - misClass_fSFG
acuracia_fSFG

# Matriz de confusão
data.s$pred_class_fSFG <- ifelse(data.s$pred_fSFG >= optimal_fSFG, "Star-forming", "Quiescent") # Converter probabilidades em classes
conf_matrix_fSFG       <- confusionMatrix(factor(data.s$pred_class_fSFG), factor(data.s$SF_char))
conf_matrix_fSFG

# Extraindo métricas de desempenho para fSFG
acc_fSFG         <- conf_matrix_fSFG$overall["Accuracy"]
sensitivity_fSFG <- conf_matrix_fSFG$byClass["Sensitivity"]
specificity_fSFG <- conf_matrix_fSFG$byClass["Specificity"]

# Curva ROC para fSFG
roc_fSFG <- roc(data.s$SF, data.s$pred_fSFG)
#plot(roc_fSFG, main="ROC Curve - fSFG", col="blue", lwd=2)
#cat("AUC para fSFG:", auc(roc_fSFG), "\n")

# Modelo para fLTG ----

data.s$LT <- ifelse(data.s$TType >= TType_lim, 1, 0)
data.s$LT <- as.factor(data.s$LT)

fit_fLTG <- glm(LT ~ logvelDisp_e + logRproj_rvir, family = binomial(link = "logit"), data = data.s)

# Resumo do modelo
summary(fit_fLTG)

# Intervalo de confiança para os coeficientes
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

# Importância das covariáveis
importance_fLTG <- as.data.frame(varImp(fit_fLTG))
importance_fLTG <- data.frame(overall = importance_fLTG$Overall, names = rownames(importance_fLTG))
importance_fLTG[order(importance_fLTG$overall, decreasing = T),]

# Predição
data.s$pred_fLTG <- predict(fit_fLTG, data.s, type = "response")
optimal_fLTG     <- optimalCutoff(data.s$LT, data.s$pred_fLTG)[1]
misClass_fLTG    <- misClassError(data.s$LT, data.s$pred_fLTG, threshold = optimal_fLTG)
acuracia_fLTG    <- 1 - misClass_fLTG
acuracia_fLTG

# Matriz de confusão
data.s$pred_class_fLTG <- ifelse(data.s$pred_fLTG >= optimal_fLTG, "Late-type", "Early-type") # Converter probabilidades em classes
conf_matrix_fLTG       <- confusionMatrix(factor(data.s$pred_class_fLTG), factor(data.s$morph_char))
conf_matrix_fLTG

# Extraindo métricas de desempenho para fSFG
acc_fLTG         <- conf_matrix_fLTG$overall["Accuracy"]
sensitivity_fLTG <- conf_matrix_fLTG$byClass["Sensitivity"]
specificity_fLTG <- conf_matrix_fLTG$byClass["Specificity"]

# Curva ROC para fSFG
roc_fLTG <- roc(data.s$LT, data.s$pred_fLTG)
#plot(roc_fLTG, main="ROC Curve - fLTG", col="blue", lwd=2)
#cat("AUC para fLTG:", auc(roc_fLTG), "\n")
# --------------------------------------------------------------------------------------------------------

## FIGURAS LOGÍSTICAS:

# Definição dos quantis
vetor_prob <- quantile(data.s$logvelDisp_e, probs = seq(0, 1, length.out = 4))

velDisp_min <- 10^unname(vetor_prob[1])
velDisp1    <- 10^unname(vetor_prob[2])
velDisp2    <- 10^unname(vetor_prob[3])
velDisp_max <- 10^unname(vetor_prob[4])

vetor_prob <- log10(c(velDisp_min, velDisp1, velDisp2, velDisp_max))

modelos <- data.frame(
  logvelDisp_e = numeric(),
  logRproj_rvir = numeric(), 
  pred_fSFG_prob = numeric(),
  prob_fSFG_upr = numeric(),
  prob_fSFG_lwr = numeric(),
  pred_fLTG_prob = numeric(),
  prob_fLTG_upr = numeric(), 
  prob_fLTG_lwr = numeric(),
  painel_logvelDisp_e = character(),
  stringsAsFactors = FALSE
)

bins_total <- data.frame(
  painel_logvelDisp_e = character(),
  meio_bin = numeric(),
  fSFG = numeric(), 
  fSFG_up = numeric(), 
  fSFG_lw = numeric(),
  fLTG = numeric(), 
  fLTG_up = numeric(), 
  fLTG_lw = numeric(),
  stringsAsFactors = FALSE
)

# Loop para calcular os modelos e bins
for (a in 1:(length(vetor_prob)-1)) {
  
  quantil_lw_logvelDisp_e <- vetor_prob[a]
  quantil_up_logvelDisp_e <- vetor_prob[a+1]
  
  painel <- data.s$logvelDisp_e >= quantil_lw_logvelDisp_e & 
    data.s$logvelDisp_e < quantil_up_logvelDisp_e 
  
  data.s_painel <- data.s[painel,]
  ngal_painel <- nrow(data.s_painel)
  
  painel_logvelDisp_e <- paste0("[",
                                round(10^quantil_lw_logvelDisp_e, digits = 0), ", ",
                                round(10^quantil_up_logvelDisp_e, digits = 0), ")")
  
  # Criando os modelos
  bins_logRproj_rvir <- seq(min(data.s_painel$logRproj_rvir), max(data.s_painel$logRproj_rvir), by = 0.01)
  
  modelo <- data.frame(logvelDisp_e = median(data.s_painel$logvelDisp_e), 
                       logRproj_rvir = bins_logRproj_rvir)
  
  pred_fSFG <- predict(fit_fSFG, type = "link", se.fit = TRUE, newdata = modelo)
  pred_fLTG <- predict(fit_fLTG, type = "link", se.fit = TRUE, newdata = modelo)
  
  modelo$pred_fSFG_prob <- fit_fSFG$family$linkinv(pred_fSFG$fit)
  modelo$prob_fSFG_upr  <- fit_fSFG$family$linkinv(pred_fSFG$fit + (critval * pred_fSFG$se.fit))
  modelo$prob_fSFG_lwr  <- fit_fSFG$family$linkinv(pred_fSFG$fit - (critval * pred_fSFG$se.fit))
  
  modelo$pred_fLTG_prob <- fit_fLTG$family$linkinv(pred_fLTG$fit)
  modelo$prob_fLTG_upr  <- fit_fLTG$family$linkinv(pred_fLTG$fit + (critval * pred_fLTG$se.fit))
  modelo$prob_fLTG_lwr  <- fit_fLTG$family$linkinv(pred_fLTG$fit - (critval * pred_fLTG$se.fit))
  
  modelo$painel_logvelDisp_e <- painel_logvelDisp_e
  modelos <- rbind(modelos, modelo)
  
  # Criando os bins
  aux <- quantile(data.s_painel$logRproj_rvir, probs = seq(0, 1, by = 0.2))
  
  for (i in 1:(length(aux)-1)) {
    
    bin <- painel & data.s$logRproj_rvir >= aux[i] & data.s$logRproj_rvir < aux[i+1]
    data.s_bin <- data.s[bin,]
    ngal_bin <- nrow(data.s_bin)
    
    if (ngal_bin > 0) {
      meio_bin <- median(c(aux[i], aux[i+1]))
      
      ngal_SF <- sum(data.s_bin$SF == 1)
      ngal_LT <- sum(data.s_bin$LT == 1)
      
      fSFG <- ngal_SF / ngal_bin
      fLTG <- ngal_LT / ngal_bin
      
      fSFG_up <- binom.confint(ngal_SF, ngal_bin, conf.level = 0.95, method = "bayes")$upper
      fSFG_lw <- binom.confint(ngal_SF, ngal_bin, conf.level = 0.95, method = "bayes")$lower
      
      fLTG_up <- binom.confint(ngal_LT, ngal_bin, conf.level = 0.95, method = "bayes")$upper
      fLTG_lw <- binom.confint(ngal_LT, ngal_bin, conf.level = 0.95, method = "bayes")$lower
      
      bins <- data.frame(
        painel_logvelDisp_e = painel_logvelDisp_e,
        meio_bin = meio_bin,
        fSFG = fSFG, fSFG_up = fSFG_up, fSFG_lw = fSFG_lw,
        fLTG = fLTG, fLTG_up = fLTG_up, fLTG_lw = fLTG_lw)
      
      bins_total <- rbind(bins_total, bins)
    }
  }
}

# Probs específicas: Radial variation of the fraction ~ delta fraction

# Criando um data frame vazio para armazenar os resultados
resultados_df <- data.frame(
  logvelDisp_e = numeric(),
  pred_fSFG_prob_min = numeric(),
  pred_fSFG_prob_max = numeric(),
  pred_fSFG_prob_abs = numeric(),
  pred_fLTG_prob_min = numeric(),
  pred_fLTG_prob_max = numeric(),
  pred_fLTG_prob_abs = numeric(),
  stringsAsFactors = FALSE
)

# Obtendo os valores únicos de logvelDisp_e
valores_logvelDisp_e <- unique(modelos$logvelDisp_e)

# Loop sobre os três painéis
for (i in seq_along(valores_logvelDisp_e)) {
  
  painel <- subset(modelos, modelos$logvelDisp_e == valores_logvelDisp_e[i])
  
  # Obtendo os índices dos valores desejados
  min_logRproj <- which.min(painel$logRproj_rvir)
  max_logRproj <- which.max(painel$logRproj_rvir)
  
  pred_fSFG_prob_min = painel$pred_fSFG_prob[min_logRproj]
  pred_fLTG_prob_min = painel$pred_fLTG_prob[min_logRproj]
  
  pred_fSFG_prob_max = painel$pred_fSFG_prob[max_logRproj]
  pred_fLTG_prob_max = painel$pred_fLTG_prob[max_logRproj]
  
  pred_fSFG_prob_abs = abs(pred_fSFG_prob_max - pred_fSFG_prob_min)
  pred_fLTG_prob_abs = abs(pred_fLTG_prob_max - pred_fLTG_prob_min)
  
  # Criando uma nova linha de resultados
  nova_linha <- data.frame(
    logvelDisp_e = valores_logvelDisp_e[i],
    
    pred_fSFG_prob_min = pred_fSFG_prob_min,
    pred_fSFG_prob_max = pred_fSFG_prob_max,
    pred_fSFG_prob_abs = pred_fSFG_prob_abs,
    
    pred_fLTG_prob_min = pred_fLTG_prob_min,
    pred_fLTG_prob_max = pred_fLTG_prob_max,
    pred_fLTG_prob_abs = pred_fLTG_prob_abs
    
  )
  
  # Adicionando a nova linha ao data frame de resultados
  resultados_df <- rbind(resultados_df, nova_linha)
}

# Diferença para fSFG:
round((resultados_df$pred_fSFG_prob_max - resultados_df$pred_fSFG_prob_min) * 100)

# Diferença para fLTG:
round((resultados_df$pred_fLTG_prob_max - resultados_df$pred_fLTG_prob_min) * 100)

# Garantindo que painel_logvelDisp_e é um fator
modelos$painel_logvelDisp_e <- factor(modelos$painel_logvelDisp_e,
                                      levels = c(paste0("[", round(10^vetor_prob)[1], ", ", round(10^vetor_prob)[2], ")"), 
                                                 paste0("[", round(10^vetor_prob)[2], ", ", round(10^vetor_prob)[3], ")"),
                                                 paste0("[", round(10^vetor_prob)[3], ", ", round(10^vetor_prob)[4], ")")))
bins_total$painel_logvelDisp_e <- factor(bins_total$painel_logvelDisp_e,
                                         levels = c(paste0("[", round(10^vetor_prob)[1], ", ", round(10^vetor_prob)[2], ")"), 
                                                    paste0("[", round(10^vetor_prob)[2], ", ", round(10^vetor_prob)[3], ")"),
                                                    paste0("[", round(10^vetor_prob)[3], ", ", round(10^vetor_prob)[4], ")")))

# Plotagem
# 1. Adiciona os rótulos dos painéis ao resultados_df com níveis ordenados
resultados_df <- resultados_df %>%
  mutate(
    painel_logvelDisp_e = case_when(
      logvelDisp_e < log10(velDisp1) ~ paste0("[", round(10^vetor_prob[1]), ", ", round(10^vetor_prob[2]), ")"),
      logvelDisp_e < log10(velDisp2) ~ paste0("[", round(10^vetor_prob[2]), ", ", round(10^vetor_prob[3]), ")"),
      TRUE                           ~ paste0("[", round(10^vetor_prob[3]), ", ", round(10^vetor_prob[4]), ")")
    ),
    painel_logvelDisp_e = as.character(painel_logvelDisp_e)
  )

# 2. Cria extremos_df com valores min/max e rótulos posicionados
extremos_df <- modelos %>%
  group_by(painel_logvelDisp_e) %>%
  summarise(
    x_min = min(logRproj_rvir),
    x_max = max(logRproj_rvir),
    .groups = "drop"
  ) %>%
  left_join(resultados_df, by = "painel_logvelDisp_e") %>%
  pivot_longer(
    cols = c(pred_fSFG_prob_min, pred_fSFG_prob_max, pred_fLTG_prob_min, pred_fLTG_prob_max),
    names_to = "tipo",
    values_to = "y"
  ) %>%
  mutate(
    modelo = ifelse(grepl("SFG", tipo), "fSFG", "fLTG"),
    x = ifelse(grepl("min", tipo), x_min, x_max),
    tipo_ponto = case_when(
      grepl("min", tipo) ~ "min",
      grepl("max", tipo) ~ "max"
    ),
    vjust_label = case_when(
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[1]), ", ", round(10^vetor_prob[2]), ")") & modelo == "fSFG" & tipo_ponto == "min" ~ -0.6,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[1]), ", ", round(10^vetor_prob[2]), ")") & modelo == "fLTG" & tipo_ponto == "min" ~  1.6,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[1]), ", ", round(10^vetor_prob[2]), ")") & modelo == "fSFG" & tipo_ponto == "max" ~ -0.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[1]), ", ", round(10^vetor_prob[2]), ")") & modelo == "fLTG" & tipo_ponto == "max" ~  1.5,
      
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[2]), ", ", round(10^vetor_prob[3]), ")") & modelo == "fSFG" & tipo_ponto == "min" ~ -0.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[2]), ", ", round(10^vetor_prob[3]), ")") & modelo == "fLTG" & tipo_ponto == "min" ~ -0.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[2]), ", ", round(10^vetor_prob[3]), ")") & modelo == "fSFG" & tipo_ponto == "max" ~ -0.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[2]), ", ", round(10^vetor_prob[3]), ")") & modelo == "fLTG" & tipo_ponto == "max" ~ -0.5,
      
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[3]), ", ", round(10^vetor_prob[4]), ")") & modelo == "fSFG" & tipo_ponto == "min" ~ -0.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[3]), ", ", round(10^vetor_prob[4]), ")") & modelo == "fLTG" & tipo_ponto == "min" ~ 1.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[3]), ", ", round(10^vetor_prob[4]), ")") & modelo == "fSFG" & tipo_ponto == "max" ~ -0.5,
      painel_logvelDisp_e == paste0("[", round(10^vetor_prob[3]), ", ", round(10^vetor_prob[4]), ")") & modelo == "fLTG" & tipo_ponto == "max" ~ 1.5,
      
      TRUE ~ 0
    ),
    x_label = case_when(
      tipo_ponto == "min" ~ x - 0.25,
      tipo_ponto == "max" ~ x + 0.3,
      TRUE ~ x
    )
  ) %>%
  mutate(
    painel_logvelDisp_e = factor(painel_logvelDisp_e,
                                 levels = c(paste0("[", round(10^vetor_prob)[1], ", ", round(10^vetor_prob)[2], ")"), 
                                            paste0("[", round(10^vetor_prob)[2], ", ", round(10^vetor_prob)[3], ")"),
                                            paste0("[", round(10^vetor_prob)[3], ", ", round(10^vetor_prob)[4], ")")))
  )

# 3. Cria padding_df para forçar expansão por painel com pontos invisíveis
limites_x <- modelos %>%
  group_by(painel_logvelDisp_e) %>%
  summarise(
    x_min = min(logRproj_rvir),
    x_max = max(logRproj_rvir),
    .groups = "drop"
  ) %>%
  mutate(
    x_min_exp = x_min - 0.3,
    x_max_exp = x_max + 0.35
  )

padding_df <- limites_x %>%
  mutate(
    y = NA,
    modelo = "fSFG"
  ) %>%
  pivot_longer(cols = c(x_min_exp, x_max_exp), names_to = NULL, values_to = "x") %>%
  select(painel_logvelDisp_e, x, y, modelo)


label_df <- data.frame(
  x = -0.65,
  y = 0.22,
  label = 'italic(IH-M["★"]~sample)',
  painel_logvelDisp_e = modelos$painel_logvelDisp_e[1]
)

# 4. Gera o gráfico final
ggplot() + 
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fSFG_prob, color = "fSFG", group = painel_logvelDisp_e), linewidth = 1) + 
  geom_ribbon(data = modelos, aes(ymin = prob_fSFG_lwr, ymax = prob_fSFG_upr, x = logRproj_rvir, fill = "fSFG", group = painel_logvelDisp_e), alpha = 0.5) + 
  
  geom_line(data = modelos, aes(x = logRproj_rvir, y = pred_fLTG_prob, color = "fLTG", group = painel_logvelDisp_e), linewidth = 1) + 
  geom_ribbon(data = modelos, aes(ymin = prob_fLTG_lwr, ymax = prob_fLTG_upr, x = logRproj_rvir, fill = "fLTG", group = painel_logvelDisp_e), alpha = 0.5) + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fSFG, color = "fSFG", group = painel_logvelDisp_e), size = 2) +
  geom_errorbar(data = bins_total, aes(ymin = fSFG_lw, ymax = fSFG_up, x = meio_bin, color = "fSFG", group = painel_logvelDisp_e), width = 0.1) + 
  
  geom_point(data = bins_total, aes(x = meio_bin, y = fLTG, color = "fLTG", group = painel_logvelDisp_e), size = 2) +
  geom_errorbar(data = bins_total, aes(ymin = fLTG_lw, ymax = fLTG_up, x = meio_bin, color = "fLTG", group = painel_logvelDisp_e), width = 0.1) + 
  
  geom_point(data = extremos_df, 
             aes(x = x, y = y, color = modelo, group = painel_logvelDisp_e),
             size = 3, shape = 21, fill = "white", stroke = 1.2) +
  
  geom_text(data = extremos_df,
            aes(x = x_label, y = y, label = percent(y, accuracy = 1), 
                color = modelo, group = painel_logvelDisp_e, 
                vjust = vjust_label),
            size = 4, show.legend = FALSE)+
  
  geom_blank(data = padding_df, aes(x = x, group = painel_logvelDisp_e)) + 
  
  facet_grid(. ~ painel_logvelDisp_e) + 
  
  ylab("Fraction") + 
  xlab(expression(log[10](R[proj]/r[vir]))) + 
  
  scale_color_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F"), 
                     name = "",
                     labels = c("fSFG" = expression(f[SFG]), "fLTG" = expression(f[LTG])),
                     guide = guide_legend(override.aes = list(linewidth = 1.5, size = 4))) + 
  scale_fill_manual(values = c("fSFG" = "#9B5DE5", "fLTG" = "#7FD23F")) + 
  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , expression(paste(sigma[e], " (km/s)")), breaks = NULL, labels = NULL)) + 
  
  geom_text(data = label_df,
            aes(x = x, y = y, label = label),
            parse = TRUE,
            size = 5,
            inherit.aes = FALSE) + 
  
  theme_Publication() + 
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "inside",
        legend.position.inside = c(0.14, 0.11),
        legend.key.size = unit(0.8, "cm"),
        legend.box.background = element_rect(color = "black", linewidth = 1)) + 
  guides(fill = "none")

# Exporta o gráfico
ggsave(path = wdfigs,
       filename = paste0("logistic_IHM.pdf"),
       device = cairo_pdf, width = width_figs, height = height_figs, 
       units = "in", dpi = 600)
