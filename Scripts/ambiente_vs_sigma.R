setwd("~/Work/Research/Astronomy/Projects/envquenching/")

library(ggplot2)
library(dplyr)
library(scales)

# Setting:
TType_lim     <- 0
critval       <- 1.96 
Ma            <- 12.3

LIMass <- read.csv("inputdata_zmax0.03_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min9.17.csv")
IHMass <- read.csv("inputdata_zmax0.1_Rlim2.5_Ma12.3_flag_good==1_MANGLE_logMstar_min10.5.csv")

# Massa do halo

LIMass$logMgroup_char <- cut(LIMass$logMgroup, breaks = c(min(LIMass$logMgroup), 14, max(LIMass$logMgroup)), 
                             include.lowest = TRUE, labels = c("Groups", "Clusters"))
IHMass$logMgroup_char <- cut(IHMass$logMgroup, breaks = c(min(IHMass$logMgroup), 14, max(IHMass$logMgroup)), include.lowest = TRUE,
                             labels = c("Groups", "Clusters"))

# Calcular os quantis por grupo
quantis <- LIMass %>%
  group_by(logMgroup_char) %>%
  summarize(
    q25 = quantile(logvelDisp_e, 0.25, na.rm = TRUE),
    q50 = quantile(logvelDisp_e, 0.50, na.rm = TRUE),
    q75 = quantile(logvelDisp_e, 0.75, na.rm = TRUE)
  ) %>%
  tidyr::pivot_longer(cols = c(q25, q50, q75), names_to = "quantil", values_to = "valor")

# Gráfico com densidade e linhas verticais nos quantis
ggplot(LIMass, aes(x = logvelDisp_e, fill = logMgroup_char)) + 
  geom_density(alpha = 0.7) + 
  facet_grid(. ~ logMgroup_char) +
  geom_vline(data = quantis, aes(xintercept = valor, colour = quantil), linetype = "dashed", linewidth = 0.8, show.legend = TRUE) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  scale_fill_brewer(palette = "Set2") + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "logvelDisp_e", y = "Density", colour = "Quantile") + 
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "bottom")

limiar <- 2

tabela <- LIMass %>%
  mutate(disp_alta = logvelDisp_e > limiar) %>%
  group_by(logMgroup_char) %>%
  summarize(
    total = n(),
    alta = sum(disp_alta, na.rm = TRUE),
    prop_alta = mean(disp_alta, na.rm = TRUE)
  )

tabela

# Se tiver poucos "sucessos", use fisher.test:
fisher.test(
  matrix(c(tabela$alta[1], tabela$total[1] - tabela$alta[1],
           tabela$alta[2], tabela$total[2] - tabela$alta[2]),
         nrow = 2, byrow = TRUE)
)

# Ou, com números maiores, um teste de proporções:
prop.test(x = tabela$alta, n = tabela$total)

# Base de dados com proporções e erro padrão
tabela_plot <- LIMass %>%
  mutate(disp_alta = logvelDisp_e > 2.0) %>%
  group_by(logMgroup_char) %>%
  summarize(
    total = n(),
    alta = sum(disp_alta, na.rm = TRUE),
    prop_alta = mean(disp_alta, na.rm = TRUE),
    se = sqrt(prop_alta * (1 - prop_alta) / total)  # erro padrão da proporção
  )

# Gráfico
ggplot(tabela_plot, aes(x = logMgroup_char, y = prop_alta, fill = logMgroup_char)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = prop_alta - 1.96 * se, ymax = prop_alta + 1.96 * se), width = 0.2) +
  labs(x = "", y = "Porportion with logvelDisp_e > 2.0") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set2") + 
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

