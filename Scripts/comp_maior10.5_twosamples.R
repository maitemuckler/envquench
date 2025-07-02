input_file <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/inputModel/GSWLC/inputdata_zmax003_logMstar_min_9.17.csv"
df         <- fread(input_file)

df003_maior105 <- subset(df, df$lgm_tot_p50 >= 10.5)
df003_maior105 <- subset(df003_maior105, df003_maior105$type == "Satellite")

# amostra zmax = 0.1

input_file <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/inputModel/GSWLC/inputdata_zmax010_logMstar_min_10.5.csv"
df  <- fread(input_file)

df010_maior105 <- subset(df, df$lgm_tot_p50 >= 10.5)
df010_maior105 <- subset(df010_maior105, df010_maior105$type == "Satellite")

set.seed(2222)
igal_unicos    <- setdiff(df010_maior105$igal, df003_maior105$igal)
df010_maior105 <- df010_maior105[which(df010_maior105$igal %in% sample(igal_unicos, size = nrow(df003_maior105))),]

table(df010_maior105_novo$igal %in% df003_maior105$igal) # não há nenhuma galáxia igual nas duas amostras

# Comparação

df003_maior105 <- df003_maior105 %>% select(logvelDisp_e, logRproj_rvir, TType_label)
df003_maior105$sample <- "zmax003"

df010_maior105 <- df010_maior105 %>% select(logvelDisp_e, logRproj_rvir, TType_label)
df010_maior105$sample <- "zmax010"

df <- rbind(df003_maior105, df010_maior105)

library(ggplot2)

ggplot(df, aes(x = logvelDisp_e, color = sample)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  theme_minimal() +
  labs(title = "ECDF das amostras", x = "Valor", y = "ECDF") +
  scale_color_manual(values = c("blue", "red"))

ks.test(
  df$logvelDisp_e[df$sample == "zmax003"],
  df$logvelDisp_e[df$sample == "zmax010"]
)

ggplot(df, aes(x = logRproj_rvir, color = sample)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  theme_minimal() +
  labs(title = "ECDF das amostras", x = "Valor", y = "ECDF") +
  scale_color_manual(values = c("blue", "red"))

ks.test(
  df$logRproj_rvir[df$sample == "zmax003"],
  df$logRproj_rvir[df$sample == "zmax010"]
)

# T-Type

tabela <- table(df$sample, df$TType_label)
chisq.test(tabela)

prop.table(tabela, margin = 1)
