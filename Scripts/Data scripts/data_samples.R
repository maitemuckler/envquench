library(data.table)
library(dplyr)

# Diretórios:
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"

# Parâmetros gerais:
Rlim <- 2.5
Ma   <- 12.3

# Arquivos:
input_clean_file  <- "clean_SDSS_DR18_Legacy_MGS_QSO_kcorr+GSWLC-X2_MGS+Simard11+DS18"

df_clean          <- fread(paste0(wddata, input_clean_file, ".csv"))
zlim_logMstar     <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")

### Sample LI-Mass -----
zmax <- 0.03

# Ler arquivo assignment:
input_assing_file <- paste0("assignment2groups_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good")
df_assign         <- fread(paste0(wddata, "Assignment2groups/", input_assing_file, ".csv"))

# Unindo tabelas
df <- merge(df_assign, df_clean, by = "igal")

## Removendo colunas repetidas 
colunas_para_remover <- grep("\\.y$", names(df), value = TRUE)
df <- df %>% select(-all_of(colunas_para_remover))

# Renomeando colunas repetidas 
colunas_para_renomear <- grep("\\.x$", names(df), value = TRUE)
colnames(df)[which(colnames(df) %in% colunas_para_renomear)] <- gsub("\\.x$", "", colunas_para_renomear)

# Criando variáveis pós-assignment ----
df$logRproj_rvir <- log10(df$Rproj_rvir)
df$type <- factor(ifelse(df$Rproj_rvir > 0, "Satellite", "Central"))

# Conversão de variáveis
lapply(df, class)
colunas_para_converter <- c("type", "groupID")
df <- df %>% mutate(across(all_of(colunas_para_converter), factor))

# Verificação de NAs e valores infinitos 
names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

# Completeza em massa estelar
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
ngals_bin     <- zlim_logMstar$Ngals_bin[which(zlim_logMstar$z_max_seq == zmax)]

df_cutmass <- subset(df, df$logMstar >= logMstar_min)
df_cutmass <- subset(df_cutmass, df_cutmass$logMstar < 10.5) # Além disso, corte em massa estelar logMstar < 10.5

# Salvando tabelas 
output_file         <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, ".csv")
output_file_cutMass <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_logMstar_min", round(logMstar_min, digits = 2),".csv")

write.csv(df, paste0(wddata, "inputModel/", output_file), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", output_file_cutMass) , row.names = F)

### Sample IH-Mass -----
zmax <- 0.10

# Ler arquivo assignment:
input_assing_file <- paste0("assignment2groups_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good")
df_assign         <- fread(paste0(wddata, "Assignment2groups/", input_assing_file, ".csv"))

# Unindo tabelas
df <- merge(df_assign, df_clean, by = "igal")

## Removendo colunas repetidas ----
colunas_para_remover <- grep("\\.y$", names(df), value = TRUE)
df <- df %>% select(-all_of(colunas_para_remover))

## Renomeando colunas repetidas ----
colunas_para_renomear <- grep("\\.x$", names(df), value = TRUE)
colnames(df)[which(colnames(df) %in% colunas_para_renomear)] <- gsub("\\.x$", "", colunas_para_renomear)

## Criando variáveis pós-assignment ----
df$logRproj_rvir <- log10(df$Rproj_rvir)
df$type <- factor(ifelse(df$Rproj_rvir > 0, "Satellite", "Central"))

## Conversão de variáveis
lapply(df, class)
colunas_para_converter <- c("type", "groupID")
df <- df %>% mutate(across(all_of(colunas_para_converter), factor))

# Verificação de NAs e valores infinitos 
names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

# Completeza em massa estelar
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
ngals_bin     <- zlim_logMstar$Ngals_bin[which(zlim_logMstar$z_max_seq == zmax)]

df_cutmass <- subset(df, df$logMstar >= logMstar_min)

# Salvando tabelas 
output_file         <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, ".csv")
output_file_cutMass <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_logMstar_min", round(logMstar_min, digits = 2),".csv")

write.csv(df, paste0(wddata, "inputModel/", output_file), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", output_file_cutMass) , row.names = F)