# Diretórios:
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"

# Parâmetros gerais:
Rlim <- 2.5
Ma   <- 12.3

# Arquivos:
zlim_logMstar     <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")

input_clean_file  <- "clean_SDSS_DR18_Legacy_MGS_QSO_kcorr+GSWLC-X2_MGS+Simard11+DS18"
df_clean          <- fread(paste0(wddata, input_clean_file, ".csv"))

### Sample LI-Mass -----
zmax <- 0.03

# Ler arquivo assignment:
input_assing_file <- paste0("assignment2groups_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma)
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


## Verificação de NAs e valores infinitos ----
names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

## Completeza em massa estelar
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
ngals_bin     <- zlim_logMstar$Ngals_bin[which(zlim_logMstar$z_max_seq == zmax)]

df_cutmass    <- subset(df, df$logMstar >= logMstar_min)
# Além disso, corte em massa estelar logMstar < 10.5
df_cutmass <- subset(df_cutmass, df_cutmass$logMstar < 10.5)


## Salvando tabelas ----
output_file         <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE.csv")
output_file_cutMass <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min",
                              round(logMstar_min, digits = 2),".csv")

write.csv(df, paste0(wddata, "inputModel/", output_file), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", output_file_cutMass) , row.names = F)

## Salvando tabela ----
write.csv(df, paste0(wddata, "Assignment2groups/clean/", output_file, ".csv"), row.names = F)

#### Para zmax = 0.10

## Definição de parâmetros ----
zmax <- 0.10
Rlim <- 2.5
Ma   <- 12.3

## Definição de arquivos ----
input_clean_file  <- "clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18"
input_assing_file <- paste0("assign_zmax", zmax, "_Ma", Ma, "_", input_clean_file, ".csv")

output_file <- paste0(input_assing_file, "__", input_clean_file)

## Leitura de dados ----
df_clean  <- fread(paste0(wddata, input_clean_file, ".csv"))
df_assign <- fread(paste0(wddata, "Assignment2groups/", input_assing_file))

## Unindo tabelas ----
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

## Conversão de variáveis ----
lapply(df, class)
colunas_para_converter <- c("Class", "bptClass", "flag_good", "flag_sed", "flag_uv", "flag_midir", "type", "SF", "TType_label", "morph_char", "groupID")
df <- df %>% mutate(across(all_of(colunas_para_converter), factor))

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

## Verificação de NAs e valores infinitos ----
names(df)[sapply(df, function(x) any(is.na(x)))]
names(df)[sapply(df, function(x) any(is.infinite(x)))]

## Salvando tabela ----
output_file <- paste0(input_assing_file, "__", input_clean_file)
write.csv(df, paste0(wddata, "Assignment2groups/clean/", output_file, ".csv"), row.names = F)
