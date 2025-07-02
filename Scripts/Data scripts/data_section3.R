# Criando inputModel files
wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"

#### Para zmax = 0.03

## Definir qual tabela assignment ----
zmax    <- 0.03
Rlim    <- 2.5
Ma      <- 12.3

## Definindo input e output files ----

input_file <- paste0("assign_zmax", zmax, "_Ma", Ma, 
                     "_clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18")

## Lendo os dados ----
df  <- fread(paste0(wddata, "Assignment2groups/clean/", input_file, ".csv"))

## Completeza em massa estelar
zlim_logMstar <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
ngals_bin     <- zlim_logMstar$Ngals_bin[which(zlim_logMstar$z_max_seq == zmax)]

df_cutmass    <- subset(df, df$logMstar >= logMstar_min)
rm(zlim_logMstar)

# Além disso, corte em massa estelar logMstar < 10.5
df_cutmass <- subset(df_cutmass, df_cutmass$logMstar < 10.5)

# Trabalhando só com as satélites agora:
df_cutmass_sat <- subset(df_cutmass, df_cutmass$type == "Satellite")

# quantas satélites e quantos grupos no final:
table(df_cutmass_sat$type)
length(unique(df_cutmass_sat$groupID))

# propriedades:
table(df_cutmass_sat$SF)
table(df_cutmass_sat$TType_label)

## Verificando valores infinitos e NA ----
names(df_cutmass)[sapply(df_cutmass, function(x) any(is.na(x)))]
names(df_cutmass)[sapply(df_cutmass, function(x) any(is.infinite(x)))]

## Salvando tabelas ----
output_file         <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE.csv")
output_file_cutMass <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min",
                              round(logMstar_min, digits = 2),".csv")

write.csv(df, paste0(wddata, "inputModel/", output_file), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", output_file_cutMass) , row.names = F)

# ------------------------------------------------------------------------------------

#### Para zmax = 0.10

## Definir qual tabela assignment ----
zmax    <- 0.10
Rlim    <- 2.5
Ma      <- 12.3

## Definindo input e output files ----
input_file <- paste0("assign_zmax", zmax, "_Ma", Ma, 
                     "_clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18")

## Lendo os dados ----
df  <- fread(paste0(wddata, "Assignment2groups/clean/", input_file, ".csv"))

## Completeza em massa estelar ----
zlim_logMstar <- fread("~/Work/Research/Astronomy/Data/zlim_logMstar/logM_zlim_MAITE.csv")
logMstar_min  <- zlim_logMstar$logM[which(zlim_logMstar$z_max_seq == zmax)]
ngals_bin     <- zlim_logMstar$Ngals_bin[which(zlim_logMstar$z_max_seq == zmax)]

df_cutmass    <- subset(df, df$logMstar >= logMstar_min)
rm(zlim_logMstar)

# Trabalhando só com as satélites agora:
df_cutmass_sat <- subset(df_cutmass, df_cutmass$type == "Satellite")

# quantas satélites e quantos grupos no final:
table(df_cutmass_sat$type)
length(unique(df_cutmass_sat$groupID))

# propriedades:
table(df_cutmass_sat$SF)
table(df_cutmass_sat$TType_label)

## Verificando valores infinitos e NA ----
names(df_cutmass)[sapply(df_cutmass, function(x) any(is.na(x)))]
names(df_cutmass)[sapply(df_cutmass, function(x) any(is.infinite(x)))]

## Salvando tabelas ----
output_file         <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE.csv")
output_file_cutMass <- paste0("inputdata_zmax", zmax, "_Rlim", Rlim, "_Ma", Ma, "_flag_good==1_MANGLE_logMstar_min",
                              round(logMstar_min, digits = 2),".csv")

write.csv(df, paste0(wddata, "inputModel/", output_file), row.names = F)
write.csv(df_cutmass, paste0(wddata, "inputModel/", output_file_cutMass) , row.names = F)
