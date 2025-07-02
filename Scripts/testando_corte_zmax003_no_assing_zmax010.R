wddata <- "~/Work/Research/Astronomy/Data/environmental-quenching-data/"

# Arquivo assignment 0.03
input_file <- paste0("assignment2groups_zmax0.03_Rlim2.5_Ma12.3_flag_good__clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18")
assign_003 <- fread(paste0(wddata, "Assignment2groups/clean/", input_file, ".csv"))

# Arquivo assignment 0.10
input_file <- paste0("assignment2groups_zmax0.1_Rlim2.5_Ma12.3_flag_good__clean_SDSS_DR18_Legacy_MGS_QSO_localDensity+GSWLC+Simard11+DS18")
assign_010 <- fread(paste0(wddata, "Assignment2groups/clean/", input_file, ".csv"))

# Arquivo assignment 0.10, cortado em zmax = 0.03
assign_010_003 <- subset(assign_010_003, assign_010_003$z <= 0.03)

# Qual a diferença?
anti <- anti_join(assign_010_003, assign_003, by = "igal")

# porque 197795 não está na tabela de assigment 0.03?
which(assign_003$igal == 197795)

# ela está no 0.10
teste <- assign_010[which(assign_010$igal == 197795),]
# ela está no grupo 1 no assignment 0.10! mas o grupo 1 não faz parte do assignment 0.03!

# no assignment 0.03, o grupo 1 sai quando retiramos os grupos perto das bordas do survey
# dai, não há grupo para essa galáxia o qual a distancia dela do centro fique menor que 20 rvir. 
# rovavelmente ela está a mais de 20 rvir de qualquer outro grupo disponivel.
