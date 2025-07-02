library(data.table)
library(dplyr)

clean <- fread("~/Work/Research/Astronomy/Data/environmental-quenching-data/Assignment2groups/clean/zmax_0.1/assignment2groups_zmax0.1_Rlim2.5_Ma12.3.csv")
todos <- fread("~/Work/Research/Astronomy/Data/environmental-quenching-data/Assignment2groups/todos/zmax_0.1/assignment2groups_zmax0.1_Rlim2.5_Ma12.3.csv")


# Ver linhas de df1 que NÃO estão em df2
anti <- anti_join(todos, clean, by = "igal")

# uma galáxia que não está no clean
anti[5,]

# grupo dela
anti$groupID[1]

# procurar esse grupo no clean
grupo_clean <- clean[which(clean$groupID == anti$groupID[10]),]
grupo_todos <- todos[which(todos$groupID == anti$groupID[10]),]

# para uma galáxia que está nos dois, mudou algo?
igals <- grupo_todos$igal[which(grupo_todos$igal %in% grupo_clean$igal)]

for (i in 1:length(igals)) {
  print(identical(grupo_todos[which(grupo_todos$igal == igals[i]),], grupo_clean[which(grupo_clean$igal == igals[i]),]))
}
