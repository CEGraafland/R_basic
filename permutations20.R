###########################################################################################
# Permutations HC 6 tot en met 25
###########################################################################################
rm(list = ls())
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")

initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
runs <- 20

datapermutations20 <- list()
for (i in 1:runs){
  datapermutations20[[i]] <- initialdata[,sample(ncol(initialdata))]
}
save(datapermutations20, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations20.rda")

permutations20 <- list()
for (j in 1:length(datapermutations20)){
  indsame <- c()
  for (i in 1:ncol(initialdata)){
    int <- which(colnames(initialdata) == colnames(datapermutations20[[j]][i]))
    indsame[i] <- int
  }
  permutations20[[j]] <- indsame
}

save(permutations20, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations20.rda")
