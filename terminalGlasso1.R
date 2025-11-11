rm(list = ls())
library(glasso)
library(igraph)
library(transformeR)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)

covmatrix<- cov(data, method = "spearman")
cormatrix <- cor(data, method = "spearman")
covRMSmatrix <- cov(dataRMS, method ="spearman")

rho <- c(0.01,0.05,0.063,0.07,0.09,0.12,0.13,0.165,0.17,0.20,0.225,0.25)
rhochar <- c(01,05,063,07,09,12,13,165,17,20,225,25)

rhoplus <- c(c(0.02,0.03,0.04,0.27,0.28,0.30,0.32,0.34,0.38,0.42,0.46,0.50),seq(0.52,0.98,0.02))
rhocharplus <- c("02","03","04","27","28","30","32","34","38","42","46","50","52", "54", "56", "58", "60", "62","64", "66" ,"68", "70", "72",
                 "74","76", "78", "80", "82", "84", "86", "88", "90", "92", "94","96", "98")
# names <- c()
# for(i in 1:length(rho)){
#   names[i] <- paste0("glasso",rhochar[i])
# }
# namesRMS <- c()
# for(i in 1:length(rho)){
#   names[i] <- paste0("glasso",rho[i])
# }
# 
###########################################################################
# normal glasso
###########################################################################
i <- length(rho)
for(i in 1:length(rho)){
  p <- glasso(cormatrix, rho = rho[i])
  filename <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso",rhochar[i],".rda")
  assign(paste0("glasso",rhochar[i]), p)
  save(list = paste0("glasso",rhochar[i]), file = filename)
}

# extra terms
i <- length(rhoplus)
for(i in 1:length(rhoplus)){
  p <- glasso(cormatrix, rho = rhoplus[i])
  filename <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/cornormal/glasso",rhocharplus[i],".rda")
  assign(paste0("glasso",rhocharplus[i]), p)
  save(list = paste0("glasso",rhocharplus[i]), file = filename)
}

################################################################################
# different aproaches
################################################################################


for(i in 1:length(rho)){
  p <- glasso(covmatrix, rho = rho[i])
  filename <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glassoCOV",rho[i],".rda")
  assign(paste0("glassoCOV",rho[i]), p)
  save(list = paste0("glassoCOV",rho[i]), file = filename)
}

for(i in 1:length(rho)){
  p <- glasso(covRMSmatrix, rho = rho[i])
  filename <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glassoCOVrms_",rho[i],".rda")
  assign(paste0("glassoCOVrms_",rho[i]), p)
  save(list = paste0("glassoCOVrms_",rho[i]), file = filename)
}
