#####################################################################################
# pcstable with nets.
#####################################################################################
rm(list = ls())
# install.packages(pkgs = "/home/catharina/Documents/Installations/modified_4.3.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")

initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
gammas <- c(25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0.75,0.5,0.25,0.1)
runs <- 5

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  
  assign(paste0("pcstable2_withnet_",i),iteration_pcstable_2_withnet(data = data, gammas = gammas))
  save(list = paste0("pcstable2_withnet_",i),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_withnet_",i,".rda"))
  
}

#####################################################################################
# gs with nets.
#####################################################################################
for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  
  assign(paste0("gs_withnet_",i),iteration_gs_2_withnet(data = data, gammas = gammas))
  save(list = paste0("gs_withnet_",i),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gs_withnet_",i,".rda"))
  
}