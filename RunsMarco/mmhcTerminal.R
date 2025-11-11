#####################################################################################
# Hill Climbing and Tabu Search. 
#####################################################################################
rm(list = ls())
install.packages(pkgs = "/home/catharina/Documents/Installations/modified_4.3.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")


initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
gammas <- c(5:1)
gammas
runs <- 5
bnlearn

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

for (i in c(3,4,5)){
  data <- datapermutations[[i]]
  
  # assign(paste0("tabu_",i),iteration_tabu(data, 8000))
  # save(list = paste0("tabu_",i),
  #       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,".rda"))
  # 
  # assign(paste0("hc_",i),iteration_hc(data, 8000))
  # save(list = paste0("hc_",i),
  #       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,".rda"))
  # 
  # # assign(paste0("pc_",i),iteration_pc(data = data, correlation = NULL, alpha = alphas))
  # # save(list = paste0("pc_",i), 
  # #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_",i,".rda"))
  # 
  # assign(paste0("mmhc_",i),iteration_mmhc_split(data = data, iterations = 2000))
  # save(list = paste0("mmhc_",i),
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,".rda"))
  
  assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  save(list = paste0("pcstable_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
}

for (i in c(1,4,5)){
  data <- datapermutations[[i]]
  
  # assign(paste0("tabu_",i),iteration_tabu(data, 8000))
  # save(list = paste0("tabu_",i),
  #       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,".rda"))
  # 
  # assign(paste0("hc_",i),iteration_hc(data, 8000))
  # save(list = paste0("hc_",i),
  #       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,".rda"))
  # 
  # # assign(paste0("pc_",i),iteration_pc(data = data, correlation = NULL, alpha = alphas))
  # # save(list = paste0("pc_",i), 
  # #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_",i,".rda"))
  # 
  assign(paste0("mmhc_",i),iteration_mmhc_split(data = data, iterations = 2000))
  save(list = paste0("mmhc_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,".rda"))
  
  # assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstable_",i), 
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
}

runs
