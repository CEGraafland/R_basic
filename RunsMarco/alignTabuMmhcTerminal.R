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
runs <- 5

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2.rda")
mmhcs <- list(mmhc_1$networks[[3]], 
              mmhc_2$networks[[3]], 
              mmhc_3$networks[[3]], 
              mmhc_4$networks[[3]], 
              mmhc_5$networks[[3]])
mmhcspure <- list(mmhc_1, 
                  mmhc_2, 
                  mmhc_3, 
                  mmhc_4, 
                  mmhc_5)
arcsmmhcs <- sapply(mmhcs, narcs)

arcsgs <- c(305,454,307,457,422)
gammasgs <- c(10,6,10,6,7)

tabus <- list(0,
              tabu_2$networks[[1]], 
              tabu_3$networks[[1]],
              tabu_4$networks[[1]],
              tabu_5$networks[[1]])


extraits <- c(0,
              narcs(mmhc_2$networks[[3]])- narcs(tabu_2$networks[[1]]),
              narcs(mmhc_3$networks[[3]])- narcs(tabu_3$networks[[1]]),
              narcs(mmhc_4$networks[[3]])- narcs(tabu_4$networks[[1]]),
              0)

for (i in c(2,3,4)){
  data <- datapermutations[[i]]
  
  assign(paste0("tabu_mmhc_",i),iteration_tabu_sel(data, iterations = extraits[i], select = arcsmmhcs[i], start = tabus[[i]]))
  save(list = paste0("tabu_mmhc_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_",i,".rda"))
  
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
  # 
  # assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstable_",i), 
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
}

for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  
  assign(paste0("tabu_gs_simple",i),tabu(x = data, max.iter = arcsgs[i]))
  save(list = paste0("tabu_gs_simple",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/tabu_gs_simple",i,".rda"))
  
  # assign(paste0("hc_",i),iteration_hc(data, 8000))
  # save(list = paste0("hc_",i),
  #       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,".rda"))
  # 
  # # assign(paste0("pc_",i),iteration_pc(data = data, correlation = NULL, alpha = alphas))
  # # save(list = paste0("pc_",i), 
  # #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_",i,".rda"))
  # 
  
  assign(paste0("mmhc_gs_simple",i),iteration_mmhc_split_simple(data = data, iterations = arcsgs[i], begin = mmhcspure[[i]]$begin))
  save(list = paste0("mmhc_gs_simple",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/mmhc_gs_simple",i,".rda"))
  # 
  # assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstable_",i), 
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
 
  assign(paste0("gs_simple",i),gs(x = data, test = "bic-gt", B = gammasgs[i]))
  save(list = paste0("gs_simple",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/gs_simple",i,".rda"))
  
}



load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple5.rda")

tabusimples <- list(tabu_mmhc_simple1,
              tabu_mmhc_simple2, 
              tabu_mmhc_simple3,
              tabu_mmhc_simple4,
              tabu_mmhc_simple5)

for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  
  assign(paste0("tabu_simple_1600_",i),tabu(x = data, max.iter = (1600 - narcs(tabusimples[[i]])), start = tabusimples[[i]]))
  save(list = paste0("tabu_simple_1600_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_",i,".rda"))
  

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
  # 
  # assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstable_",i), 
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
}
runs
