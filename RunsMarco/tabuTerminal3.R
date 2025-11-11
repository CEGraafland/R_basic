#####################################################################################
# Hill Climbing and Tabu Search. 
#####################################################################################
rm(list = ls())
# install.packages(pkgs = "/home/catharina/Documents/Installations/modified_4.3.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")


initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
gammas <- c(25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0.75,0.5,0.25,0.1)
runs <- 5

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    assign(paste0("tabu_",i,"_eBIC_g",gammas[j]),gamma_tabu(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("tabu_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,"_eBIC_g",gammas[j],".rda"))
    
    assign(paste0("mmhc_",i,"_eBIC_g",gammas[j]),gamma_mmhc_split(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("mmhc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g",gammas[j],".rda"))
  
    assign(paste0("hc_",i,"_eBIC_g",gammas[j]),gamma_hc(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("hc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,"_eBIC_g",gammas[j],".rda"))
  }
  
  # assign(paste0("hc_",i),iteration_hc(data, 8000))
  # save(list = paste0("hc_",i),
  #       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,".rda"))
  # 
  # # assign(paste0("pc_",i),iteration_pc(data = data, correlation = NULL, alpha = alphas))
  # # save(list = paste0("pc_",i), 
  # #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_",i,".rda"))
  # 
  # assign(paste0("mmhc_",i,"_eBIC_g",gammas[j]),iteration_mmhc_split(data = data, iterations = 2000))
  # save(list = paste0("mmhc_",i,"_eBIC_g",gammas[j]),
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g",gammas[j],".rda"))
  
  # assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstable_",i), 
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
}

for (i in c(2)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    # assign(paste0("tabu_sel",i,"_eBIC_g",gammas[j]),iteration_tabu_sel(data, iterations = 2000, gamma = gammas[j]))
    # save(list = paste0("tabu_sel",i,"_eBIC_g",gammas[j]),
    #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_sel",i,"_eBIC_g",gammas[j],".rda"))
    # 

    assign(paste0("mmhc_",i,"_eBIC_g",gammas[j]),gamma_mmhc_split(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("mmhc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g",gammas[j],".rda"))
    
    }
  
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
  
  # assign(paste0("pcstable_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstable_",i), 
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_",i,".rda"))
  
}
for (i in c(3)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    assign(paste0("tabu_",i,"_eBIC_g",gammas[j]),gamma_tabu(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("tabu_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,"_eBIC_g",gammas[j],".rda"))
    
    assign(paste0("mmhc_",i,"_eBIC_g",gammas[j]),gamma_mmhc_split(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("mmhc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g",gammas[j],".rda"))
  }
}

for (i in c(4)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    assign(paste0("tabu_",i,"_eBIC_g",gammas[j]),gamma_tabu(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("tabu_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,"_eBIC_g",gammas[j],".rda"))
    
    assign(paste0("mmhc_",i,"_eBIC_g",gammas[j]),gamma_mmhc_split(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("mmhc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g",gammas[j],".rda"))
  }
}

for (i in c(2,5)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    assign(paste0("tabu_",i,"_eBIC_g",gammas[j]),gamma_tabu(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("tabu_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,"_eBIC_g",gammas[j],".rda"))
    
    assign(paste0("mmhc_",i,"_eBIC_g",gammas[j]),gamma_mmhc_split(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("mmhc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g",gammas[j],".rda"))
  }
}

for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    assign(paste0("hc_",i,"_eBIC_g",gammas[j]),gamma_hc(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("hc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,"_eBIC_g",gammas[j],".rda"))
    
  }
}
runs
