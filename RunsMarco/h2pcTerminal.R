############################################################################################
# h2pc terminal. 
############################################################################################
rm(list = ls())
# install.packages(pkgs = "/home/catharina/Documents/Installations/modified2.tar.gz", repos = NULL,type ='source')
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

############################################################################################
# First of all: gamma runs.
############################################################################################
# run 1,5
for (i in c(1,5)){
  data <- datapermutations[[i]]
  for( j in 1:length(gammas)){
    
    assign(paste0("h2pc_",i,"_eBIC_g",gammas[j]),gamma_h2pc_split(data = data, iterations = 5000, gamma = gammas[j]))
    save(list = paste0("h2pc_",i,"_eBIC_g",gammas[j]),
         file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_",i,"_eBIC_g",gammas[j],".rda"))
    

  }
}

# run 2,4
for (i in c(2,4)){
    data <- datapermutations[[i]]
    for( j in 1:length(gammas)){
      
      assign(paste0("h2pc_",i,"_eBIC_g",gammas[j]),gamma_h2pc_split(data = data, iterations = 5000, gamma = gammas[j]))
      save(list = paste0("h2pc_",i,"_eBIC_g",gammas[j]),
           file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_",i,"_eBIC_g",gammas[j],".rda"))
      
      
    }
}
# run 3
    for (i in c(3)){
      data <- datapermutations[[i]]
      for( j in 1:length(gammas)){
        
        assign(paste0("h2pc_",i,"_eBIC_g",gammas[j]),gamma_h2pc_split(data = data, iterations = 5000, gamma = gammas[j]))
        save(list = paste0("h2pc_",i,"_eBIC_g",gammas[j]),
             file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_",i,"_eBIC_g",gammas[j],".rda"))
        
        
      }
}
############################################################################################
# Then gamma = 0 run with h2pc_iteration_split
############################################################################################
for (i in c(1,2,3,4,5)){
  data <- datapermutations[[i]]
  assign(paste0("h2pc_",i,"_eBIC_g0"),iteration_h2pc_split(data = data, iterations = 5000))
  save(list = paste0("h2pc_",i,"_eBIC_g0"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_",i,"_eBIC_g0.rda"))
}

############################################################################################
# Then gamma = 0 run with h2pc_iteration_split 2 parts to check
############################################################################################
data <- gaussian.test
for (i in c(1)){
  data <- datapermutations[[i]]
  hpcbegin <- iteration_h2pc_split_part1(data = data)
  assign(paste0("h2pc_",i,"_eBIC_g0_p1"), hpcbegin)
  save(list = paste0("h2pc_",i,"_eBIC_g0_p1"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_",i,"_eBIC_g0_p1.rda"))
}
# cuidado con hpcbegin (from which permutations comes?)
for (i in c(1)){
  data <- datapermutations[[i]]
  assign(paste0("h2pc_",i,"_eBIC_g0_p2"), iteration_h2pc_split_part2(data = data, iterations = 2000, start = hpcbegin))
  save(list = paste0("h2pc_",i,"_eBIC_g0_p2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_",i,"_eBIC_g0_p2.rda"))
}

