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
gammas <- (40:5)
gammas2 <- (4:1)
runs <- 5

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

for (i in c(1,2,3,4,5)){
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

  assign(paste0("pcstable2_",i),iteration_pcstable_2(data = data, gammas = gammas))
  save(list = paste0("pcstable2_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_",i,".rda"))

}

for (i in c(1,2,3,4,5)){
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
  
  # assign(paste0("pcstablem_",i),iteration_pcstable(data = data, gammas = gammas))
  # save(list = paste0("pcstablem_",i),
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstablem_",i,".rda"))
  
  assign(paste0("gslow_",i),iteration_GS_2(data = data, gammas = gammas2))
  save(list = paste0("gslow_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gslow_",i,".rda"))
}

# i<- 2
# iteration_GS(data, c(5))
# iteration_GS_2(data, c(5))
# iteration_pcstable_2(data, c(100))
# 
# if (pdag2dag(getGraph(mat))$success == TRUE){1}
# test2 <- gs(data, test = "bic-gt", B = 20 )
# mat <- amat(test2)
# mat2 <- amat2dag(mat)
# all.equal(mat, mat2)
# cextend(test2)
# test <- gs(data, test = "bic-gt", B = 5 )
# cextend(test)
# test <- pc.stable(data, test = "bic-gt", B = 15)
# ext <- cextend(test)
# mat <- amat(test)
# amat2dag(mat)
# isValidGraph(mat, type = "pdag", verbose = TRUE)
# pdag2dag(getGraph(mat))
# addBgKnowledge(mat, verbose = TRUE)
# acyclic(ext)
# directed(ext)
# igraph <-igraph.from.graphNEL(as.graphNEL(ext))
# is.dag(igraph)
# 
# pdag2dag(mat)
