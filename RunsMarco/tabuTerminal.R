#####################################################################################
# Hill Climbing and Tabu Search. 
#####################################################################################
rm(list = ls())
install.packages(pkgs = "/home/catharina/Documents/Installations/modified_bnlearn.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions.R")


initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
alphas <- c(2:10 %o% 10^(-112:-2))
runs <- 5

datapermutations <- list()
datapermutations[[1]] <- initialdata
for (i in 2:runs){
  datapermutations[[i]] <- initialdata[,sample(ncol(initialdata))]
}
save(datapermutations, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

permutations <- list()
for (j in 1:length(datapermutations)){
  indsame <- c()
  for (i in 1:ncol(datapermutations[[1]])){
    int <- which(colnames(datapermutations[[1]]) == colnames(datapermutations[[j]][i]))
    indsame[i] <- int
  }
  permutations[[j]] <- indsame
}

save(permutations, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

for (i in 2:runs){
  data <- datapermutations[[i]]
  
  assign(paste0("tabu_",i),iteration_tabu(data, 8000))
  save(list = paste0("tabu_",i),
        file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,".rda"))

  assign(paste0("hc_",i),iteration_hc(data, 8000))
  save(list = paste0("hc_",i),
        file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,".rda"))
  
  assign(paste0("pc_",i),iteration_pc(data = data, correlation = NULL, alpha = alphas))
  save(list = paste0("pc_",i), 
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_",i,".rda"))
  
  assign(paste0("mmhc_",i),iteration_mmhc_split(data = data, iterations = 2000))
  save(list = paste0("mmhc_",i), 
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,".rda"))

}

runs

