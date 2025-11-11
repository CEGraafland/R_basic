#####################################################################################
# H2PC
#####################################################################################
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
runs <- 5
h2pc_1_eBIC_g0 <- h2pc(initialdata, whitelist = NULL, blacklist = NULL, restrict.args = list(),
     maximize.args = list(), debug = FALSE)
h2pc_1
h2pc_1_eBIC_g0 <- rsmax2(initialdata, whitelist = NULL, blacklist = NULL, restrict = "hpc",
       maximize = "hc", restrict.args = list(test = "bic-gt"), maximize.args = list(score = "bic-g"), debug = FALSE)

h2pc_1_eBIC_g0_1 <- rsmax2(initialdata, whitelist = NULL, blacklist = NULL, restrict = "hpc",
                        maximize = "hc", restrict.args = list(test = "bic-gt", B = 0), maximize.args = list(score = "bic-g"), debug = FALSE)

gamma <- 10
hpc0.01 <- hpc(initialdata,test = "bic-gt", B = 0.01,debug = TRUE)
hpc(initialdata,test = "bic-gt", B = 0,debug = TRUE)
hpc(initialdata,test = "bic-gt", B = 0.01,debug = TRUE)
h2pc25 <- gamma_h2pc_split(initialdata, 5000, gamma = 25)
h2pc25
data <- gaussian.test
iterations <- 2000
start <- NULL
gaussian0 <- iteration_h2pc_split(data = gaussian.test, iterations = 2000, start = NULL)
gaussian0mmhc <- iteration_mmhc_split(gaussian.test, 2000, start = NULL)
gaussian0mmhc
mmhc_5_eBIC_g25
mmhc_1_eBIC_g0.1
mmhc_1_eBIC_g1
mmhc_1_eBIC_g0

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")


i <- 1
data <- datapermutations[[i]]


for (i in c(5)){
  data <- datapermutations[[i]]
  
  assign(paste0("tabu_",i),iteration_tabu(data, 8000))
  save(list = paste0("tabu_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,".rda"))
  
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

for (i in c(1,5)){
  data <- datapermutations[[i]]
  
  # assign(paste0("tabu_",i),iteration_tabu(data, 8000))
  # save(list = paste0("tabu_",i),
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,".rda"))
  
  assign(paste0("hc_",i),iteration_hc(data, 8000))
  save(list = paste0("hc_",i),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,".rda"))
  
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
