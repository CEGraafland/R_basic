##############################################################################
# Learning final networks with HC and TabuSearch for 5 datapermutations. 
# Use of old package bnlearn, because of Score function
# Use of imedeate learning in stead of intermedeate, because of reliable ntests.
##############################################################################
# install RIGHT VERSION of BNLEARN
# install.packages("bnlearn")
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")



# startlist <- list(hc_2,hc_3,hc_4,hc_5)
# grepp <- function(x){y <- x$networks[[length(x$networks)]];return(y)}
# startnetworks <- lapply(startlist, grepp)


for (i in 1:2){
  data <- datapermutations[[i]]
  assign(paste0("hc_",i,"_eBIC_g0"), gamma_hc(data = data, iterations = Inf, gamma = 0, start = NULL))
  save(list = paste0("hc_",i,"_eBIC_g0"),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,"_eBIC_g0.rda"))
}

for (i in 3:5){
  data <- datapermutations[[i]]
  assign(paste0("hc_",i,"_eBIC_g0"), gamma_hc(data = data, iterations = Inf, gamma = 0, start = NULL))
  save(list = paste0("hc_",i,"_eBIC_g0"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,"_eBIC_g0.rda"))
}


for (i in 1:2){
  data <- datapermutations[[i]]
  assign(paste0("tabu_",i,"_eBIC_g0"), gamma_tabu(data = data, iterations = Inf, gamma = 0, start = NULL))
  save(list = paste0("tabu_",i,"_eBIC_g0"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,"_eBIC_g0.rda"))
}

for (i in 3:5){
  data <- datapermutations[[i]]
  assign(paste0("tabu_",i,"_eBIC_g0"), gamma_tabu(data = data, iterations = Inf, gamma = 0, start = NULL))
  save(list = paste0("tabu_",i,"_eBIC_g0"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_",i,"_eBIC_g0.rda"))
}

# TESTEN
# tabu_1_eBIC_g0$networkdata
# 
# hc_1_eBIC_g0$networks
# hc_1_eBIC_g0$networkdata
# hc_edges_loglik_10d_100i
# x <- hc_2$networkdata
# tail(x,6)
# score(x = hc_2_eBIC_g0$networks$hc_10d_g0, data = datapermutations[[2]])
# ntests(hc_2_eBIC_g0$networks$hc_10d_g0)
# 
# hc_edges_loglik_10d_100i
# d <- gamma_hc(gaussian.test, iterations = Inf, gamma = 0, start = NULL)
# a <- hc(gaussian.test)
# ntests(a)
# nparams(a, gaussian.test)
# score(a, gaussian.test)
# b <- gamma_hc(gaussian.test, iterations = 2, gamma = 0, start = NULL)
# c <- gamma_hc(gaussian.test, iterations = Inf, gamma = 0, start = b$networks$hc_10d_g0)
# b$networkdata
# c$networkdata
# d$networkdata
# 
# for(i in c(1,2)){print(i)}
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")
# ntests(hc_edges_loglik_10d_200i$networks)
# e  <-gamma_hc(datapermutations[[1]], iterations = 200, gamma = 0, start = NULL)
# e$networkdata
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_1.rda")
# tabu_1_eBIC_g0$networkdata
# tabu_1$networkdata[73,]
# 1306 - 648
# nparams(, x = tabu_1_eBIC_g0$networks$tabu_10d_g0, effective = FALSE)
# nparams(tabu_1$networks$tabu_10d_73i, datapermutations[[1]])
# 1369 - 721
# tabu_1$networks$tabu_10d_253i