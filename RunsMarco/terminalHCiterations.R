#####################################################################################
# Hill Climbing iterations of all 25 permutations
#####################################################################################
rm(list = ls())
install.packages("bnlearn")
library(bnlearn)
library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations20.rda")


####################################################################################
# different terminals, permutations 1 - 5
####################################################################################
for(k in c(1,5)){
  data <- datapermutations[[k]]
  start <- NULL
  steps <- 100
  last <- 10000
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hc",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("hc",k,"_",i,"_",j,"i"), 
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",k,"/hc",k,"_",i,"_",j,"i.rda"))
    
    start <- berekening
    }
}


for(k in c(2)){
  data <- datapermutations[[k]]
  start <- NULL
  steps <- 100
  last <- 10000
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hc",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("hc",k,"_",i,"_",j,"i"), 
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",k,"/hc",k,"_",i,"_",j,"i.rda"))
    
    start <- berekening
  }
}


for(k in c(3)){
  data <- datapermutations[[k]]
  start <- NULL
  steps <- 100
  last <- 10000
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hc",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("hc",k,"_",i,"_",j,"i"), 
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",k,"/hc",k,"_",i,"_",j,"i.rda"))
    
    start <- berekening
  }
}

for(k in c(4)){
  data <- datapermutations[[k]]
  start <- NULL
  steps <- 100
  last <- 10000
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hc",k,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("hc",k,"_",i,"_",j,"i"), 
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",k,"/hc",k,"_",i,"_",j,"i.rda"))
    
    start <- berekening
  }
}

########################################################################################################
# Permutations 6 tot en met 25 Stopgezet na 10 (begin 11)
########################################################################################################

for(k in 1:20){
  permnumber <- k + 5
  data <- datapermutations20[[k]]
  start <- NULL
  steps <- 100
  last <- 9900
  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hc",permnumber,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("hc",permnumber,"_",i,"_",j,"i"), 
         file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",permnumber,"/hc",permnumber,"_",i,"_",j,"i.rda"))
    
    start <- berekening
  }
}

##############################################################################
# Only middles 11 tot en met 25 (k = 6:20)
##############################################################################
for(k in 6:20){
  permnumber <- k + 5
  data <- datapermutations20[[k]]
  begin <- 1100
  start <- NULL
  steps <- 100
  last <- 1900
  for (m in (begin/steps):(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hc",permnumber,"_",i,"_",j,"i"), berekening)
    
    save(list = paste0("hc",permnumber,"_",i,"_",j,"i"), 
         file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",permnumber,"/hc",permnumber,"_",i,"_",j,"i.rda"))
    
    start <- berekening
  }
}