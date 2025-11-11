#####################################################################################
# Sensitivity Analysis Bayesian Networks
#####################################################################################
rm(list = ls())
library(transformeR)
library(magrittr)
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
#############################################################################################
# Load data
#######################################################################################
# mmhc
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_1.rda")
# Tabu simples
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple5.rda")
# Tabuslarge
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_5.rda")
#permutations
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

mmhcs <- list(mmhc_1$networks[[3]], 
              mmhc_2$networks[[3]], 
              mmhc_3$networks[[3]], 
              mmhc_4$networks[[3]], 
              mmhc_5$networks[[3]])

tabus <- list(tabu_mmhc_simple1,
              tabu_mmhc_simple2, 
              tabu_mmhc_simple3,
              tabu_mmhc_simple4,
              tabu_mmhc_simple5)

tabuslarge <- list(tabu_simple_1600_1,
                   tabu_simple_1600_2,
                   tabu_simple_1600_3,
                   tabu_simple_1600_4,
                   tabu_simple_1600_5)




#########################################################################################
# General 
#########################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

tabufits <- list()
mmhcfits <- list()
tabulargefits<- list()

for (i in 1:5){
tabufits[[i]] <- bn.fit(x = tabus[[i]], data = dataframes[[i]])
mmhcfits[[i]] <- bn.fit(x = mmhcs[[i]], data = dataframes[[i]])
tabulargefits[[i]] <- bn.fit(x = tabuslarge[[i]], data = dataframes[[i]])
}

for (i in c(1,2,3,4,5)){
  
  # assign(paste0("prop_tabu",i,"_",narcs(tabus[[i]]),"_V205_equal2"),
  #        PropagationExactGeneralPerm(baysnet = tabufits[[i]],
  #                                    nodesEvents = 1:648,
  #                                    valueEvent = ">= 1",
  #                                    nodesEvidence = c(205),
  #                                    valueEvidence = c(2),
  #                                    perm = permutations[[i]]))
  # save(list = paste0("prop_tabu",i,"_",narcs(tabus[[i]]),"_V205_equal2"),
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop/prop_tabu",i,"_",narcs(tabus[[i]]),"_V205_equal2.rda"))
  # 
  # assign(paste0("prop_mmhc",i,"_",narcs(mmhcs[[i]]),"_V205_equal2"), 
  #        PropagationExactGeneralPerm(baysnet = mmhcfits[[i]], 
  #                                    nodesEvents = 1:648,
  #                                    valueEvent = ">= 1",
  #                                    nodesEvidence = c(205),
  #                                    valueEvidence = c(2),
  #                                    perm = permutations[[i]]))
  # save(list = paste0("prop_mmhc",i,"_",narcs(mmhcs[[i]]),"_V205_equal2"),
  #      file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop/prop_mmhc",i,"_",narcs(mmhcs[[i]]),"_V205_equal2.rda"))

  assign(paste0("prop_tabu_1600_",i,"_",narcs(tabuslarge[[i]]),"_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = tabulargefits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[i]]))
  save(list = paste0("prop_tabu_1600_",i,"_",narcs(tabuslarge[[i]]),"_V81_equal2"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop/prop_tabu_1600_",i,"_",narcs(tabuslarge[[i]]),"_V81_equal2.rda"))
  
  
  }
  

##############################################################################################
# specific tabu g propagation
##############################################################################################
##########################################################################
# Load tabu eBIC
##########################################################################
tabus <- list()
for(j in c(1,2,5)){
  pattern <- paste0("tabu_",j,"_eBIC_g")
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(), networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1$networks[[i]] <- variablepos$networks[[1]]
    names[i] <- names(variablepos$networks)
    data.frames[i,] <- as.matrix(variablepos$networkdata)
  }
  colnames(data.frames) <- names(variablepos$networkdata)
  data.frames <- as.data.frame(data.frames)
  variablelisttabu1$networkdata<- data.frames
  variablelisttabu1$networkdata
  
  names(variablelisttabu1$networks) <- names
  assign(paste0("tabu_",j,"_eBIC"),variablelisttabu1)
  tabus[[j]] <- variablelisttabu1
}


#permutations
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")


####################################################################################
# Specific tabu propagation
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)


tabufits <- list()
for (i in c(1,2)){
  tabufits[[i]] <- lapply(tabus[[i]]$networks, bn.fit, data = dataframes[[i]])
}

for (i in c(1,2)){
  for(j in 1:length(tabufits[[1]])) { 
    assign(paste0("prop_",i,"_",names[j],"_V81_equal2"),
           PropagationExactGeneralPerm(baysnet = tabufits[[i]][[j]],
                                       nodesEvents = 1:648,
                                       valueEvent = ">= 1",
                                       nodesEvidence = c(81),
                                       valueEvidence = c(2),
                                       perm = permutations[[i]]))
    save(list = paste0("prop_",i,"_",names[j],"_V81V_equal2"),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_",i,"_",names[j],"_V81_equal2.rda"))
    
  }
}

##########################################################################################
# add propagation 5
##########################################################################################
tabufits <- list()
for (i in c(1,2,5)){
  tabufits[[i]] <- lapply(tabus[[i]]$networks, bn.fit, data = dataframes[[i]])
}

for (i in c(5)){
  for(j in 1:length(tabufits[[1]])) { 
    assign(paste0("prop_",i,"_",names[j],"_V81_equal2"),
           PropagationExactGeneralPerm(baysnet = tabufits[[i]][[j]],
                                       nodesEvents = 1:648,
                                       valueEvent = ">= 1",
                                       nodesEvidence = c(81),
                                       valueEvidence = c(2),
                                       perm = permutations[[i]]))
    save(list = paste0("prop_",i,"_",names[j],"_V81_equal2"),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_",i,"_",names[j],"_V81_equal2.rda"))
    
  }
}

####################################################################################
# Specific tabu propagation V280 V81 Extra evidence
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)


tabufits <- list()
for (i in c(1,2)){
  tabufits[[i]] <- lapply(tabus[[i]]$networks, bn.fit, data = dataframes[[i]])
}


for (i in c(1,2)){
  for(j in c(3,4,5)) { 
    assign(paste0("prop_",i,"_",names[j],"_V81V280_equal22"),
           PropagationExactGeneralPerm(baysnet = tabufits[[i]][[j]],
                                       nodesEvents = 1:648,
                                       valueEvent = ">= 1",
                                       nodesEvidence = c(81,280),
                                       valueEvidence = c(2,2),
                                       perm = permutations[[i]]))
    save(list = paste0("prop_",i,"_",names[j],"_V81V280_equal22"),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_",i,"_",names[j],"_V81V280_equal22.rda"))
    
  }
}
####################################################################################
# load HC iteration data and make list
####################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")
x <- seq(200,8600,200)
y <- seq(400,8800,200)
i <- 200
iterationnetworks <- list()


for(i in 1:length(x)){
  iterationnetworks[[i]] <- eval(parse(text = paste0("hc_edges_loglik_10d_",x[i],"_",y[i],"i$networks")))
}

#permutations
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

####################################################################################
# Specific HC iteration propagation V280 V81 extra evidence
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

for (i in 1:15){

    assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
           PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                       nodesEvents = 1:648,
                                       valueEvent = ">= 1",
                                       nodesEvidence = c(81,280),
                                       valueEvidence = c(2,2),
                                       perm = permutations[[1]]))
    save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
         file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/prop_hc_",x[i],"_",y[i],"i_V81V280_equal22.rda"))
    
}
############################################################################
# double V81V280 negative
#############################################################################
for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,280),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V280neg/propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22.rda"))
  
}

#################################################################################
# Single evidence. V81 postive
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos/prop_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V81
#################################################################################

for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81neg/propneg_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
}