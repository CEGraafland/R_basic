#####################################################################################
# Sensitivity Analysis Bayesian Networks
#####################################################################################
rm(list = ls())
library(transformeR)
library(magrittr)
library(bnlearn)
# library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
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
# h2pc eBIC 0.1
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_1_eBIC_g0.1.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_2_eBIC_g0.1.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_3_eBIC_g0.1.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_4_eBIC_g0.1.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/h2pc/h2pc_5_eBIC_g0.1.rda")


#permutations
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

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

h2pcs0.1 <- list(h2pc_1_eBIC_g0.1,
              h2pc_2_eBIC_g0.1,
              h2pc_3_eBIC_g0.1,
              h2pc_4_eBIC_g0.1,
              h2pc_5_eBIC_g0.1)
h2pcs0.1 <- lapply(h2pcs0.1, function(x) x$networks)
h2pcs0.1 <- lapply(h2pcs0.1, function(x) x[[1]])


#########################################################################################
# Tabu General 
#########################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

tabufits <- list()
mmhcfits <- list()
tabulargefits<- list()
h2pcfits <- list()

for (i in 1:5){
tabufits[[i]] <- bn.fit(x = tabus[[i]], data = dataframes[[i]])
mmhcfits[[i]] <- bn.fit(x = mmhcs[[i]], data = dataframes[[i]])
tabulargefits[[i]] <- bn.fit(x = tabuslarge[[i]], data = dataframes[[i]])
}
for (i in 1:5){
  h2pcfits[[i]] <- bn.fit(x = h2pcs0.1[[i]], data = dataframes[[i]])
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
  
#,2,3,4,5
for (i in c(2,3,4,5)){

  assign(paste0("prop_h2pc_g0.1_",i,"_",narcs(h2pcs0.1[[i]]),"_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = h2pcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[i]]))
  save(list = paste0("prop_h2pc_g0.1_",i,"_",narcs(h2pcs0.1[[i]]),"_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC/prop_h2pc_g0.1_",i,"_",narcs(h2pcs0.1[[i]]),"_V81_equal2.rda"))
  
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
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")
x <- seq(200,8600,200)
y <- seq(400,8800,200)
i <- 200
iterationnetworks <- list()


for(i in 1:length(x)){
  iterationnetworks[[i]] <- eval(parse(text = paste0("hc_edges_loglik_10d_",x[i],"_",y[i],"i$networks")))
}

#permutations
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

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

# last network (DONE)
x[length(x)]
y[length(y)]
length(hcfits)
for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,280),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_",x[i],"_",y[i],"i_V81V280_equal22.rda"))
  
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

# last network negative (NOT YET)
for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,280),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22.rda"))
  
}
#################################################################################
# Single evidence. V81 postive (+ +)
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
# last network positive
for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V81 (+ -)
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

# last network negative
for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
}

#################################################################################
# Negative Niño. V81 = -2 Positive (- +)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81pos/prop_hc_",x[i],"_",y[i],"i_V81_equalmin2.rda"))
  
}



for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81pos/prop_hc_",x[i],"_",y[i],"i_V81_equalmin2.rda"))
  
}
#################################################################################
# Negative Niño. V81 = -2 Negative (- -)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81neg/propneg_hc_",x[i],"_",y[i],"i_V81_equalmin2.rda"))
  
}



for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81neg/propneg_hc_",x[i],"_",y[i],"i_V81_equalmin2.rda"))
  
}

####################################################################################
# Specific HC iteration propagation V280 V424 extra evidence
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

for (i in 6:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(280,424),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280V424pos/prop_hc_",x[i],"_",y[i],"i_V280V424_equal22.rda"))
  
}
############################################################################
# double V280 V424 negative
#############################################################################
for (i in 6:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(280,424),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280V424neg/propneg_hc_",x[i],"_",y[i],"i_V280V424_equal22.rda"))
  
}
####################################################################################
# Specific HC iteration propagation V280 V424 extra evidence
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

for (i in 6:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(280,424),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280V424pos/prop_hc_",x[i],"_",y[i],"i_V280V424_equal22.rda"))
  
}
############################################################################
# double V280 V424 negative
#############################################################################
for (i in 6:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(280,424),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V280V424_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280V424neg/propneg_hc_",x[i],"_",y[i],"i_V280V424_equal22.rda"))
  
}


####################################################################################
# Specific HC iteration propagation V81 V171 extra evidence 
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

for (i in 6:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171pos/prop_hc_",x[i],"_",y[i],"i_V81V171_equal22.rda"))
  
}
####################################################################################
# Specific HC iteration propagation V81 V171 extra evidence WARM POOL 
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])
y[c(8,13,43)]
for (i in c(8,13,43)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal20"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal20"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171WP/prop_hc_",x[i],"_",y[i],"i_V81V171_equal20.rda"))
  
}

####################################################################################
# V81 V171 extra evidence. Negative.  WARM POOL 
####################################################################################
y[c(8,13,43)]
for (i in c(8,13,43)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V171_equal20"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V171_equal20"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171WPneg/propneg_hc_",x[i],"_",y[i],"i_V81V171_equal20.rda"))
  
}

####################################################################################
# Specific HC iteration propagation V28 V171 extra evidence WARM POOL 2
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])
y[c(8,13,43)]
for (i in c(8,13,43)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V28V171_equal20"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(28,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V28V171_equal20"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28V171WP2pos/prop_hc_",x[i],"_",y[i],"i_V28V171_equal20.rda"))
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V28V171_equal20"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(28,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V28V171_equal20"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28V171WP2neg/propneg_hc_",x[i],"_",y[i],"i_V28V171_equal20.rda"))
}

####################################################################################
# V28 V171 extra evidence. Negative.  WARM POOL 2
####################################################################################
y[c(8,13,43)]
for (i in c(8,13,43)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V28V171_equal20"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(28,171),
                                     valueEvidence = c(2,0),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V28V171_equal20"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28V171WP2neg/propneg_hc_",x[i],"_",y[i],"i_V28V171_equal20.rda"))
  
}

####################################################################################
# Specific HC iteration propagation V81 V171 extra evidence COLD TONG forzada: (1,2)
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

y[c(8,13,43)]
for (i in c(8,13,43)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal12"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(1,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal12"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171CT/prop_hc_",x[i],"_",y[i],"i_V81V171_equal12.rda"))
  
}

####################################################################################
# Single evidence. V171 pos (Simulating COLD TONG )
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

y[c(8,13,43)]
for (i in c(8,13,43)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V171_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(171),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V171_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V171pos/prop_hc_",x[i],"_",y[i],"i_V171_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V171 (Simulating COLD TONG)
#################################################################################

for (i in c(8,13,43)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V171_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(171),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V171_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V171neg/propneg_hc_",x[i],"_",y[i],"i_V171_equal2.rda"))
  
}

####################################################################################
# Specific HC iteration propagation V81 V171 extra evidence POS NEG
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

for (i in 6:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,171),
                                     valueEvidence = c(2,-2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171posneg/prop_hc_",x[i],"_",y[i],"i_V81V171_equal22.rda"))
  
}

#################################################################################
# EVIDENCIA V205 / V227 / V568 (V81)
#################################################################################
####################################################################################
# Specific HC iteration propagation V205 V227 extra evidence (positive)
####################################################################################
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
initialdatapermutations <- lapply(permutations, function(x){initialdata[x]})
dataframes <- lapply(initialdatapermutations, as.data.frame)

nedgesnetworks <- as.character(sapply(iterationnetworks, narcs))
nedgesnetworks
hcfits <- lapply(iterationnetworks, bn.fit, data = dataframes[[1]])

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V205V227_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(205,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V205V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V205V227pos/prop_hc_",x[i],"_",y[i],"i_V205V227_equal22.rda"))
  
}
###################################################################################
# Specific HC iteration propagation V227 V568 extra evidence (positive)
###################################################################################
for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V227V568_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227,568),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V227V568_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227V568pos/prop_hc_",x[i],"_",y[i],"i_V227V568_equal22.rda"))
  
}

###################################################################################
# Specific HC iteration propagation V81 V227 extra evidence (positive)
###################################################################################
for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V227pos/prop_hc_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
  
}

for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
  
}
###################################################################################
# Specific HC iteration propagation V81 V227 extra evidence (negative)
###################################################################################
for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V227neg/propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
  
}

for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
  
}

#######################################################################################
# Specific HC iteration propagation min V81 plus V227 extra evidence (positive) (- + +)
#######################################################################################
for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equalmin2plus2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(-2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equalmin2plus2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81plusV227pos/prop_hc_",x[i],"_",y[i],"i_V81V227_equalmin2plus2.rda"))
  
}

for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equalmin15plus15"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(-1.5,1.5),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equalmin15plus15"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81plusV227pos/prop_hc_",x[i],"_",y[i],"i_V81V227_equalmin15plus15.rda"))
  
}
#######################################################################################
# Specific HC iteration propagation min V81 plus V227 extra evidence (negative) (- + -)
#######################################################################################
for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equalmin2plus2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(-2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equalmin2plus2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81plusV227neg/propneg_hc_",x[i],"_",y[i],"i_V81V227_equalmin2plus2.rda"))
  
}

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equalmin15plus15"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(-1.5,1.5),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equalmin15plus15"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81plusV227neg/propneg_hc_",x[i],"_",y[i],"i_V81V227_equalmin15plus15.rda"))
  
}



#################################################################################
# Single evidence. V205 postive
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V205_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(205),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V205_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V205pos/prop_hc_",x[i],"_",y[i],"i_V205_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V205
#################################################################################

for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V205_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(205),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V205_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V205neg/propneg_hc_",x[i],"_",y[i],"i_V205_equal2.rda"))
  
}
#################################################################################
# Single evidence. V227 postive HERE
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227pos/prop_hc_",x[i],"_",y[i],"i_V227_equal2.rda"))
  
}

x[length(x)]
y[length(y)]
for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_",x[i],"_",y[i],"i_V227_equal2.rda"))
  
}
#################################################################################
# Negative evidence. V227 HERE
#################################################################################

for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227neg/propneg_hc_",x[i],"_",y[i],"i_V227_equal2.rda"))
  
}

for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_",x[i],"_",y[i],"i_V227_equal2.rda"))
  
}
#################################################################################
# Single evidence. V568 postive (+ +)
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V568_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(568),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V568_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V568pos/prop_hc_",x[i],"_",y[i],"i_V568_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V568 (+ -)
#################################################################################

for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V568_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(568),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V568_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V568neg/propneg_hc_",x[i],"_",y[i],"i_V568_equal2.rda"))
  
}

#################################################################################
# Negative V568 = -2 Positive (- +)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V568_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(568),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V568_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV568pos/prop_hc_",x[i],"_",y[i],"i_V568_equalmin2.rda"))
  
}


#################################################################################
# Negative V568 = -2 Negative (- -)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V568_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(568),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V568_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV568neg/propneg_hc_",x[i],"_",y[i],"i_V568_equalmin2.rda"))
  
}

#################################################################################
# Single evidence. V550 postive (+ +)
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V550_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(550),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V550_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V550pos/prop_hc_",x[i],"_",y[i],"i_V550_equal2.rda"))
  
}


#################################################################################
# Negative evidence. V550 (+ -)
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V550_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(550),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V550_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V550neg/propneg_hc_",x[i],"_",y[i],"i_V550_equal2.rda"))
  
}
#################################################################################
# Negative V550 = -2 Positive (- +)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V550_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(550),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V550_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV550pos/prop_hc_",x[i],"_",y[i],"i_V550_equalmin2.rda"))
  
}


#################################################################################
# Negative V550 = -2 Negative (- -)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V550_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(550),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V550_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV550neg/propneg_hc_",x[i],"_",y[i],"i_V550_equalmin2.rda"))
  
}

#################################################################################
# Single evidence. V280 postive
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V280_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(280),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V280_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280pos/prop_hc_",x[i],"_",y[i],"i_V280_equal2.rda"))
  
}

# last network positive
for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V280_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(280),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V280_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_",x[i],"_",y[i],"i_V280_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V280
#################################################################################

for (i in 1:15){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V280_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(280),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V280_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280neg/propneg_hc_",x[i],"_",y[i],"i_V280_equal2.rda"))
  
}
# last network negative 
for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V280_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(280),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V280_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_",x[i],"_",y[i],"i_V280_equal2.rda"))
  
}
###################################################################################
# Extra for comunities 1800 only evidence propagation in 1800. i = 8
################################################################################### 
i <-8
x[8]
y[8]
#################################################################################
# Single evidence. V424 postive
#################################################################################

for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V424_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(424),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V424_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V424pos/prop_hc_",x[i],"_",y[i],"i_V424_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V424
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V424_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(424),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V424_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V424neg/propneg_hc_",x[i],"_",y[i],"i_V424_equal2.rda"))
  
}
#################################################################################
# Single evidence. V459 postive
#################################################################################

for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V459_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(459),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V459_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V459pos/prop_hc_",x[i],"_",y[i],"i_V459_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V459
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V459_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(459),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V459_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V459neg/propneg_hc_",x[i],"_",y[i],"i_V459_equal2.rda"))
  
}

# #################################################################################
# # Single evidence. V568postive (AL GEDAAN)
# #################################################################################
# 
# for (i in c(8)){
#   
#   assign(paste0("prop_hc_",x[i],"_",y[i],"i_V568_equal2"),
#          PropagationExactGeneralPerm(baysnet = hcfits[[i]],
#                                      nodesEvents = 1:648,
#                                      valueEvent = ">= 1",
#                                      nodesEvidence = c(568),
#                                      valueEvidence = c(2),
#                                      perm = permutations[[1]]))
#   save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V568_equal2"),
#        file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V568pos/prop_hc_",x[i],"_",y[i],"i_V568_equal2.rda"))
#   
# }

#################################################################################
# Negative evidence. V568
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V568_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(568),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V568_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V568neg/propneg_hc_",x[i],"_",y[i],"i_V568_equal2.rda"))
  
}

#################################################################################
# Single evidence. V532 postive
#################################################################################

for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V532_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(532),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V532_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V532pos/prop_hc_",x[i],"_",y[i],"i_V532_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V532
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V532_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(532),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V532_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V532neg/propneg_hc_",x[i],"_",y[i],"i_V532_equal2.rda"))
  
}
#################################################################################
# Negative V532 = -2 Positive (- +)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V532_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(532),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V532_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV532pos/prop_hc_",x[i],"_",y[i],"i_V532_equalmin2.rda"))
  
}


#################################################################################
# Negative V532 = -2 Negative (- -)
#################################################################################
x[8]
y[8]
for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V532_equalmin2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(532),
                                     valueEvidence = c(-2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V532_equalmin2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV532neg/propneg_hc_",x[i],"_",y[i],"i_V532_equalmin2.rda"))
  
}
#################################################################################
# Single evidence. V85 postive
#################################################################################

for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V85_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(85),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V85_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V85pos/prop_hc_",x[i],"_",y[i],"i_V85_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V85
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V85_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(85),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V85_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V85neg/propneg_hc_",x[i],"_",y[i],"i_V85_equal2.rda"))
  
}
#################################################################################
# Single evidence. V28 postive
#################################################################################

for (i in c(8)){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V28_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(28),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V28_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28pos/prop_hc_",x[i],"_",y[i],"i_V28_equal2.rda"))
  
}

#################################################################################
# Negative evidence. V28
#################################################################################

for (i in c(8)){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V28_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(28),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V28_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28neg/propneg_hc_",x[i],"_",y[i],"i_V28_equal2.rda"))
  
}
#################################################################################
# Single evidence. V81 postive
#################################################################################

for (i in 1:15){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V640_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos/prop_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
}
# last network positive
for (i in c(length(x))){
  
  assign(paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
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

# last network negative
for (i in c(length(x))){
  
  assign(paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
}

