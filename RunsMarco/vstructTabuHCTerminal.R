library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(igraph)
library(transformeR)
library(gridExtra)
library(visualizeR)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")

rm(list = ls())

for(j in c(1,2,3,4,5)){
  pattern <- paste0("tabu_",j,"_eBIC_g")
  filestabu1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
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
}

##########################################################################
# Load HC eBIC
##########################################################################

for(j in c(1,2,3,4,5)){
  pattern <- paste0("hc_",j,"_eBIC_g")
  filestabu1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
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
  assign(paste0("hc_",j,"_eBIC"),variablelisttabu1)
}


vstructsprop <- function(dag){
  # Select graph
  g <- dag
  
  # Give names (To avoid problems)
  names <- c()
  for(i in 1:648) names <- append(names,paste0("V",i))
  # obtain vstructs immoral (with BNLEARN)
  voefimmor <- vstructs(g, moral = FALSE)
  nvoefimmor <- nrow(voefimmor)
  # Which of those vstructs leads to dseperation?
  if (!nvoefimmor == 0){
    dsepvimmor <- logical()
    for(i in 1:nrow(voefimmor)){
      answer <- bnlearn::dsep(g,voefimmor[i,1],voefimmor[i,3])
      dsepvimmor[i] <- answer
    }
    blockedim <- which(dsepvimmor)
    # Matrix with Dseperating Vstructs
    permv <- voefimmor[blockedim,]
  } else {permv <- blockedim <- dsepvimmor <- 0}
  ###########################################################################################
  # Now calculate all paths of length 2 in the graph 
  degarray <- array(dim = c(length(nodes(g))))
  twopathsarray <- array(dim = c(length(nodes(g))))
  
  # per node 
  deg <- c()
  twopaths <- c()
  for(i in 1:length(nodes(g))){
    node <- nodes(g)[i]
    deg[i] <- bnlearn::degree(g,node)
    if(deg[i]>1){
      twopaths[i] <- ncol(combn(deg[i],2))
    } else{twopaths[i] <- 0}
  }
  degarray <- deg
  twopathsarray <- twopaths
  
  # All possible two paths in the graph are the sum of possible twopaths per node. 
  possible2pathsmatrix <-  sum(twopathsarray)
  
  # relative d seperating vstructs
  if(length(blockedim) ==0){
    nvstructs <- 0
  } else if(length(blockedim) ==1){
    nvstructs <- 1
  } else {nvstructs <- nrow(permv)}
  
  relnumbervstructs <- nvstructs/possible2pathsmatrix
  
  return(c(nvoefimmor, nvstructs,relnumbervstructs))
}


ejes_gamma_vstructs <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- (nets$networkdata$gamma)*2
  miss <- !is.na(xeje)
  vs <- lapply(network,vstructsprop)
  yeje <- sapply(vs, function(x) x[1])
  return(list(xeje,yeje))
}

tabueBICs <- list(tabu_1_eBIC,tabu_2_eBIC,tabu_3_eBIC,tabu_4_eBIC,tabu_5_eBIC)
ejestabuBIC <- list()
for (i in 1:5){
  ej <- ejes_gamma_vstructs(tabueBICs[[i]])
  assign(paste0("ejestabuBIC_vstruct_",i), ejes_gamma_vstructs(tabueBICs[[i]]))
  save(list = paste0("ejestabuBIC_vstruct_",i), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dsepvstructseBIC/ejestabuBIC_vstruct_",i,".rda"))
  ejestabuBIC[[i]] <- ej
}

save(ejestabuBIC, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dsepvstructseBIC/ejestabuBIC_vstruct.rda")


HCeBICs <- list(hc_1_eBIC, hc_2_eBIC, hc_3_eBIC,hc_4_eBIC,hc_5_eBIC)
ejeshcBIC <- list()
for (i in 1:5){
  ej <- ejes_gamma_vstructs(HCeBICs[[i]])
  assign(paste0("ejeshcBIC_vstruct_",i), ejes_gamma_vstructs(HCeBICs[[i]]))
  save(list = paste0("ejeshcBIC_vstruct_",i), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dsepvstructseBIC/ejeshcBIC_vstruct_",i,".rda"))
  ejeshcBIC[[i]] <- ej
}

save(ejeshcBIC, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dsepvstructseBIC/ejeshcBIC_vstruct.rda")
