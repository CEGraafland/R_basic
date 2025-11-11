##################################################################################
# Load Bayesian Network perm 1
# Create gridGraphsBN
##################################################################################
library(bnlearn)
library(igraph)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
# IFCA: 
# filesinmap <- list.files("/media/catharina/ubuntu/outofoffice/hc_perm1", full.names = T)
# filesnames <- list.files("/media/catharina/ubuntu/outofoffice/hc_perm1")
# MACBOOK:
# filesinmap <- list.files("/Volumes/ubuntu/outofoffice/hc_perm1", full.names = T)
# filesnames <- list.files("/Volumes/ubuntu/outofoffice/hc_perm1")
# OCEANO
# filesinmap <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1", full.names = T)
# filesnames <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1")
# OCEANO LOCAL
filesinmap <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1", full.names = T)
filesnames <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1")
filesnames <- gsub(".rda", "", filesnames)

networklist <- list()
names <- c()

for (i in 1:length(filesinmap)){
  variablepos <- get(load(filesinmap[i]))
  networklist[[i]] <- variablepos
}
names(networklist) <- filesnames
rm(list = filesnames)



# Convert the graphs to undirected igraphs
graphsNEL <- lapply(networklist, as.graphNEL) 
igraphsdir <- lapply(graphsNEL, igraph.from.graphNEL)
igraphsske<- lapply(igraphsdir, as.undirected)
edgelists <- lapply(igraphsske, E)
nedgeslists <- sapply(edgelists, length)
nedgeslists
indsame <- c()
# Sort the graphs from small to big
for (i in 1:length(nedgeslists)){
  int <- which(nedgeslists == sort(nedgeslists)[i])
  indsame[i] <- as.vector(int[1])
}
igraphsskesort <- igraphsske[indsame]
perm1sort <- networklist[indsame]

# Create the graphObject as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(igraphsskesort))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- igraphsskesort[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsskesort[[i]])
}

gridGraphsBN <- graphObjects
names(gridGraphsBN) <- names(networklist)[indsame]
GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})
gridGraphsBN$

 save(perm1sort, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
 save(gridGraphsBN, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBN.rda")

#######################################################################################################
# Nu voor 2 tot en emt 5 
#######################################################################################################

for(p in 2:5){
grab <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",p)
filesinmap <- list.files(grab, full.names = T)
filesnames <- list.files(grab)
filesnames <- gsub(".rda", "", filesnames)

networklist <- list()
names <- c()

for (i in 1:length(filesinmap)){
  variablepos <- get(load(filesinmap[i]))
  networklist[[i]] <- variablepos
}
names(networklist) <- filesnames
rm(list = filesnames)
length(networklist)


# Convert the graphs to undirected igraphs
graphsNEL <- lapply(networklist, as.graphNEL) 
igraphsdir <- lapply(graphsNEL, igraph.from.graphNEL)
igraphsske<- lapply(igraphsdir, as.undirected)
edgelists <- lapply(igraphsske, E)
nedgeslists <- sapply(edgelists, length)
length(unique(nedgeslists))
indsame <- c()
# Sort the graphs from small to big

for (i in 1:length(unique(nedgeslists))){
  int <- which(nedgeslists == sort(nedgeslists)[i])
  if(length(int)>1){
  indsame[i] <- as.vector(int[2])} else {indsame[i] <- int}
}
length(indsame)
igraphsskesort <- igraphsske[indsame]

assign(paste0("perm",p,"sort"), networklist[indsame])

# Create the graphObjects as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(igraphsskesort))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- igraphsskesort[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsskesort[[i]])
}

names(graphObjects) <- names(networklist)[indsame]
assign(paste0("gridGraphsBNperm",p),graphObjects)
length(gridGraphsBNperm2)

save(list = paste0("perm",p,"sort"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",p,"sort.rda"))
save(list = paste0("gridGraphsBNperm",p), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBNperm",p,".rda"))
}

######################################################################
# Nu voor 6 tot en met 10
######################################################################
 
 for(p in 6:10){
   grab <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",p)
   filesinmap <- list.files(grab, full.names = T)
   filesnames <- list.files(grab)
   filesnames <- gsub(".rda", "", filesnames)
   
   networklist <- list()
   names <- c()
   
   for (i in 1:length(filesinmap)){
     variablepos <- get(load(filesinmap[i]))
     networklist[[i]] <- variablepos
   }
   names(networklist) <- filesnames
   rm(list = filesnames)
   length(networklist)
   
   
   # Convert the graphs to undirected igraphs
   graphsNEL <- lapply(networklist, as.graphNEL) 
   igraphsdir <- lapply(graphsNEL, igraph.from.graphNEL)
   igraphsske<- lapply(igraphsdir, as.undirected)
   edgelists <- lapply(igraphsske, E)
   nedgeslists <- sapply(edgelists, length)
   length(unique(nedgeslists))
   indsame <- c()
   # Sort the graphs from small to big
   
   for (i in 1:length(unique(nedgeslists))){
     int <- which(nedgeslists == sort(nedgeslists)[i])
     if(length(int)>1){
       indsame[i] <- as.vector(int[2])} else {indsame[i] <- int}
   }
   length(indsame)
   igraphsskesort <- igraphsske[indsame]
   
   assign(paste0("perm",p,"sort"), networklist[indsame])
   
   # Create the graphObjects as would have been obtained by graph_from_Grid
   graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
   graphObjects <- rep(list(graphObject),length(igraphsskesort))
   for (i in 1:length(graphObjects)){
     graphObjects[[i]]$graph <- igraphsskesort[[i]]
     graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsskesort[[i]])
   }
   
   names(graphObjects) <- names(networklist)[indsame]
   assign(paste0("gridGraphsBNperm",p),graphObjects)
   
   save(list = paste0("perm",p,"sort"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",p,"sort.rda"))
   save(list = paste0("gridGraphsBNperm",p), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBNperm",p,".rda"))
 }
#############################################################################
# Nu speciale 11 tot en met 25 (Hebben er minder)
#############################################################################
 for(p in 11:25){
   grab <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",p)
   filesinmap <- list.files(grab, full.names = T)
   filesnames <- list.files(grab)
   filesnames <- gsub(".rda", "", filesnames)
   
   networklist <- list()
   names <- c()
   
   for (i in 1:length(filesinmap)){
     variablepos <- get(load(filesinmap[i]))
     networklist[[i]] <- variablepos
   }
   names(networklist) <- filesnames
   rm(list = filesnames)
   length(networklist)
   
   
   # Convert the graphs to undirected igraphs
   graphsNEL <- lapply(networklist, as.graphNEL) 
   igraphsdir <- lapply(graphsNEL, igraph.from.graphNEL)
   igraphsske<- lapply(igraphsdir, as.undirected)
   edgelists <- lapply(igraphsske, E)
   nedgeslists <- sapply(edgelists, length)
   length(unique(nedgeslists))
   indsame <- c()
   # Sort the graphs from small to big
   
   for (i in 1:length(unique(nedgeslists))){
     int <- which(nedgeslists == sort(nedgeslists)[i])
     if(length(int)>1){
       indsame[i] <- as.vector(int[2])} else {indsame[i] <- int}
   }
   length(indsame)
   igraphsskesort <- igraphsske[indsame]
   
   assign(paste0("perm",p,"sort"), networklist[indsame])
   
   # Create the graphObjects as would have been obtained by graph_from_Grid
   graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
   graphObjects <- rep(list(graphObject),length(igraphsskesort))
   for (i in 1:length(graphObjects)){
     graphObjects[[i]]$graph <- igraphsskesort[[i]]
     graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsskesort[[i]])
   }
   
   names(graphObjects) <- names(networklist)[indsame]
   assign(paste0("gridGraphsBNperm",p),graphObjects)
   
   save(list = paste0("perm",p,"sort"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",p,"sort.rda"))
   save(list = paste0("gridGraphsBNperm",p), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBNperm",p,".rda"))
 }
 