##################################################################
# Maps Erdos-renyi, lattice, bn, cn
# assortativity global
# assortativity degree
# diameter of giant component
##################################################################
setwd("~/data/Untitled/Trabajo/R_practice/R/")
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")
setwd("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")


# Load CN_constructionandMEasuresFunctions
# Load Basic NEtwork functions
# Load tas ncep 10d
# OCEANO:
rm(list = ls())
library(bnlearn)
source("Functions/BasicNetworkFunctions.R")
source("Functions/CN_ConstructionandMeasuresFunctions.R")
source("Functions/propagationFunctions.R")
load("../Data/tas_ncep_10d.rda")
load("../Data/gridGraphs.rda")

# MACBOOK: 
# load("/Volumes/ubuntu/Trabajo/R_practice/Data/tas_ncep_10d.rda")
# load(file = "/Volumes/ubuntu/Trabajo/R_practice/Data/gridGraphs.rda")
# load libraries PaperScript
# IFCA: 
# load("/media/catharina/ubuntu/Trabajo/R_practice/Data/tas_ncep_10d.rda")
# load("/media/catharina/ubuntu/Trabajo/R_practice/Data/gridGraphs.rda")
# load libraries PaperScript
#################################################################
# Made 100 correlation graphs
# Saved in ubuntu
# loaded above
# renamed
#################################################################
# MAke list correlation networks and load
# taus <- seq(from = 0, to = 1, length.out = 101)
# gridGraphs <- lapply(taus,graph_from_Grid, grid = tas_ncep_10d)
# save(gridGraphs, file = "/Volumes/ubuntu/Trabajo/R_practice/Data/gridGraphs.rda")
gridGraphsCN <- gridGraphs
rm(gridGraphs)
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- lapply(edgelistsCN, length)
names(gridGraphsCN) <- as.character(numberofedgesCN)

##################################################################################
# Load Bayesian Network perm 1
# Create gridGraphsBN
##################################################################################
# IFCA: 
# filesinmap <- list.files("/media/catharina/ubuntu/outofoffice/hc_perm1", full.names = T)
# filesnames <- list.files("/media/catharina/ubuntu/outofoffice/hc_perm1")
# MACBOOK:
# filesinmap <- list.files("/Volumes/ubuntu/outofoffice/hc_perm1", full.names = T)
# filesnames <- list.files("/Volumes/ubuntu/outofoffice/hc_perm1")
# OCEANO
filesinmap <- list.files("../Data/Struct_learn/hciterations/perm1", full.names = T)
filesnames <- list.files("../Data/Struct_learn/hciterations/perm1")



filesnames <- gsub(".rda", "", filesnames)
filesnames

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

# save(perm1sort, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
# save(gridGraphsBN, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBN.rda")
#################################################################
# Create random erdös-Renyi graph
# Check other possibility in BNLEARN !
#################################################################
edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- lapply(edgelistsBN, length)
part2 <- numberofedgesCN[0:35]
part2 <- rev(part2)

edgesRenyi <- sapply(edgelistsBN, length)
edgesRenyi <- c(numberofedgesBN,part2)

renyigraphs <- lapply(edgesRenyi, sample_gnm, n = 648)
renyiObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
renyiObjects <- rep(list(renyiObject),length(renyigraphs))
for (i in 1:length(renyiObjects)){
  renyiObjects[[i]]$graph <- renyigraphs[[i]]
  renyiObjects[[i]]$adjacency <- as_adjacency_matrix(renyigraphs[[i]])
}

gridGraphsRenyi <- renyiObjects
names(gridGraphsRenyi) <- as.character(edgesRenyi)

# #################################################################
# # Create random regular graph
# # Check other possibility in BNLEARN !
# #################################################################
# edgesRegular <- sapply(edgelistsBN, length)
# degreesRegular <- sapply(edgesRegular, function(x){(2*x)/648})
# degreesRegular <- ceiling(degreesRegular)
# degreesRegular <- 1:25
# regulargraphs <- lapply(degreesRegular, sample_k_regular, no.of.nodes = 648)
# regulargraphs[[1]]
# regularObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
# regularObjects <- rep(list(regularObject),length(regulargraphs))
# for (i in 1:length(regularObjects)){
#   regularObjects[[i]]$graph <- regulargraphs[[i]]
#   regularObjects[[i]]$adjacency <- as_adjacency_matrix(regulargraphs[[i]])
# }
# 
# gridGraphsRegular <- regularObjects
# names(gridGraphsRegular) <- as.character(edgesRegular)
#################################################################
# Create deterministic lattices graph
# Check other possibility in BNLEARN !
#################################################################
g <- graph.lattice( c(18,36), circular = FALSE )
plot(g)
latticesgraphs <- lapply(1:25, connect.neighborhood, graph = g)
# # lapply(latticesgraphs, layout_on_sphere)
# layout(matrix(1:4, nrow=2, byrow=TRUE))
# sapply(graphsLattices[1:4], plot.igraph, vertex.label=NA)
edgelistsLattices <- lapply(latticesgraphs, E)
numberofedgesLattices <- lapply(edgelistsLattices, length)
names(latticesgraphs) <- as.character(numberofedgesLattices)

latticesObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
latticesObjects <- rep(list(latticesObject),length(latticesgraphs))
for (i in 1:length(latticesObjects)){
  latticesObjects[[i]]$graph <- latticesgraphs[[i]]
  latticesObjects[[i]]$adjacency <- as_adjacency_matrix(latticesgraphs[[i]])
}

gridGraphsLattices <- latticesObjects
names(gridGraphsLattices) <- as.character(numberofedgesLattices)

dev.off()

#################################################################
# Create small world graph
# Check other possibility in BNLEARN !
#################################################################
smallworldgraphs <- lapply(latticesgraphs, rewire, with = each_edge(prob = 0.05))


edgelistsSmallworlds <- lapply(smallworldgraphs, E)
numberofedgesSmallworlds <- lapply(smallworldgraphs, length)
  
names(smallworldgraphs) <- as.character(numberofedgesSmallworlds)

smallworldsObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
smallworldsObjects <- rep(list(smallworldsObject),length(smallworldgraphs))
for (i in 1:length(smallworldsObjects)){
  smallworldsObjects[[i]]$graph <- smallworldgraphs[[i]]
  smallworldsObjects[[i]]$adjacency <- as_adjacency_matrix(smallworldgraphs[[i]])
}

gridGraphsSmallworlds <- smallworldsObjects
names(gridGraphsSmallworlds) <- as.character(numberofedgesSmallworlds)

dev.off()

#################################################################
# Plots: igraph always allocate at different posititions!
#################################################################
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})
GraphsRenyi <- lapply(gridGraphsRenyi, function(x){x$graph})
# GraphsRegular <- lapply(gridGraphsRegular, function(x){x$graph})
GraphsLattices <- lapply(gridGraphsLattices,function(x){x$graph})
GraphsSmallworlds <- lapply(gridGraphsSmallworlds, function(x){x$graph})
names(GraphsCN) <- names(gridGraphsCN)
names(GraphsBN) <- names(gridGraphsBN)
names(GraphsRenyi) <- names(gridGraphsRenyi)
# names(GraphsRegular) <- names(gridGraphsRegular)
names(GraphsLattices) <- names(gridGraphsLattices)
names(GraphsSmallworlds) <- names(gridGraphsSmallworlds)

par(mfrow = c(2,2))
load("../Data/lisworld.rda")

plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),meteodag = gridGraphsRenyi$`1398`$graph, lis = TRUE)
plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),meteodag = gridGraphsLattices$`3620`$graph, lis = TRUE)
plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),meteodag = gridGraphsBN$hc1_1600_1700i$graph, lis = TRUE)
plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),meteodag = gridGraphsCN$`1512`$graph, lis = TRUE)
# plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),meteodag = gridGraphsSmallworlds$, lis = TRUE)

#PAra la thesis tambien
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/graphvariety.pdf")
pdf(plotname, height = 7, width = 9)
par(mfrow = c(2,2))
mar = c(0,0,0,0)
plot_long_distances(data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),dag = gridGraphsRenyi$`1497`$graph, minimdist = 10000, title = "B")
plot_long_distances(data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),dag = gridGraphsLattices$`3620`$graph, minimdist = 10000, title = "B")
plot_long_distances(data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),dag = gridGraphsBN$hc1_1400_1500i$graph, minimdist = 10000, title = "B")
plot_long_distances(data.dag  = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),dag = gridGraphsCN$`1512`$graph, minimdist = 10000, title = "B")
dev.off()

# Para la thesis:
plot_long_distances(data.dag  = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),dag = gridGraphsCN$`3239`$graph, minimdist = 10000, title = "B")
#################################################################
# Calculate measures
# global
# local
# degree
#################################################################
names(gridGraphsLattices)
names <- c("CN","BN",'Renyi',"Lattices")
gridGraphlists <- list(gridGraphsCN,gridGraphsBN,gridGraphsRenyi,gridGraphsLattices)


# without weights

for(i in 1:length(names)){
  name <- names[i]
  gridGraphlist <- gridGraphlists[[i]]
  locallist <- lapply(gridGraphlist, graph2measure.local)
  globallist <- mapply(graph2measure.global, graphObj = gridGraphlist, localmeasureObj = locallist, SIMPLIFY = FALSE)
  assign(paste0("globmeas",name),globallist)
  assign(paste0("globalasslist",name),lapply(globallist, function(x){x$globalass}))
  assign(paste0("globalcluslist",name),lapply(globallist, function(x){x$globalclus}))
  assign(paste0("globalbarratcluslist",name),lapply(globallist, function(x){x$globalbarratclus}))
  assign(paste0("globalbarratweightcluslist",name),lapply(globallist, function(x){x$globalbarratweightclus}))
  assign(paste0("globaledenslist",name), lapply(globallist, function(x){x$edens}))
  assign(paste0("globaldiamLCClist",name), lapply(globallist,function(x){x$diameterLCC}))
  assign(paste0("globaldiamGlist",name), lapply(globallist,function(x){x$diameterG}))
  assign(paste0("meandistLCClist",name), lapply(globallist,function(x){x$meandistLCC}))
  assign(paste0("meandistGlist",name), lapply(globallist,function(x){x$meandistG}))
  assign(paste0("meanavpathlist",name), lapply(globallist,function(x){x$meanavpathlength}))
  globallist <- NULL
  
  degreelist <- mapply(graph2measure.degree,graphObj = gridGraphlist, localmeasureObj = locallist, SIMPLIFY = FALSE)
  assign(paste0("degreemeas",name),degreelist)
  assign(paste0("degreeclustering",name),lapply(degreelist,function(x){x$degreeclustering}))
  assign(paste0("degreedist",name), lapply(degreelist, function(x){x$degrdist}))
  assign(paste0("meanstrength",name), lapply(degreelist, function(x){x$meanstrength}))
  assign(paste0("meandistance",name), lapply(degreelist, function(x){x$meandistance}))
  degreelist <- NULL
}





####################################################################
# set names
####################################################################
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(GraphsCN) <- as.character(numberofedgesCN)

edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- sapply(edgelistsBN, length)
names(GraphsBN) <- as.character(numberofedgesBN)

edgelistsRenyi <- lapply(GraphsRenyi, E)
numberofedgesRenyi <- lapply(edgelistsRenyi, length)
names(GraphsRenyi) <- as.character(numberofedgesRenyi)

edgelistsLattices <- lapply(latticesgraphs, E)
numberofedgesLattices <- lapply(edgelistsLattices, length)
names(latticesgraphs) <- as.character(numberofedgesLattices)
######################################################################
# global assortivity
######################################################################

# BN CN CN
# plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalass.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalass.pdf")
pdf(plotname)
par(mfrow = c(1,2))
# plot(numberofedgesBN,globalasslistBN, col = "red", xlab = "edges", ylab = "global assortativity")
# # plot(numberofedgesRegular,globalasslistRegular) doet het niet
# points(numberofedgesRenyi,globalasslistRenyi, col = "blue")
# points(numberofedgesLattices,globalasslistLattices, col = "green")
# m1 <- seq(0,16,4)
# text(numberofedgesBN[m1],globalasslistBN[m1], pos = 4,labels = names(GraphsBN)[m1])
plot(numberofedgesCN,globalasslistCN, xlab = "edges", ylab = "global assortativity")
points(numberofedgesBN,globalasslistBN, col = "red")
points(numberofedgesRenyi,globalasslistRenyi, col = "blue")
points(numberofedgesLattices,globalasslistLattices, col = "green")
m2 <- seq(0,60,10)
# text(numberofedgesCN[m2],globalasslistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
plot(numberofedgesCN,globalasslistCN, xlim = c(0,10000), xlab = "edges", ylab = "global assortativity")
points(numberofedgesBN,globalasslistBN, col = "red",xlab = "edges", ylab = "global assortativity")
points(numberofedgesRenyi,globalasslistRenyi, col = "blue")
points(numberofedgesLattices,globalasslistLattices, col = "green")

# text(numberofedgesCN[m2],globalasslistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
dev.off()


###########################################################################
# diameter
###########################################################################
# plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/diameter.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/diameter.pdf")
pdf(plotname)
par(mfrow = c(1,2))
# plot(numberofedgesBN,globaldiamLCClistBN, col = "indianred", xlab = "edges", ylab = "diameter D")
# points(numberofedgesBN,globaldiamGlistBN, col = "red")
# points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
# points(numberofedgesRenyi,globaldiamGlistRenyi, col = "blue")
# points(numberofedgesLattices,globaldiamLCClistLattices, col = "lightgreen")
# points(numberofedgesLattices,globaldiamGlistLattices,col = "green")
# text(numberofedgesLattices,globaldiamGlistLattices, pos = 4,labels = names(GraphsLattices))
# m1 <- seq(0,16,4)
# text(numberofedgesBN[m1],globaldiamLCClistBN[m1], pos = 4, labels = names(GraphsBN)[m1])

plot(numberofedgesCN,globaldiamLCClistCN, col = "grey", xlab = "edges", ylab = "diameter D")
points(numberofedgesCN,globaldiamGlistCN, col = "black")
points(numberofedgesBN,globaldiamGlistBN, col = "red")
# points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,globaldiamGlistRenyi, col = "blue")
points(numberofedgesLattices,globaldiamLCClistLattices, col = "lightgreen")
points(numberofedgesLattices,globaldiamGlistLattices,col = "green")
m2 <- seq(0,60,10)
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])

plot(numberofedgesCN,globaldiamLCClistCN, col = "grey", xlab = "edges", ylab = "diameter D",xlim = c(0,60000) )
points(numberofedgesCN,globaldiamGlistCN, col = "black")
points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,globaldiamGlistRenyi, col = "blue")
points(numberofedgesBN,globaldiamLCClistBN, col = "indianred", xlab = "edges", ylab = "diameter D")
points(numberofedgesBN,globaldiamGlistBN, col = "red")
points(numberofedgesLattices,globaldiamGlistLattices,col = "green")
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
dev.off()

#################################################################################################
# Mean absolute distance (always unweighted)
#################################################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/meandistance.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/meandistance.pdf")
pdf(plotname)
par(mfrow = c(1,2))
# plot(numberofedgesBN,meandistLCClistBN, col = "indianred", xlab = "edges", ylab = "meandistance D")
# points(numberofedgesCN,meandistGlistCN, col = "black")
# points(numberofedgesBN,meandistGlistBN, col = "red")
# points(numberofedgesRenyi,meandistLCClistRenyi, col = "lightblue")
# points(numberofedgesRenyi,meandistGlistRenyi, col = "blue")
# points(numberofedgesLattices,meandistLCClistLattices, col = "lightgreen")
# points(numberofedgesLattices,meandistGlistLattices,col = "green")
# text(numberofedgesLattices,meandistGlistLattices, pos = 4,labels = names(GraphsLattices))
# m1 <- seq(0,16,4)
# text(numberofedgesBN[m1],globaldiamLCClistBN[m1], pos = 4, labels = names(GraphsBN)[m1])

plot(numberofedgesCN,meandistLCClistCN, col = "grey", xlab = "edges", ylab = "meandistance L(|A|)")
points(numberofedgesCN,meandistGlistCN, col = "black")
points(numberofedgesBN,meandistGlistBN, col = "red")
points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,meandistGlistRenyi, col = "blue")
points(numberofedgesLattices,meandistLCClistLattices, col = "lightgreen")
points(numberofedgesLattices,meandistGlistLattices,col = "green")
m2 <- seq(0,60,10)
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])

plot(numberofedgesCN,meandistLCClistCN, col = "grey", xlab = "edges", ylab = "meandistance L(|A|)",xlim = c(0,10000) )
points(numberofedgesCN,meandistGlistCN, col = "black")
points(numberofedgesRenyi,meandistLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,meandistGlistRenyi, col = "blue")
points(numberofedgesBN,meandistLCClistBN, col = "indianred", xlab = "edges", ylab = "meandistance L(|A|)")
points(numberofedgesBN,meandistGlistBN, col = "red")
points(numberofedgesLattices,meandistGlistLattices,col = "green")
# text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
dev.off()


#################################################################################################
# Mean av path length absolute  = mean distance (absolute) for unweighted. 
#################################################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/avpathlength.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/avpathlength.pdf")
pdf(plotname)
par(mfrow = c(1,2))
# plot(numberofedgesBN,meanavpathlistBN, col = "indianred", xlab = "edges", ylab = "meanavpathlength L(|A|)", xlim = c(0,8000), ylim = c(0,20))
# points(numberofedgesBN,meandistGlistBN, col = "red")
# points(numberofedgesCN,meandistGlistCN, col = "black")
# points(numberofedgesRenyi,meanavpathlistRenyi, col = "lightblue")
# points(numberofedgesRenyi,meandistGlistRenyi, col = "blue")
# points(numberofedgesLattices,meanavpathlistLattices, col = "lightgreen")
# points(numberofedgesLattices,meandistGlistLattices,col = "green")
# text(numberofedgesLattices,meandistGlistLattices, pos = 4,labels = names(GraphsLattices))
# m1 <- seq(0,16,4)
# text(numberofedgesBN[m1],globaldiamLCClistBN[m1], pos = 4, labels = names(GraphsBN)[m1])

plot(numberofedgesCN,meanavpathlistCN, col = "grey", xlab = "edges", ylab = "meanavpathlength L(|A|)", ylim = c(0,25))
points(numberofedgesCN,meandistGlistCN, col = "black")
points(numberofedgesBN,meanavpathlistBN, col = "red")
# points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "blue")
points(numberofedgesLattices,meanavpathlistLattices, col = "lightgreen")
points(numberofedgesLattices,meandistGlistLattices,col = "green")
m2 <- seq(0,60,10)
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])

plot(numberofedgesCN,meanavpathlistCN, col = "grey", xlab = "edges", ylab = "meanavpathlength L(|A|)",xlim = c(0,10000), ylim = c(0,20))
points(numberofedgesCN,meanavpathlistCN, col = "black")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "lightblue")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "blue")
points(numberofedgesBN,meanavpathlistBN, col = "indianred", xlab = "edges", ylab = "meanavpathlength L(|A|)")
points(numberofedgesBN,meanavpathlistBN, col = "red")
points(numberofedgesLattices,meanavpathlistLattices,col = "green")
# text(numberofedgesCN[m2],meanavpathlistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
dev.off()



###################################################################
# Hier global clustering (TRIANGLES)
####################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalclust_tri_unw.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalclust_tri_unw.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalclust_3diff.pdf")
pdf(plotname)
par(mfrow = c(3,2))
par(mfrow = c(1,2))
# plot(numberofedgesBN,globalcluslistBN, ylim = c(0,0.25), xlab = "edges", ylab = "global clustering C", col = "red")
# points(numberofedgesRenyi,globalcluslistRenyi, col = "lightblue")
# points(numberofedgesLattices,globalcluslistLattices,col = "green")
# text(numberofedgesBN[m1],globalcluslistBN[m1], labels = names(GraphsBN)[m1], pos = 4)

plot(numberofedgesCN,globalcluslistCN, xlab = "edges", ylab = "global clustering C",  ylim = c(0,1))
points(numberofedgesRenyi,globalcluslistRenyi, col = "lightblue")
points(numberofedgesLattices,globalcluslistLattices,col = "green")
points(numberofedgesBN,globalcluslistBN, col = "red")


plot(numberofedgesCN,globalcluslistCN, xlim = c(0,10000), ylim = c(0,1), xlab = "edges", ylab = "global clustering C")
points(numberofedgesBN,globalcluslistBN, col = "red")
points(numberofedgesRenyi,globalcluslistRenyi, col = "lightblue")
points(numberofedgesLattices,globalcluslistLattices,col = "green")
text(numberofedgesCN[m2],globalcluslistCN[m2], pos = 3,labels = names(GraphsCN)[m2])
dev.off()


###################################################################
# Hier global clustering (BARRAT)
####################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalclust_bar_unw.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/globalclust_bar_unw.pdf")
pdf(plotname)
par(mfrow = c(2,2))
# plot(numberofedgesBN,globalbarratcluslistBN, ylim = c(0,0.25), xlab = "edges", ylab = "global barrat clustering C", col = "red")
# points(numberofedgesRenyi,globalbarratcluslistRenyi, col = "lightblue")
# points(numberofedgesLattices,globalbarratcluslistLattices,col = "green")
# text(numberofedgesBN[m1],globalbarratcluslistBN[m1], labels = names(GraphsBN)[m1], pos = 4)

plot(numberofedgesCN,globalbarratcluslistCN, xlab = "edges", ylab = "global clustering C", ylim = c(0,1))
points(numberofedgesRenyi,globalbarratcluslistRenyi, col = "lightblue")
points(numberofedgesLattices,globalcluslistLattices,col = "green")
points(numberofedgesBN,globalbarratcluslistBN, col = "red")
text(numberofedgesCN,globaldiamLCClistCN, pos = 4,labels = names(GraphsCN))

plot(numberofedgesCN,globalbarratcluslistCN, xlim = c(0,10000), ylim = c(0,1), xlab = "edges", ylab = "global clustering C")
points(numberofedgesBN,globalbarratcluslistBN, col = "red")
points(numberofedgesRenyi,globalbarratcluslistRenyi, col = "lightblue")
points(numberofedgesLattices,globalbarratcluslistLattices,col = "green")
points(numberofedgesBN,globalbarratcluslistBN, col = "red")
text(numberofedgesLattices,globalbarratcluslistLattices, pos = 3,labels = names(GraphsLattices))
dev.off()

plot(numberofedgesCN,globalbarratweightcluslistCN, xlab = "edges", ylab = "global clustering C", ylim = c(0,1))
points(numberofedgesRenyi,globalbarratweightcluslistRenyi, col = "lightblue")
points(numberofedges,globalbarratweightcluslistLattices,col = "green")
points(numberofedgesBN,globalbarratweightcluslistBN, col = "red")
text(numberofedgesCN,globaldiamLCClistCN, pos = 4,labels = names(GraphsCN))

plot(numberofedgesCN,globalbarratweightcluslistCN, xlim = c(0,10000),  ylim = c(0,1), xlab = "edges", ylab = "global clustering C")
points(numberofedgesBN,globalbarratweightcluslistBN, col = "red")
points(numberofedgesRenyi,globalbarratweightcluslistRenyi, col = "lightblue")
points(numberofedgesLattices,globalbarratweightcluslistLattices,col = "green")

text(numberofedgesLattices,globalbarratcluslistLattices, pos = 3,labels = names(GraphsLattices))

all.equal(globalbarratweightcluslistBN,globalbarratcluslistBN)
all.equal(globalbarratweightcluslistCN,globalbarratcluslistCN)
dev.off()
##################################################################
# How many friends do your friends have? "local assortativity"
# How many friends on average do people with k friends have? "k degree assortativity"
##################################################################
knnBN <- lapply(GraphsBN, knn)
names(knnBN) <- as.character(numberofedgesBN)
knnkBN <- lapply(knnBN, function(x){x$knnk})

knnCN <- lapply(GraphsCN, knn)
names(knnCN) <- names(GraphsCN)
knnkCN <- lapply(knnCN, function(x){x$knnk})

knnRenyi <- lapply(GraphsRenyi, knn)
names(knnRenyi) <- as.character(numberofedgesRenyi)
knnkRenyi <- lapply(knnRenyi, function(x){x$knnk})

knnLattices<- lapply(GraphsLattices, knn)
names(knnLattices) <- as.character(numberofedgesLattices)
knnkLattices <- lapply(knnLattices, function(x){x$knnk})
#####################################################################
# plots knnk
#####################################################################
plotseqBN <- seq(6,length(knnkBN),length.out = 8)
lastBN <- plotseqBN[length(plotseqBN)]
restBN <- plotseqBN[1:(length(plotseqBN)-1)]
colrestBN <- rainbow(length(restBN))
numberofedgesBN[plotseqBN]

plotseqRenyi <- plotseqBN
lastRenyi <- lastBN
restRenyi <- restBN
colrestRenyi <- rainbow(length(restRenyi))

plotseqLattices <- seq(1,length(knnkLattices),length.out = 5)
plotseqLattices <- seq(1,5,1)
lastLattices <- plotseqLattices[length(plotseqLattices)]
restLattices <- plotseqLattices[1:(length(plotseqLattices)-1)]
colrestLattices <- rainbow(length(restLattices))
numberofedgesLattices[plotseqLattices]

#plotseqCN <- seq(1,length(knnkCN),length.out = 5)
# plotseqCN <- seq(1,20,1)
plotseqCN <- seq(25,85,length.out = 10)
# plotseqCN <- seq(8,38,length.out = 5)
plotseqCN <- round(plotseqCN)
plotseqCN <- rev(plotseqCN)
lastCN <- plotseqCN[length(plotseqCN)]
restCN <- plotseqCN[1:(length(plotseqCN)-1)]
colrestCN <- rainbow(length(restCN))
numberofedgesCN[plotseqCN]

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/degreeassunw.pdf")
pdf(plotname)
drawdegree = FALSE
par(mfrow = c(2,2))

plot(1:length(knnkRenyi[[lastRenyi]]),knnkRenyi[[lastRenyi]], 
     xlim = c(0,100),ylim = c(0,50), 
     xlab = "degree k", ylab = "k_nn,k",
     main = "Random Renyi")
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}
legend(x = "topright",legend = names(knnkRenyi)[plotseqRenyi], col = c(colrestRenyi,"black"), pch = 1, cex = 0.5 )

plot(1:length(knnkLattices[[lastLattices]]),knnkLattices[[lastLattices]],
     xlim = c(0,100),ylim = c(0,50), 
     xlab = "degree k", ylab = "k_nn,k",
     main = "Regular Lattices")
abline(0,1)
for (i in 1:length(restLattices)){
  j <- restLattices[i]
  points(1:length(knnkLattices[[j]]),knnkLattices[[j]],col = colrestLattices[i])
  if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
  text(x= length(knnkLattices[[j]]), y=mean(knnkLattices[[j]],na.rm = TRUE), pos=4, labels= names(knnLattices)[j])
}
legend(x = "topright",legend = names(knnkLattices)[plotseqLattices], col = colrestLattices, pch = 1, cex = 0.5 )
# plot(1:length(knnkLattices[[lastLattices]]),knnkLattices[[lastLattices]],xlim = c(0,700),ylim = c(0,550), xlab = "degree", ylab = "k_nn,k")
# for (i in 1:length(restLattices)){
#   j <- restLattices[i]
#   points(1:length(knnkLattices[[j]]),knnkLattices[[j]],col = colrestLattices[i])
#   if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
#   text(x= length(knnkLattices[[j]]), y=mean(knnkLattices[[j]],na.rm = TRUE), pos=4, labels= names(knnLattices)[j])
# }


plot(1:length(knnkBN[[lastBN]]),knnkBN[[lastBN]], 
     xlim = c(0,100),ylim = c(0,50), 
     xlab = "degree k", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(knnkBN[[j]]),knnkBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  text(x= length(knnkBN[[j]]), y=mean(knnkBN[[j]],na.rm = TRUE), pos=4, labels= names(knnBN)[j])
}
legend(x = "topright",legend = names(knnkBN)[plotseqBN], col = c(colrestBN, "black"), pch = 1, cex = 0.5 )

# plot(1:length(knnkCN[[lastCN]]),knnkCN[[lastCN]], xlim = c(0,800),ylim = c(0,800), xlab = "degree", ylab = "k_nn,k")
# text(x= length(knnkCN[[lastCN]]), y=mean(knnkCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(knnCN)[lastCN])
# abline(0,1)
# for (i in 1:length(restCN)){
#   j <- restCN[i]
#   points(1:length(knnkCN[[j]]),knnkCN[[j]],col = colrestCN[i])
#   if(drawdegree == TRUE){lines(1:length(degreedistCN[[j]]),degreedistCN[[j]]*648,col = colrestCN[i])}
#   text(x= length(knnkCN[[j]]), y=mean(knnkCN[[j]],na.rm = TRUE), pos=4, labels= names(knnCN)[j])
# }

plot(1:length(knnkCN[[lastCN]]),knnkCN[[lastCN]], 
     # xlim = c(0,150),
     # ylim = c(0,100), 
     xlab = "degree k", ylab = "k_nn,k",
     main = "Correlation Networks")
abline(0,1)
# text(x= length(knnkCN[[lastCN]])+30, y=max(knnkCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(knnCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(knnkCN[[j]]),knnkCN[[j]],col = colrestCN[i])
  if(drawdegree == TRUE){lines(1:length(degreedistCN[[j]]),degreedistCN[[j]]*648,col = colrestCN[i])}
  # text(x= length(knnkCN[[j]])+50, y=max(knnkCN[[j]],na.rm = TRUE), pos=4, labels= names(knnCN)[j])
}
legend(x = "bottomright",legend = names(knnkCN)[plotseqCN], col = c(colrestCN,"black"), pch = 1, cex = 0.5 )
dev.off()

#####################################################################################
# degree clustering
# How do vertices of degree k cluster? 
#####################################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/degreeclustunw.pdf")
pdf(plotname)
par(mfrow = c(2,2))

drawdegree = FALSE


plot(1:length(degreeclusteringRenyi[[lastRenyi]]),degreeclusteringRenyi[[lastRenyi]], 
     # xlim = c(0,50),
     ylim = c(0,0.6), 
      xlab = "degree", ylab = "C(k) Barrat",
     main = "Random Renyi")
if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[lastRenyi]]),degreedistRenyi[[lastRenyi]]*648)}
# text(x= length(degreeclusteringRenyi[[lastRenyi]]), y=mean(degreeclusteringRenyi[[lastRenyi]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[lastRenyi])
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(degreeclusteringRenyi[[j]]),degreeclusteringRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  # text(x= length(degreeclusteringRenyi[[j]]), y=mean(degreeclusteringRenyi[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[j])
}
legend(x = "topright",legend = names(degreeclusteringRenyi)[plotseqRenyi], col = c(colrestRenyi,"black"), pch = 1, cex = 0.5 )

plot(1:length(degreeclusteringLattices[[lastLattices]]),degreeclusteringLattices[[lastLattices]],
     # xlim = c(0,50),
     ylim = c(0,1), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Regular Lattices")
if(drawdegree == TRUE){lines(1:length(degreedistLattices[[lastLattices]]),degreedistLattices[[lastLattices]]*648)}
for (i in 1:length(restLattices)){
  j <- restLattices[i]
  points(1:length(degreeclusteringLattices[[j]]),degreeclusteringLattices[[j]],col = colrestLattices[i])
  if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
  # text(x= length(degreeclusteringLattices[[j]]), y=mean(degreeclusteringLattices[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringLattices)[j])
}
legend(x = "topright",legend = names(degreeclusteringLattices)[plotseqLattices], col = c(colrestLattices,"black"), pch = 1, cex = 0.5 )


plot(1:length(degreeclusteringBN[[lastBN]]),degreeclusteringBN[[lastBN]], 
     xlim = c(0,100),ylim = c(0,0.6), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Bayesian Networks")
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
# text(x= length(degreeclusteringBN[[lastBN]]), y=mean(degreeclusteringBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(GraphsBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(degreeclusteringBN[[j]]),degreeclusteringBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  # text(x= length(degreeclusteringBN[[j]]), y=mean(degreeclusteringBN[[j]],na.rm = TRUE), pos=4, labels= names(GraphsBN)[j])
}
legend(x = "topright",legend = names(knnBN)[plotseqBN], col = c(colrestBN,"black"), pch = 1, cex = 0.5 )


plot(1:length(degreeclusteringCN[[lastCN]]),degreeclusteringCN[[lastCN]], 
     xlim = c(0,120),ylim = c(0,1), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Correlation Network")
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(degreeclusteringCN[[j]]),degreeclusteringCN[[j]],col = colrestCN[i])
  
  # text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}
legend(x = "topright",legend = names(degreeclusteringCN)[plotseqCN], col = c(colrestCN,"black"), pch = 1, cex = 0.5 )
dev.off()

#####################################################################################
# strength per degree
#####################################################################################
plot(1:length(meanstrengthBN[[lastBN]]),meanstrengthBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(meanstrengthBN[[lastBN]]), y=mean(meanstrengthBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(numberofedgesBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(meanstrengthBN[[j]]),meanstrengthBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthBN[[j]]),meanstrengthBN[[j]]*648,col = colrestBN[i])}
  text(x= length(meanstrengthBN[[j]]), y=mean(meanstrengthBN[[j]],na.rm = TRUE), pos=4, labels= names(meanstrengthBN)[j])
}






plot(1:length(meanstrengthCN[[lastCN]]),meanstrengthCN[[lastCN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meanstrengthCN[[lastCN]]), y=mean(meanstrengthCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]]*648,col = colrestCN[i])}
  text(x= length(meanstrengthCN[[j]]), y=mean(meanstrengthCN[[j]],na.rm = TRUE), pos=4, labels= names(meanstrengthCN)[j])
}
###########################################################################
# mean distances (= mean strength per degree)
###########################################################################

plot(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(meandistanceBN[[lastBN]]), y=mean(meandistanceBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(numberofedgesBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(meandistanceBN[[j]]),meandistanceBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meandistanceBN[[j]]),meandistanceBN[[j]]*648,col = colrestBN[i])}
  text(x= length(meandistanceBN[[j]]), y=mean(meandistanceBN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceBN)[j])
}



plot(1:length(meandistanceCN[[lastCN]]),meandistanceCN[[lastCN]], 
     # xlim = c(0,120),ylim = c(0,12000), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
points(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]],col = "red") 
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meandistanceCN[[lastCN]]), y=mean(meandistanceCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meandistanceCN[[j]]),meandistanceCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]]*648,col = colrestCN[i])}
  text(x= length(meandistanceCN[[j]]), y=mean(meandistanceCN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceCN)[j])
}




###########################################################################
# HC Graph departing from PCStable graph 
###########################################################################
data <- TimeCoordsAnom_from_Grid_std(tas_ncep_10d)
dir <- cextend(pc_1_eBIC_g5)
class(dir)
combi <- hc(x = as.data.frame(data), start = dir)
############################################################################
# WEIGHTED GRAPHS
############################################################################
 gridGraphwCN1 <- lapply(gridGraphsCN[seq(60,100,1)], set_distances)
 gridGraphwCN2 <- lapply(gridGraphsCN[seq(34,59,1)], set_distances)
 gridGraphwCN3 <- lapply(gridGraphsCN[seq(25,33,1)], set_distances)
 gridGraphwCN4 <- lapply(gridGraphsCN[seq(2,24,1)],set_distances)
 gridGraphwCN5 <- set_distances(gridGraphsCN[[1]])
# gridGraphwBN <- lapply(gridGraphsBN, set_distances)
 gridGraphwRenyi <- lapply(gridGraphsRenyi, set_distances)
 gridGraphwLattices <- lapply(gridGraphsLattices[seq(1,4,1)], set_distances)
 save(gridGraphwCN1,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN1.rda")
 save(gridGraphwCN2,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN2.rda")
 save(gridGraphwCN3,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN3.rda")
 save(gridGraphwCN4,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN4.rda")
# save(gridGraphwBN,file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwBN.rda")
# save(gridGraphwRenyi,file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwRenyi.rda")
 save(gridGraphwRenyi,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwRenyi.rda")
 save(gridGraphwLattices,file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwLattices.rda")



load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN2.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwBN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwRenyi.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwLattices.rda")

names <- c("CN","BN",'Renyi',"Lattices")
gridGraphlists <- list(gridGraphwCN2,gridGraphwBN,gridGraphwRenyi,gridGraphwLattices)


# with weights
for(i in 1:length(names)){
  name <- names[i]
  gridGraphlist <- gridGraphlists[[i]]
  globallist <- lapply(gridGraphlist, graph2measure.global)
  assign(paste0("globmeas",name),globallist)
  assign(paste0("globalasslist",name),lapply(globallist, function(x){x$globalass}))
  assign(paste0("globalcluslist",name),lapply(globallist, function(x){x$globalclus}))
  assign(paste0("globaledenslist",name), lapply(globallist, function(x){x$edens}))
  assign(paste0("globaldiamLCClist",name), lapply(globallist,function(x){x$diameterLCC}))
  assign(paste0("globaldiamGlist",name), lapply(globallist,function(x){x$diameterG}))
  globallist <- NULL
  
  degreelist <- lapply(gridGraphlist,graph2measure.degree)
  assign(paste0("degreemeas",name),degreelist)
  assign(paste0("degreedist",name), lapply(degreelist, function(x){x$degrdist}))
  degreelist <- NULL
}
