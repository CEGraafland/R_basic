##############################################################################
# local measures on Bayesian networks and correlation networks
##############################################################################
rm(list = ls())
library(gridExtra)
library(gridGraphics)
library(grid)
library("bnlearn")
library(raster)
################################################
# For MAC:
# source("/Users/lisettegraafland/Desktop/mnt/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
############################################################################################
# Correlation networks: 
# Complex measures with lapply and measure2clim
# used in Resumen 2
############################################################################################
# taus that belong to gridgraphs:
taus <- seq(from = 0, to = 1, length.out = 101)
seqCN2 <- 34:59
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
gridGraphsCN <- gridGraphs[seqCN2]
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBN.rda")
# weighteds.
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN2.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwBN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwRenyi.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwLattices.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBNtoCNperm1.rda")


names <- c("CN","BN", 
           "BNtoCN", 
           "CN2w", "BNw")

gridGraphlists <- list(gridGraphsCN,gridGraphsBN,
                       gridGraphsBNtoCN, 
                       gridGraphwCN2,gridGraphwBN)

for(i in 1:length(names)){
  name <- names[i]
  gridGraphlist <- gridGraphlists[[i]]
  locallist <- lapply(gridGraphlist, graph2measure.local)
  
  assign(paste0("degreemeas",name),locallist)

  assign(paste0("betweennesses",name),lapply(locallist, measure2clim, what = "betweenness", ref.grid = tas_ncep_10d))
  assign(paste0("awconnectivities",name), lapply(locallist, measure2clim, what = "awconnectivity", ref.grid = tas_ncep_10d))
  assign(paste0("degrees",name),lapply(locallist, measure2clim, what = "degree", ref.grid = tas_ncep_10d))
  assign(paste0("strengths",name),lapply(locallist, measure2clim, what = "strength", ref.grid = tas_ncep_10d))
  assign(paste0("closenesses",name), lapply(locallist, measure2clim, what = "closeness", ref.grid = tas_ncep_10d))
  assign(paste0("localclusterings",name), lapply(locallist, measure2clim, what = "localclustering", ref.grid = tas_ncep_10d))
  assign(paste0("barratclusterings",name), lapply(locallist, measure2clim, what = "barratclustering", ref.grid = tas_ncep_10d))
  assign(paste0("barratweightclusterings",name), lapply(locallist, measure2clim, what = "barratweightclustering", ref.grid = tas_ncep_10d))
  locallist <- NULL
}

####################################################################
# set names
####################################################################
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(GraphsCN) <- as.character(numberofedgesCN)
edge_attr(gridGraphsCN[[1]]$graph)

GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})
edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- sapply(edgelistsBN, length)
names(GraphsBN) <- as.character(numberofedgesBN)

GraphsBNtoCN  <- lapply(gridGraphsBNtoCN, function(x){x$graph})
edgelistsBNtoCN <- lapply(GraphsBNtoCN, E)
numberofedgesBNtoCN <- lapply(edgelistsBNtoCN, length)

GraphsCN2w <- lapply(gridGraphwCN2, function(x){x$graph}) 
edgelistsCN2w <- lapply(GraphsCN2w, E)
numberofedgesCN2w <- sapply(edgelistsCN2w, length)
names(GraphsCN2w) <- as.character(numberofedgesCN2w)
edge_attr(gridGraphwCN2[[1]]$graph)

GraphsBNw <- lapply(gridGraphwBN, function(x){x$graph})
edgelistsBNw <- lapply(GraphsBNw, E)
numberofedgesBNw <- sapply(edgelistsBNw, length)
names(GraphsBNw) <- as.character(numberofedgesBNw)
# edgelistsRenyi <- lapply(GraphsRenyi, E)
# numberofedgesRenyi <- lapply(edgelistsRenyi, length)
# names(GraphsRenyi) <- as.character(numberofedgesRenyi)
# 
# edgelistsLattices <- lapply(latticesgraphs, E)
# numberofedgesLattices <- lapply(edgelistsLattices, length)
# names(latticesgraphs) <- as.character(numberofedgesLattices)


length(whatCN)
length(whatBN)
whatBN <- betweennessesBN
multisBN = list()
for (i in 1:length(whatBN)){
  plot <- plotClimatology(whatBN[[i]], main = list(paste0("BN: ",numberofedgesBN[i]), cex = 0.5),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisBN[[i]] <- plot
}

whatCN <- barratweightclusteringsCN
multisCN= list()
for (i in 1:length(whatCN)){
  plot <- plotClimatology(whatCN[[i]], main = list(paste0("CN: ",numberofedgesCN[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisCN[[i]] <- plot
}

whatBNtoCN <- barratweightclusteringsBNtoCN
multisBNtoCN= list()
for (i in 1:length(whatBNtoCN)){
  plot <- plotClimatology(whatBNtoCN[[i]], main = list(paste0("BNtoCN: ",numberofedgesBNtoCN[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisBNtoCN[[i]] <- plot
}


whatCN2w <- barratweightclusteringsCN2w
multisCN2w= list()
for (i in 1:length(whatCN2w)){
  plot <- plotClimatology(whatCN2w[[i]], main = list(paste0("CN2w: ",numberofedgesCN2w[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisCN2w[[i]] <- plot
}

whatBNw <- barratweightclusteringsBNw
multisBNw = list()
for (i in 1:length(whatBNw)){
  plot <- plotClimatology(whatBNw[[i]], main = list(paste0("BNw: ",numberofedgesBNw[i]), cex = 0.5),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisBNw[[i]] <- plot
}


selectCN <- multisCN[seq(1,length(multisCN),5)]
n <- length(selectCN)
nCol <- floor(sqrt(n))
measurename <- attr(whatCN[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectCN, ncol=nCol,top = paste0(measurename)))
dev.off()

selectBN <- multisBN[seq(6,length(multisBN),length.out = 10)]
n <- length(selectBN)
nCol <- floor(sqrt(n))
measurename <- attr(whatBN[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectBN, ncol=nCol,top = paste0(measurename)))
dev.off()

selectBNtoCN <- multisBNtoCN[1:10]
n <- length(selectBNtoCN)
nCol <- floor(sqrt(n))
measurename <- attr(whatCN[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectBNtoCN, ncol=nCol,top = paste0(measurename)))
dev.off()

selectCN2w <- multisCN2w[seq(1,length(multisCN2w),5)]
n <- length(selectCN2w)
nCol <- floor(sqrt(n))
measurename <- attr(whatCN2w[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectCN2w, ncol=nCol,top = paste0(measurename)))
dev.off()


selectBNw <- multisBNw[seq(6,length(multisBNw),length.out = 10)]
n <- length(selectBNw)
nCol <- floor(sqrt(n))
measurename <- attr(whatBNw[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectBNw, ncol=nCol,top = paste0(measurename)))
dev.off()
####################################################################################
#Change plotClimatology in spatialPlot
####################################################################################
length(whatCN)
length(whatBN)
center <- 180
whatBN <- localclusteringsBN
multisBN = list()
for (i in 1:length(whatBN)){
  plot <- spatialPlot(whatBN[[i]], backdrop.theme = "coastline", 
                      lonCenter = center, 
                      main = list(paste0("BN: ",numberofedgesBN[i]), cex = 0.5),
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisBN[[i]] <- plot
}

whatCN <- localclusteringsCN
multisCN= list()
for (i in 1:length(whatCN)){
  plot <- spatialPlot(whatCN[[i]], backdrop.theme = "coastline",
                      lonCenter = center,
                      main = list(paste0("CN: ",numberofedgesCN[i]), cex = 0.5, pos = 0.25),
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisCN[[i]] <- plot
}

whatBNtoCN <- localclusteringsBNtoCN
multisBNtoCN= list()
for (i in 1:length(whatBNtoCN)){
  plot <- spatialPlot(whatBNtoCN[[i]], backdrop.theme = "coastline",
                      lonCenter = center,
                      main = list(paste0("BNtoCN: ",numberofedgesBNtoCN[i]), cex = 0.5, pos = 0.25),
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisBNtoCN[[i]] <- plot
}


whatCN2w <- localclusteringsCN2w
multisCN2w= list()
for (i in 1:length(whatCN2w)){
  plot <- spatialPlot(whatCN2w[[i]], backdrop.theme = "coastline",
                      lonCenter = center,
                      main = list(paste0("CN2w: ",numberofedgesCN2w[i]), cex = 0.5, pos = 0.25),
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisCN2w[[i]] <- plot
}

whatBNw <- betweennessesBNw
multisBNw = list()
for (i in 1:length(whatBNw)){
  plot <- spatialPlot(whatBNw[[i]], backdrop.theme = "coastline", 
                      lonCenter = center,
                      main = list(paste0("BNw: ",numberofedgesBNw[i]), cex = 0.5),
                      # region = TRUE,
                      # col.regions = colRedsamp,
                      rev.colors = TRUE,
                      set.max = 20,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisBNw[[i]] <- plot
}


selectCN <- multisCN[seq(1,length(multisCN),5)]
n <- length(selectCN)
nCol <- floor(sqrt(n))
measurename <- attr(whatCN[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectCN, ncol=nCol,top = paste0(measurename)))
dev.off()

selectBN <- multisBN[seq(6,length(multisBN),length.out = 10)]
n <- length(selectBN)
nCol <- floor(sqrt(n))
measurename <- attr(whatBN[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectBN, ncol=nCol,top = paste0(measurename)))
dev.off()

selectBNtoCN <- multisBNtoCN[1:10]
n <- length(selectBNtoCN)
nCol <- floor(sqrt(n))
measurename <- attr(whatCN[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectBNtoCN, ncol=nCol,top = paste0(measurename)))
dev.off()

selectCN2w <- multisCN2w[seq(1,length(multisCN2w),5)]
n <- length(selectCN2w)
nCol <- floor(sqrt(n))
measurename <- attr(whatCN2w[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectCN2w, ncol=nCol,top = paste0(measurename)))
dev.off()


selectBNw <- multisBNw[seq(6,length(multisBNw),length.out = 10)]
n <- length(selectBNw)
nCol <- floor(sqrt(n))
measurename <- attr(whatBNw[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(selectBNw, ncol=nCol,top = paste0(measurename)))
dev.off()
####################################################################
# Hier verder gaan.
####################################################################
list  <- append(multisBN, multisCN)


measurename <- attr(whatCN[[1]]$Data,"climatology:fun")
measurename

l = 2:6
m = length(list)/2 + l
n<- c(l,m)

do.call(grid.arrange, c(list[n], nrow = length(multisCN[n])/2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

l = 7:10
m = length(list)/2 + l
n<- c(l,m)

do.call(grid.arrange, c(list[n], nrow = length(multis[n])/2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))






##############################################################
# Dit is alleen maar gekopieerd en geplakt vanuit localMeasure. NIET NODIG
##########################################################




#############################################################
# Now, get information global scalar measures. 
# (Possible to insert this in graph2measure later)
#############################################################
grid <- tas_ncep_10d

numberofedgesCN <- c()
transCN <- c()
assdegreeCN <- c()

for (i in 1:length(tau)){
  graphObj <- graph_from_Grid(grid, th = tau[i])
  numberofedgesCN[i] <- length(E(graphObj$graph))
  transCN[i] <- transitivity(graphObj$graph, type = "global")
  assdegreeCN[i] <- assortativity_degree(graphObj$graph, directed = FALSE)
}

plot(numberofedgesCN,transCN, xlab = "|E|", ylab = "C", main = "Transitivity CN", col = "black")
text(numberofedgesCN, transCN, labels = c(numberofedgesCN), cex= 0.7, pos = 3)
plot(numberofedgesCN,assdegreeCN, xlab = "|E|", ylab = "r", main = "Degree assortativity CN", col = "black")
text(numberofedgesCN, assdegreeCN, labels = c(numberofedgesCN), cex= 0.7, pos = 3)
##########################################################################################
# Bayesian network: measures working exmample *
##########################################################################################
# TCA_from_Grid only produces time.coords.matrix (data_coords) and VertexCoords

TimeCoordsAnom_from_Grid_aslist <- function(grid) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% localScaling() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <- do.call("bindGrid.time", seas.list) %>% redim(drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  
  out <- list("data_coords" = time.coords.matrix,
              "VertexCoords" = ref.coords)
  
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}

exampleB <- hc_edges_loglik_10d_1200_1400i$networks
igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
igraphske<- as.undirected(igraphdir)
# create graph_from_Grid object with TimeCoordsAnom_from_Grid 
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
#check
attr(graphObject, "VertexCoords")
E(graphObject$graph)
9*36+5
10*36+14
323 + 36
class(graphObject$graph)
fitted <- bn.fit(exampleB, as.data.frame(graphObject$data_coords))
mutilate <- mutilated(fitted, evidence = list(V356 = 1))
# cpquery(fitted, event = (V374 >= 0.5), evidence = (V375 >= 0.5))
# cpquery(fitted, event = (V356 >= 0.5), evidence = (V375 >= 0.8))
# cpquery(fitted, event = (V356 >= 0.5), evidence = ((V375 >= 0.8)&(V375 <= 0.8)))
# cpquery(fitted, event = (V356 >= 0.5), evidence = (V375 <= 0.8))

mutilate$V374
fitted$V374

# apply measure2clim and plotClimatology
# on graphObject
measuresBay <- graph2measure(graphObject)
clim.awcon.Bay10d <- measure2clim(measuresBay, what = "awconnectivity", ref.grid = tas_ncep_10d)
clim.betw.Bay10d<- measure2clim(measuresBay, what = "betweenness", ref.grid = tas_ncep_10d)
clim.close.Bay10d <- measure2clim(measuresBay, what = "closeness", ref.grid = tas_ncep_10d)
clima.lclus.Bay10d <- measure2clim(measuresBay, what = "localclustering", ref.grid = tas_ncep_10d)

plotClimatology(clim.betw.Bay10d, backdrop.theme = "coastline")
plotClimatology(clim.close.Bay10d, backdrop.theme = "coastline")
plotClimatology(clim.awcon.Bay10d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.Bay10d,backdrop.theme = "coastline")

#######################################################################################
# Working example measures complex networks for list of bayesian networks * 
# Measures are applicated on non rescaled data (only anomaly)
# (propagation is applicated on rescaled data)
# TMS_as_list is therefore adequate: Graph from non rescaled data in list with non rescaled
# data from function.
#######################################################################################
# Load data in PaperScript.R
tabu_list <- list(tabu_1_eBIC_g0, tabu_1_eBIC_g0.1, tabu_1_eBIC_g0.25, tabu_1_eBIC_g0.5, tabu_1_eBIC_g0.75, tabu_1_eBIC_g1)
networks <- lapply(tabu_list, function(m) m[["networks"]][[1]])
nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")

hc_list <- list(hc_edges_loglik_10d_200_400i,
                hc_edges_loglik_10d_600_800i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1000_1200i,
                hc_edges_loglik_10d_1200_1400i,
                hc_edges_loglik_10d_1400_1600i,
                hc_edges_loglik_10d_1600_1800i,
                hc_edges_loglik_10d_1800_2000i,
                hc_edges_loglik_10d_2400_2600i,
                hc_edges_loglik_10d_3000_3200i,
                hc_edges_loglik_10d_4800_5000i,
                hc_edges_loglik_10d_8600_8695i)
# select only the networks to calculate nedges
networks <- lapply(hc_list, function(m) m[["networks"]])
nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]

# Convert the graphs to undirected igraphs
igraphsdir <- lapply(networks, as.graphNEL) 
igraphsdir <- lapply(igraphsdir, igraph.from.graphNEL)
igraphsske<- lapply(igraphsdir, as.undirected)
# Create the graphObject as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(networks))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- igraphsske[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsske[[i]])
}
# apply measure2clim and plotClimatology
# on list with graphObjects
measureslistBay <- lapply(graphObjects,graph2measure)
degrees <- lapply(measureslistBay, measure2clim,what = "degree", ref.grid = tas_ncep_10d)
closenesses <- lapply(measureslistBay, measure2clim,what = "closeness", ref.grid = tas_ncep_10d)
betweenness <- lapply(measureslistBay, measure2clim,what = "betweenness", ref.grid = tas_ncep_10d)
awconnectivities <- lapply(measureslistBay, measure2clim,what = "awconnectivity", ref.grid = tas_ncep_10d)
localclusterings <- lapply(measureslistBay, measure2clim,what = "localclustering", ref.grid = tas_ncep_10d)
# Choose measure to plot 
what <- awconnectivities
multis = list()
for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.6),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5))
                          # ,
                          #  at = seq(0,0.035,0.005)
                          #   ,
                          #    set.max = 0.03
  )
  multis[[i]] <- plot
}

n <- length(multis)
nCol <- floor(sqrt(n))
measurename <- attr(what[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(multis, ncol=nCol, top = paste0(measurename)))
dev.off()

#############################################################
# Now, get information global scalar measures. 
# (Possible to insert this in graph2measure later)
#############################################################
coefficientsBN <- c()
numberofedgesBN <- c()
transBN <- c()
assdegreeBN <- c()

for (i in 1:length(networks)){
  exampleB <- networks[[i]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedgesBN[i] <- length(E(igraphske))
  transBN[i] <- transitivity(igraphske, type = "global")
  assdegreeBN[i] <- assortativity_degree(igraphske, directed = FALSE)
}


plot(numberofedgesBN[4:length(numberofedgesBN)],coefficientsBN[4:length(numberofedgesBN)])
plot(numberofedgesBN,coefficientsBN, col = "blue", xlab = "|E|", ylab = "slope", main = "coefficients BN")
text(numberofedgesBN, coefficientsBN, labels = c(numberofedgesBN), cex= 0.7, pos = 3)
plot(numberofedgesBN,transBN, xlab = "|E|", ylab = "C", main = "Transitivity BN", col = "blue")
text(numberofedgesBN, transBN, labels = c(numberofedgesBN), cex= 0.7, pos = 3)
plot(numberofedgesBN,assdegreeBN, xlab = "|E|", ylab = "r", main = "Degree assortativity BN", col = "blue")
text(numberofedgesBN, assdegreeBN, labels = c(numberofedgesBN), cex= 0.7, pos = 3)
#######################################################################################
# Working example measures complex networks for list of precision matrix networks
#######################################################################################
# Choose precision matrix 
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
tau_list <- c(0.1,0.06,0.045,0.026)
precmatrix <- glasso01$wi # beta = 0.01

graph_from_precisionmatrix <- function(precmatrix, tau){
  # convert precision matrix to partial variances matrix
  f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
  Vecf <- Vectorize(f,vectorize.args = c('r','c'))
  PartVar <- outer(1:nrow(precmatrix),1:ncol(precmatrix),Vecf,precmatrix)
  # create adj.matrix 
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[adj.matrix <= tau ] <- 0
  adj.matrix[abs(adj.matrix) > tau ] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  nedgesprecgraph <- E(precgraph)
  return(precgraph)
}

# select only the networks to calculate nedges
networks <- lapply(X = tau_list, FUN = graph_from_precisionmatrix, precmatrix = glasso01$wi)
networks
edgesnetworks <- sapply(networks, E)
nedgesnetworksPN <- sapply(edgesnetworks,length)
nedgesnetworksPN
firstel <- nedgesnetworksPN[1]
lastel <- nedgesnetworksPN[length(nedgesnetworksPN)]

# Create the graphObject as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(networks))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- networks[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(networks[[i]])
}
# apply measure2clim and plotClimatology
# on list with graphObjects
measureslistPrec <- lapply(graphObjects,graph2measure)
degreesPrec <- lapply(measureslistPrec, measure2clim,what = "degree", ref.grid = tas_ncep_10d)
closenessesPrec <- lapply(measureslistPrec, measure2clim,what = "closeness", ref.grid = tas_ncep_10d)
betweennessPrec <- lapply(measureslistPrec, measure2clim,what = "betweenness", ref.grid = tas_ncep_10d)
awconnectivitiesPrec <- lapply(measureslistPrec, measure2clim,what = "awconnectivity", ref.grid = tas_ncep_10d)
localclusteringsPrec <- lapply(measureslistPrec, measure2clim,what = "localclustering", ref.grid = tas_ncep_10d)
# Choose measure to plot 
what <- awconnectivitiesPrec
multis = list()
for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("PN: ",nedgesnetworksPN[i]), cex = 0.5),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.25)))
  multis[[i]] <- plot
}

n <- length(multis)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(multis, ncol=nCol))
dev.off()


###############################################################################################
# combining measures bayesian and complex networks and precision networks
###############################################################################################
what <- awconnectivities
multis = list()
for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.5),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multis[[i]] <- plot
}

whatCM <- awconnectivitiesCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- plotClimatology(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisCM[[i]] <- plot
}

whatPN <- closenessesPrec
multisPN = list()
for (i in 1:length(whatPN)){
  plot <- plotClimatology(whatPN[[i]], main = list(paste0("PN: ",nedgesnetworksPN[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multisPN[[i]] <- plot
}

list  <- append(multis, multisCM)
listBN_CM_PN <- append(list, multisPN)

measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename

l = 2:6
m = length(list)/2 + l
n<- c(l,m)

do.call(grid.arrange, c(list[n], nrow = length(multis[n])/2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

l = 7:10
m = length(list)/2 + l
n<- c(l,m)

do.call(grid.arrange, c(list[n], nrow = length(multis[n])/2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))


do.call(grid.arrange, c(listBN_CM_PN, nrow = length(multis), ncol = 3, as.table = FALSE, top = paste0(measurename) ))
do.call(marrangeGrob, c(list, nrow = length(multis), ncol = 2))
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/comparemeasures/",measurename,firstel,"_",lastel,"_",length(list)/2,".pdf")
plotname
#pdf(plotname, height = 10, width = 4) # para 8
pdf(plotname)
do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
dev.off()

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/comparemeasures/",measurename,firstel,"_",lastel,"_",length(list)/2,"_BNCNPN.pdf")
plotname
#pdf(plotname, height = 10, width = 4) # para 8
pdf(plotname)
do.call(grid.arrange, c(listBN_CM_PN, nrow = length(multis), ncol = 3, as.table = FALSE, top = paste0(measurename) ))
dev.off()


###############################################################################################
# combining measures: BUT WITH equal colorkey: bayesian and complex networks and precision networks
# Localclusterings
###############################################################################################
what <- localclusterings
multis = list()
seq <- 2:2
a <- 1
at <- seq(0,0.25,0.025)

for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5))
                          # ,
                          # at = at,
                          # set.max = a
  )
  multis[[i]] <- plot
}

whatCM <- localclusteringsCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- plotClimatology(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5))
                          # ,
                          # at = at,
                          # set.max = a
  )
  
  multisCM[[i]] <- plot
}

whatPN <- localclusteringsPrec
multisPN = list()
for (i in length()){
  plot <- plotClimatology(whatPN[[i]], main = list(paste0("PN: ",nedgesnetworksPN[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)),
                          at = at,
                          set.max = a)
  multisPN[[i]] <- plot
}


listBN_CM <- list(multis[[2]],multis[[7]], multisCM[[8]], multisCM[[10]])
measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename


#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
##################################################################
# Closenesses
##################################################################
what <- closenesses
multis = list()
seq <- 2:2
at1 <- seq(0,2.5*10^-5,0.25*10^-5)
at2 <- seq(0,0.00065,0.00005)
at
a <- 0.00039


for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)),
                          #at = at2,
                          set.max = a)
  multis[[i]] <- plot
}

whatCM <- closenessesCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- plotClimatology(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  #at = at2,
  #set.max = a)
  
  multisCM[[i]] <- plot
}
listBN_CM <- list(multis[[2]],multis[[3]], multis[[6]], multisCM[[7]],multisCM[[9]], multisCM[[10]])
measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename


#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 3, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

##################################################################
# Betweenness
##################################################################
what <- betweenness
multis = list()
seq <- 2:2
at1 <- seq(0,12,12/16)
at2 <- seq(0,0.00065,0.00005)
at
a <- 12


for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5))
                          # ,
                          # at = at1
                          # ,
                          # set.max = a
  )
  multis[[i]] <- plot
}

whatCM <- betweennessCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- plotClimatology(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5))
                          # ,
                          # at = at1
                          # ,
                          # set.max = a
  )
  
  multisCM[[i]] <- plot
}
listBN_CM <- list(multis[[2]],multis[[5]], multisCM[[7]], multisCM[[10]])
measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename


#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

##################################################################
# area weighted connectivity
##################################################################
what <- awconnectivities
multis = list()
seq <- 2:2
at1 <- seq(0,0.018,0.018/16)
at2 <- seq(0,0.00065,0.00005)
at
a <- 1


for (i in 1:length(what)){
  plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.6))
                          # ,
                          # at = at1
                          #set.max = a)
  )
  multis[[i]] <- plot
}

whatCM <- awconnectivitiesCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- plotClimatology(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.6))
                          # ,
                          # at = at1
                          #set.max = a)
  )
  
  multisCM[[i]] <- plot
}
listBN_CM <- list(multis[[4]],multis[[5]], multisCM[[3]], multisCM[[10]])
listBN_CM <- list(multis[[4]],multisCM[[9]])
measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename


#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE))