##############################################################################
# Complex measures on Bayesian networks and correlation networks
##############################################################################
setwd("~/data/Untitled/Trabajo/R_practice/R/")
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")
rm(list = ls())
library(gridExtra)
library(gridGraphics)
library(grid)
library("bnlearn")
library("visualizeR")
################################################
# For MAC:
# source("/Users/lisettegraafland/Desktop/mnt/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
################################################
source("Functions/BasicNetworkFunctions.R")
source("Functions/HillclimbingFunctions2.R")
source("Functions/CN_ConstructionandMeasuresFunctions.R")
load("../Data/tas_ncep_10d.rda")
# MAC:
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
##############################################################################
# Complex measures con correlation model
# and powerlaw on correlation model
# Used in Resumen 1. 
# Saved in catharina/Documents/PlotsResumen
##############################################################################
#Plaatjes Complex Networks 10d
graph_from_Grid
i <- 1
for (i in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)){
tau = i
graph10d <- graph_from_Grid(tas_ncep_10d, th = tau, subind = NULL, method = "pearson")
measures10d <- graph2measure(graph10d)

numberofedges <- length(E(graph10d$graph))
edgedensity <- round(measures10d$edens, digits = 4)

clim.awcon.10d <- measure2clim(measures10d, what = "awconnectivity", ref.grid = tas_ncep_10d)
clim.betwe.10d <- measure2clim(measures10d, what = "betweenness", ref.grid = tas_ncep_10d)
clim.close.10d <- measure2clim(measures10d, what = "closeness", ref.grid = tas_ncep_10d)
clim.lclus.10d <- measure2clim(measures10d, what = "localclustering", ref.grid = tas_ncep_10d)

assign(paste0("awctau",tau),clim.awcon.10d)
assign(paste0("betwtau",tau),clim.betwe.10d)
assign(paste0("closetau",tau),clim.close.10d)
assign(paste0("locclustau",tau),clim.lclus.10d)

c <- spatialPlot(clim.awcon.10d, backdrop.theme = "coastline", main = paste0("tau = ",tau," edges = ", numberofedges,"\nedge density = ", edgedensity))
d <- spatialPlot(clim.betwe.10d, backdrop.theme = "coastline", main = paste0("tau = ",tau," edges = ", numberofedges,"\nedge density = ", edgedensity))#, at = seq(0,6,0.5))
e <- spatialPlot(clim.close.10d, backdrop.theme = "coastline", main = paste0("tau = ",tau," edges = ", numberofedges,"\nedge density = ", edgedensity))#, at = seq(2.38e-6,2.53e-6,0.01e-6))
f <- spatialPlot(clim.lclus.10d, backdrop.theme = "coastline", main = paste0("tau = ",tau," edges = ", numberofedges,"\nedge density = ", edgedensity))

assign(paste0("plotawctau",tau),c)
assign(paste0("plotbetwtau",tau),d)
assign(paste0("plotclosetau",tau),e)
assign(paste0("plotlocclustau",tau),f)

}

awcs = list(awctau0.1,awctau0.2,awctau0.3,awctau0.4,awctau0.5,awctau0.6,awctau0.7,awctau0.8)
betws = list(betwtau0.8,betwtau0.7,betwtau0.6,betwtau0.5,betwtau0.4,betwtau0.3,betwtau0.2,betwtau0.1)
closes = list(closetau0.8,closetau0.7,closetau0.6,closetau0.5,closetau0.4,closetau0.3,closetau0.2,closetau0.1)
loccluss = list(locclustau0.8,locclustau0.7,locclustau0.6,locclustau0.5,locclustau0.4,locclustau0.3,locclustau0.2,locclustau0.1)


textawc <- textGrob("Area weighted connectivity")
textbetw <- textGrob("Betweenness Centrality")
textclose <- textGrob("Closeness centrality")
textlocclus <- textGrob("Local Clustering coefficient")

grid.arrange(plotawctau0.8,plotawctau0.7,plotawctau0.6,plotawctau0.5,plotawctau0.4,plotawctau0.3,plotawctau0.2,plotawctau0.1,textawc,nrow = 5,ncol = 2,heights = c(1,1,1,1,0.2))
grid.arrange(plotbetwtau0.8,plotbetwtau0.7,plotbetwtau0.6,plotbetwtau0.5,plotbetwtau0.4,plotbetwtau0.3,plotbetwtau0.2,plotbetwtau0.1,textbetw,nrow = 5,ncol = 2,heights = c(1,1,1,1,0.2))
grid.arrange(plotclosetau0.8,plotclosetau0.7,plotclosetau0.6,plotclosetau0.5,plotclosetau0.4,plotclosetau0.3,plotclosetau0.2,plotclosetau0.1,textclose,nrow = 5,ncol = 2,heights = c(1,1,1,1,0.2))
grid.arrange(plotlocclustau0.8,plotlocclustau0.7,plotlocclustau0.6,plotlocclustau0.5,plotlocclustau0.4,plotlocclustau0.3,plotlocclustau0.2,plotlocclustau0.1,textlocclus,nrow = 5,ncol = 2,heights = c(1,1,1,1,0.2))


#visualize degree distribution (ddist) with respect to the tresholds
plotname <- "/home/catharina/Documents/PlotsResumen/Powerlaw10dtau01_08.pdf"
pdf(plotname, paper = 'US')
plplotlist(tas_ncep_10d, tau = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))
dev.off()

############################################################################################
# Correlation networks: 
# Complex measures with lapply and measure2clim
# used in Resumen 2
############################################################################################
# 1 3 6 7 +4
tau <- c(0.87,
         0.8,
         0.75,
         0.7,
         0.66,
         0.62,
         0.52,
         0.49,
         0.41,
         0.345,
         0.32,
         0.3,
         0.2,
         0.1,
         0.05,
         0.0)
tau <- c(0.60,0.55,0.34,0.31)
tau <- c(0.5,0.42,0.35,0.31)
graphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = NULL, method = 'pearson')
measures10d <- lapply(graphs10d, graph2measure)

graphssolo  <- lapply(graphs10d, function(m) m$graph)
graphssolo
nedgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(nedgesnetworksCM, length)
nedgesnetworksCM

degreesCM <- lapply(measures10d, measure2clim, what = "degree", ref.grid = tas_ncep_10d)
closenessesCM <- lapply(measures10d, measure2clim, what = "closeness", ref.grid = tas_ncep_10d)
betweennessCM <- lapply(measures10d, measure2clim, what = "betweenness", ref.grid = tas_ncep_10d)
awconnectivitiesCM <- lapply(measures10d, measure2clim, what = "awconnectivity", ref.grid = tas_ncep_10d)
localclusteringsCM <- lapply(measures10d, measure2clim, what = "localclustering", ref.grid = tas_ncep_10d)

whatCM <- degreesCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.6, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5))
                          #at = seq(0,0.03,0.005),
                          #set.max = 0.03
                          )
  multisCM[[i]] <- plot
}

n <- length(multisCM)
nCol <- floor(sqrt(n))
measurename <- attr(whatCM[[1]]$Data,"climatology:fun")
do.call("grid.arrange", c(multisCM, ncol=nCol,top = paste0(measurename)))
dev.off()
#############################################################
# Now, get information global scalar measures. 
# (Possible to insert this in graph2measure later)
#############################################################
grid <- tas_ncep_10d

numberofedgesCN <- c()
transCN <- c()
assdegreeCN <- c()

for (i in 1:length(tau)){
  graphObj <- graph_from_Grid(grid, th = tau[i], method = 'pearson')
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

# TimeCoordsAnom_from_Grid_aslist <- function(grid) {
#   seas <- getSeason(grid)
#   coords <- getCoordinates(grid)
#   x <- coords$x
#   y <- coords$y
#   ref.coords <- expand.grid(y, x)[2:1]
#   names(ref.coords) <- c("x", "y")
#   ref.dates <- getRefDates(grid)
#   seas.list <- lapply(1:length(seas), function(i) {
#     subsetGrid(grid, season = seas[i]) %>% localScaling() %>% redim(drop = TRUE)
#   })
#   grid <- NULL
#   aux <- do.call("bindGrid.time", seas.list) %>% redim(drop = TRUE)
#   seas.list <- NULL
#   time.coords.matrix <- array3Dto2Dmat(aux$Data)
#   
#   out <- list("data_coords" = time.coords.matrix,
#               "VertexCoords" = ref.coords)
#   
#   
#   attr(out, "Xcoords") <- x
#   attr(out, "Ycoords") <- y
#   attr(out, "ref.dates") <- ref.dates
#   attr(out, "VertexCoords") <- ref.coords
#   return(out)
# }

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

spatialPlot(clim.betw.Bay10d, backdrop.theme = "coastline")
spatialPlot(clim.close.Bay10d, backdrop.theme = "coastline")
spatialPlot(clim.awcon.Bay10d,backdrop.theme = "coastline")
spatialPlot(clima.lclus.Bay10d,backdrop.theme = "coastline")

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

load("../Data/networks_hcnetworks10d.rda")

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

what <- closenesses
multis = list()
for (i in 1:length(what)){
  plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.6),
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
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
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
what <- localclusteringsPrec
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
  plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.5),
                          backdrop.theme = "coastline", 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  multis[[i]] <- plot
}

whatCM <- awconnectivitiesCM
multisCM = list()
for (i in 1:length(whatCM)){
  plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.5, pos = 0.25),
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
  plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
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
  plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
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
  plot <- spatialPlot(whatPN[[i]], main = list(paste0("PN: ",nedgesnetworksPN[i]), cex = 0.5, pos = 0.25),
                          backdrop.theme = "coastline",
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)),
                          at = at,
                          set.max = a)
  multisPN[[i]] <- plot
}

#listBN_CM <- list(multis[[2]],multis[[3]], multisCM[[2]], multisCM[[3]])
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
at2 <- seq(0.00020,0.00050,0.000033)
at
a <- 0.00050

display.brewer.all()
for (i in 1:length(what)){
  #plot <- spatialPlot(grid = what[[i]], backdrop.theme = "coastline", set.min = NULL,
   #                   set.max = a, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
    #                  color.theme = "Purples",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                          backdrop.theme = "coastline", rev.colors = TRUE,
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
                          # ,
                          #   at = at2,
                          #  set.max = a)
  multis[[i]] <- plot
}

whatCM <- closenessesCM
multisCM = list()
for (i in 1:length(whatCM)){
 # plot <- spatialPlot(grid = whatCM[[i]], backdrop.theme = "coastline", set.min = NULL,
  #                    set.max = 0.016, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8),
   #                   color.theme = "Purples",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
                          backdrop.theme = "coastline", rev.colors = TRUE,
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
                          # ,
                          # at = at2,
                          # set.max = a)
  
  multisCM[[i]] <- plot
}
listBN_CM <- list(multis[[2]],multis[[3]], multis[[6]], multisCM[[7]],multisCM[[9]], multisCM[[10]])

listBN_CM <- list(multis[[2]],multis[[3]], multis[[6]], multisCM[[1]],multisCM[[2]], multisCM[[3]])
measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename


#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 3, ncol = 2, as.table = FALSE, top = paste0(measurename) ))


cl.bn <-6
movawc1<-movingmedias(measureslistBay[[cl.bn]]$closeness)
movawc1matrixmean <- apply(movawc1, MARGIN = 3, FUN = mean)

memClim <- quantity2clim(movawc1matrixmean, "awc mov", tas_ncep_10d)
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average BN: ",nedgesnetworks[cl.bn]), cex = 0.8), rev.colors = TRUE,colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")

plot1

cl.cn <- 3
movawc1.cm<-movingmedias(measures10d[[cl.cn]]$closeness)
movawc1matrixmean.cm <- apply(movawc1.cm, MARGIN = 3, FUN = mean)

memClim.cm <- quantity2clim(movawc1matrixmean.cm, "awc mov", tas_ncep_10d)
plot2 <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[cl.cn]), cex = 0.8), rev.colors = TRUE,colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot2

plotname <- paste0("../plots/comparemeasures_thesis/",measurename,"_raw.pdf")
plotname
#pdf(plotname, height = 10, width = 4) # para 8
pdf(plotname)
listBN_CM <- list(multis[[cl.bn]],plot1,multisCM[[cl.cn]],plot2)
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
dev.off()

listBN <- list(multis[[cl.bn]],plot1)
do.call(grid.arrange, c(listBN, nrow = 2, ncol = 1, as.table = FALSE, top = paste0(measurename) ))






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
  plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
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
  plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
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
# Betweenness for borrador other colorscale and for Thesis!
##################################################################
what <- betweenness
multis = list()
seq <- 2:2
at1 <- seq(0,12,1)
at2 <- seq(0,0.00065,0.00005)
at
a <- 12

cb <- colorRampPalette(brewer.pal(9, "Reds"))(18)

for (i in 1:length(what)){
  plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                          backdrop.theme = "coastline", lonCenter = 180, color.theme = "Reds", set.max = 10,
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,12,2)
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
  plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
                          backdrop.theme = "coastline", lonCenter = 0, color.theme = "Reds", set.max = 10,
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,10,1)
                          # ,
                          # at = at1
                          # ,
                          # set.max = a
  )
  
  multisCM[[i]] <- plot
}
blank <- grid.rect(gp=gpar(col="white"))
listBN_CM <- list(multis[[3]],multis[[5]],multis[[7]],blank, multisCM[[1]], multisCM[[2]],multisCM[[3]],multisCM[[4]])
measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename

# FOR THESIS
multisCM[[3]]


cl.cn <- 3
movbet1.cm<-movingmedias(measures10d[[cl.cn]]$betweenness)
movbet1matrixmean.cm <- apply(movbet1.cm, MARGIN = 3, FUN = mean)

memClim.cm <- quantity2clim(movbet1matrixmean.cm, "bet mov", tas_ncep_10d)
plot2 <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                     set.max = 10, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[cl.cn]), cex = 0.8), color.theme = "Reds",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot2


plotname <- paste0("../plots/comparemeasures_thesis/",measurename,"_raw.pdf")
pdf(plotname)
grid.arrange(multisCM[[cl.cn]],plot2,nrow =1)
dev.off()


grid.arrange(multisCM[[cl.cn]],plot2,nrow =1)

#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 4, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
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
  plot <- spatialPlot(grid = what[[i]], backdrop.theme = "coastline", set.min = NULL,
              set.max = 0.016, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
              color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  #plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
   #                       backdrop.theme = "coastline", 
    #                      colorkey = list(width = 0.6, lables = list(cex = 0.6))
     #                     # ,
      #                    # at = at1
   #                       #set.max = a)
  #)
  multis[[i]] <- plot
}

whatCM <- awconnectivitiesCM
multisCM = list()
for (i in 1:length(whatCM)){
  
  plot <- spatialPlot(grid = whatCM[[i]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = 0.020, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8),
                      color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
 # plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
  #                        backdrop.theme = "coastline",
   #                       colorkey = list(width = 0.6, lables = list(cex = 0.6))
    #                      # ,
     ##                    #set.max = a)
  #)
  
  multisCM[[i]] <- plot
}

measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename

# For resumen 2:
listBN_CM <- list(multis[[4]],multis[[5]], multisCM[[3]], multisCM[[10]])
listBN_CM <- list(multis[[4]],multisCM[[9]])



#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE))
        
#########################
# Voor methods thesis
#########################

awc.bn <-7
movawc1<-movingmedias(measureslistBay[[awc.bn]]$awconnectivity)
movawc1matrixmean <- apply(movawc1, MARGIN = 3, FUN = mean)

memClim <- quantity2clim(movawc1matrixmean, "awc mov", tas_ncep_10d)
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average BN: ",nedgesnetworks[awc.bn]), cex = 0.8), color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1

awc.cn <- 2
movawc1.cm<-movingmedias(measures10d[[awc.cn]]$awconnectivity)
movawc1matrixmean.cm <- apply(movawc1.cm, MARGIN = 3, FUN = mean)

memClim.cm <- quantity2clim(movawc1matrixmean.cm, "awc mov", tas_ncep_10d)
plot2 <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[awc.cn]), cex = 0.8), color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot2

plotname <- paste0("../plots/comparemeasures_thesis/",measurename,"_raw.pdf")
pdf(plotname)
listBN_CM <- list(multis[[awc.bn]],plot1,multisCM[[awc.cn]],plot2)
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

dev.off()


##################################################################
# Local Clustering for theiss
##################################################################
what <- localclusterings
multis = list()
seq <- 2:2
at1 <- seq(0,0.018,0.018/16)
at2 <- seq(0,0.00065,0.00005)
at
a <- 1


for (i in 1:length(what)){
  plot <- spatialPlot(grid = what[[i]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                      color.theme = "Blues",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  #plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
  #                       backdrop.theme = "coastline", 
  #                      colorkey = list(width = 0.6, lables = list(cex = 0.6))
  #                     # ,
  #                    # at = at1
  #                       #set.max = a)
  #)
  multis[[i]] <- plot
}

whatCM <- localclusteringsCM
multisCM = list()
for (i in 1:length(whatCM)){
  
  plot <- spatialPlot(grid = whatCM[[i]], backdrop.theme = "coastline", set.min = NULL, at = seq(0,1,0.1),
                      set.max = NULL, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8),
                      color.theme = "Blues",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  # plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
  #                        backdrop.theme = "coastline",
  #                       colorkey = list(width = 0.6, lables = list(cex = 0.6))
  #                      # ,
  ##                    #set.max = a)
  #)
  
  multisCM[[i]] <- plot
}

measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename

# For resumen 2:
listBN_CM <- list(multis[[4]],multis[[5]], multisCM[[3]], multisCM[[10]])
listBN_CM <- list(multis[[4]],multisCM[[9]])



#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE))

#########################
# Voor methods thesis
#########################
measureslistBay[[6]]$edens
graphObjects[[6]]$graph

lcl.bn<-7
movawc1<-movingmedias(measureslistBay[[lcl.bn]]$localclustering)
movawc1matrixmean <- apply(movawc1, MARGIN = 3, FUN = mean)

memClim <- quantity2clim(movawc1matrixmean, "awc mov", tas_ncep_10d)
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL, 
                     set.max = NULL, lonCenter =0, main =  list(paste0("Neighbor's average BN: ",nedgesnetworks[lcl.bn]),cex =0.8), color.theme = "Blues",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1

lcl.cn <-2
movawc1.cm<-movingmedias(measures10d[[lcl.cn]]$localclustering)
movawc1matrixmean.cm <- apply(movawc1.cm, MARGIN = 3, FUN = mean)

memClim.cm <- quantity2clim(movawc1matrixmean.cm, "awc mov", tas_ncep_10d)
plot2 <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[lcl.cn]), cex = 0.8), color.theme = "Blues",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot2

listBN_CM <- list(multis[[lcl.bn]],plot1,multisCM[[lcl.cn]],plot2)
plotname <- paste0("../plots/comparemeasures_thesis/",measurename,"_raw.pdf")
pdf(plotname)

do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
dev.off()

multisCM[[2]]

##################################################################
# Decrees for theiss
##################################################################
what <- degrees
multis = list()
seq <- 2:2
at1 <- seq(0,0.018,0.018/16)
at2 <- seq(0,0.00065,0.00005)
at
a <- 1


for (i in 1:length(what)){
  plot <- spatialPlot(grid = what[[i]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = 9, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
                      color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  #plot <- spatialPlot(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.8),
  #                       backdrop.theme = "coastline", 
  #                      colorkey = list(width = 0.6, lables = list(cex = 0.6))
  #                     # ,
  #                    # at = at1
  #                       #set.max = a)
  #)
  multis[[i]] <- plot
}

whatCM <- degreesCM
multisCM = list()
for (i in 1:length(whatCM)){
  
  plot <- spatialPlot(grid = whatCM[[i]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8),
                      color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
  
  # plot <- spatialPlot(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.8, pos = 0.25),
  #                        backdrop.theme = "coastline",
  #                       colorkey = list(width = 0.6, lables = list(cex = 0.6))
  #                      # ,
  ##                    #set.max = a)
  #)
  
  multisCM[[i]] <- plot
}

measurename <- attr(what[[1]]$Data,"climatology:fun")
measurename

# For resumen 2:
listBN_CM <- list(multis[[4]],multis[[5]], multisCM[[3]], multisCM[[10]])
listBN_CM <- list(multis[[4]],multisCM[[9]])



#do.call(grid.arrange, c(list, nrow = length(multis), ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE, top = paste0(measurename) ))

do.call(grid.arrange, c(listBN_CM, nrow = 1, ncol = 2, as.table = FALSE))

#########################
# Voor methods thesis
#########################
measureslistBay[[6]]$edens
graphObjects[[6]]$graph

d.bn <- 4
movawc1<-movingmedias(measureslistBay[[d.bn]]$degree)
movawc1matrixmean <- apply(movawc1, MARGIN = 3, FUN = mean)

memClim <- quantity2clim(movawc1matrixmean, "awc mov", tas_ncep_10d)
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main =  list(paste0("Neighbor's average BN: ",nedgesnetworks[d.bn]),cex =0.8), color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1

d.cn <- 2
movawc1.cm<-movingmedias(measures10d[[d.cn]]$degree)
movawc1matrixmean.cm <- apply(movawc1.cm, MARGIN = 3, FUN = mean)

memClim.cm <- quantity2clim(movawc1matrixmean.cm, "awc mov", tas_ncep_10d)
plot2 <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[d.cn]), cex = 0.8), color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot2

plotname <- paste0("../plots/comparemeasures_thesis/",measurename,"_raw.pdf")
pdf(plotname)

listBN_CM <- list(multis[[d.bn]],plot1,multisCM[[d.cn]],plot2)
do.call(grid.arrange, c(listBN_CM, nrow = 2, ncol = 2, as.table = FALSE, top = paste0(measurename) ))
dev.off()
