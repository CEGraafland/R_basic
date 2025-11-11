library(igraph)

library(transformeR)
library(devtools)
#library(mopa)
#data(wrld)
library(sp)
library(magrittr)
library(abind)

#introduce grid (NCEP_Reanalysis)

graph_from_Grid <- function(grid, th = 0.8) {
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
  
  # Correlation matrix
  cor.matrix <- cor(time.coords.matrix, method = "spearman")
  abs.cor.matrix <- abs(cor.matrix)
  adj.matrix <- abs.cor.matrix
  # Adjacency matrix
  diag(adj.matrix) <- 0
  adj.matrix[adj.matrix <= th ] <- 0
  adj.matrix[adj.matrix > th ] <- 1
  # Graph
  graph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  
  out <- list("graph" = graph,
              "data_coords" = time.coords.matrix,
              "correlation" = cor.matrix,
              "VertexCoords" = ref.coords ,
              "adjacency" = adj.matrix)
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  return(out)
}
grid =tas_ncep_10d
th = 0.4
subind = ind1981_2010
subind

localScaling(tas_ncep_10d, scale = TRUE)

graph_from_Grid <- function(grid,th = 0.8, subind = NULL) {
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
  time.coords.matrix[1,]
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]}
  
  # Correlation matrix
  cor.matrix <- cor(time.coords.matrix, method = "spearman")
  
  # Graph
  graph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  
  out <- list("graph" = graph,
              "data_coords" = time.coords.matrix,
              "correlation" = cor.matrix,
              "VertexCoords" = ref.coords ,
              "adjacency" = adj.matrix)
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  return(out)
}

#Test graph_from_Grid
test5d <- graph_from_Grid(tas_ncep_5d2, th = .68)
test10d <- graph_from_Grid(tas_ncep_10d, th = .3)
igraph::ed

str(test5d$adjacency)
str(test5d$VertexCoords$y)
str(test5d$VertexCoords$x)
str(test5d$correlation)
str(test5d$adjacency)

graph2measure <- function(graphObj) {
  
  # edge density
  edens <- edge_density(graphObj$graph, loops = TRUE)
  # degree distrib
  ddist <- degree_distribution(graphObj$graph)
  # betweenness
  betw <- log1p(betweenness(graphObj$graph))
  # betweenness edges
  edge.betw <- edge.betweenness(graphObj$graph)
  # closeness edges
  close <- closeness(graphObj$graph, mode ="all")
  # area weighted connectivity
  sumArea <- sum(cos(graphObj$VertexCoords$y/(180)*pi))
  awconnectivity <- as.vector(cos(graphObj$VertexCoords$y/(180)*pi)%*%graphObj$adjacency) / sumArea
  # transitivity (local clustering)
  localclus <- igraph::transitivity(graphObj$graph, type = "local")
  
  out <- list("edens" = edens,
              "ddist" = ddist, 
              "betweenness" = betw, 
              "closeness" = close, 
              "awconnectivity" = awconnectivity,
              "localclustering" = localclus)
  attr(out, "Xcoords") <- attr(graphObj, "Xcoords")
  attr(out, "Ycoords") <- attr(graphObj, "Ycoords")
  attr(out, "ref.dates") <- attr(graphObj, "ref.dates")
  return(out)
}

#Test graph2measure
measuretest5d <- graph2measure(test5d)
measuretest10d<- graph2measure(test10d)

for (a in c(0.2,0.3,0.4,)){
print(paste0("edgedensity =", graph2measure(graph_from_Grid(tas_ncep_10d, th = a))$edens))}

measuretest5d$edens
str(measuretest5d$ddist)

str(measuretest5d$awconnectivity)
str(measuretest5d$closeness)
str(measuretest5d$betweenness)
str(measuretest5d$localclustering)

measuretest10d$awconnectivity
test10d$graph

sumArea <- sum(cos(graphObj$VertexCoords$y/(180)*pi))
awconnectivity <- as.vector(cos(graphObj$VertexCoords$y/(180)*pi)%*%graphObj$adjacency) / sumArea

measure2clim <- function(measureObj, what = c("betweenness","closeness","awconnectivity","localclustering"), ref.grid) {
  if (what == "betweenness") {
    mat <- matrix(measureObj$betweenness, nrow = 1)
  }
  if (what == "closeness") {
    mat <- matrix(measureObj$closeness, nrow = 1)
  }
  if (what == "awconnectivity") {
    mat <- matrix(measureObj$awconnectivity, nrow = 1)
  }
  if (what == "localclustering") {
    mat <- matrix(measureObj$localclustering, nrow = 1)  
  }
  ref.grid$Data <- mat2Dto3Darray(mat, x = attr(measureObj, "Xcoords"), y = attr(measureObj, "Ycoords"))
  attr(ref.grid$Data, "climatology:fun") <- what
  return(ref.grid)
}

#Climatology of measures
clim.awcon.test5d <- measure2clim(measuretest5d, what = "awconnectivity", ref.grid = tas_ncep_5d2)
clim.betw.test5d <- measure2clim(measuretest5d, what = "betweenness", ref.grid = tas_ncep_5d2)
clim.close.test5d <- measure2clim(measuretest5d, what = "closeness", ref.grid = tas_ncep_5d2)
clima.lclus.test5d <- measure2clim(measuretest5d, what = "localclustering", ref.grid = tas_ncep_5d2)

#plotClimatology of Climatology of measures
plotClimatology(clim.betw.test5d, backdrop.theme = "coastline", at = seq(0,13,.5), col.regions = topo.colors(69))
plotClimatology(clim.betw.test5d, backdrop.theme = "coastline", at = seq(0,15,0.5), col.regions = colores2(69))
plotClimatology(clim.close.test5d, backdrop.theme = "coastline",at = seq(0.0000025,0.000003,0.00000005))
plotClimatology(clim.awcon.test5d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.test5d,backdrop.theme = "coastline")

# Makes complex graph with graph_from_Grid.
# Plots degreedistribution of graph for different cor. tresholds.
# Entries:  Grid, correlation treshold
# Outcome:  two plots:  log degreedistr vs log degree k . 
#                       degreedistr vs degree k  

plplot <- function(grid, tau){
  graphObj <- graph_from_Grid(grid, th = tau)
  list_graphObj <- graph2measure(graphObj)
  
  par(mfrow = c(1, 2))
  
  graph1 <- plot(1:length(list_graphObj$ddist),
                 list_graphObj$ddist, 
                 main = paste0("degree_distribution\n tau = ", tau), 
                 type = "l", 
                 xlab = "degree K" , 
                 ylab = "p(K)")
  graph2 <- plot(log(1:length(list_graphObj$ddist)),
                 log(list_graphObj$ddist), main = paste0("degree_distribution loglog  \ntau = ", tau), 
                 type = "l", 
                 xlab = "log (degree K)" , 
                 ylab = "log (p(K))")
}

# Same function as plplot but with capacidad to 
# entrie entire list of tau's.
# Outcome: Organized plots of different tau's.

grid <- tas_ncep_10d
th <- 0.5

plplotlist <- function(grid, tau){
  par(mar=c(4,4,2,2))
  par(mfrow = c(length(tau)/2, 4))
  
  for (i in 1:length(tau)){
    graphObj <- graph_from_Grid(grid, th = tau[i])
    list_graphObj <- graph2measure(graphObj)
    numberofedges <- length(E(graphObj$graph))
   plot(1:length(list_graphObj$ddist),
         list_graphObj$ddist, 
         main = paste0("degree_distribution\n tau = ", tau[i]," |E| = ",numberofedges), 
         type = "l", 
         xlab = "degree K" , 
         ylab = "p(K)")
   plot(log(1:length(list_graphObj$ddist)),
        log(list_graphObj$ddist), main = paste0("degree_distribution loglog\n tau = ", tau[i]," |E| = ",numberofedges), 
        type = "l", 
        xlab = "log (degree K)" , 
        ylab = "log (p(K))")
  }
}

load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_10d2.rda")

plplotlist(tas_ncep_10d,c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))
plplotlist(tas_ncep_5d2,c(0.1,0.15,0.2,0.3,0.4,0.6,0.7,0.8))
dev.off()


# Show powerlaw plots quickly in RSTudio 
for( a in c(0.2,0.4,0.6,0.7,0.8)){
  par(mfrow = c(5, 1))
  plplot(tas_ncep_10d, a)
}
for( a in c(0.2,0.4,0.6,0.7,0.8)){
  plplot(tas_ncep_5d2, a)
}
# create powerlaw plots
for(i in c(0.25,0.3,0.35,0.45,0.5,0.55,0.65)){
  tau<-i
  plplotname <- paste0("/Users/lisettegraafland/Documents/R_practice/plots/powerlaw/powlaw10da",tau,".pdf")
  pdf(plplotname)
  plplot(tas_ncep_10d2,tau)
  dev.off()
}
