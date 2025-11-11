###########################################################################################
# Basic complex network functions
###########################################################################################
library(igraph)
library(transformeR)
library(devtools)
# library(mopa)
# data(wrld)
library(sp)
library(magrittr)
library(abind)
###########################################################################################
# Get complex graph from grid.
###########################################################################################
graph_from_Grid <- function(grid,th = 0.8, subind = NULL, method = c("spearman", "pearson")) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  
  
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]}
  
  # Correlation matrix
  cor.matrix <- cor(time.coords.matrix, method = method)
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

graph_from_cov <- function(grid, covar,th = 0.8) {
  
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  
  
  cor.matrix <- cov2cor(covar)
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

graph_from_cor <- function(grid, cormat,th = 0.8) {
  
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  
  
  cor.matrix <- cormat
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
################################################################################################
# Get complex measure for each gridpoint in grid
# different measures: edge density / degree distribution / betweenness centrality / closeness
# centrality / Area weighted connectivitiy / local clustering coefficient
################################################################################################
graph2measure.clos <- function(graphObj, quan) {
  
  # closeness edges
  close <- closeness(graphObj$graph, mode ="all")
 
  
  out <- list("closeness" = close)
  attr(out, "Xcoords") <- attr(graphObj, "Xcoords")
  attr(out, "Ycoords") <- attr(graphObj, "Ycoords")
  attr(out, "ref.dates") <- attr(graphObj, "ref.dates")
  return(out)
}

graph2measure <- function(graphObj) {
  
  # edge density
  edens <- edge_density(graphObj$graph, loops = TRUE)
  # degree distrib
  ddist <- degree_distribution(graphObj$graph)
  # degree distrib
  degree <- igraph::degree(graphObj$graph)
  # strength
  if(!is.null(edge_attr(graphObj$graph))){
    strength <- igraph::strength(graphObj$graph)
    } else {strength <- NA}

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
              "localclustering" = localclus,
              "degree" = degree)
  attr(out, "Xcoords") <- attr(graphObj, "Xcoords")
  attr(out, "Ycoords") <- attr(graphObj, "Ycoords")
  attr(out, "ref.dates") <- attr(graphObj, "ref.dates")
  return(out)
}

# graphObj <- gridGraphsCN$`855`
# graphObj <- gridGraphwCN2$`8637`
# graphObj <- gridGraphsBN$hc1_8600_8700i
 # localmeasObj <- gridGraphwCN2$`8637`
 # # graphObj <- gridGraphsCN[[1]]
# 
# edge_attr(gridGraphwBN$hc1_1300_1400i$graph)
# edge_attr(gridGraphsBN$hc1_1300_1400i$graph)
# edge_attr(gridGraphs[[1]]$graph)
# edge_attr(gridGraphwCN2$`8637`$graph)

graph2measure.local <- function(graphObj) {
  
  # degree distrib
  degree <- igraph::degree(graphObj$graph)
  
  # strength
  if(!is.null(edge_attr(graphObj$graph))){
    strength <- igraph::strength(graphObj$graph)
  } else {strength <- NA}
  
  # betweenness
  betw <- log1p(betweenness(graphObj$graph))
  
  # closeness edges
  close <- closeness(graphObj$graph, mode ="all", normalized = TRUE)
  # area weighted connectivity
  sumArea <- sum(cos(graphObj$VertexCoords$y/(180)*pi))
  awconnectivity <- as.vector(cos(graphObj$VertexCoords$y/(180)*pi)%*%graphObj$adjacency) / sumArea
  # transitivity (local clustering)
  localclus <- igraph::transitivity(graphObj$graph, type = "local")
  barratclus <- igraph::transitivity(graphObj$graph, type = "barrat", weights = NA)
  # for weighted:
  if(!is.null(edge_attr(graphObj$graph))){
  barratweightclus <- igraph::transitivity(graphObj$graph, type = "weighted", weights =  E(graphObj$graph)$weight)
  } else {barratweightclus <- NA}
  
  # Average length of shortest path (with weights)
  avpathlength <- (1/close)
  
  
  out <- list("degree" = degree,
              "strength" = strength,
              "betweenness" = betw, 
              "closeness" = close, 
              "awconnectivity" = awconnectivity,
              "localclustering" = localclus,
              "barratclustering"= barratclus,
              "barratweightclustering" = barratweightclus,
              "avpathlength" = avpathlength
              )
  attr(out, "Xcoords") <- attr(graphObj, "Xcoords")
  attr(out, "Ycoords") <- attr(graphObj, "Ycoords")
  attr(out, "ref.dates") <- attr(graphObj, "ref.dates")
  return(out)
}

# graphObj <- gridGraphsCN$`855`
# localmeasureObj <- graph2measure.local(graphObj)
# graphObj <- gridGraphwCN2$`8637`
# localmeasureObj <- graph2measure.local(gridGraphwCN2$`8637`)
graph2measure.global <- function(graphObj, localmeasureObj = NULL) {
  
  # edge density
  edens <- edge_density(graphObj$graph, loops = TRUE)
  # diameter
  diameterLCC <- diameter(graphObj$graph, directed = FALSE, unconnected = TRUE, weights = NULL)
  diameterG <- diameter(graphObj$graph, directed = FALSE, unconnected = FALSE, weights = NULL)
  # mean distance
  meandistLCC <- mean_distance(graphObj$graph, directed = FALSE, unconnected = TRUE)
  meandistG <- mean_distance(graphObj$graph, directed = FALSE, unconnected = FALSE)
  # mean average path length
  if(!is.null(localmeasureObj)){
    apl<- localmeasureObj$avpathlength
    meanavpathlength <- sum(apl, na.rm = TRUE)/length(apl[!is.na(apl)])
  } else {meanavpathlength <- NA}
  # transitivity (global clustering)
  globalclus <- igraph::transitivity(graphObj$graph, type = "global")
  # transitivity 2 (globalbarratclustering)
  if(!is.null(localmeasureObj)){
    lbc<- localmeasureObj$barratclustering
    globalbarratclus <- sum(lbc, na.rm = TRUE)/length(lbc[!is.na(lbc)])
  } else {globalbarratclus <- NA}
  # transtiivity 3 (globalweightedbarratclustering)
  if(!is.null(localmeasureObj)&!is.null(edge_attr(graphObj$graph))){
    lwbc <- localmeasureObj$barratweightclustering
    globalbarratweightclus <- sum(lwbc, na.rm = TRUE)/length(lwbc[!is.na(lwbc)])
  } else {globalbarratweightclus <- NA}
  # global assortativity
  globalass <- igraph::assortativity_degree(graphObj$graph, directed = "false")
  
  out <- list("edens" = edens,
              "diameterLCC" = diameterLCC,
              "diameterG" = diameterG,
              "meandistLCC" = meandistLCC,
              "meandistG" = meandistG,
              "meanavpathlength" = meanavpathlength,
              "globalclus" = globalclus,
              "globalbarratclus" = globalbarratclus,
              "globalbarratweightclus"= globalbarratweightclus,
              "globalass" = globalass)
  attr(out, "Xcoords") <- attr(graphObj, "Xcoords")
  attr(out, "Ycoords") <- attr(graphObj, "Ycoords")
  attr(out, "ref.dates") <- attr(graphObj, "ref.dates")
  return(out)
}
# test <- gridGraphsBN[[10]]
# test$graph
# graphObj <- test
# testlocalmeasure <- graph2measure.local(test)
# testdegreemeasure <- graph2measure.degree(test,testlocalmeasure)
# testdegreemeasure$degrdist
# testdegreemeasure$degreeclustering
# graphObj <- gridGraphsBN[[4]] 
# localmeasureObj <- graph2measure.local(graphObj)
####################################################
# meandistance per degree = meanstrength per degree / degree
# 
######################################################

graph2measure.degree <- function(graphObj, localmeasureObj) {
  
  # degree distrib
  degree <- igraph::degree(graphObj$graph)
  
  # strength
  if(!is.null(edge_attr(graphObj$graph))){
    strength <- igraph::strength(graphObj$graph)
    distances <- edge_attr(graphObj$graph, name = "weight")
    maxdegree <- max(degree)
    meanstrength <- numeric(maxdegree+1)
    meandistance <- numeric(maxdegree+1)
    for( k in 0:maxdegree){
      kaas <- which(degree == k, arr.ind = TRUE)
      s_k <- sum(strength[kaas])/length(kaas)
      d_k <- sum(strength[kaas])/(length(kaas)*k)
      meanstrength[k+1] <- s_k
      meandistance[k+1] <- d_k
    } 
  } else {strength <- NA
          meanstrength <- NA
          meandistance <- NA}
  
  # degree distrib
  ddist <- degree_distribution(graphObj$graph)
  
  # degree barrat / weighted clustering
  if(is.null(edge_attr(graphObj$graph))){
  maxdegree <- max(degree)
  degreeclus <- numeric(maxdegree+1)
  for( k in 0:maxdegree){
    kaas <- which(degree == k, arr.ind = TRUE)
    c_k <- sum(localmeasureObj$barratclustering[kaas])/length(kaas)
    degreeclus[k+1] <- c_k
    } 
  } else {
    maxdegree <- max(degree)
    degreeclus <- numeric(maxdegree+1)
    for( k in 0:maxdegree){
      kaas <- which(degree == k, arr.ind = TRUE)
      c_k <- sum(localmeasureObj$barratweightclustering[kaas])/length(kaas)
      degreeclus[k+1] <- c_k
    } 
  }
  
  
  out <- list("degrdist" = ddist,
              "degreeclustering" = degreeclus,
              "meanstrength" = meanstrength,
              "meandistance" = meandistance)
  attr(out, "Xcoords") <- attr(graphObj, "Xcoords")
  attr(out, "Ycoords") <- attr(graphObj, "Ycoords")
  attr(out, "ref.dates") <- attr(graphObj, "ref.dates")
  return(out)
}


#######################################################################################
# Convert local measure vector to climatology object for visualizing.
#######################################################################################

measure2clim <- function(measureObj, what = c("betweenness","closeness","awconnectivity","localclustering","degree","strength", "barratclustering",
                                              "barratweightclustering", "strength"), ref.grid) {
  if (what == "betweenness") {
    mat <- matrix(measureObj$betweenness, nrow = 1)
  }
  if (what == "degree") {
    mat <- matrix(measureObj$degree, nrow = 1)
  }
  if (what == "strength") {
    mat <- matrix(measureObj$strength, nrow = 1)
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
  
  if (what == "barratclustering") {
    mat <- matrix(measureObj$barratclustering, nrow = 1)  
  }
  
  if (what == "barratweightclustering") {
    mat <- matrix(measureObj$barratweightclustering, nrow = 1)  
  }
  
  ref.grid$Data <- mat2Dto3Darray(mat, x = attr(measureObj, "Xcoords"), y = attr(measureObj, "Ycoords"))
  attr(ref.grid$Data, "climatology:fun") <- what
  return(ref.grid)
}

#######################################################################################
# Convert measure/quantity/probabilty vector to climatology object for visualizing.
#######################################################################################
quantity2clim <- function(quantity, what, ref.grid, backperm = NULL) {
  if(!is.null(backperm)){quantity <- quantity[backperm]}
  mat <- matrix(quantity, nrow = 1)  
  ref.grid$Data <- mat2Dto3Darray(mat, x = ref.grid$xyCoords$x , y = ref.grid$xyCoords$y)
  attr(ref.grid$Data, "climatology:fun") <- what
  return(ref.grid)
}
# str(tas_ncep_10d)
# ref.grid <- tas_ncep_10d
#######################################################################################
# Plot of degree distributions
# 
#######################################################################################
plplot <- function(grid, tau,method = c("pearson")){
  graphObj <- graph_from_Grid(grid, th = tau,method = method)
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



plplotlist <- function(grid, tau){
  par(mar=c(4,4,2,2))
  par(mfrow = c(length(tau)/2, 4))
  
  for (i in 1:length(tau)){
    graphObj <- graph_from_Grid(grid, th = tau[i],method = 'pearson')
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
#########################################################################################
# logarithmic bining of degree distribution. 
#########################################################################################
plplotlist_logbin <- function(grid, tau){
  
  ggplots = list()
  for (i in 1:length(tau)){
    graphObj <- graph_from_Grid(grid, th = tau[i])
    numberofedges <- length(E(graphObj$graph))
    
    degrees <- igraph::degree(graphObj$graph)
    k_max <- max(degrees)
    nbins <- ceiling(log2(k_max))
    binlims <- 2^(1:nbins)
    binlims <- append(binlims, 0, after = 0)
    
    histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
    histdegrees2 <- histdegrees
    histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
    
    # Convert to adequate dataframe for ggplot.
    # Plot on log2 log2 scale. 
    datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
    names(datahist) <- c("breaks","counts")
    plot <- ggplot(datahist, aes(breaks,counts)) + 
      scale_x_continuous(trans = "log2") + 
      scale_y_continuous(trans = "log2", labels = fmt_dcimals(2)) + 
      geom_point() +
      ggtitle(paste0("CN: ",numberofedges," tau = ",tau[i])) +
      xlab("degree k") +
      ylab("frequency k /\nbinwidth")
    ggplots[[i]] <- plot
  }
  
  return(ggplots)
}
