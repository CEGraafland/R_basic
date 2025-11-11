library(igraph)
library(transformeR)

## Calculating Adj matrix from grid without rescaling -----------------------------
adj.matrix.function <- function(grid){
time.coords.matrix <- array3Dto2Dmat(grid$Data)

cor.matrix <- cor(time.coords.matrix,method="spearman")

adj.matrix <- cor.matrix
diag(adj.matrix) <- 0
adj.matrix[adj.matrix <= .8 ] <- 0
adj.matrix[adj.matrix > .8 ] <- 1

return(adj.matrix)
}

adj.matrix.function(tas_ncep_5d_Iberia)

## Calculating Adj matrix from grid with rescaling (manual) -----------------------------------
adj.matrix.function2 <- function(grid){
  
  months <- list()
  for(i in 1:12){
    sea <- subsetGrid(grid, season = i)
    months[[i]] <- climatology(sea)
  }
  
  grid2 <- redim(bindGrid.time(months), drop = T)
  
  matref <- array3Dto2Dmat(grid2$Data)
  mat <- array3Dto2Dmat(grid$Data)
  ind <- rep(1:12,30)
  ind
  
  outmat <- mat 
  
  for (i in 1:12){
    month <- which(ind == i)
    outmat[month,] <- mat[month,] - matref[i,]
  }
  
  
  cor.matrix <- cor(outmat,method="spearman")
  
  adj.matrix2 <- cor.matrix
  diag(adj.matrix2) <- 0
  adj.matrix2[adj.matrix2 <= .8 ] <- 0
  adj.matrix2[adj.matrix2 > .8 ] <- 1
  
  return(adj.matrix2)
}

adj.matrix.function2(tas_ncep_5d_Iberia)

adj.matrix.function3 <- function(grid){
  Rescale <-rescaleGrid(grid)
  time.coords.matrix <- array3Dto2Dmat(Rescale$Data)
  
  cor.matrix <- cor(time.coords.matrix,method="spearman")
  
  adj.matrix3 <- cor.matrix
  diag(adj.matrix3) <- 0
  adj.matrix3[adj.matrix3 <= .8 ] <- 0
  adj.matrix3[adj.matrix3 > .8 ] <- 1
  
  return(adj.matrix3)
}

adj.matrix.function3(tas_ncep_5d_Iberia)

## adjacency matrix function with localscale: ------------------------
adj.matrix.function4 <- function(grid){
  Rescale2 <-localScaling(grid, time.frame = "monthly")
  Rescale <-redim(Rescale2, drop = T)
  time.coords.matrix <- array3Dto2Dmat(Rescale$Data)
  
  cor.matrix4 <- cor(time.coords.matrix,method="spearman")
  
  adj.matrix4 <- cor.matrix4
  diag(adj.matrix4) <- 0
  adj.matrix4[adj.matrix4 <= .8 ] <- 0
  adj.matrix4[adj.matrix4 > .8 ] <- 1
  
  return(adj.matrix4)
}



load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_5d.rda",verbose = TRUE)
adj.matrix.function3(tas_ncep_5d)
tas.adj.matrix2 <- adj.matrix.function2(tas_ncep_5d)
tas.adj.matrix2.nais0 <-tas.adj.matrix2
tas.adj.matrix2.nais0[is.na(tas.adj.matrix2.nais0)] <- 0
adj.matrix.function(tas_ncep_5d)
tas.adj.matrix4 <-adj.matrix.function4(tas_ncep_5d)

## NETWORK ANALYSE
library(igraph)
# seems like rescalegrid does not change a lot:
#No rescaling


load("/Users/lisettegraafland/Documents/R_practice/Data/coord_5d.rda")
load("/Users/lisettegraafland/Documents/R_practice/Data/tas_adj_matrix.rda")

#First impression
tas.graph2<- graph_from_adjacency_matrix(tas.adj.matrix2, mode = "undirected")
plot(tas.graph2, layout = coord.5d)
tas.graph4<- graph_from_adjacency_matrix(tas.adj.matrix4, mode = "undirected")
plot(tas.graph4, layout = coord.5d)
tas.graph2.nais0 <- graph_from_adjacency_matrix(tas.adj.matrix2.nais0, mode = "undirected")
#edges are iguql with or without NA treatment:
identical(degree_distribution(tas.graph2), degree_distribution(tas.graph2.nais0))

# edge density
edge_density(tas.graph2, loops = TRUE)
edge_density(tas.graph4, loops = TRUE)

# degree distribution
degree_distribution(tas.graph4)
plot(degree_distribution(tas.graph4))

#betweenness
betw.tas.graph2 <- betweenness(tas.graph2)
#betweenness edges
edge.betw.tas.graph2 <- edge.betweenness(tas.graph2)

#convert spatial object for example Betweenneess BEKIJKEN COORDINATEN SYSTEEM
library(devtools)
install_github('SantanderMetGroup/mopa')
library(mopa)
data(wrld)
library(sp)

betw.matrix <- cbind(coord.5d,betw.tas.graph2)
points <- SpatialPoints(betw.matrix[,1:2])
betw.sp <- SpatialPixelsDataFrame(points,data.frame(betw.matrix[,3]))
spplot(betw.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

betw.matrix.yx <- cbind(coord.5dyx,betw.tas.graph2)
pointsyx <- SpatialPoints(betw.matrix.yx[,c(2,1)])
betw.sp.yx <- SpatialPixelsDataFrame(pointsyx,data.frame(betw.matrix.yx[,3]))
spplot(betw.sp.yx, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

#convert spatial object for example density BEKIJKEN COORDINATEN SYSTEEM   

#Find NA values
# in adjacency matrix:
which(is.na(tas.adj.matrix2), arr.ind=TRUE)
# in Dataset 
str(tas_ncep_5d$Data)
NAlocation <- which(is.na(tas_ncep_5d$Data), arr.ind=TRUE)
unique(NAlocation[,2])
length(unique(NAlocation[,2])) #in every latitude
unique(NAlocation[,3]) 
length(unique(NAlocation[,3])) # in every longitude
unique(NAlocation[,1]) 
length(unique(NAlocation[,1])) #in every month

numberNA <-sum(is.na(tas_ncep_5d$Data))
length.tas <- length(tas_ncep_5d$Data)
length.tas/numberNA

colSums(is.na(tas_ncep_5d$Data))[1:13,]
colSums(is.na(tas_ncep_5d$Data))[14:26,]
colSums(is.na(tas_ncep_5d$Data))[27:37,]



tas_ncep_5d.clim <- climatology(tas_ncep_5d)

plotClimatology(tas_ncep_5d.clim, backdrop.theme = "countries")
