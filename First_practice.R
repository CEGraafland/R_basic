library(transformeR)
load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_5d_Iberia.rda",verbose = TRUE)


months <- list()
for(i in 1:12){
sea <- subsetGrid(tas_ncep_5d_Iberia, season = i)
months[[i]] <- climatology(sea)
}

tas2 <- redim(bindGrid.time(months), drop = T)

matref <- array3Dto2Dmat(tas2$Data)
mat <- array3Dto2Dmat(tas_ncep_5d_Iberia$Data)
ind <- rep(1:12,30)
ind
outmat <- mat

for (i in 1:12){
  month <- which(ind == i)
  outmat[month,] <- mat[month,] - matref[i,]
}



 #Create annual means
tas_annual_5d <- aggregateGrid(tas_ncep_5d, aggr.d = list(FUN = "mean", na.rm = TRUE),
aggr.m = list(FUN = "mean", na.rm = TRUE),
aggr.y = list(FUN = "mean", na.rm = TRUE))
                                                                 
#create a grid structure that fits with the output correlation matrix
                                                                 
corgrid <- redim(climatology(obs), drop = T)
                                                                 
#extract coordinates from one of the grids (assuming that both grids have the same spatial grid)
                                                                 
x.coord <- getCoordinates(tas_ncep_5d)$x
y.coord <- getCoordinates(tas_ncep_5d)$y
                                                                 
#create a coordinate index for the loop + coordinate vector
                                                                 
coords <- expand.grid(1:length(x.coord), 1:length(y.coord))
coord.5d <- expand.grid( x.coord , y.coord )
coord.5d <- as.matrix(coord.5d)
save(coord.5d, file = "/Users/lisettegraafland/Documents/R_practice/Data/coord_5d.rda")
                                                                 
#Convert 3D grid of Data [time,lat,lon] to 2D matrix [time, grid-Point]
tas.time.coords.matrix <- array3Dto2Dmat(obs$Data)
save(tas.time.coords.matrix, file = "/Users/lisettegraafland/Documents/R_practice/Data/tas_time_coords_matrix.rda")
                                                                 
#Convert 3D grid of annual means [time,lat,lon] to 2D matrix [time, grid-Point]
tas.annual.time.coords.matrix <- array3Dto2Dmat(tas_annual_5d$Data)
                                                                 
#Matrix annual means for each month (x12)
tas.annual.time.coords.matrix[rep(seq_len(nrow(tas.annual.time.coords.matrix)), each=12),]
                                                                 
#anomaly time series.
tas.anomaly.time.coords.matrix <- tas.time.coords.matrix - tas.annual.time.coords.matrix[rep(seq_len(nrow(tas.annual.time.coords.matrix)), each=12),]
#Normalize to unit variance
                                                                 
 #calculate correlation between the timeseries of the coordinates
tas.cor.matrix <- cor(tas.time.coords.matrix,method="spearman")
                                                                 
                                                                 
#calculate correlation between the anomaly timeseries of the coordinates
tas.anomaly.cor.matrix <- cor(tas.anomaly.time.coords.matrix,method="spearman")
                                                                 
#create Adjacency matrix \theta = 0.8 
tas.adj.matrix <- tas.anomaly.cor.matrix
diag(tas.adj.matrix) <- 0
tas.adj.matrix[tas.adj.matrix <= .8 ] <- 0
tas.adj.matrix[tas.adj.matrix > .8 ] <- 1
save(tas.adj.matrix, file = "/Users/lisettegraafland/Documents/R_practice/Data/tas_adj_matrix.rda")
                                                                 
                                                                 