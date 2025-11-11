#Calculating Adj matrix from grid without rescaling

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

#Calculating Adj matrix from grid with rescaling (al mano)
adj.matrix.function2 <- function(grid){
  time.coords.matrix <- array3Dto2Dmat(grid$Data)
  
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

load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_5d.rda",verbose = TRUE)
adj.matrix.function3(tas_ncep_5d)
adj.matrix.function2(tas_ncep_5d)
adj.matrix.function(tas_ncep_5d)
identical(adj.matrix.function3(tas_ncep_5d),
          adj.matrix.function(tas_ncep_5d))
