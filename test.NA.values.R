#TEST NA VALUES
library(transformeR)
library(igraph)
adj.matrix.function2 <- function(grid){
  
  months <- list()
  for(i in 1:12){
    sea <- subsetGrid(tas_ncep_5d, season = i)
    months[[i]] <- climatology(sea)
  }
  
  grid2 <- redim(bindGrid.time(months), drop = T)
  
  matref <- array3Dto2Dmat(tas_ncep_5d$Data)
  
mat <- array3Dto2Dmat(tas_ncep_5d$Data)
  ind <- rep(1:12,30)
  ind
  
  outmat <- mat 
  
  for (i in 1:12){
    month <- which(ind == i)
    outmat[month,] <- mat[month,] - matref[i,]
  }
  
  
  cor.matrix.naoc <- cor(outmat, use = "na.or.complete", method="spearman")
  adj.matrix.naoc <- cor.matrix.naoc
  diag(adj.matrix.naoc) <- 0
  adj.matrix.naoc[adj.matrix.naoc <= .8 ] <- 0
  adj.matrix.naoc[adj.matrix.naoc > .8 ] <- 1
  adj.matrix.naoc[is.na(adj.matrix.naoc)] <- 0
  
  cor.matrix.pac <- cor(outmat,use = "pairwise.complete.obs", method = "spearman")
  adj.matrix2 <- cor.matrix
  diag(adj.matrix2) <- 0
  adj.matrix2[adj.matrix2 <= .8 ] <- 0
  adj.matrix2[adj.matrix2 > .8 ] <- 1
  
  return(adj.matrix2)
}

str(cor.matrix.naoc)
str(adj.matrix.naoc)
identical(cor.matrix.naoc,cor.matrix.pac)
#list rows of data dat having missing values
if complete.cases(tas_ncep_5d$Data)
str(tas_ncep_5d$Data)

complete.cases(tas.adj.matrix2)
?complete.cases
tas_ncep_5d$Data[1,,]
tas_ncep_5d$Data[!complete.cases(tas_ncep_5d$Data),,]
?

graph_from_adjacency_matrix(adj.matrix.naoc)



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

#of gewoon in het plaatje kijken
tas_ncep_5d.clim <- climatology(tas_ncep_5d)

plotClimatology(tas_ncep_5d.clim, backdrop.theme = "countries")

tas_ncep_10d.clim <- climatology(tas_ncep_10d)

plotClimatology(tas_ncep_10d.clim, backdrop.theme = "countries")
