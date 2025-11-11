library(bnlearn)
library(mvtnorm)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
###################################################################################
# MORAL GRAPHS OF BAYESIAN NETWORKS: 
# QUESTION IS: HOW TO FIND LOGLIKELIHOODS??
###################################################################################

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm1sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm5sort.rda")
moral2 <- moral(iterationnetworks[[1]])
moral1 <- moral(perm1sort$hc1_8600_8700i)
admat1 <- amat(moral1)
zeros <- which(admat1 == 0, arr.ind = TRUE)

dataCOV <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cormatrix <- cor(data, method = "spearman")
covmatrix<- cov(data, method = "spearman")

test <- glasso(covmatrix, rho = 0, zero = zeros)
littlestest <- dmvnorm(dataCOV, sigma = test$w, log = TRUE)
sum(littlestest)


morals <- lapply(iterationnetworks, moral)
nedgesmorals <- sapply(morals, narcs)

morals <- lapply(perm1sort, moral)
admats <- lapply(morals, amat)
preczeros <- list()
for (i in 1:length(admats)){
  preczeros[[i]] <- which(admats[[i]] == 0, arr.ind = TRUE)
}

precmatslasso <- list()
for (i in 1:length(admats)){
precmatslasso[[i]] <- glasso(cormatrix, thr = 1*exp(-8),rho = 0.01, zero = preczeros[[i]], penalize.diagonal = FALSE)
}

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
M <- colMeans(data)
logliksGLfromBAY <- c()
nedgesGLfromBAY <- c()

precmatslasso[[80]]
precmatslasso[[16]]

lapply(precmatslasso, function(x)x$w)

for (i in 1:length(precmatslasso)){
  littles <- dmvnorm(data,  sigma = precmatslasso[[i]]$w, log = TRUE)
  logliksGLfromBAY[i] <- sum(littles)
}

logliksGLfromBAY <- sort(logliksGLfromBAY, decreasing = TRUE)
plot(rho2, logliksGL)

for (i in 1:length(precmatslasso)){
  m <- precmatslasso[[i]]$wi
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) != 0] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(precgraph))
  nedgesGLfromBAY[i] <- numberofedges
}
plot(nedgesGL, rho2)
plot(nedgesGLfromBAY,logliksGLfromBAY, col = "green")
points(nedgesGLfromBAY,logliksGLfromBAY, col = "blue")
plot(nedgesGLfromBAY, logliks)
