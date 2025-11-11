#######################################################################################
# Obtain loglikelihoods Bayesian networks
#######################################################################################
rm(list = ls())
library(glasso)
library(igraph)
library(transformeR)
library(visualizeR)
library(bnlearn)
library(mvtnorm)
library(corpcor)
library(gridExtra)
library(condMVNorm)
library(RColorBrewer)
# source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
# source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
# source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")

source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("../R/Functions/propagationFunctions.R")
load("../Data/tas_ncep_10d.rda")


f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
#####################################################################################
# Load data HC eBIC 
#######################################################################################
for(j in c(1,2,3,4,5)){
  pattern <- paste0("hc_",j,"_eBIC_g")
  filestabu1 <- list.files("../Data/Struct_learn/HC", full.names = T, pattern = pattern)
  filestabu1names <- list.files("../Data/Struct_learn/HC", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1$networks[[i]] <- variablepos$networks[[1]]
    names[i] <- names(variablepos$networks)
    data.frames[i,] <- as.matrix(variablepos$networkdata)
  }
  colnames(data.frames) <- names(variablepos$networkdata)
  data.frames <- as.data.frame(data.frames)
  variablelisttabu1$networkdata<- data.frames
  variablelisttabu1$networkdata
  
  names(variablelisttabu1$networks) <- names
  assign(paste0("hc_",j,"_eBIC"),variablelisttabu1)
}
####################################################################################
# load HC iteration data and make list
####################################################################################
load("../Data/networks_hcnetworks10d.rda")
#load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda")
load(file ="../Data/hc_edges_loglik_10d_200i.rda")
x <- seq(200,8600,200)
y <- seq(400,8800,200)
i <- 200
iterationnetworks <- list()

iterationnetworks[[1]] <- hc_edges_loglik_10d_200i$networks
for(i in 1:length(x)){
  iterationnetworks[[i+1]] <- eval(parse(text = paste0("hc_edges_loglik_10d_",x[i],"_",y[i],"i$networks")))
}
iterationnetworks
####################################################################################
# load HC iteration data permutation 1 and make list
####################################################################################
for(j in c(1)){
  pattern <- paste0("hc1_")
  filestabu1 <- list.files("../Data/Struct_learn/hciterations/perm1", full.names = T, pattern = pattern)
  filestabu1names <- list.files("../Data/Struct_learn/hciterations/perm1", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list()
  names <- c()
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1[[i]] <- variablepos
  }

iterationnetworks <- variablelisttabu1
}
x <- seq(0,8600,100)
y <- seq(100,8700,100)
i <- 100

iterationnetworks <- list()

iterationnetworks[[1]] <- hc_edges_loglik_10d_200i$networks
for(i in 1:length(x)){
  iterationnetworks[[i+1]] <- eval(parse(text = paste0("hc_edges_loglik_10d_",x[i],"_",y[i],"i$networks")))
}

####################################################################################
# hc iteration analyse
####################################################################################3
dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
iterationfits <- lapply(iterationnetworks, bn.fit, data = dataRMS)

logliksIT <- c()
for (i in 1:length(iterationfits)){
  logliksIT[i] <- logLik(iterationfits[[i]], data = dataRMS)
}
logliksIT
nedgesIT <- sapply(iterationnetworks, narcs)

data.frame(nedges = nedgesIT, logliks = logliksIT)
plot(nedgesIT,logliksIT, col = "blue", xlim = c(0,3000))
####################################################################################
# hc eBIC analyse
####################################################################################3
dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
fits <- lapply(hc_1_eBIC$networks, bn.fit, data = dataRMS)

loglikseBIC <- c()
for (i in 1:length(fits)){
  loglikseBIC[i] <- logLik(fits[[i]], data = dataRMS)
}
nedgeseBIC <- sapply(hc_1_eBIC$networks, narcs)
sort(nedgeseBIC)
data.frame(nedges = nedgeseBIC, logliks = loglikseBIC)
points(nedgeseBIC,loglikseBIC, col = "red")
###################################################################################
# glasso different rhos without tresholds PRECISION
# Loglik of glasso with rho
# Nedges of precision network (wi) without treshold
###################################################################################
rho <- c(0.01,0.05,0.063,0.07,0.09,0.12,0.13,0.165,0.17,0.20,0.225,0.25)
rhoplus <- c(c(0.02,0.03,0.04,0.27,0.28,0.30,0.32,0.34,0.38,0.42,0.46,0.50),seq(0.52,0.98,0.02))
rho2 <- c(rho,rhoplus)
rho2 <- sort(rho2)
length(rho2)
filesglasso <- list.files("../Data/glasso/cornormal", full.names = T)
length(filesglasso)
glassolist <- list()
names <- c()
temp.space <- new.env()

for (i in 1:length(filesglasso)){
  temp.space <- new.env()
  variablepos <- get(load(filesglasso[i], temp.space),envir = temp.space)
  names[i] <- ls(envir = temp.space)
  glassolist[[i]] <- variablepos
  rm(temp.space)
}
names(glassolist) <- names
glassolist

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
M <- colMeans(data)
logliksGL <- c()
nedgesGL <- c()

for (i in 1:length(glassolist)){
littles <- dmvnorm(data, mean = M, sigma = glassolist[[i]]$w, log = TRUE)
logliksGL[i] <- sum(littles)
}

logliksGL <- sort(logliksGL, decreasing = TRUE)
plot(rho2, logliksGL)

for (i in 1:length(glassolist)){
  m <- glassolist[[i]]$wi
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) != 0] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(precgraph))
  nedgesGL[i] <- numberofedges
}
plot(nedgesGL, rho2)
plot(nedgesGL,logliksGL, col = "green")
points(nedgesGL,logliksGL, col = "green")
###################################################################################
# glasso Covariance of RMS 
# different rhos without tresholds 
# Loglik of glasso with rho
# Nedges of precision network (wi) without treshold
###################################################################################
rho <- c(0.01,0.05,0.063,0.07,0.09,0.12,0.13,0.165,0.17,0.20,0.225,0.25)
filesglasso <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/covRMS", full.names = T)

glassolistCOVRMS <- list()
names <- c()
temp.space <- new.env()

for (i in 1:length(filesglasso)){
  temp.space <- new.env()
  variablepos <- get(load(filesglasso[i], temp.space),envir = temp.space)
  names[i] <- ls(envir = temp.space)
  glassolistCOVRMS[[i]] <- variablepos
  rm(temp.space)
}
names(glassolistCOVRMS) <- names
####################################################################################
# glasso COVRMS op data gewoon / RMS (TRUE/FALSE)!
#####################################################################################
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)
M <- colMeans(data)
logliksGL_COVRMS <- c()
nedgesGL_COVRMS <- c()

for (i in 1:length(glassolistCOVRMS)){
  littles <- dmvnorm(data, mean = M, sigma = glassolistCOVRMS[[i]]$w, log = TRUE)
  logliksGL_COVRMS[i] <- sum(littles)
}

logliksGL_COVRMS <- sort(logliksGL_COVRMS, decreasing = TRUE)
plot(rho, logliksGL_COVRMS)

#################
# glasso COVRMS op dataRMS correlation.
#################
glassolistCOVRMS_W <- lapply(glassolistCOVRMS, function (x){x$w})
glassolistCOVRMS_cor <- lapply(glassolistCOVRMS_W, cov2cor)

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
M <- colMeans(data)
logliksGL_COVRMS <- c()
nedgesGL_COVRMS <- c()

for (i in 1:length(glassolistCOVRMS_cor)){
  littles <- dmvnorm(data, mean = M, sigma = glassolistCOVRMS_cor[[i]], log = TRUE)
  logliksGL_COVRMS[i] <- sum(littles)
}

logliksGL_COVRMS <- sort(logliksGL_COVRMS, decreasing = TRUE)
plot(rho, logliksGL_COVRMS)

for (i in 1:length(glassolistCOVRMS)){
  m <- glassolistCOVRMS[[i]]$wi
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) != 0] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(precgraph))
  nedgesGL_COVRMS[i] <- numberofedges
}
plot(nedgesGL_COVRMS, rho)
plot(nedgesGL_COVRMS,logliksGL_COVRMS, col = "green", ylim = c(-200000000,100000000))
points(nedgesGL_COVRMS,logliksGL_COVRMS, col = "green")
###################################################################################
# glasso different rhos without tresholds CORRELATION
# Loglik of glasso with rho (= logliksGL)
# nedges of correlation matrix (w) without treshold
###################################################################################
rho <- c(0.01,0.05,0.063,0.07,0.09,0.12,0.13,0.165,0.17,0.20,0.225,0.25)

nedgesGL_C <- c()

for (i in 1:length(glassolist)){
  m <- glassolist[[i]]$w
  adj.matrix <- m
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) != 0] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(precgraph))
  nedgesGL_C[i] <- numberofedges
}
plot(rho2, nedgesGL_C)
plot(nedgesGL_C,logliksGL, col = "darkgreen")
points(nedgesGL_C,logliksGL, col = "darkgreen")
648*648
# RESULT:
# all covariancematrices have equal number of zeros -> equal number in networks

# k <- glassolist[[6]]$w
# m <- glassolist[[1]]$w
# length(m[m != 0])
# length(k[k != 0])
###################################################################################
# glasso different rhos with tresholds CORRELATION NOT YET
###################################################################################
# rho <- c(0.01,0.05,0.063,0.07,0.09,0.12,0.13,0.165,0.17,0.20,0.225,0.25)
# 
# logliksGL_C <- c()
# nedgesGL_C <- c()
# 
# for (i in 1:length(glassolist)){
#   littles <- dmvnorm(data, mean = M, sigma = glassolist[[i]]$w, log = TRUE)
#   logliksGL_C[i] <- sum(littles)
# }
# 
# logliksGL_C <- sort(logliksG, decreasing = TRUE)
# plot(rho, logliksGL_C)
# 
# for (i in 1:length(glassolist)){
#   m <- glassolist[[i]]$w
#   adj.matrix <- m
#   diag(adj.matrix) <- 0
#   adj.matrix[abs(adj.matrix) != 0] <- 1
#   precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
#   numberofedges <- length(E(precgraph))
#   nedgesGL_C[i] <- numberofedges
# }
# plot(rho, nedgesGL_C)
# plot(nedgesGL_C,logliksGL_C, col = "darkgreen")

###################################################################################
# glasso 1 rho (0.01) with different treshold PRECISION
###################################################################################
# load glasso 1
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
# set different tresholds.
th <-  c(0.11,0.067,0.055,0.0353)

nedgesGL1_T <- c()
for (i in 1:length(th)){
  m <- glasso01$wi
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) <= th[i] ] <- 0
  adj.matrix[abs(adj.matrix) > th[i] ] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(precgraph))
  nedgesGL1_T[i] <- numberofedges
  
}
plot(th, nedgesGL1_T)

# hoort bij initial loglikelihood van glass01:
littles <- dmvnorm(data, mean = M, sigma = glasso01$w, log = TRUE)
logliksGL001 <- sum(littles)
logliksGL001

# Hier moet komen: logliks van glasso inverse. 
for (i in 1:length(th)){
  m <- glasso01$wi
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) <= th[i] ] <- 0
  adj.matrix[abs(adj.matrix) > th[i] ] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(precgraph))
  nedgesGL1_T[i] <- numberofedges
  
  
  sigma <- m
  sigma[abs(sigma) <= th[1] ] <- 0
  diag(sigma) <- diag(m)
  inv <- solve(sigma)
  sigma*inv
  # NEED: inverse of glasso!!!!! with GLASSO
  # sigma <- m
  # sigma[abs(sigma) <= th[i] ] <- 0
  # diag(sigma) <- diag(m)
  # glasso(sigma, rho = 0, penalize.diagonal = FALSE)
  nnsigma <- make.positive.definite(inv)
  density <- mvtnorm::dmvnorm(dataRMS, sigma = nnsigma,log = TRUE)
  logliksGL_CORR_T[i] <- sum(density)
  
}
plot(nedgesGL1_T,logliksGL_T, col = "orange")
plot(th, nedgesGL1_T)
###################################################################################
# glasso 1 rho (0.01) with different treshold CORRELATION
###################################################################################
# load glasso 1
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
# set different tresholds.
tau <- c(0.87,
         0.8,
         0.7,
         0.62,
         0.52,
         0.49,
         0.41,
         0.345,
         0.3,
         0.2,
         0.1,
         0.05,
         0)

nedgesGL_CORR_T <- c()
logliksGL_CORR_T <- c()

i <- 1
for (i in 1:length(tau)){
  m <- glasso01$w
  adj.matrix <- m
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) <= tau[i] ] <- 0
  adj.matrix[abs(adj.matrix) > tau[i] ] <- 1
  corgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(corgraph))
  nedgesGL_CORR_T[i] <- numberofedges
  
  sigma <- m
  sigma[abs(sigma) <= tau[i] ] <- 0
  diag(sigma) <- diag(m)
  is.positive.definite(sigma)
  nnsigma <- make.positive.definite(sigma
                                    #, tol = 1e-3
                                    )
  is.positive.definite(nnsigma)

  # nnsigma <- glasso(sigma, rho = 0.1)
  density <- mvtnorm::dmvnorm(dataRMS, sigma = nnsigma,log = TRUE)
  logliksGL_CORR_T[i] <- sum(density)
}



library(LaplacesDemon)
M <- colMeans(data)
x <- dmvnp(c(1,2,3), c(0,1,2), diag(3))

density <- mvtnorm::dmvnorm(dataRMS, sigma = glasso01$w,log = TRUE)
sum(density)
density <- mvtnorm::dmvnorm(dataRMS, sigma = sigma,log = TRUE)
sum(density)

plot(nedgesGL_CORR_T,logliksGL_CORR_T, col = "orange")
plot(tau, nedgesGL_CORR_T)
plot(tau, logliksGL_CORR_T, ylim = c(-300000 ,-100000))


###################################################################################
# glasso 1 rho (0.01) with different treshold CORRELATION Usage of Mvtnorm with precision.
###################################################################################
# load glasso 1
library(LaplacesDemon)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
M <- colMeans(data)
# set different tresholds.
tau <- c(0.87,
         0.8,
         0.7,
         0.62,
         0.52,
         0.49,
         0.41,
         0.345,
         0.3,
         0.2,
         0.1,
         0.05,
         0)

nedgesGL_CORR_T <- c()
logliksGL_CORR_T <- c()

i <- 1
for (i in 1:length(tau)){
  m <- glasso01$w
  adj.matrix <- m
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) <= tau[i] ] <- 0
  adj.matrix[abs(adj.matrix) > tau[i] ] <- 1
  corgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  numberofedges <- length(E(corgraph))
  nedgesGL_CORR_T[i] <- numberofedges
  
  sigma <- m
  sigma[abs(sigma) <= tau[i] ] <- 0
  diag(sigma) <- diag(m)
  is.positive.definite(sigma)
  nnsigma <- make.positive.definite(sigma
                                    #, tol = 1e-3
  )
  is.positive.definite(glasso01$wi)
  
  # nnsigma <- glasso(sigma, rho = 0.1)
  density <- dmvnp(dataRMS,mu = M, Omega = glasso01$wi)
  logliksGL_CORR_T[i] <- sum(density)
}



density <- mvtnorm::dmvnorm(dataRMS, sigma = glasso01$w,log = TRUE)
sum(density)
density <- mvtnorm::dmvnorm(dataRMS, sigma = sigma,log = TRUE)
sum(density)

plot(nedgesGL_CORR_T,logliksGL_CORR_T, col = "orange")
plot(tau, nedgesGL_CORR_T)
plot(tau, logliksGL_CORR_T, ylim = c(-300000 ,-100000))


###################################################################################
# glasso different rhos with different tresholds
###################################################################################
# load glasso 1
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso01.rda")
# set different tresholds.
rho <- c(0.01,0.05,0.063,0.07,0.09,0.12,0.13,0.165,0.17,0.20,0.225,0.25)
th <-  c(0.15, 0.11,0.067,0.055,0.0353)

# utility function
f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
#
nedgesG <- matrix(data = NA, nrow = length(rho), ncol = length(th))
nedgesG
j <- 1

for(j in 1:length(rho)){
  for (i in 1:length(th)){
    m <- glassolist[[j]]$wi
    PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
    adj.matrix <- PartVar
    diag(adj.matrix) <- 0
    adj.matrix[abs(adj.matrix) <= th[i] ] <- 0
    adj.matrix[abs(adj.matrix) > th[i] ] <- 1
    precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
    numberofedges <- length(E(precgraph))
    nedgesG[j,i] <- numberofedges
  
  }
}

plot(th, nedgesG[1,])

for (i in 2:length(rho)){
     lines(th,nedgesG[i,], col = rainbow(length(rho))[i])
  }

nedgesG


##################################################################################
# Complex Network treshold
##################################################################################
tau <- c(0.87,
  0.8,
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
  0.06,
  0.05,
  0.0375,
  0.025,
  0)

tau <- seq(0.00,0.99,0.01)
tau <- c(tau,0.026)
# tau <- c(0.87,
#          0.8,
#          0.75,
#          0.7,
#          0.66,
#          0.62,
#          0.52,
#          0.49,
#          0.41,
#          0.345,
#          0.32,
#          0.3,
#          0.25,
#          0.2,
#          0.1,
#          0.15,
#          0.05,
#          0)
length(tau)
tauchar <- as.character(tau)
labelsCM <- c()

for(i in 1:length(tauchar)){
  labelsCM[i] <- paste0("= ",tauchar[i])
}

labelsCM
graphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, method = "spearman",subind = NULL)

# measures10d <- lapply(graphs10d, graph2measure)

graphssolo  <- lapply(graphs10d, function(m) m$graph)
graphssolo
nedgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(nedgesnetworksCM, length)
nedgesnetworksCM

corssolo <- lapply(graphs10d, function(m) m$correlation)
all.equal(corssolo[[1]],corssolo[[2]])
corcutsolo <- list()
for (i in 1:length(corssolo)){
  corcutsolo[[i]] <- corssolo[[i]]
  corcutsolo[[i]][abs(corcutsolo[[i]]) <= tau[i] ] <- 0
}

logliksCM <- c()
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cormatrix <- cor(data, method = "spearman")
for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolo[[i]], tol = 0.06)
  density <- mvtnorm::dmvnorm(dataRMS, sigma = nncor,log = TRUE)
  logliksCM[i] <- sum(density)
}
logliksCM

################################################################################
# Combine plots loglik vs number of edges networks
################################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Loglikelihood/logliks_CNall.pdf")
pdf(plotname)
dev.off()
# Nearest positive definite correlation matrix, different tresholds
plot(nedgesnetworksCM, logliksCM, xlim = c(0,210000),ylim = c(-360000,-100000), xlab = "number of edges", ylab = bquote("log P(X|G,"~Theta~")"))
for (i in 1:length(nedgesnetworksCM)){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
}
abline(h = logliksCM[1], lty = 2)
# Glasso obtained correlation matrix, different rhos, without tresholds
# EDGES OF PRECISION NETWORK
points(nedgesGL,logliksGL, col = "green")
# Bayesian network: HC iteration obtained parameters
points(nedgesIT,logliksIT, col = "blue")
# Bayesian network: eBIC obtained parameters
points(nedgeseBIC,loglikseBIC, col = "red")
# markov network: by bayesian network
points(nedgesmorals, logliksIT, col = "purple")
legend("bottomright", c("BN","CN"),bty = "n", pch = 1, col =c("blue","black", cex = 0.6))
legend("bottomright", c("BN: HC BIC_0", "BN: HC eBIC", "CN", "MN: precision", "MN: moralBN"),bty = "n", pch = 1, col =c("blue","red","black","green","purple"), cex = 0.6)
# Glasso estimated correlation matrix 
points(nedgesGL_C,logliksGL, col = "darkgreen")
# points(nedgesGL_CORR_T,logliksGL_CORR_T, col = "orange")
# #
# points(nedgesGL_C,logliksGL, col = "darkgreen")
# #
# points(nedgesGL_C,logliksGL_COVRMS, col = "pink")

# Nearest positive definite correlation matrix, different tresholds
plot(nedgesnetworksCM, logliksCM, 
     xlim = c(0,8000),ylim = c(-360000,-100000), 
     xlab = "number of edges", ylab = bquote("log P(X|G,"~Theta~")"))
platpoints <- c(which(logliksCM == maxplat),which(logliksCM == lastplat))
for (i in 1:length(nedgesnetworksCM)){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
}
label0 <- bquote(tau ~ "=" ~ .(0))
abline(h = logliksCM[length(logliksCM)], lty = 2)
text(8000,-130000, labels = label0, cex = 0.6, col = "blue")
text(8000,-140000, labels = round(logliksCM[length(logliksCM)]), cex = 0.6)
# Glasso obtained correlation matrix, different rhos, without tresholds
# EDGES OF PRECISION NETWORK
points(nedgesGL,logliksGL, col = "green")
# for (i in 1:length(nedgesGL)){
#   value <- rhochar[i]
#   label <- bquote(tau ~ "=" ~ .(value))
#   text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
# }
# iteration obtained parameters
points(nedgesIT,logliksIT, col = "blue")
# eBIC obtained parameters
points(nedgeseBIC,loglikseBIC, col = "red")
# morals pink
points(nedgesmorals, logliksIT, col = "purple")

legend("topleft", c("BN: HC BIC_0", "BN: HC eBIC", "CN", "MN: precision", "MN: moralBN"),bty = "n", pch = 1, col =c("blue","red","black","green","purple"), cex = 0.6)
#points(nedgesGL_CORR_T,logliksGL_CORR_T, col = "orange")

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Loglikelihood/logliks_8500.pdf")
pdf(plotname)
dev.off()
plot(nedgesnetworksCM, logliksCM, 
     xlim = c(0,8500),ylim = c(-360000,-100000), 
     xlab = "number of edges", ylab = bquote("log P(X|G,"~Theta~")"))

firstplat <- logliksCM[53:100][which(logliksCM[53:100] == max(logliksCM[53:100]))]
nedgesnetworksCM[53:100][which(logliksCM[53:100] == max(logliksCM[53:100]))]
which(logliksCM[58:31] == max(logliksCM[58:31]))
logliksCM[58:31]
maxplat <- logliksCM[58:31][which(logliksCM[58:31] == max(logliksCM[58:31]))]
lastplat <- logliksCM[30:1][10]
platpoints <- c(which(logliksCM == maxplat),which(logliksCM == firstplat),which(logliksCM == lastplat))

for (i in platpoints){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
}
label0 <- bquote(tau ~ "=" ~ .(0))
abline(h = logliksCM[1], lty = 2)
text(8000,-130000, labels = label0, cex = 0.6, col = "blue")
text(8000,-140000, labels = round(logliksCM[1]), cex = 0.6)

points(nedgesIT,logliksIT, col = "blue")

legend("topleft", c("BN: HC BIC_0", "CN"),bty = "n", pch = 1, col =c("blue","black"), cex = 0.6)
#points(nedgesGL_CORR_T,logliksGL_CORR_T, col = "orange")


#############################################################################
# Event logliks 
#############################################################################
logliksCMevent <- array(data = NA, dim = c(nrow(dataRMS),length(tau)))
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
nedgesnetworksCM[35]
cormatrix <- cor(data, method = "spearman")

for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolo[[i]], tol = 0.06)
  density <- mvtnorm::dmvnorm(dataRMS, sigma = nncor,log = TRUE)
  logliksCMevent[,i] <- density
}

dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
dataRMSframes <- apply(dataRMS, MARGIN = 1, FUN = function(x) as.data.frame(x))
dataRMSframes <- lapply(dataRMSframes, t)
dataRMSframes <- lapply(dataRMSframes, as.data.frame)


logliksITevent <- array(data = NA, dim = c(nrow(dataRMS),length(iterationfits)))
logliksITevent

for (i in 1:length(iterationfits)){
  logliksITevent[,i] <- sapply(dataRMSframes, FUN = logLik, object = iterationfits[[i]])
}

logliksITevent[,1]
nedgesIT <- sapply(iterationnetworks, narcs)
nedgesnetworksCM[4]
data.frame(nedges = nedgesIT, logliks = logliksIT)

bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(3)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/Loglik_datarealizations.pdf")
pdf(plotname)

plot(logliksCMevent[,1],col = blackcn[3],ylim = c(-1000,500), xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), main = "Loglik data realizations given models")
points(logliksITevent[,9], col = bluebn[1], pch = 2)
points(logliksCMevent[,35],col = blackcn[2],pch = 3)
points(logliksCMevent[,56],col = blackcn[1], pch = 4)
points(logliksITevent[,9]-logliksCMevent[,1], col = "green")
legend("topleft", c(paste0("BN ",nedgesIT[9]), paste0("CN ",nedgesnetworksCM[1]), paste0("CN ",nedgesnetworksCM[35]), paste0("CN ",nedgesnetworksCM[56]), paste0("BN ",nedgesIT[9],"- CN ",nedgesnetworksCM[1])),bty = "n", pch = c(2,1,3,4,1), col =c(bluebn[1],blackcn[3],blackcn[2],blackcn[1], "green"), cex = 0.6)

dev.off()

orderCNs <- matrix(data = NA, nrow = nrow(dataRMS), ncol = length(tau))
for(i in 1:length(tau)){
  orderCNs[,i] <- order(logliksCMevent[,i])
}
orderBNs <- matrix(data = NA, nrow = nrow(dataRMS), ncol = length(iterationfits))
for(i in 1:length(iterationfits)){
  orderBNs[,i] <- order(logliksITevent[,i])
}

orderCNs[,1]
# all.equal(orderBNs[,18],orderBN)
# orderCM <- order(logliksCMevent[,35])
# orderBN <- order(logliksITevent[,18])
 orderdif <- order(logliksITevent[,9]-logliksCMevent[,101])
 (logliksITevent[,9]-logliksCMevent[,1])[orderdif]
##########################################################################
# Evidenc propagation correlation network
##########################################################################
# All networks.
# ONE evidence: V81
#########################################################################
evolutionCN <- list()
i <- 1
for (i in (1:length(tau))[seq(1,100,10)]){
  possingleCNtau <- propagationCorr(corcutsolo[[i]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(459,280), valueEvidence = c(2,2)) 
  evolutionCN[[i]] <- possingleCNtau
  }

nedgesplot <- nedgesnetworksCM
logliksplot <- logliksCM
props <- evolutionCN[seq(1,100,10)]
Center1 <- 180
plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center1, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )),rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center1, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center1, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)),rev.colors = TRUE)
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper2/figures/EvolucionV81CN.pdf")
pdf(plotname, height = 7, width = 7)
do.call("grid.arrange", c(all, ncol = 4, nrow = 5))
dev.off()

##########################################################################
# All networks.
# Double evidence: V81 V280
#########################################################################
evolutionCN <- list()
i <- 1
for (i in (1:length(tau))[seq(1,100,10)]){
  possingleCNtau <- propagationCorr(corcutsolo[[i]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2)) 
  evolutionCN[[i]] <- possingleCNtau
}

nedgesplot <- nedgesnetworksCM
logliksplot <- logliksCM
props <- evolutionCN[seq(1,100,10)]
Center1 <- 180
plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center1, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )),rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center1, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center1, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)),rev.colors = TRUE)
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper2/figures/EvolucionV81CN.pdf")
pdf(plotname, height = 7, width = 7)
do.call("grid.arrange", c(all, ncol = 4, nrow = 5))
dev.off()




###############################################################################
# Propagation part Integrate functions: Positive function: DOUBLE EVIDENCE 
# Bayesian network
###############################################################################
# Load HC propagation iterations
pattern <- paste0("prop_")
fileshciteration <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration", full.names = T, pattern = pattern)
fileshciteration
hcproplist <- list()
names <- c()
temp.space <- new.env()

for (i in 1:length(fileshciteration)){
  temp.space <- new.env()
  variablepos <- get(load(fileshciteration[i], temp.space),envir = temp.space)
  names[i] <- ls(envir = temp.space)
  hcproplist[[i]] <- variablepos
  rm(temp.space)
}
names(hcproplist) <- names

x <- seq(200,3000,200)
y <- seq(400,3200,200)

hcproplistord <- list()
namesord <- c()
for(i in 1:length(x)){
  hcproplistord[[i]] <- eval(parse(text =paste0("hcproplist$prop_hc_",x[i],"_",y[i],"i_V81V280_equal22")))
  namesord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81V280_equal22")
    }
names(hcproplistord) <- namesord

nedgesord <- sapply(iterationnetworks[2:length(iterationnetworks)], narcs) # Begins with 400 !
logliksord <- logliksIT[2:length(logliksIT)]
props <- hcproplistord

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center2 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center2, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )),rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center2, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05),rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center2, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n |E| = ",nedgesord[i]," log = ",round(logliksord[i])),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)),rev.colors = TRUE)
  
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, ncol = 4, nrow = 4))

###############################################################################
# Propagation part Integrate functions: NEGATIVE FUNCTION: Double evidence
# Bayesian network
###############################################################################
# Load HC propagation iterations
pattern <- paste0("propneg_")
fileshciterationneg <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V280neg", full.names = T, pattern = pattern)
fileshciterationneg
hcpropneglist <- list()
namesneg <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationneg)){
  temp.space <- new.env()
  variableneg <- get(load(fileshciterationneg[i], temp.space),envir = temp.space)
  namesneg[i] <- ls(envir = temp.space)
  hcpropneglist[[i]] <- variableneg
  rm(temp.space)
}
names(hcpropneglist) <- namesneg

x <- seq(200,2800,200)
y <- seq(400,3000,200)

hcpropneglistord <- list()
namesnegord <- c()
for(i in 1:length(x)){
  hcpropneglistord[[i]] <- eval(parse(text =paste0("hcpropneglist$propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22")))
  namesnegord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22")
}
names(hcpropneglistord) <- namesnegord

nedgesord <- sapply(iterationnetworks, narcs)
nedgesord

props <- hcpropneglistord
Center3 <- 180
plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center3, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center3, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center3, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n |E| = ",nedgesord[i]," log = ",round(logliksIT[i])),cex = 0.5), at = seq(0.05,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  b
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, ncol = 4, nrow = 4))

###############################################################################
# Propagation part Integrate functions: positive FUNCTION: single evidence
# Bayesian network
###############################################################################
# Load HC propagation iterations
pattern <- paste0("prop_")
fileshciterationpos1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos", full.names = T, pattern = pattern)
fileshciterationpos1
hcproppos1list <- list()
namespos1 <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationpos1 )){
  temp.space <- new.env()
  variablepos1 <- get(load(fileshciterationpos1 [i], temp.space),envir = temp.space)
  namespos1[i] <- ls(envir = temp.space)
  hcproppos1list[[i]] <- variablepos1
  rm(temp.space)
}
names(hcproppos1list) <- namespos1

x <- seq(200,3000,200)
y <- seq(400,3200,200)

hcproppos1listord <- list()
namespos1ord <- c()
for(i in 1:length(x)){
  hcproppos1listord[[i]] <- eval(parse(text =paste0("hcproppos1list$prop_hc_",x[i],"_",y[i],"i_V81_equal2")))
  namespos1ord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2")
}
names(hcproppos1listord) <- namespos1ord

nedgesord <- sapply(iterationnetworks, narcs)
nedgesord




props <- hcproppos1listord
 Center4 <- 180


plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center4, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )),rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center4, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05),rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center4, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n |E| = ",nedgesord[i]," log = ",round(logliksIT[i])),cex = 0.5), at = seq(0.01,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)),rev.colors = TRUE)
  b
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, ncol = 4, nrow = 4))

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper2/figures/EvolucionV81BN.pdf")
pdf(plotname, height = 7, width = 7)
do.call("grid.arrange", c(all, ncol = 4, nrow = 4))
dev.off()
###############################################################################
# Propagation part Integrate functions: negative FUNCTION: single evidence
# Bayesian network
###############################################################################
# Load HC propagation iterations
pattern <- paste0("propneg_")
fileshciterationneg1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81neg", full.names = T, pattern = pattern)
fileshciterationneg1
hcpropneg1list <- list()
namesneg1 <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationneg1 )){
  temp.space <- new.env()
  variableneg1 <- get(load(fileshciterationneg1 [i], temp.space),envir = temp.space)
  namesneg1[i] <- ls(envir = temp.space)
  hcpropneg1list[[i]] <- variableneg1
  rm(temp.space)
}
names(hcpropneg1list) <- namesneg1

x <- seq(200,3000,200)
y <- seq(400,3200,200)

hcpropneg1listord <- list()
namesneg1ord <- c()
for(i in 1:length(x)){
  hcpropneg1listord[[i]] <- eval(parse(text =paste0("hcpropneg1list$propneg_hc_",x[i],"_",y[i],"i_V81_equal2")))
  namesneg1ord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2")
}
names(hcpropneg1listord) <- namesneg1ord

nedgesord <- sapply(iterationnetworks, narcs)
nedgesord


props <- hcpropneg1listord
Center5 <- 180
plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center5, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center5, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center5, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n |E| = ",nedgesord[i]," log = ",round(logliksIT[i])),cex = 0.5), at = seq(0.01,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  b
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, ncol = 4, nrow = 4))
##############################################################################
# Compare from Correlation network FULL and Bayesian network equivalent.
##############################################################################
numberCM <- 1
numberBN <- 12

posdoubleCNFull <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
posdoubleBNFull <- hcproplistord[[numberBN]]
negdoubleCNFull <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
negdoubleBNFull <- hcpropneglistord[[numberBN]]
possingleCNFull <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
possingleBNFull <- hcproppos1listord[[numberBN]]
negsingleCNFull <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negsingleBNFull <- hcpropneg1listord[[numberBN]]

attr(posdoubleCNFull, "sign") <- attr(posdoubleBNFull, "sign") <- attr(possingleCNFull, "sign") <-  attr(possingleBNFull, "sign") <- "pos"
attr(negdoubleCNFull, "sign") <- attr(negdoubleBNFull, "sign") <- attr(negsingleCNFull, "sign") <-  attr(negsingleBNFull, "sign") <- "neg"


props <- list(possingleBNFull, possingleCNFull,
              negsingleBNFull,negsingleCNFull,
              posdoubleBNFull,posdoubleCNFull,
              negdoubleBNFull,negdoubleCNFull)

logCM <- logliksCM[numberCM]
logBN <- logliksIT[[numberBN+1]]
logliksplot <- c(rep(c(logBN,logCM),4))
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBN+1]
nedgesplot <- c(rep(c(nBM,nCM),4))
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center6 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline",lonCenter = Center6, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center6, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center6, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.00,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
all

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper2/figures/probabilitiesCompFull.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/probabilitiesCNFullBNeq.pdf")
pdf(plotname, height = 7, width = 5)
do.call("grid.arrange", c(all, ncol=2))
dev.off()
##############################################################################
# Compare from CM and BN 1800
##############################################################################
numberCM <- length(corcutsolo)-1
numberBN <- 9

posdoubleCN14 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
posdoubleBN <- hcproplistord[[numberBN]]
negdoubleCN14 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
negdoubleBN <- hcpropneglistord[[numberBN]]
possingleCN14 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
possingleBN <- hcproppos1listord[[numberBN]]
negsingleCN14 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negsingleBN <- hcpropneg1listord[[numberBN]]

props <- list(possingleBN, possingleCN14,negsingleBN,negsingleCN14,posdoubleBN,posdoubleCN14,negdoubleBN,negdoubleCN14)
logCM <- logliksCM[numberCM]
logBN <- logliksIT[[numberBN+1]]
logliksplot <- c(rep(c(logBN,logCM),4))
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBN+1]
nedgesplot <- c(rep(c(nBM,nCM),4))
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
all

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper2/figures/probabilitiesComp1800.pdf")
pdf(plotname, height = 7, width = 5)
do.call("grid.arrange", c(all, ncol=2))
dev.off()

##############################################################################
# Compare from CM and BN FULL double evidence
##############################################################################
numberCM <- length(corcutsolo)
numberBN <- 13

test <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
test2 <- hcproplistord[[numberBN]]

props <- list(test,test2)
logliksplot <- c(logliksCM[numberCM],logliksIT[[numberBN]])
nedgesplot <- c(nedgesnetworksCM[numberCM] ,nedgesIT[numberBN])
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center8 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center8, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center8, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center8, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all))
##############################################################################
# Compare from CM and BN FULL Single evidence
##############################################################################
# Load HC propagation iterations
pattern <- paste0("prop_")
fileshciterationpos1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos", full.names = T, pattern = pattern)
fileshciterationpos1
hcproppos1list <- list()
namespos1 <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationpos1 )){
  temp.space <- new.env()
  variablepos1 <- get(load(fileshciterationpos1 [i], temp.space),envir = temp.space)
  namespos1[i] <- ls(envir = temp.space)
  hcproppos1list[[i]] <- variablepos1
  rm(temp.space)
}
names(hcproppos1list) <- namespos1

x <- seq(200,3000,200)
y <- seq(400,3200,200)

hcproppos1listord <- list()
namespos1ord <- c()
for(i in 1:length(x)){
  hcproppos1listord[[i]] <- eval(parse(text =paste0("hcproppos1list$prop_hc_",x[i],"_",y[i],"i_V81_equal2")))
  namespos1ord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2")
}
names(hcproppos1listord) <- namespos1ord

nedgesord <- sapply(iterationnetworks, narcs)
nedgesord






numberCM <- length(corcutsolo)-5
numberBN <- 6

test <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
test <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227,187), valueEvidence = c(2,-1))
test2 <- hcproppos1listord[[numberBN]]

props <- list(test,test2)
logliksplot <- c(logliksCM[numberCM],logliksIT[[numberBN]])
nedgesplot <- c(nedgesnetworksCM[numberCM] ,nedgesIT[numberBN])
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center9 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.01,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all))
##############################################################################
# Compare from CM GLASSO and BN
##############################################################################
g <- glassolist$glasso09$w
numberBN <- 13
numberGL <- 1


posdoubleCN <- propagationCorr(g, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
posdoubleBN <- hcproplistord[[numberBN]]
negdoubleCN <- propagationnegCorr(g, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
negdoubleBN <- hcpropneglistord[[numberBN]]
possingleCN <- propagationCorr(g, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
possingleBN <- hcproppos1listord[[numberBN]]
negsingleCN <- propagationnegCorr(g, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negsingleBN <- hcpropneg1listord[[numberBN]]

props <- list(posdoubleCN,possingleCN,negdoubleCN,negsingleCN,posdoubleBN,possingleBN,negdoubleBN,negsingleBN)
logCM <- logliksGL[numberGL]
logBN <- logliksIT[[numberBN]]
logliksplot <- c(rep(logCM,4),rep(logBN,4))
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBN]
nedgesplot <- c(rep(nCM,4), rep(nBM,4))
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

Center10 <- 180
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center10, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center10, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center10, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all,nrow = 2))

##############################################################################
# Glasso COV RMS COVARIANCES ARE VERY BIG: 
##############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/covRMS/glassoCOVrms_0.01.rda")
g <- glassoCOVrms_0.01$w
gcor <- cov2cor(g)
all.equal(gcor,glasso01$w)
glasso2pos <- propagationCorr(g, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
glasso2neg <- propagationnegCorr(g, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
glasso1pos <- propagationCorr(g, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
glasso1neg <- propagationnegCorr(g, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))

glasso2posc <- propagationCorr(gcor, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
glasso2negc <- propagationnegCorr(gcor, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
glasso1posc <- propagationCorr(gcor, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
glasso1negc <- propagationnegCorr(gcor, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))

glasso2posc <- propagationCorr(gcor, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
glasso2negc <- propagationnegCorr(gcor, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
glasso1posc <- propagationCorr(gcor, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
glasso1negc <- propagationnegCorr(gcor, nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))



props <- list(glasso2pos,glasso2neg,glasso1pos,glasso1neg)
props <- list(glasso2posc,glasso2negc,glasso1posc,glasso1negc)
logCM <- logliksGL[numberGL]
logBN <- logliksIT[[numberBN]]
logliksplot <- c(rep(logCM,4),rep(logBN,4))
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBN]
nedgesplot <- c(rep(nCM,4), rep(nBM,4))
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center11 <- NULL

for(i in 1:length(props)){
  prop <- props[[i]]
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center11, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center11, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center11, main = list(paste0(attr(prop_dif$Data,"climatology:fun")),cex = 0.5), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all,nrow = 2))

coef(iterationfits[[1]]$V81)
names(coef(iterationfits[[1]]$V63))
coef(iterationfits[[1]])
