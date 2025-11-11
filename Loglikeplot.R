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
# source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
# source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
# source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")

source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")


f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
#####################################################################################
# Load data HC eBIC 
#######################################################################################
for(j in c(1,2,3,4,5)){
  pattern <- paste0("hc_",j,"_eBIC_g")
  filestabu1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC", pattern = pattern)
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
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")
#load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda")
load(file ="/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/hc_edges_loglik_10d_200i.rda")
x <- seq(200,8600,200)
y <- seq(400,8800,200)
i <- 200
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

##################################################################################
# Complex Network treshold
##################################################################################
tau <- seq(0,1,0.01)

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
graphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = NULL)
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
# Nearest positive definite correlation matrix, different tresholds
plot(nedgesnetworksCM, logliksCM, xlim = c(0,100000),ylim = c(-360000,-100000), xlab = "number of edges", ylab = bquote("log(P(X|(G,"~Theta~")"))
tautexts <- seq(10,60,10)
for (i in tautexts){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
}
abline(h = logliksCM[1], lty = 2)
# Glasso obtained correlation matrix, different rhos, without tresholds
# EDGES OF PRECISION NETWORK
# points(nedgesGL,logliksGL, col = "green")
# Bayesian network: HC iteration obtained parameters
points(nedgesIT,logliksIT, col = "blue")
# Bayesian network: eBIC obtained parameters
points(nedgeseBIC,loglikseBIC, col = "red")
# markov network: by bayesian network
# points(nedgesmorals, logliksIT, col = "purple")
legend("bottomright", c("BN: HC BIC_0", "BN: HC eBIC", "CN", "MN: precision", "MN: moralBN"),bty = "n", pch = 1, col =c("blue","red","black","green","purple"), cex = 0.6)
legend("bottomright", c("BN: HC BIC_0", "BN: HC eBIC", "CN"),bty = "n", pch = 1, col =c("blue","red","black"), cex = 0.6)
#
# points(nedgesGL_CORR_T,logliksGL_CORR_T, col = "orange")
# #
# points(nedgesGL_C,logliksGL, col = "darkgreen")
# #
# points(nedgesGL_C,logliksGL_COVRMS, col = "pink")

# Nearest positive definite correlation matrix, different tresholds
plot(nedgesnetworksCM, logliksCM, xlim = c(0,8000),ylim = c(-360000,-100000), xlab = "number of edges", ylab = bquote("log(P(X|(G,"~Theta~")"))
for (i in 1:length(nedgesnetworksCM)){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
}
label0 <- bquote(tau ~ "=" ~ .(0))
abline(h = logliksCM[1], lty = 2)
text(8000,-130000, labels = label0, cex = 0.6, col = "blue")
text(8000,-140000, labels = round(logliksCM[1]), cex = 0.6)
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


nedgesnetworksCM[54]
