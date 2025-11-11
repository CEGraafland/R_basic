######################################################################################
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
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")


f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
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
plot(nedgesIT,logliksIT, col = "blue", xlim = c(0,8000))
########################################################################################
# to sparseBN package
########################################################################################
NELS <- lapply(iterationfits,as.graphNEL)
edgelistssparse <- lapply(NELS,as.edgeList)
sparsedata <- sparsebnData(dataRMS, type = "continuous")
COVS <- lapply(edgelistssparse, get.covariance, data = sparsedata)
mats <- lapply(COVS, as.matrix)

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
M <- colMeans(data)
logliksSP <- c()
nedgesSP<- c()

for (i in 1:length(NELS)){
  littles <- dmvnorm(data, mean = M, sigma = mats[[i]], log = TRUE)
  logliksSP[i] <- sum(littles)
}


testfit <- iterationfits[[1]]


# 
# NEL <- as.graphNEL(testfit)
# edgelist <- as.edgeList(NEL)
# sparsedata <- sparsebnData(dataRMS, type = "continuous")
# cov <- get.covariance(edgelist,sparsedata)
# mat <- as.matrix(cov)
# edgelist[3]



plot(nedgesIT,logliksIT, col = "blue", xlim = c(0,8000))
points(nedgesIT,logliksSP)
all.equal(logliksSP,logliksIT)
#########################################################################
# plogp 
##########################################################################
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
M <- colMeans(data)
plogpsSP <- c()
nedgesSP<- c()

for (i in 1:length(NELS)){
  littles <- dmvnorm(data, mean = M, sigma = mats[[i]], log = FALSE)*dmvnorm(data, mean = M, sigma = mats[[i]], log = TRUE)
  plogpsSP[i] <- -sum(littles)
}

plot(nedgesIT,log(plogpsSP), ylim = c(-1000,0))
points(nedgesnetworksCM,log(plogpsCM))



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

plogpsCM <- c()
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cormatrix <- cor(data, method = "spearman")
for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolo[[i]], tol = 0.06)
  plogp <- mvtnorm::dmvnorm(dataRMS, sigma = nncor,log = FALSE)*mvtnorm::dmvnorm(dataRMS, sigma = nncor,log = TRUE)
  plogpsCM[i] <- -sum(plogp)
}
plogpsCM 

################################################################################
# Combine plots loglik vs number of edges networks
################################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Loglikelihood/plogp_CNBNbig.pdf")
pdf(plotname)
dev.off()
# Nearest positive definite correlation matrix, different tresholds
plot(nedgesnetworksCM, log(plogpsCM), xlim = c(0,210000), ylim = c(-800,0), xlab = "number of edges", ylab = bquote("log (-P(X|G,"~Theta~") log P(X|G,"~Theta~"))"))

points(nedgesIT,log(plogpsSP), col = "blue")

legend("bottomright", c("BN","CN"),bty = "n", pch = 1, col =c("blue","black", cex = 0.6))


# Nearest positive definite correlation matrix, different tresholds
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Loglikelihood/plogp_CNBNsmall.pdf")
pdf(plotname)
plot(nedgesnetworksCM, log(plogpsCM), 
     xlim = c(0,8000),ylim = c(-800,-200), 
     xlab = "number of edges", ylab = bquote("log (-P(X|G,"~Theta~") log P(X|G,"~Theta~"))"))
# platpoints <- c(which(logliksCM == maxplat),which(logliksCM == lastplat))
# for (i in 1:length(nedgesnetworksCM)){
#   value <- tauchar[i]
#   label <- bquote(tau ~ "=" ~ .(value))
#   text(nedgesnetworksCM[i],logliksCM[i],label= label,col='blue', cex = 0.6, pos = 3)
# }
# label0 <- bquote(tau ~ "=" ~ .(0))
# abline(h = logliksCM[length(logliksCM)], lty = 2)
# text(8000,-130000, labels = label0, cex = 0.6, col = "blue")
# text(8000,-140000, labels = round(logliksCM[length(logliksCM)]), cex = 0.6)
# Glasso obtained correlation matrix, different rhos, without tresholds

points(nedgesIT,log(plogpsSP), col = "blue")
legend("bottomright", c("BN","CN"),bty = "n", pch = 1, col =c("blue","black", cex = 0.6))
dev.off()
