# install.packages("sparsebn")
library(sparsebn)
library(mvtnorm)
library(igraph)
library(bnlearn)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
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


# Convert the graphs to undirected igraphs
igraphs1 <- lapply(mats, graph_from_adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
edgelists1 <- lapply(igraphs1, E)
nedgeslists1 <- sapply(edgelists1, length)


# Create the graphObject as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(igraphs1))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- igraphs1[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphs1[[i]])
}

plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid(tas_ncep_10d),graphObjects[[4]]$graph)
plot_long_distances(graphObjects[[2]]$graph, as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d)), minimdist = 10000)
dev.off()


par(mfrow = c(2,2))
test4 <- graph_from_cov(grid = tas_ncep_10d, covar = mats[[4]], th= 0.4)
test16 <- graph_from_cov(grid = tas_ncep_10d, covar = mats[[16]], th= 0.4)
test21 <- graph_from_cov(grid = tas_ncep_10d, covar = mats[[21]], th= 0.4)
test44<- graph_from_cov(grid = tas_ncep_10d, covar = mats[[44]], th= 0.4)
plot_long_distances(test4$graph, as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d)), minimdist = 10000)
plot.Meteodag(data,test4$graph)
plot.Meteodag(data,test16$graph)
plot.Meteodag(data,test21$graph)
plot.Meteodag(data,test44$graph)
length(E(test44$graph))

GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- lapply(edgelistsCN, length)
plot.Meteodag(data, gridGraphsCN[[11]]$graph)
##########################################################################################
# Convert bn permutation 1 hc to covariance matrix
##########################################################################################
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")

dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
hcperm1fits <- lapply(perm1sort, bn.fit, data = dataRMS)
bnarcs <- lapply(perm1sort,narcs)

GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})
edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- lapply(edgelistsBN, length)

NELShc1 <- lapply(hcperm1fits,as.graphNEL)
edgelistssparsehc1 <- lapply(NELShc1,as.edgeList)
sparsedata <- sparsebnData(dataRMS, type = "continuous")
COVShc1 <- lapply(edgelistssparsehc1, get.covariance, data = sparsedata)
mats <- lapply(COVShc1, as.matrix)

gridGraphsBNtoCN <- lapply(mats, graph_from_cov, grid = tas_ncep_10d, th= 0)
gridGraphsBNtoCN
save(gridGraphsBNtoCN, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBNtoCNperm1.rda")



#########################################################################
# plogp (see plogp_data_given_models.R)
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