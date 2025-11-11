#############################################################################
# Precision matrix with Glasso method
#############################################################################
library(huge)
library(glasso)
library(igraph)
library(mopa)
library(igraph)
library(gridExtra)

data(wrld)
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
#############################################################################
# Simple emperical cormatrix = S is used. and than with blockdescent way glasso estimate. 
#############################################################################


# # precisionr05 <- glasso(cormatrix, rho = 0.5)
# # precisionr04 <- glasso(cormatrix, rho = 0.4)
# # precisionr025 <- glasso(cormatrix, rho = 0.25) 677 opnieuw
#  precisionr0225 <- glasso(cormatrix, rho = 0.225) # 847
# # precisionr022 <- glasso(cormatrix, rho = 0.22) 894
# # precisionr020<-  glasso(cormatrix, rho = 0.20) 1122
# # precisionr0175 <- glasso(cormatrix, rho = 0.175) 1451
# # precisionr017 <- glasso(cormatrix, rho = 0.17) 1556
# precisionr0165 <- glasso(cormatrix, rho = 0.165) # 1643
# precisionr013 <- glasso(cormatrix, rho = 0.13) # 2584
# #precisionr011 <- glasso(cormatrix, rho = 0.11) # 3432
# #precisionr01 <- glasso(cormatrix, rho = 0.1) # 3981 
# #precisionr009 <- glasso(cormatrix, rho = 0.09) #4686
# #precisionr0075 <- glasso(cormatrix, rho = 0.075) # 6136
# #precisionr007<- glasso(cormatrix, rho = 0.07) #6765
# #precisionr0065 <- glasso(cormatrix, rho = 0.065) #7476
# precisionr0063 <- glasso(cormatrix, rho = 0.063) #7789
# #precisionr006 <- glasso(cormatrix, rho = 0.06) #8239
# #precisionr005 <- glasso(cormatrix, rho = 0.05, trace = TRUE) # 10204
# # precision <- glasso(cormatrix, rho = 0.01) 36356 
(length(glasso01$wi[glasso01$wi != 0])-648)/2
(length(which(glasso01$wi != 0)-648))/2
648*648

glasso01$wi
# save(precision, precisionr005,precisionr006,precisionr0065,precisionr009, precisionr0075,precisionr007, precisionr01,precisionr011,precisionr013,precisionr0165,precisionr017,precisionr0175, precisionr020,precisionr022, precisionr0225,precisionr025, precisionr04, precisionr05, file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/HC10dprecision.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/HC10dprecision.rda")
##################################################################
# Change names
##################################################################
# glasso09 <- glasso9
# save(glasso09, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/glasso/glasso09.rda")

# which(precision$wi != 0)
# adj.matrix <- precisionr0063$wi
# m <- matrix(data = 1:9, nrow =3, ncol = 3)
# m
#-m[2,3]/sqrt(m[2,2]*m[3,3])
f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
# f(2,3,m)

###############################################################################################
# Compare graph from covvariance matrix and precision matrix. 
###############################################################################################
m <- precisionr05$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
diag(adj.matrix) <- 0
adj.matrix[adj.matrix != 0 ] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
dev.off()
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d),precgraph)

m <- precisionr05$w
CoVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
CoVar
adj.matrix <- CoVar
diag(adj.matrix) <- 0
adj.matrix[abs(adj.matrix) >= 0.05 ] <- 1
adj.matrix[abs(adj.matrix) < 0.05] <- 0
covvargraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(covvargraph)
dev.off()
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d),covvargraph)


################################################################################################
# Figure out teleconnections in precision matrix:
# Makes use of plot.meteodag in Bayesian Graph functions
# 
# TO DO: PLOT LONG DISTANCES OVER GRAPH.
################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/HC10dprecision.rda")
m <- precision$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
diag(adj.matrix) <- 0
adj.matrix[adj.matrix < .09 & adj.matrix > 0.08] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
dev.off()
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d),precgraph)
plot_long_distances(precgraph,TimeCoordsAnom_from_Grid(tas_ncep_10d), 10000)

adj.matrix <- PartVar
################################################################################################
# Plot correlation matrix with respect to one node (All dependencies)
################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/ncep/tas_ncep.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")

library('RColorBrewer')
display.brewer.all()
col.blue <- rev(brewer.pal(9,"Blues"))
col.blue
col.red <- brewer.pal(9,"Reds")
col.b <- c(col.blue,"white","white",col.red)
length(col.b)

plotCorFixedV <- function(node, cormat, refgrid, RMS = TRUE, at = seq(0.05,1,0.05)){
  if (RMS == FALSE){
    dataRMS <- TimeCoordsAnom_from_Grid(refgrid)
    }
  else {dataRMS <- TimeCoordsAnom_from_Grid_rms(refgrid)}
  respcor<- cormat[node,]

  matrespcor <- matrix(respcor, nrow = 1)
  clim_matresp <- refgrid
  clim_matresp$Data <- mat2Dto3Darray(matrespcor, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  attr(clim_matresp$Data, "climatology:fun") <- paste0("corr with respect to V",node)
  plotClimatology(clim_matresp,backdrop.theme = "coastline", main = list(paste0("corr(VX,V",node,")")),at = at, col.regions= col.b)
}

plotParFixedV <- function(node, precmat, refgrid, RMS = TRUE, at = seq(0.05,1,0.05)){
  if (RMS == FALSE){
    dataRMS <- TimeCoordsAnom_from_Grid(refgrid)
  }
  else {dataRMS <- TimeCoordsAnom_from_Grid_rms(refgrid)}
 
  m <- precmat
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  respPart<- PartVar[node,]
  
  matrespPart <- matrix(respPart, nrow = 1)
  clim_matresp <- refgrid
  clim_matresp$Data <- mat2Dto3Darray(matrespPart, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  attr(clim_matresp$Data, "climatology:fun") <- paste0("Partcor with respect to V",node)
  plotClimatology(clim_matresp,backdrop.theme = "coastline", main = list(paste0("partcorr(VX,V",node,")")),at = at, col.regions= col.b)
}



# cormatrix2.5 <- cor(TimeCoordsAnom_from_Grid_rms(tas_ncep), method = "spearman")
 cormatrix10d <- cor(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d), method = "spearman")
plotCorFixedV(280, cormatrix10d, tas_ncep_10d, RMS = FALSE, at = seq(-1,1,0.1))
plotCorFixedV(81, cormatrix10d, tas_ncep_10d, RMS = TRUE, at = seq(-1,1,0.1))
plotCorFixedV(441, cormatrix10d, tas_ncep_10d, RMS = TRUE, at = seq(-1,1,0.1))
plotcorFixedV(171, cormatrix2.5, tas_ncep, RMS = TRUE)
plotParFixedV(606, glasso01$wi, tas_ncep_10d, RMS = FALSE, at = seq(-0.2,0.2,0.05))
0.0001
blub <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
attr(blub, "VertexCoords")[81,]
31.73/2
################################################################################################
# Perform probability query correlation matrix
################################################################################################
library(mvtnorm)
library(transformeR)
library(magrittr)
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
str(data)
cormatrix <- cor(data, method = "spearman")
corRMSmatrix <- cor(dataRMS, method = "spearman")
covRMSmatrix <- cov(dataRMS, method ="spearman")



############################################### LOGLIK IS ZERO!! 
################################ For positive definite covpos returns value:
covpos <- make.positive.definite(covRMSmatrix)
R <- dmvnorm(dataRMS,sigma = covpos, log = T)
sum(R)
rank.condition(cormatrix)

m <- length(attr(data, "VertexCoords")$x)
sigma <- cormatrix 
upper <- rep(Inf, m)
lower2 <- rep(-Inf,m)
lower2[c(81)] <-1
prob2 <- pmvnorm(mean=rep(0, m), sigma, lower= lower2, upper= upper)
prob2

probs1 <- c()
i <- 1
for(i in 1:648){
lower <- rep(-Inf, m)
lower[c(81,i)] <- 1
probs1[i] <- pmvnorm(mean=rep(0, m), sigma, lower= lower, upper= upper)
}


probswithout <- c()
i <- 1
for(i in 1:648){
  lower <- rep(-Inf, m)
  lower[c(i)] <- 1
  probswithout[i] <- pmvnorm(mean=rep(0, m), sigma, lower= lower, upper= upper)
}

probscorwith <- probs1/prob2
probscorwith[117]
matcorwith <- matrix(probscorwith, nrow = 1)
matcorwithout <- matrix(probswithout, nrow = 1)
matcordif <- matcorwith - matcorwithout 
prop_corwith <- tas_ncep_10d
prop_corwithout <- tas_ncep_10d
prop_cordif <- tas_ncep_10d

prop_corwith$Data <- mat2Dto3Darray(matcorwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_corwithout$Data <- mat2Dto3Darray(matcorwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_cordif$Data <- mat2Dto3Darray(matcordif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
attr(prop_corwith$Data, "climatology:fun") <- "with evidence"
attr(prop_corwithout$Data, "climatology:fun") <- "without evidence"
attr(prop_cordif$Data, "climatology:fun") <- "with - without evidence"
plotClimatology(prop_corwith,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ")")),at = seq(0,1,0.05))
plotClimatology(prop_corwithout,backdrop.theme = "coastline", main = list(paste0("P(V >=1)")),at = seq(0,1,0.05))
plotClimatology(prop_cordif,backdrop.theme = "coastline", main = list(paste0("P(V >=1)")),at = seq(0.05,1,0.05))

################################################################################################
# Perform probability query correlation matrix package condMVnorm
###############################################################################################
library(condMVNorm)

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
str(data)
cormatrix <- cor(data, method = "spearman")
covmatrix <- cov(data, method = "spearman")
m <- length(attr(data, "VertexCoords")$x)
sigma <- cormatrix 

# Step 1: Obtaining conditional independence mean and covariance of {X1....X648}/{X81}|X81
# Warning: sigma is not positive definite. 
# Check if Precision matrix correlation matrix is positive definite. 
dependent <- 1:648
dependent <- dependent[-81]
test <- condMVN(mean = rep(0, m), sigma = cormatrix, dependent.ind = dependent, given.ind = c(81), X.given = 2, check.sigma= FALSE)
test$condMean
test$condVar

# Step 2: Perform query
# OR FIRST OUTINTEGRATE OTHER VARIABLES.
upper <- rep(Inf, m-1)
lower2 <- rep(-Inf,m-1)
lower2[c(116)] <-1
step2 <- pmvnorm(lower= lower2, upper= upper,mean=test$condMean, sigma = test$condVar)
step2
# or directly:
pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
         dependent.ind = dependent, given.ind = c(81), X.given = 1,
         check.sigma= FALSE)
pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
         dependent.ind = dependent, given.ind = c(81), X.given = 3,
         check.sigma= FALSE)
pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
         dependent.ind = dependent, given.ind = c(81), X.given = 2,
         check.sigma= FALSE)


# in loop:

Cprobwiths <- c()
Cprobwithouts <- c()
dependent <- 1:648
dependent <- dependent[-81]

upper <- rep(Inf, m-1)
upperw <- rep(Inf, m)

for(i in 1:647){
  lower2 <- rep(-Inf,m-1)
  lower2[c(i)] <-1
  Cprobwiths[i] <- pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
                            dependent.ind = dependent, given.ind = c(81), X.given = 1,
                            check.sigma= FALSE)
}
for(i in 1:648){
  lower2w <- rep(-Inf,m)
  lower2w[c(i)] <-1
  Cprobwithouts[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), corr = cormatrix)
}
length(Cprobwiths)
Cprobwiths <- append(Cprobwiths, c(1), after = 80)

Cprobwithouts

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
str(data)
cormatrix <- cor(data, method = "spearman")
all.equal(cormatrix, precision$w)
m

propagationCorr <- function(cormatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- cormatrix
  # nodesEvents <- 1:648
  # valueEvent <- 1
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  m <- nrow(cormatrix)
  names <- c()
  names
  for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V >", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
  }
  
  
  if (length(nodesEvidence) > 1){
    proves <- c()
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0("V",nodesEvidence[j]," = ", valueEvidence[j])
    }
    
    text <- "list("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],",")
    }
    text <- paste0(text,proves[length(nodesEvidence)],")")
    str2 <- text
    
    
    probname <- paste0("P(V >", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  upper <- rep(Inf, (m-length(nodesEvidence)))
  upperw <- rep(Inf, m)
  
  for(i in 1:(m-length(nodesEvidence))) {
    lower2 <- rep(-Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    lower2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
                              dependent.ind = dependent, given.ind = nodesEvidence, X.given = valueEvidence,
                              check.sigma= FALSE)
  }
  with[nodesEvidence] <- 1
  
  for(i in 1:m){
    lower2w <- rep(-Inf,m)
    lower2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), sigma = cormatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V >", valueEvent,")")
  df <- data.frame(names = names[nodesEvents], with = with, without = without)
  return(df)
}

propagationCov <- function(covmatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- cormatrix
  # nodesEvents <- 1:648
  # valueEvent <- 1
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  m <- nrow(covmatrix)
  names <- c()
  names
  for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V >", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
  }
  
  
  if (length(nodesEvidence) > 1){
    proves <- c()
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0("V",nodesEvidence[j]," = ", valueEvidence[j])
    }
    
    text <- "list("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],",")
    }
    text <- paste0(text,proves[length(nodesEvidence)],")")
    str2 <- text
    
    
    probname <- paste0("P(V >", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  upper <- rep(Inf, (m-length(nodesEvidence)))
  upperw <- rep(Inf, m)
  
  for(i in 1:(m-length(nodesEvidence))) {
    lower2 <- rep(-Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    lower2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower2, upper=upper, mean = rep(0, m), covmatrix,
                        dependent.ind = dependent, given.ind = nodesEvidence, X.given = valueEvidence,
                        check.sigma= FALSE)
  }
  with[nodesEvidence] <- 1
  
  for(i in 1:m){
    lower2w <- rep(-Inf,m)
    lower2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), sigma = covmatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V >", valueEvent,")")
  df <- data.frame(names = names[nodesEvents], with = with, without = without)
  return(df)
}


# nodesEvents <- 1:648
# valueEvent <- 1
# nodesEvidence <- c(81,280)
# valueEvidence <- c(1,1)
test <- propagationCorr(corcutsolo[[8]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
test <- propagationCorr(corRMSmatrix, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
testCOV <- propagationCov(covRMSmatrix, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
test2 <- propagationCorr(nnsigma, nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
test3 <- propHC1600_Vpos_V81V280_equal22
test3 <- prop_hc_200_400_V81V280_equal22
test3 <- prop_1_tabu_10d_g0_V81V280_equal22
dataRMS
testCOV$without
test2$without

props<- list(test,test2,test3)
props <- list(test3)
length(props)


source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
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
  
  a <- plotClimatology(prop_without,backdrop.theme = "coastline",main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- plotClimatology(prop_with,backdrop.theme = "coastline", main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- plotClimatology(prop_dif,backdrop.theme = "coastline", main = list(attr(prop_dif$Data,"climatology:fun"),cex = 0.5), at = seq(0.05,0.5,0.01), set.max = 0.5, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  b
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotwithouts, plotwiths, plotdifferences)
all
do.call("grid.arrange", c(all, ncol=3, nrow=3))
################################################################################################
# Perform probability query correlation matrix package tmvtnorm
###############################################################################################
library(matrixcalc)
library(tmvtnorm)
library(corpcor)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
dfRMS <- as.data.frame(dataRMS) 

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
str(data)
cormatrix <- cor(data, method = "spearman")
m <- length(attr(data, "VertexCoords")$x)
sigma <- cormatrix 

lower3 <- rep(-Inf,m)
lower3[c(81)] <-1

is.positive.definite(precision$w)
is.positive.definite(cormatrix)
det(precision$w)
det(cormatrix)


nncor <- make.positive.definite(cormatrix, tol = 1e-1)
det(nncor)

nn<- precision$w
nn <- make.positive.definite(nn,tol=1e-1)
det(nn)

truncated <- mtmvnorm(mean = rep(0, nrow(nn)), 
         sigma = nncor, 
         lower = lower3, 
         upper = rep(Inf, length = m), 
         doComputeVariance=FALSE,
         pmvnorm.algorithm=GenzBretz())

Tprobs <- c()
for(i in 1:648){
  lower4 <- rep(-Inf,m)
  lower4[c(i)] <-1
  Tprobs[i] <- pmvnorm(lower = lower4, upper = rep(Inf, length = m), mean = truncated$tmean, sigma =truncated$tvar)
}

Tprobswithout <- c()
for(i in 1:648){
  lower4 <- rep(-Inf,m)
  lower4[c(i)] <-1
  Tprobswithout[i] <- pmvnorm(lower = lower4, upper = rep(Inf, length = m), mean = rep(0, nrow(nn)) , sigma = nn)
}

probscorwith <- Tprobs
probscorwithout<- Tprobswithout
matcorwith <- matrix(probscorwith, nrow = 1)
matcorwithout <- matrix(probscorwithout, nrow = 1)
matcordif <- matcorwith - matcorwithout 
prop_corwith <- tas_ncep_10d
prop_corwithout <- tas_ncep_10d
prop_cordif <- tas_ncep_10d

prop_corwith$Data <- mat2Dto3Darray(matcorwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_corwithout$Data <- mat2Dto3Darray(matcorwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_cordif$Data <- mat2Dto3Darray(matcordif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
attr(prop_corwith$Data, "climatology:fun") <- "with evidence"
attr(prop_corwithout$Data, "climatology:fun") <- "without evidence"
attr(prop_cordif$Data, "climatology:fun") <- "with - without evidence"
plotClimatology(prop_corwith,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ")")),at = seq(0,1,0.05))
plotClimatology(prop_corwithout,backdrop.theme = "coastline", main = list(paste0("P(V >=1)")),at = seq(0,1,0.05))
plotClimatology(prop_cordif,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ") -P(V >=1)")),at = seq(0.05,0.5,0.01), set.max = 0.5)
################################################################################################
# Perform probability query precision matrix
################################################################################################
# library("mvtnorm")
# data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
# cormatrix <- cor(data, method = "spearman")
# precisionRMS <- glasso(cormatrix, rho = 0.01)
# precisioncovRMS <- glasso(cormatrix, rho = 0.01, wi.init=precisionRMS$wi)
# precisionRMS$wi[1,1]
# 
# m <- 648
# sigma <- precisioncovRMS$w
# lower <- rep(-Inf, m)
# lower[c(81,99)] <- 1
# upper <- rep(Inf, m)
# prob1 <- pmvnorm(mean=rep(0, m), sigma, lower= lower, upper= upper)
# lower2 <- rep(-Inf,m)
# lower2[c(81)] <-1
# prob2 <- pmvnorm(mean=rep(0, m), sigma, lower= lower2, upper= upper)
# prob1/prob2
# 



################################################################################################
# or obtain precision graph with fixed rho and treshold tau as in CN: 
# Case with th = 0.1 and precision (rho = 0.01) is case of Zerenner et al.
################################################################################################
th = 0.0353
# Adjacency matrix
f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
# f(2,3,m)
m <- glasso01$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
# adj.matrix <- abs(precision$wi) # beta = 0.01
diag(adj.matrix) <- 0
adj.matrix[abs(adj.matrix) <= th ] <- 0
adj.matrix[abs(adj.matrix) > th ] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
sum(degree(precgraph))/648
dev.off()
plot.Meteodag(data,precgraph)


################################################################################################
# Create measuress precision graph
################################################################################################
igraphske<- precgraph
# create graph_from_Grid object with TimeCoordsAnom_from_Grid 
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
#check
measuresPrec <- graph2measure(graphObject)
clim.awcon.Prec10d <- measure2clim(measuresPrec, what = "awconnectivity", ref.grid = tas_ncep_10d)
clim.degr.Prec10d <- measure2clim(measuresPrec, what = "degree", ref.grid = tas_ncep_10d)
clim.betw.Prec10d<- measure2clim(measuresPrec, what = "betweenness", ref.grid = tas_ncep_10d)
clim.close.Prec10d <- measure2clim(measuresPrec, what = "closeness", ref.grid = tas_ncep_10d)
clima.lclus.Prec10d <- measure2clim(measuresPrec, what = "localclustering", ref.grid = tas_ncep_10d)
plotClimatology(clim.betw.Prec10d, backdrop.theme = "coastline")
plotClimatology(clim.close.Prec10d, backdrop.theme = "coastline")
plotClimatology(clim.awcon.Prec10d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.Prec10d,backdrop.theme = "coastline")
plotClimatology(clim.degr.Prec10d,backdrop.theme = "coastline")
# Compare with bayesian 1400 
plotClimatology(clim.betw.Bay10d, backdrop.theme = "coastline")
plotClimatology(clim.close.Bay10d, backdrop.theme = "coastline")
plotClimatology(clim.awcon.Bay10d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.Bay10d,backdrop.theme = "coastline")

##################################################################################
# Investigate locations of vertices on worldmap
##################################################################################
example <- hc_edges_loglik_10d_1400_1600i$networks
str(example)
igraphdir <- igraph.from.graphNEL(as.graphNEL(example))
igraphske<- as.undirected(igraphdir)
# create graph_from_Grid object with TimeCoordsAnom_from_Grid
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
check
plot.Meteodag(graphObject,igraphdir)
dev.off()
# Make bn 
fitted <- bn.fit(example, as.data.frame(graphObject$data_coords))
# mutilate <- mutilated(fitted, evidence = list("V81" > 0.5))
# dagfitted <- bn.net(fitted)
# dagfitted
# dagfittedigraph <- igraph.from.graphNEL(as.graphNEL(dagfitted))
# dagfittedigraph
# dagmutilateigraph <- igraph.from.graphNEL(as.graphNEL(mutilate))
# dagmutilateigraph
# V(dagfittedigraph)[c("V338","V229")]$color <- "purple"
# V(dagfittedigraph)[V(dagfittedigraph)$name != "V338"]$color <- "green"
# V(dagfittedigraph)$color <- "yellow"
# plot
# V(dagfittedigraph)!="V338"
# V(dagfittedigraph)["V338"]$size <- 300
# V(dagfittedigraph)["V338"]$color
# V(dagfittedigraph)["V229"]$color
# dagmutilate$nodes$V338



meteodag <- dagfittedigraph
class(meteodag)

vertex_attr(meteodag, "color", index = V(meteodag))
vertex_attr(meteodag, "name", index = V(meteodag))

if (class(meteodag) == "igraph") 
  igraphDegree <- meteodag

if (class(meteodag) == "graphNEL") 
  igraphDegree <- igraph.from.graphNEL(meteodag)

if (class(meteodag) == "bn") 
  igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))

x <- attr(time.coords, "Xcoords", exact = FALSE)
y <- attr(time.coords, "Ycoords", exact = FALSE)

points <- expand.grid(y, x)[2:1]
points

plot(wrld)
plot.igraph(meteodag, 
            vertex.label = vertex_attr(meteodag, "name", index = V(meteodag)),
            vertex.label.cex = 0.5,
            vertex.label.degree = -pi/2,
            vertex.color = vertex_attr(meteodag, "color", index = V(meteodag)),
            vertex.size = 100,
            edge.arrow.size = 0,
            layout = as.matrix(points), add = TRUE, rescale = FALSE)




pdf("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Nodesgraph.pdf",
    width = 10, height = 6 )
dev.off()



##################################################################################
# Correlation matrixes 
# Indices come from Nino script
##################################################################################
library(corrplot)
# Select data for correlations
data <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
datanino <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = ind1981_2010)
dataColdbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = coldmonths)
dataWarmbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = warmmonths)
# make different correlation matrices
cormatrix <- cor(data, method = "spearman")
cormatrixnino <- cor(datanino, method = "spearman")
cormatrixwarm <- cor(dataWarmbasin, method ="spearman")
# Give variable names (Vertex Names)
# names as in Nino script
class(cormatrix)
names
rownames(cormatrix) <- colnames(cormatrix) <- names
rownames(cormatrixnino) <- colnames(cormatrixnino) <- names
rownames(cormatrixwarm) <- colnames(cormatrixwarm) <- names
# Plot some scatterplot
par("mar")
par(mar = c(1,1,1,1))
par(mfrow = c(1,3))
plot(cormatrix)
plot(cormatrixnino)
plot(cormatrixwarm)
# Corrplot from library corrplot
# Corrplots: "heatmaps" 
par("mar")
par(mar = c(1,1,1,1))
par(mfrow = c(1,3))
corrplot(cormatrix, method = "color", type = "full", is.corr = FALSE)
corrplot(cormatrixnino, method = "color", type = "full", is.corr = FALSE)
corrplot(cormatrixwarm, method = "color", type = "full", is.corr = FALSE)


pairs(cormatrix[1:5,1:5]) 



dev.off()


########################################################################################
# Correlation matrix with ggplot
########################################################################################
# Convert correlation matrices to ggplot format
library(reshape2)
library(ggplot2)
str(cormatrix)
melted_cormat <- melt(cormatrix, as.is = TRUE) # GOED: met variable namen!!! NU GROOT UITPRINTEN
str(melted_cormat)
melted_cormatabs <- melt(abs(cormatrix), as.is = TRUE)
melted_cormatnino <- melt(cormatrixnino, as.is = TRUE)
melted_cormatabsnino <- melt(abs(cormatrixnino), as.is = TRUE)
melted_cormatwarm <- melt(cormatrixwarm, as.is = TRUE)
melted_cormatabswarm <- melt(abs(cormatrixwarm), as.is = TRUE)

melted_cormatdif <- melted_cormat
melted_cormatdif[,3] <- melted_cormat[,3]-melted_cormatnino[,3]
melted_cormatdifabs <- melted_cormatabs
melted_cormatdifabs[,3] <- abs(melted_cormat[,3]-melted_cormatnino[,3])
melted_cormatdifabs

# Plot correlation matrices with geom_raster() 
# and highlight / select variables with scale _ discrete breaks /limits.
# For example Nino basin or various regions
somevar <- c(81,286,297,205,604)
indNino
somevar <- indNino
allvar <- as.vector(seq(1:ncol(data)))
offvar <- setdiff(allvar,somevar)

# normal correlation matrices
plotcor <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill= value))+ 
  geom_raster()+
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor all data\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])        
plotcor
plotcornino <- ggplot(data = melted_cormatnino, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster()+
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor nino exact\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrixnino)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrixnino)[somevar]) 
plotcornino
plotcorwarm <- ggplot(data = melted_cormatwarm, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor warm basin\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrixwarm)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrixwarm)[somevar]) 

grid.arrange(plotcor,plotcornino, plotcorwarm)

# Absolute correlation matrices 
plotcorabs <- ggplot(data = melted_cormatabs, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = c("white","orange","red"), 
                       #midpoint = 0, 
                       limit = c(0,.8), space = "Lab", 
                       name="|cor| all data\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])   

plotcorabsnino <- ggplot(data = melted_cormatabsnino, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = c("white","orange","red"), 
                       #midpoint = 0, 
                       limit = c(0,.8), space = "Lab", 
                       name="|cor| nino exact\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])   

plotcorabswarm <- ggplot(data = melted_cormatabswarm, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = c("white","orange","red"), 
                       #midpoint = 0, 
                       limit = c(0,.8), space = "Lab", 
                       name="|cor| warm basin\nSpearman") + 
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])   
grid.arrange(plotcorabs,plotcorabsnino, plotcorabswarm)

# differencias between all data and nino data or warm months.
plotcordif <- ggplot(data = melted_cormatdif, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-0.5,0.5), space = "Lab", 
                       name="Spearman\n cor all - cor nino") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar]) 

plotcordifabs <- ggplot(data = melted_cormatdifabs, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("white","yellow","green","blue"), 
                       #midpoint = 0, 
                       limit = c(0,0.5), space = "Lab", 
                       name="Spearman\n |cor all - cor nino|") + 
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar]) 

grid.arrange(plotcordifabs, plotcordif)


####################################################################################
# Try to order by magnitude correlation. (NOT YET: VARIABLE NAMES UNCLEAR)
#####################################################################################
order.AOE <- corrMatOrder(cormatrix, order = "AOE")
corrAOE <- cormatrix[order.AOE,order.AOE]
rownames(corrAOE)
colnames(corrAOE)
melted_cormatAOE <- melt(corrAOE, as.is = TRUE)
melted_cormatAOE

orderedcordifabs <- melted_cormatAOE
rownames(corrAOE)[somevar,]
rownames(corrAOE)[527]
indvector <- c()
for (i in 1:length(somevar)){
indvector <- append(indvector, which(rownames(corrAOE) == names[somevar][i]))
       }
rownames(corrAOE)[indvector]
plotorderedcordifabs<- ggplot(data = melted_cormatAOE, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  scale_x_discrete(breaks = rownames(corrAOE)[indvector]) +
  scale_y_discrete(breaks = rownames(corrAOE)[indvector])
plotorderedcordifabs
ord <- hclust(melted_cormatdifabs, method = "complete" )
dist(melted_cormatdifabs[,3])
melted_cormatdifabs
########################################################################################
# Correlation matrix with ggcorrplot: Cluster option
########################################################################################
library("ggcorrplot")
hclust(cormatrix)
names <- c()
names
for(i in 1:648) names <- append(names,paste0("V",i))
names
matrix
rownames(cormatrix) <- names
colnames(cormatrix) <- names
rownames(cormatrixnino) <- names
colnames(cormatrixnino) <- names
rownames(cormatrixwarm) <- names
colnames(cormatrixwarm) <- names

cormatrix
?hclust
x <- hclust(as.dist(cormatrix), method ="complete")
ggcorrplot(cormatrix, hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6) + scale_x_discrete(limits = c(81,286,315,205,604))
ggcorrplot(cormatrixnino,hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6) # HOW ORDER?
plotcorabsO <- ggcorrplot(abs(cormatrix), outline.color = NA, tl.cex = 2, tl.srt = pi/6)
plotcorabsninoO <- ggcorrplot(abs(cormatrixnino), outline.color = NA, tl.cex = 2, tl.srt = pi/6)
ggcorrplot(cormatrixnino - cormatrix, hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6)
plotcorabsdifO <- ggcorrplot(abs(cormatrixnino - cormatrix), hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6)

grid.arrange(plotcorabsO,plotcorabsninoO,plotcorabsdifO)
###################################################################################
# Clustering from data correlation matrix.
###################################################################################
x <- hclust(as.dist(cormatrix), method ="single")
x$order
x$labels
cormatrix
cormatrixPer <- cormatrix[x$order,x$order]
melt(cormatrixPer, as.is = TRUE)
plotcorPer <- ggplot(data = melt(cormatrixPer), aes(x=Var1, y=Var2, fill= value))+ 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor all data\nSpearman") +
  scale_x_discrete(breaks = c("V81","V286","V315","V205","V604")) +
  scale_y_discrete(breaks = c("V81","V286","V315","V205","V604"))

plotcorPer
####################################################################################
# Clustering by mean temperature ?
####################################################################################
data <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plotClimatology(climatology(tas_ncep_10d))
tas_ncep_10d$Data
time.coords.matrix <- array3Dto2Dmat(tas_ncep_10d$Data)
colMeans(time.coords.matrix)

datanino <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = ind1981_2010)
dataColdbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = coldmonths)
dataWarmbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = warmmonths)
# Different correlation matrices
str(data)

cormatrix <- cor(data, method = "spearman")
cormatrixnino <- cor(datanino, method = "spearman")
cormatrixwarm <- cor(dataWarmbasin, method ="spearman")


####################################################################################
# Investigate detection long range arcs.
####################################################################################
example <- hc_edges_loglik_10d_800_1000i$networks
igraphdir <- igraph.from.graphNEL(as.graphNEL(example))
igraphske<- as.undirected(igraphdir)
is_connected(igraphdir)
compdir <- component_distribution(igraphdir)
plot(compdir)
components(igraphdir)
subcomponent(igraphdir, "V244", mode = "all")



##############################################################################
# Try cpquery and mutilate on gaussian.test
#############################################################################
data(gaussian.test)
dag <- hc(gaussian.test)
fitted <- bn.fit(dag, gaussian.test)
# the result should be around 0.04.
cpquery(fitted,
        event = ((A >= 0) & (A <= 1)) & ((B >= 0) & (B <= 3)),
        evidence = (C + D < 10))
muti <- mutilated(fitted, evidence = list(C = 10))
bn.net(muti)
bn.net(fitted)
hamming(muti,fitted)
as.bn(fitted)
##############################################################################
# loglikelihood testing
###############################################################################
library("LAM")
library("MASS")
library("corpcor")
glasso25$loglikç
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)

cormatrix <- cor(data, method = "spearman")
nncor <- make.positive.definite(cormatrix, tol = 1e-1)
det(nncor)
covmatrix <- cov(data, method = "spearman")
nncov <- make.positive.definite(covmatrix, tol = 1e-4)
det(nncov) 

density <- dmvnorm(data, sigma = nncor,log = TRUE)
 sum(density)
 density2 <- dmvnorm(data, sigma = glasso01$w, log = TRUE)
 sum(density2)

 
 
 
M <- colMeans(data)
M
covmatrix <- cov(data, method =  "spearman")
suff <- suff_stat_NA_pattern(data)
LAM::loglike_mvnorm( M=M , S=cormatrix , mu=M, Sigma=precisionr020$w , n = 360 , lambda = 0 )
LAM::loglike_mvnorm_NA_pattern(suff, mu = M, Sigma = glasso01$w)
dmvnorm(x=c(0,0), mean=c(1,1))
glasso225
norm(e, type = "1")
glasso05$loglik

extra <- sum(abs(0.20*precisionr020$w))*180
precisionr020$loglik - extra



# the same;
hc_edges_loglik_10d_1400_1600i$networkdata$logliks
tabu_1_eBIC_g0.5$networkdata
glasso
critfun
bic
2809 -1513
1296/2
test <- hc(gaussian.test, score = "bic-g", k = 0.5)
test
score(test,gaussian.test, type = "bic-g", k = 0)
nparams(test,gaussian.test)
test


# new
gaussian()

library(mvtnorm)
library(bnlearn)
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
littles <- dmvnorm(data, mean = M, sigma = glasso05$w, log = TRUE)
sum(littles)
logLik(glm)
AIC
stats:::logLik.glm
stats:::AIC
AIC.logLik

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
fit <- bn.fit(hc_edges_loglik_10d_1000_1200i$networks, data = as.data.frame(data))
fit <- bn.fit(tabu_1_eBIC_g0.1$networks$tabu_10d_g0.1, data = as.data.frame(data))
logLik(fit,as.data.frame(data))
logLik.mlm
bnlearn::sigma(fit)
lm
fitted(fit)






vcov(summary(fit))

data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
data2 <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)
S <- cor(data, method = "spearman")
S2 <- cov(data2, method = "spearman")
all.equal(S,S2)
AIC
M <- colMeans(data)
loglike_mvnorm(M = M, S = S, mu = M, Sigma = glasso01$w, n = 180,log = TRUE)
