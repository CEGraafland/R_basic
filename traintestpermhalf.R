#######################################################################################
# Propagation train half
#######################################################################################
rm(list = ls())
# traintest.
library(mvtnorm)
library(corpcor)
library(gridExtra)
library(condMVNorm)
library(transformeR)
library(visualizeR)
library(magrittr)
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(RColorBrewer)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
#########################################################################################
# Load data train
#########################################################################################
#######################################################################################
# Half om half data.
#######################################################################################
initialdata <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
indTRAINhalf <-1:180
indTESThalf <- (1:360)[-indTRAINhalf]


# make data TEST and TRAIN
learndata <- initialdata[indTRAINhalf,]
attributes(learndata) <- attributes(initialdata)[2:5]
attr(learndata, "dim") <- c(180,648)
testdata <- initialdata[indTESThalf,]
attributes(testdata) <- attributes(initialdata)[2:5]
attr(testdata, "dim") <- c(180,648)
##########################################################################################
# load data train
filestrain1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permhalftrain", full.names = T)
filestrain1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permhalftrain")
filestrain1names <- gsub(".rda", "", filestrain1names)
filestrain1names

variablelisttrain1 <- list(networks = list())
names <- c()

for (i in 1:length(filestrain1)){
  variablepos <- get(load(filestrain1[i]))
  variablelisttrain1$networks[[i]] <- variablepos
}

names(variablelisttrain1$networks) <- filestrain1names
assign(paste0("hc_1_trainhalf"),variablelisttrain1)
nedgestrain1 <- sapply(hc_1_trainhalf$networks, narcs)
sorted <- sort(nedgestrain1, index.return = TRUE)
hc_1_trainhalf$networks <- hc_1_trainhalf$networks[sorted$ix]
sapply(hc_1_trainhalf$networks, narcs)
##############################################################################################
# load data test
#############################################################################################
filestest1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permhalftest", full.names = T)
filestest1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permhalftest")
filestest1names <- gsub(".rda", "", filestest1names)
filestest1names

variablelisttest1 <- list(networks = list())
names <- c()

for (i in 1:length(filestest1)){
  variablepos <- get(load(filestest1[i]))
  variablelisttest1$networks[[i]] <- variablepos
}

names(variablelisttest1$networks) <- filestest1names
assign(paste0("hc_1_testhalf"),variablelisttest1)
nedgestest1 <- sapply(hc_1_testhalf$networks, narcs)
sorted <- sort(nedgestest1, index.return = TRUE)
hc_1_testhalf$networks <- hc_1_testhalf$networks[sorted$ix]
sapply(hc_1_testhalf$networks, narcs)



###########################################################################
# Make fits on structure. 
###########################################################################
dataRMStrain <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTRAINhalf))
dataRMStest <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTESThalf))
fitstrain <- lapply(hc_1_trainhalf$networks, bn.fit, data = dataRMStrain)
fitstest <- lapply(hc_1_testhalf$networks, bn.fit, data = dataRMStest)

loglikstrain <- c()
loglikstest <- c()
loglikstesttest <- c()
for (i in 1:length(fitstrain)){
  loglikstrain[i] <- logLik(fitstrain[[i]], data = dataRMStrain)
  loglikstesttest[i] <- logLik(fitstest[[i]], data = dataRMStest)
  loglikstest[i] <- logLik(fitstrain[[i]], data = dataRMStest)
}
nedgestrain1 <- sapply(hc_1_trainhalf$networks, narcs)
nedgestest1 <- sapply(hc_1_testhalf$networks,narcs)
x <- seq(0,8800,100)
y <- seq(100,8900,100)

dataRMStrainframes <- apply(dataRMStrain, MARGIN = 1, FUN = function(x) as.data.frame(x))
dataRMStrainframes <- lapply(dataRMStrainframes, t)
dataRMStrainframes <- lapply(dataRMStrainframes, as.data.frame)
dataRMStestframes <- apply(dataRMStest, MARGIN = 1, FUN = function(x) as.data.frame(x))
dataRMStestframes <- lapply(dataRMStestframes, t)
dataRMStestframes <- lapply(dataRMStestframes, as.data.frame)

loglikstrain1BNevent <- array(data = NA, dim = c(nrow(dataRMStrain),length(fitstrain)))
loglikstest1BNevent <- array(data = NA, dim = c(nrow(dataRMStest),length(fitstrain)))

for (i in 1:length(fitstrain)){
  loglikstrain1BNevent[,i] <- sapply(dataRMStrainframes, FUN = logLik, object = fitstrain[[i]])
  loglikstest1BNevent[,i] <- sapply(dataRMStestframes, FUN = logLik, object = fitstrain[[i]])
}



databestBN <- which(loglikstest == max(loglikstest))
nedgestest1[databestBN]

###################################################################################
# Make fit on structure CN
###################################################################################
tau<- seq(0, 0.99, 0.01)
tau <- c(tau,0.026)
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
         0.026,
         0.025,
         0)

length(tau)
tauchar <- as.character(tau)
labelsCM <- c()

for(i in 1:length(tauchar)){
  labelsCM[i] <- paste0("= ",tauchar[i])
}

labelsCM
traingraphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = indTRAINhalf)
testgraphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = indTESThalf)
# measures10d <- lapply(graphs10d, graph2measure)

graphssolo  <- lapply(traingraphs10d, function(m) m$graph)
nedgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(nedgesnetworksCM, length)

graphssolotrain <- lapply(traingraphs10d, function(m) m$graph)
nedgesnetworkstrainCM <- lapply(graphssolotrain, E)
nedgesnetworkstrainCM <- sapply(nedgesnetworkstrainCM, length)
graphssolotest <- lapply(testgraphs10d, function(m) m$graph)
nedgesnetworkstestCM <- lapply(graphssolotest, E)
nedgesnetworkstestCM <- sapply(nedgesnetworkstestCM, length)


corssolo <- lapply(traingraphs10d, function(m) m$correlation)
corssolotrain <- lapply(traingraphs10d, function(m) m$correlation)
corssolotest <- lapply(testgraphs10d, function(m) m$correlation)


corcutsolo <- list()
corcutsolotrain<- list()
corcutsolotest <- list()

for (i in 1:length(corssolo)){
  corcutsolo[[i]] <- corssolo[[i]]
  corcutsolo[[i]][abs(corcutsolo[[i]]) <= tau[i] ] <- 0
  corcutsolotrain[[i]] <- corssolotrain[[i]]
  corcutsolotrain[[i]][abs(corcutsolotrain[[i]]) <= tau[i] ] <- 0
  corcutsolotest[[i]] <- corssolotest[[i]]
  corcutsolotest[[i]][abs(corcutsolotest[[i]]) <= tau[i] ] <- 0
}

loglikstrain1CN <- c()
loglikstest1CN <- c()
loglikstesttest1CN <- c()

dataRMStrain <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTRAINhalf))
dataRMStest <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTESThalf))

for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolo[[i]], tol = 0.06)
  nncortest <- make.positive.definite(corcutsolotest[[i]], tol = 0.06)
  densitytrain <- mvtnorm::dmvnorm(dataRMStrain, sigma = nncor,log = TRUE)
  loglikstrain1CN[i] <- sum(densitytrain)
  densitytest<- mvtnorm::dmvnorm(dataRMStest, sigma = nncor,log = TRUE)
  loglikstest1CN[i] <- sum(densitytest)
  densitytesttest<- mvtnorm::dmvnorm(dataRMStest, sigma = nncortest,log = TRUE)
  loglikstesttest1CN[i] <- sum(densitytesttest)
}
loglikstrain1CN
loglikstest1CN
loglikstesttest1CN

plot(nedgesnetworksCM,loglikstrain1CN, col = "black")
points(nedgesnetworksCM,loglikstesttest1CN, col = "brown")
points(nedgesnetworksCM,loglikstest1CN, col = "grey")

databest <- which(loglikstest1CN == max(loglikstest1CN))
nedgesnetworksCM[databest]


loglikstrain1CNevent <- array(data = NA, dim = c(nrow(dataRMStrain),length(tau)))
loglikstest1CNevent <- array(data = NA, dim = c(nrow(dataRMStest),length(tau)))
loglikstesttest1CNevent <- array(data = NA, dim = c(nrow(dataRMStest),length(tau)))


for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolotrain[[i]], tol = 0.06)
  densitytrain <- mvtnorm::dmvnorm(dataRMStrain, sigma = nncor,log = TRUE)
  loglikstrain1CNevent[,i] <- densitytrain
  nncortest <- make.positive.definite(corcutsolotest[[i]], tol = 0.06)
  densitytest<- mvtnorm::dmvnorm(dataRMStest, sigma = nncor,log = TRUE)
  loglikstest1CNevent[,i] <- densitytest
}
#############################################################################
# Plot events logliks train and test
#############################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/Loglik_datarealizations_traintest15Y.pdf")
pdf(plotname)
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Loglikelihood/Loglik_datarealizations_traintest15Y.pdf")
pdf(plotname)

plot(indTRAINhalf, loglikstrain1CNevent[,1],col = "black", xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), xlim = c(0,360),ylim = c(-1300,0),main = "Loglik train and test data realizations given trainmodels\n Period: First 15 Year", pch = 1)
points(indTESThalf, loglikstest1CNevent[,1],col = "black", pch = 2)
points(indTRAINhalf, loglikstrain1BNevent[,18], col = "lightblue", pch = 1)
points(indTESThalf, loglikstest1BNevent[,18], col = "lightblue",pch = 2)

legend("topleft", c(paste0("CN train ",nedgesnetworksCM[1]), paste0("CN test ",nedgesnetworkstestCM[1]),paste0("BN train: ",nedgestrain1[18]),paste0("BN test: ",nedgestest1[18])),bty = "n", pch = c(1,2,1,2), col =c("black","black","blue","blue"), cex = 0.6)

dev.off()



############################################################################################
# HC TRAIN propagation V81 V227 extra evidence (positive)
############################################################################################
proprange <- 18  
for (i in proprange){
  
  assign(paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
  
  assign(paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
}

############################################################################################
# HC TRAIN propagation V81 V227 extra evidence (negative)
############################################################################################
proprange <- 18
x[proprange]
y[proprange]
for (i in proprange){
  
  assign(paste0("propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
  
  assign(paste0("propneg_hc_permhalftest_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_permhalftest_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/propneg_hc_permhalftest_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
}

#################################################################################
# Single evidence. V227 postive
#################################################################################
proprange <- 18
x[proprange]
y[proprange]
for (i in proprange){
  
  assign(paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_",x[i],"_",y[i],"i_V227_equal2.rda"))
  
  assign(paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_",x[i],"_",y[i],"i_V227_equal2.rda"))
}

#################################################################################
# Single evidence. V280 postive
#################################################################################
proprange <- 18
for (i in proprange){
  
  assign(paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V280_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(280),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V280_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_",x[i],"_",y[i],"i_V280_equal2.rda"))
  
  assign(paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V280_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(280),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V280_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_",x[i],"_",y[i],"i_V280_equal2.rda"))
}

#################################################################################
# Single evidence. V227 neg
#################################################################################
proprange <- 18
x[proprange]
y[proprange]
for (i in proprange){
  
  assign(paste0("propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V227_equal2.rda"))
  
  assign(paste0("propneg_hc_permhalftest_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_permhalftest_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/propneg_hc_permhalftest_",x[i],"_",y[i],"i_V227_equal2.rda"))
}


#################################################################################
# Single evidence. V81 postive
#################################################################################
proprange <- 18
for (i in proprange){
  
  assign(paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftrain_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
  assign(paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_permhalftest_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_",x[i],"_",y[i],"i_V81_equal2.rda"))
}

#################################################################################
# Single evidence. V81 neg
#################################################################################
proprange <- 18
x[proprange]
y[proprange]
for (i in proprange){
  
  assign(paste0("propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/propneg_hc_permhalftrain_",x[i],"_",y[i],"i_V81_equal2.rda"))
  
  assign(paste0("propneg_hc_permhalftest_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_permhalftest_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/propneg_hc_permhalftest_",x[i],"_",y[i],"i_V81_equal2.rda"))
}




###############################################################################
# Visualize V81 MAKE CORCUTSOLO TRAIN and TEST: HALF
###############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_1600_1700i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_1600_1700i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_1700_1800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_1700_1800i_V81_equal2.rda")
numberCM <- 1
numberBN <- 18

propcortrain1 <- propagationCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
propcortest1 <- propagationCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
diffcor1 <- propcortrain1
diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
diffcor1$without <- 0

prophctrain1 <- prop_hc_permhalftrain_1700_1800i_V81_equal2
prophctest1 <- prop_hc_permhalftest_1700_1800i_V81_equal2
diffhc <- prophctrain1 
diffhc$with <- abs(prophctrain1$with - prophctest1$with)
diffhc$without <- 0



props <- list(propcortrain1,propcortest1,diffcor1, prophctrain1,prophctest1,diffhc)

logliksplot <- c(loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA)
nedgesplot <- c(nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA)
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
  cb <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # if(attr(prop, "sign") == "pos"){cb <- col.r}
  # if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.875,0.01), region = TRUE, col.regions= cb,set.max = 0.875, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, nrow = 2))

prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))

do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))

###############################################################################
# Visualize V81 Negative MAKE CORCUTSOLO TRAIN and TEST: HALF
###############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/propneg_hc_permhalftrain_1700_1800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/propneg_hc_permhalftest_1700_1800i_V81_equal2.rda")
numberCM <- 1
numberBN <- 18

propcortrain1 <- propagationnegCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
propcortest1 <- propagationnegCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
diffcor1 <- propcortrain1
diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
diffcor1$without <- 0

prophctrain1 <- propneg_hc_permhalftrain_1700_1800i_V81_equal2
prophctest1 <- propneg_hc_permhalftest_1700_1800i_V81_equal2
diffhc <- prophctrain1 
diffhc$with <- abs(prophctrain1$with - prophctest1$with)
diffhc$without <- 0

props <- list(propcortrain1,propcortest1,diffcor1, prophctrain1,prophctest1,diffhc)

logliksplot <- c(loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA)
nedgesplot <- c(nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA)
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
  cb <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  # col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # if(attr(prop, "sign") == "pos"){cb <- col.r}
  # if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.875,0.01), region = TRUE, col.regions= cb,set.max = 0.875, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, nrow = 2))

prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))

do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))

###############################################################################
# Visualize V227 MAKE CORCUTSOLO TRAIN and TEST: HALF
###############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_1700_1800i_V227_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_1700_1800i_V227_equal2.rda")
numberCM <- 1
numberBN <- 18

propcortrain1 <- propagationCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
propcortest1 <- propagationCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
diffcor1 <- propcortrain1
diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
diffcor1$without <- 0

prophctrain1 <- prop_hc_permhalftrain_1700_1800i_V227_equal2
prophctest1 <- prop_hc_permhalftest_1700_1800i_V227_equal2
prop_hc_permhalftest_1700_1800i_V227_equal2$without[227]
prop_hc_permhalftrain_1700_1800i_V227_equal2$without[227]
prop_hc_1600_1800i_V227_equal2$without
propcortrain1$without

diffhc <- prophctrain1 
diffhc$with <- abs(prophctrain1$with - prophctest1$with)
diffhc$without <- 0

props <- list(propcortrain1,propcortest1,diffcor1, prophctrain1,prophctest1,diffhc)

logliksplot <- c(loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA)
nedgesplot <- c(nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA)
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
  cb <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # if(attr(prop, "sign") == "pos"){cb <- col.r}
  # if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.875,0.01), region = TRUE, col.regions= cb,set.max = 0.875, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, nrow = 2))

prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))

do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))

###############################################################################
# Visualize V227 Negative MAKE CORCUTSOLO TRAIN and TEST: HALF
###############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/propneg_hc_permhalftrain_1700_1800i_V227_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/propneg_hc_permhalftest_1700_1800i_V227_equal2.rda")
numberCM <- 1
numberBN <- 18

propcortrain1 <- propagationnegCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(227), valueEvidence = c(2))
propcortest1 <- propagationnegCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(227), valueEvidence = c(2))
diffcor1 <- propcortrain1
diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
diffcor1$without <- 0

prophctrain1 <- propneg_hc_permhalftrain_1700_1800i_V227_equal2
prophctest1 <- propneg_hc_permhalftest_1700_1800i_V227_equal2
diffhc <- prophctrain1 
diffhc$with <- abs(prophctrain1$with - prophctest1$with)
diffhc$without <- 0

props <- list(propcortrain1,propcortest1,diffcor1, prophctrain1,prophctest1,diffhc)

logliksplot <- c(loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA)
nedgesplot <- c(nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA)
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
  cb <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  # col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # if(attr(prop, "sign") == "pos"){cb <- col.r}
  # if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.875,0.01), region = TRUE, col.regions= cb,set.max = 0.875, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, nrow = 2))

prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))

do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))

###############################################################################
# Visualize V81V227 MAKE CORCUTSOLO TRAIN and TEST. 
###############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_1600_1700i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_1600_1700i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/prop_hc_permhalftrain_1700_1800i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/prop_hc_permhalftest_1700_1800i_V81V227_equal22.rda")
numberCM <- 1
numberCM2 <- 5
numberBN <- 18

propcortrain1 <- propagationCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
propcortest1 <- propagationCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
diffcor1 <- propcortrain1
diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
diffcor1$without <- 0
# propcortrain2 <- propagationCorr(corcutsolotrain[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
# propcortest2 <- propagationCorr(corcutsolotest[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
# diffcor2 <- propcortrain2
# diffcor2$with <- abs(propcortrain2$with - propcortest2$with)
# diffcor2$without <- 0




prophctrain1 <- prop_hc_permhalftrain_1700_1800i_V81V227_equal22
prophctest1 <- prop_hc_permhalftest_1700_1800i_V81V227_equal22
diffhc <- prophctrain1 
diffhc$with <- abs(prophctrain1$with - prophctest1$with)
diffhc$without <- 0

props <- list(prophctrain1,prophctest1,diffhc, propcortrain1,propcortest1,diffcor1)
              # ,propcortrain2,propcortest2,diffcor2)

logliksplot <- c(loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA,loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain1CN[[numberCM2]],loglikstesttest1CN[[numberCM2]],NA)
nedgesplot <- c(nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA,nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgesnetworkstrainCM[[numberCM2]],nedgesnetworkstestCM[[numberCM2]],NA)
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
  cb <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # if(attr(prop, "sign") == "pos"){cb <- col.r}
  # if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all

do.call("grid.arrange", c(all, nrow = 2, top = "train grow, test grow, abs(train-test)"))

prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n BN: loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n CN: loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# prop_dif_plus3 <- quantity2clim(abs(propcortrain2$with - propcortrain2$without - propcortest2$with + propcortest2$without), paste0(attr(propcortrain2$with, "probability"),"-", attr(propcortrain2$without, "probability")), tas_ncep_10d)
# dif3 <- spatialPlot(prop_dif_plus3,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n CN: loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))
# all.equal(test5$with-test5$without,test8$with-test8$without)
# all.equal(test1$with-test1$without,test4$with-test4$without)  
do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))

###############################################################################
# Visualize V81V227 negative MAKE CORCUTSOLO TRAIN and TEST. 
###############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/trainhalf/propneg_hc_permhalftrain_1700_1800i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/testhalf/propneg_hc_permhalftest_1700_1800i_V81V227_equal22.rda")
numberCM <- 1
numberBN <- 18

propcortrain1 <- propagationnegCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
propcortest1 <- propagationnegCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
diffcor1 <- propcortrain1
diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
diffcor1$without <- 0
# propcortrain2 <- propagationCorr(corcutsolotrain[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
# propcortest2 <- propagationCorr(corcutsolotest[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
# diffcor2 <- propcortrain2
# diffcor2$with <- abs(propcortrain2$with - propcortest2$with)
# diffcor2$without <- 0




prophctrain1 <- propneg_hc_permhalftrain_1700_1800i_V81V227_equal22
prophctest1 <- propneg_hc_permhalftest_1700_1800i_V81V227_equal22
diffhc <- prophctrain1 
diffhc$with <- abs(prophctrain1$with - prophctest1$with)
diffhc$without <- 0

props <- list(prophctrain1,prophctest1,diffhc, propcortrain1,propcortest1,diffcor1)
# ,propcortrain2,propcortest2,diffcor2)

logliksplot <- c(loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA,loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain1CN[[numberCM2]],loglikstesttest1CN[[numberCM2]],NA)
nedgesplot <- c(nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA,nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgesnetworkstrainCM[[numberCM2]],nedgesnetworkstestCM[[numberCM2]],NA)
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
  cb <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  # col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  # if(attr(prop, "sign") == "pos"){cb <- col.r}
  # if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all

do.call("grid.arrange", c(all, nrow = 2, top = "train grow, test grow, abs(train-test)"))

prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n BN: loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n CN: loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# prop_dif_plus3 <- quantity2clim(abs(propcortrain2$with - propcortrain2$without - propcortest2$with + propcortest2$without), paste0(attr(propcortrain2$with, "probability"),"-", attr(propcortrain2$without, "probability")), tas_ncep_10d)
# dif3 <- spatialPlot(prop_dif_plus3,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n CN: loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))
# all.equal(test5$with-test5$without,test8$with-test8$without)
# all.equal(test1$with-test1$without,test4$with-test4$without)  
do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))
        