# traintest.
rm(list = ls())
library(bnlearn)
library(transformeR)
library(magrittr)
library(mvtnorm)
library(corpcor)
library(gridExtra)
library(condMVNorm)
# library(pcalg)
# library(igraph)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
##############################################################################################
# TERMINAL: Trains structure learning
##############################################################################################
initialdata <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
# indTRAIN <- sample(1:360, 180, replace = FALSE, prob = NULL)
# indTRAIN1 <- sort(indTRAIN)
# save(indTRAIN1, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
indTEST1 <- (1:360)[-indTRAIN1]

# make data TEST and TRAIN
learndata <- initialdata[indTRAIN1,]
attributes(learndata) <- attributes(initialdata)[2:5]
attr(learndata, "dim") <- c(180,648)
testdata <- initialdata[indTEST1,]
attributes(testdata) <- attributes(initialdata)[2:5]
attr(testdata, "dim") <- c(180,648)

# # learning TRAIN with unpermutated data
#   data <- as.data.frame(learndata)
#   start <- NULL
#   steps <- 100
#   last <- 10000
# 
#   for (m in 0:(last/steps)) {
#     i <- m*steps
#     j <- i+steps
#     berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
#     assign(paste0("hcTrain1_",i,"_",j,"i"), berekening)
# 
#     save(list = paste0("hcTrain1_",i,"_",j,"i"),
#          file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1train/hcTrain1_",i,"_",j,"i.rda"))
# 
#     start <- berekening
#   }

# learning TEST with unpermutated data
  data <- as.data.frame(testdata)
  start <- NULL
  steps <- 100
  last <- 10000

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hcTest1_",i,"_",j,"i"), berekening)

    save(list = paste0("hcTest1_",i,"_",j,"i"),
         file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1test/hcTest1_",i,"_",j,"i.rda"))

    start <- berekening
  }




##################################################################################################
# Correlation network with same data
##################################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
  tau<- seq(0, 0.99, 0.01)
  tau <- c(tau,0.026)
  tau <- seq(0,0.9,0.025)
tau <- seq(0.00,0.99,0.01)
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
traingraphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = indTRAIN1)
testgraphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = indTEST1)
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

dataRMStrain <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTRAIN1))
dataRMStest <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTEST1))

for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolotrain[[i]], tol = 0.06)
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


##################################################################################################
# Propagation
# Loglik with traindata and testdata.
##################################################################################################
####################################################################################
# hc train1 analyse
####################################################################################
filestrain1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1train", full.names = T)
filestrain1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1train")
filestrain1names <- gsub(".rda", "", filestrain1names)
filestrain1names

variablelisttrain1 <- list(networks = list())
names <- c()

for (i in 1:length(filestrain1)){
  variablepos <- get(load(filestrain1[i]))
  variablelisttrain1$networks[[i]] <- variablepos
}

names(variablelisttrain1$networks) <- filestrain1names
assign(paste0("hc_1_train1"),variablelisttrain1)
nedgestrain1 <- sapply(hc_1_train1$networks, narcs)
sorted <- sort(nedgestrain1, index.return = TRUE)
hc_1_train1$networks <- hc_1_train1$networks[sorted$ix]
nedgestrain1 <- sapply(hc_1_train1$networks, narcs)
##############################################################################################
# load data test
#############################################################################################
filestest1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1test", full.names = T)
filestest1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1test")
filestest1names <- gsub(".rda", "", filestest1names)
filestest1names

variablelisttest1 <- list(networks = list())
names <- c()

for (i in 1:length(filestest1)){
  variablepos <- get(load(filestest1[i]))
  variablelisttest1$networks[[i]] <- variablepos
}

names(variablelisttest1$networks) <- filestest1names
assign(paste0("hc_1_test1"),variablelisttest1)
nedgestest1 <- sapply(hc_1_test1$networks, narcs)
sorted <- sort(nedgestest1, index.return = TRUE)
hc_1_test1$networks <- hc_1_test1$networks[sorted$ix]
nedgestest1 <- sapply(hc_1_test1$networks, narcs)


####################################################################################################
dataRMStrain <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTRAIN1))
dataRMStest <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTEST1))
fitstrain <- lapply(hc_1_train1$networks, bn.fit, data = dataRMStrain)
fitstest <- lapply(hc_1_test1$networks, bn.fit, data = dataRMStest)
  
loglikstrain <- c()
loglikstest <- c()
loglikstesttest <- c()
for (i in 1:length(fitstrain)){
  loglikstrain[i] <- logLik(fitstrain[[i]], data = dataRMStrain)
  loglikstest[i] <- logLik(fitstrain[[i]], data = dataRMStest)
  loglikstesttest[i] <- logLik(fitstest[[i]], data = dataRMStest)
}

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


points(nedgestrain1,loglikstrain, col = "green")
points(nedgestest1[1:length(loglikstesttest)],loglikstest, col = "blue")
points(nedgestrain1,loglikstest, col = "red")


which(loglikstest == max(loglikstest))
bndatabest <- which(loglikstest == max(loglikstest))
nedgestrain1[bndatabest]
############################################################################################
# check average markov blanket size of networks.
############################################################################################
names <- c()
names
for(i in 1:648) names <- append(names,paste0("V",i))

avmb <- c()
for(j in 1:length(hc_1_train1$networks)){
  all <- c()
  for(k in 1:length(names)){
    name <- names[k]
    all[k] <- length(mb(hc_1_train1$networks[[j]], node = name))
  }
  avmb[j] <- mean(all)
}

avparents <- c()
for(j in 1:length(hc_1_train1$networks)){
  all <- c()
  for(k in 1:length(names)){
    name <- names[k]
    all[k] <- length(parents(hc_1_train1$networks[[j]], node = name))
  }
  avparents[j] <- mean(all)
}

hc_1_train1$networks$hcTrain1_1100_1200i

dev.off()

plot(nedgestrain1,loglikstrain)
par(new = T)
plot(nedgestrain1, avmb)
par(new=T)
plot(nedgestrain1,c(0,nedgestrain1[2:length(nedgestrain1)]-nedgestrain1[1:(length(nedgestrain1)-1)]))
##########################################################################################
# Figures for in Resumen 5
##########################################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/Traintestsmall.pdf")
pdf(plotname)
plot(nedgestrain1,loglikstrain, col = "red", xlab = "number of edges", ylab = bquote("log P(X|G,"~Theta~")"),main = "Traintest 1980-2010 50%")
points(nedgestrain1,loglikstest, col = "orange")
points(nedgesnetworksCM,loglikstrain1CN, col = "black")
points(nedgesnetworksCM,loglikstest1CN, col = "grey")
taupoints <- databest
# taupoints <- 1:length(nedgesnetworksCM)
for (i in taupoints){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  label2 <- paste0("|E|=",nedgesnetworksCM[i])
  text(nedgesnetworksCM[i],loglikstrain1CN[i],label= label,col='blue', cex = 0.6, pos = 1)
  text(nedgesnetworksCM[i],loglikstrain1CN[i],label= label2,col='blue', cex = 0.6, pos = 3)
}

for (i in bndatabest){
  labelbn <- paste0("|E|=",nedgestest1[i])
  text(nedgestest1[i],loglikstest[i],label= labelbn,col='blue', cex = 0.6, pos = 3)
}

legend("topleft", c("BN: Train", "BN: test", "CN: Train", "CN: test"), bty = "n", pch = 1, col =c("red","orange", "black", "grey"), cex = 0.6)
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/Traintestmed.pdf")
pdf(plotname)
plot(nedgestrain1,loglikstrain, col = "red", xlim = c(0,95000),ylim = c(-170000,-40000),xlab = "number of edges", ylab = bquote("log P(X|G,"~Theta~")"),main = "Traintest 1980-2010 50%")
points(nedgestrain1,loglikstest, col = "orange")
points(nedgesnetworksCM,loglikstrain1CN, col = "black")
points(nedgesnetworksCM,loglikstest1CN, col = "grey")
for (i in 1:length(nedgesnetworksCM)){
  value <- tauchar[i]
  label <- bquote(tau ~ "=" ~ .(value))
  text(nedgesnetworksCM[i],loglikstrain1CN[i],label= label,col='blue', cex = 0.6, pos = 3)
}
legend("topleft", c("BN: Train", "BN: test", "CN: Train", "CN: test"), bty = "n", pch = 1, col =c("red","orange", "black", "grey"), cex = 0.6)
dev.off()





plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/Traintestbig.pdf")
pdf(plotname)
plot(nedgesnetworksCM,loglikstrain1CN, col = "black",ylim = c(-200000,-40000),xlab = "number of edges", ylab = bquote("log P(X|G,"~Theta~")"), main = "Traintest 1980-2010 50%")
points(nedgesnetworksCM,loglikstest1CN, col = "grey")
points(nedgestrain1,loglikstrain, col = "red")
points(nedgestrain1,loglikstest, col = "orange")
legend("topleft", c("BN: Train", "BN: test", "CN: Train", "CN: test"),bty = "n", pch = 1, col =c("red","orange", "black", "grey"), cex = 0.6)
dev.off()


#############################################################################
# Event logliks train and test
#############################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/Loglik_datarealizations_traintestR.pdf")
pdf(plotname)
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Loglikelihood/Loglik_datarealizations_traintestR.pdf")
pdf(plotname)

plot(indTRAIN1, loglikstrain1CNevent[,1],col = blackcn[3], xlab = "D_i", ylab = bquote("log P(D_i|G,"~Theta~")"), ylim = c(-1300,0),main = "Loglik train and test data realizations given trainmodels\n Period: Random 15 Year", pch = 1)
points(indTEST1, loglikstest1CNevent[,1],col = blackcn[3], pch = 2)
points(indTRAIN1, loglikstrain1BNevent[,18], col = bluebn[1], pch = 1)
points(indTEST1, loglikstest1BNevent[,18], col = bluebn[1],pch = 2)

legend("topleft", c(paste0("CN train ",nedgesnetworksCM[1]), paste0("CN test ",nedgesnetworkstestCM[1]),paste0("BN train: ",nedgestrain1[18]),paste0("BN test: ",nedgestest1[18])),bty = "n", pch = c(1,2,1,2), col =c(blackcn[3],blackcn[2],bluebn[1],bluebn[1]), cex = 0.6)

dev.off()







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

# learning TRAIN with unpermutated data
  data <- as.data.frame(learndata)
  start <- NULL
  steps <- 100
  last <- 9900

  for (m in 0:(last/steps)) {
    i <- m*steps
    j <- i+steps
    berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
    assign(paste0("hcTrainhalf_",i,"_",j,"i"), berekening)

    save(list = paste0("hcTrainhalf_",i,"_",j,"i"),
         file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permhalftrain/hcTrainhalf_",i,"_",j,"i.rda"))

    start <- berekening
  }

# learning TEST with unpermutated data
data <- as.data.frame(testdata)
start <- NULL
steps <- 100
last <- 9900

for (m in 0:(last/steps)) {
  i <- m*steps
  j <- i+steps
  berekening <- hc(data, max.iter = steps, score = "bic-g",start = start)
  assign(paste0("hcTesthalf_",i,"_",j,"i"), berekening)
  
  save(list = paste0("hcTesthalf_",i,"_",j,"i"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permhalftest/hcTesthalf_",i,"_",j,"i.rda"))
  
  start <- berekening
}




