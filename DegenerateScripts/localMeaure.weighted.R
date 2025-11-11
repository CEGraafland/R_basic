# OCEANO LOCAL:
rm(list = ls())
library(bnlearn)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")




load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN2.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwBN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwRenyi.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwLattices.rda")


gridGraphsCN <- gridGraphwCN2
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- lapply(edgelistsCN, length)
names(gridGraphsCN) <- as.character(numberofedgesCN)


gridGraphsBN <- gridGraphwBN
GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})
edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- sapply(edgelistsBN, length)
names(gridGraphsBN) <- as.character(numberofedgesBN)

gridGraphsRenyi <- gridGraphwRenyi
GraphsRenyi <- lapply(gridGraphsRenyi, function(x){x$graph})
edgelistsRenyi <- lapply(GraphsRenyi, E)
numberofedgesRenyi <- lapply(edgelistsRenyi, length)
names(gridGraphsRenyi) <- as.character(numberofedgesRenyi)

gridGraphsLattices <- gridGraphwLattices
GraphsLattices <- lapply(gridGraphsLattices,function(x){x$graph})
edgelistsLattices <- lapply(GraphsLattices, E)
numberofedgesLattices <- lapply(edgelistsLattices, length)
names(gridGraphsLattices) <- as.character(numberofedgesLattices)

names <- c("CN","BN",'Renyi',"Lattices")
gridGraphlists <- list(gridGraphwCN2,gridGraphwBN,gridGraphwRenyi,gridGraphwLattices)

edge_attr(gridGraphsRenyi$`100`$graph)

# with weights

for(i in 1:length(names)){
  name <- names[i]
  gridGraphlist <- gridGraphlists[[i]]
  globallist <- lapply(gridGraphlist, graph2measure.global)
  assign(paste0("globmeas",name),globallist)
  assign(paste0("globalasslist",name),lapply(globallist, function(x){x$globalass}))
  assign(paste0("globalcluslist",name),lapply(globallist, function(x){x$globalclus}))
  assign(paste0("globaledenslist",name), lapply(globallist, function(x){x$edens}))
  assign(paste0("globaldiamLCClist",name), lapply(globallist,function(x){x$diameterLCC}))
  assign(paste0("globaldiamGlist",name), lapply(globallist,function(x){x$diameterG}))
  globallist <- NULL
  
  locallist <- lapply(gridGraphlist, graph2measure.local)
  degreelist <- mapply(graph2measure.degree,graphObj = gridGraphlist, localmeasureObj = locallist, SIMPLIFY = FALSE)
  assign(paste0("degreemeas",name),degreelist)
  assign(paste0("degreeclustering",name),lapply(degreelist,function(x){x$degreeclustering}))
  assign(paste0("degreedist",name), lapply(degreelist, function(x){x$degrdist}))
  assign(paste0("meanstrength",name), lapply(degreelist, function(x){x$meanstrength}))
  assign(paste0("meandistance",name), lapply(degreelist, function(x){x$meandistance}))
  degreelist <- NULL
}


##################################################################
# How many friends do your friends have? "local assortativity"
# How many friends on average do people with k friends have? "k degree assortativity"
##################################################################
knnBN <- lapply(GraphsBN, knn)
names(knnBN) <- as.character(numberofedgesBN)
knnkBN <- lapply(knnBN, function(x){x$knnk})

knnCN <- lapply(GraphsCN, knn)
names(knnCN) <- names(GraphsCN)
knnkCN <- lapply(knnCN, function(x){x$knnk})

knnRenyi <- lapply(GraphsRenyi, knn)
names(knnRenyi) <- as.character(numberofedgesRenyi)
knnkRenyi <- lapply(knnRenyi, function(x){x$knnk})

knnLattices<- lapply(GraphsLattices, knn)
names(knnLattices) <- as.character(numberofedgesLattices)
knnkLattices <- lapply(knnLattices, function(x){x$knnk})



plotseqBN <- seq(6,length(knnkBN),length.out = 6)
lastBN <- plotseqBN[length(plotseqBN)]
restBN <- plotseqBN[1:(length(plotseqBN)-1)]
colrestBN <- rainbow(length(restBN))

plotseqRenyi <- plotseqBN
lastRenyi <- lastBN
restRenyi <- restBN
colrestRenyi <- rainbow(length(restRenyi))

plotseqLattices <- seq(1,length(knnkLattices),1)
lastLattices <- plotseqLattices[length(plotseqLattices)]
restLattices <- plotseqLattices[1:(length(plotseqLattices)-1)]
colrestLattices <- rainbow(length(restLattices))

plotseqCN <- seq(1,length(knnkCN),length.out = 6)
plotseqCN <- round(plotseqCN)
plotseqCN <- rev(plotseqCN)
lastCN <- plotseqCN[length(plotseqCN)]
restCN <- plotseqCN[1:(length(plotseqCN)-1)]
colrestCN <- rainbow(length(restCN))
numberofedgesCN[plotseqCN]

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/degreeassweight.pdf")
pdf(plotname)
drawdegree = TRUE
par(mfrow = c(2,2))


plot(1:length(knnkRenyi[[lastRenyi]]),knnkRenyi[[lastRenyi]], 
      xlim = c(0,50),ylim = c(0,30), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Random Renyi")
if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[lastRenyi]]),degreedistRenyi[[lastRenyi]]*648)}
text(x= length(knnkRenyi[[lastRenyi]]), y=mean(knnkRenyi[[lastRenyi]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[lastRenyi])
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}

plot(1:length(knnkLattices[[lastLattices]]),knnkLattices[[lastLattices]],
      xlim = c(0,50),ylim = c(0,50), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Regular Lattices")
if(drawdegree == TRUE){lines(1:length(degreedistLattices[[lastLattices]]),degreedistLattices[[lastLattices]]*648)}
abline(0,1)
for (i in 1:length(restLattices)){
  j <- restLattices[i]
  points(1:length(knnkLattices[[j]]),knnkLattices[[j]],col = colrestLattices[i])
  if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
  text(x= length(knnkLattices[[j]]), y=mean(knnkLattices[[j]],na.rm = TRUE), pos=4, labels= names(knnLattices)[j])
}


plot(1:length(knnkBN[[lastBN]]),knnkBN[[lastBN]], 
     xlim = c(0,100),ylim = c(0,60), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(knnkBN[[lastBN]]), y=mean(knnkBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(knnBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(knnkBN[[j]]),knnkBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  text(x= length(knnkBN[[j]]), y=mean(knnkBN[[j]],na.rm = TRUE), pos=4, labels= names(knnBN)[j])
}

length(knnkCN)
plot(1:length(knnkCN[[lastCN]]),knnkCN[[lastCN]], 
     xlim = c(0,100),ylim = c(0,100), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(knnkCN[[j]]),knnkCN[[j]],col = colrestCN[i])
  if(drawdegree == TRUE){lines(1:length(degreedistCN[[j]]),degreedistCN[[j]]*648,col = colrestCN[i])}
  text(x= length(knnkCN[[j]])+20, y=mean(knnkCN[[j]],na.rm = TRUE), pos=4, labels= names(knnCN)[j])
}
names(knnCN)[plotseqCN]
dev.off()

#####################################################################################
# degree clustering
# How do vertices of degree k cluster? 
#####################################################################################



drawdegree = TRUE


plot(1:length(degreeclusteringRenyi[[lastRenyi]]),degreeclusteringRenyi[[lastRenyi]], 
     # xlim = c(0,50),ylim = c(0,30), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Random Renyi")
if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[lastRenyi]]),degreedistRenyi[[lastRenyi]]*648)}
text(x= length(degreeclusteringRenyi[[lastRenyi]]), y=mean(degreeclusteringRenyi[[lastRenyi]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[lastRenyi])
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(degreeclusteringRenyi[[j]]),degreeclusteringRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  text(x= length(degreeclusteringRenyi[[j]]), y=mean(degreeclusteringRenyi[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[j])
}

plot(1:length(degreeclusteringLattices[[lastLattices]]),degreeclusteringLattices[[lastLattices]],
     # xlim = c(0,50),ylim = c(0,50), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Regular Lattices")
if(drawdegree == TRUE){lines(1:length(degreedistLattices[[lastLattices]]),degreedistLattices[[lastLattices]]*648)}
abline(0,1)
for (i in 1:length(restLattices)){
  j <- restLattices[i]
  points(1:length(degreeclusteringLattices[[j]]),degreeclusteringLattices[[j]],col = colrestLattices[i])
  if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
  text(x= length(degreeclusteringLattices[[j]]), y=mean(degreeclusteringLattices[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringLattices)[j])
}


plot(1:length(degreeclusteringBN[[lastBN]]),degreeclusteringBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(degreeclusteringBN[[lastBN]]), y=mean(degreeclusteringBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(degreeclusteringBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(degreeclusteringBN[[j]]),degreeclusteringBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  text(x= length(degreeclusteringBN[[j]]), y=mean(degreeclusteringBN[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringBN)[j])
}

plot(1:length(degreeclusteringCN[[lastCN]]),degreeclusteringCN[[lastCN]], 
     # xlim = c(0,100),ylim = c(0,50), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "CN")
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(degreeclusteringCN[[j]]),degreeclusteringCN[[j]],col = colrestCN[i])
  
  # text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}
dev.off()


#####################################################################################
# meanstrength
#####################################################################################
plot(1:length(meanstrengthBN[[lastBN]]),meanstrengthBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(meanstrengthBN[[lastBN]]), y=mean(meanstrengthBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(numberofedgesBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(meanstrengthBN[[j]]),meanstrengthBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthBN[[j]]),meanstrengthBN[[j]]*648,col = colrestBN[i])}
  text(x= length(meanstrengthBN[[j]]), y=mean(meanstrengthBN[[j]],na.rm = TRUE), pos=4, labels= names(meanstrengthBN)[j])
}

plot(1:length(meanstrengthCN[[lastCN]]),meanstrengthCN[[lastCN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meanstrengthCN[[lastCN]]), y=mean(meanstrengthCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]]*648,col = colrestCN[i])}
  text(x= length(meanstrengthCN[[j]]), y=mean(meanstrengthCN[[j]],na.rm = TRUE), pos=4, labels= names(meanstrengthCN)[j])
}

################################################################################
# meandistances
##################################################################################
plot(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(meandistanceBN[[lastBN]]), y=mean(meandistanceBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(numberofedgesBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(meandistanceBN[[j]]),meandistanceBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meandistanceBN[[j]]),meandistanceBN[[j]]*648,col = colrestBN[i])}
  text(x= length(meandistanceBN[[j]]), y=mean(meandistanceBN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceBN)[j])
}

plot(1:length(meandistanceCN[[lastCN]]),meandistanceCN[[lastCN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meandistanceCN[[lastCN]]), y=mean(meandistanceCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meandistanceCN[[j]]),meandistanceCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meandistanceCN[[j]]),meandistanceCN[[j]]*648,col = colrestCN[i])}
  text(x= length(meandistanceCN[[j]]), y=mean(meandistanceCN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceCN)[j])
}



