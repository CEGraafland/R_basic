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
  
  degreelist <- lapply(gridGraphlist,graph2measure.degree)
  assign(paste0("degreemeas",name),degreelist)
  assign(paste0("degreedist",name), lapply(degreelist, function(x){x$degrdist}))
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



plotseqBN <- seq(1,length(knnkBN),length.out = 5)
lastBN <- plotseqBN[length(plotseqBN)]
restBN <- plotseqBN[1:(length(plotseqBN)-1)]
colrestBN <- topo.colors(length(restBN))

plotseqRenyi <- plotseqBN
lastRenyi <- lastBN
restRenyi <- restBN
colrestRenyi <- cm.colors(length(restRenyi))

plotseqLattices <- seq(1,length(knnkLattices),1)
lastLattices <- plotseqLattices[length(plotseqLattices)]
restLattices <- plotseqLattices[1:(length(plotseqLattices)-1)]
colrestLattices <- cm.colors(length(restLattices))

plotseqCN <- seq(1,length(knnkCN),length.out = 5)
plotseqCN <- round(plotseqCN)
plotseqCN <- rev(plotseqCN)
lastCN <- plotseqCN[length(plotseqCN)]
restCN <- plotseqCN[1:(length(plotseqCN)-1)]
colrestCN <- rainbow(length(restCN))


drawdegree = TRUE
par(mfrow = c(1,4))


plot(1:length(knnkRenyi[[lastRenyi]]),knnkRenyi[[lastRenyi]], 
      xlim = c(0,50),ylim = c(0,30), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Random Renyi")
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
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(knnkCN[[j]]),knnkCN[[j]],col = colrestCN[i])
  if(drawdegree == TRUE){lines(1:length(degreedistCN[[j]]),degreedistCN[[j]]*648,col = colrestCN[i])}
  text(x= length(knnkCN[[j]])+20, y=mean(knnkCN[[j]],na.rm = TRUE), pos=4, labels= names(knnCN)[j])
}
names(knnCN)[plotseqCN]
dev.off()
