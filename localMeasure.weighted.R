# OCEANO LOCAL:
rm(list = ls())
library(bnlearn)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")



##############################################################################################
# WEIGHTS ARE ARCLENGTHS.
###############################################################################################
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN1.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN2.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN3.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN4.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwBN.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwRenyi.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwLattices.rda")

gridGraphwCN <- c(rev(gridGraphwCN1),rev(gridGraphwCN2),rev(gridGraphwCN3),rev(gridGraphwCN4))
GraphsCN <- lapply(gridGraphwCN, function(x){x$graph}) # s = w?
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- lapply(edgelistsCN, length)
length(numberofedgesCN)
names(gridGraphwCN) <- as.character(numberofedgesCN)



gridGraphsBN <- gridGraphwBN
GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})
edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- sapply(edgelistsBN, length)
names(gridGraphsBN) <- as.character(numberofedgesBN)
names(gridGraphwBN) <- as.character(numberofedgesBN)

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
gridGraphlists <- list(gridGraphwCN,gridGraphwBN,gridGraphwRenyi,gridGraphwLattices)
edge_attr(gridGraphwCN$`5802`$graph)
edge_attr(gridGraphs[[40]]$graph)

a <- gridGraphwCN$`5802`$graph
b <- gridGraphs[[40]]$graph

edge_attr(gridGraphsRenyi$`100`$graph)

abarratclus <- igraph::transitivity(a, type = "barrat")
bbarratclus <- igraph::transitivity(b, type = "barrat")
all.equal(abarratclus,bbarratclus)

aglobalbarratclus <- sum(abarratclus, na.rm = TRUE)/length(abarratclus[!is.na(abarratclus)])
bglobalbarratclus <- sum(bbarratclus, na.rm = TRUE)/length(bbarratclus[!is.na(bbarratclus)])

all.equal(aglobalbarratclus,bglobalbarratclus)

# for weighted:
if(!is.null(edge_attr(graphObj$graph))){
  barratweightclus <- igraph::transitivity(graphObj$graph, type = "weighted", weights =  E(graphObj$graph)$weight)


if(!is.null(localmeasureObj)){
  lbc<- localmeasureObj$barratclustering
  globalbarratclus <- sum(lbc, na.rm = TRUE)/length(lbc[!is.na(lbc)])
} else {globalbarratclus <- NA}
# transtiivity 3 (globalweightedbarratclustering)
if(!is.null(localmeasureObj)&!is.null(edge_attr(graphObj$graph))){
  lwbc <- localmeasureObj$barratweightclustering
  globalbarratweightclus <- sum(lwbc, na.rm = TRUE)/length(lwbc[!is.na(lwbc)])


# with weights
for(i in 1:length(names)){
  name <- names[i]
  gridGraphlist <- gridGraphlists[[i]]
  locallist <- lapply(gridGraphlist, graph2measure.local)
  globallist <- mapply(graph2measure.global, graphObj = gridGraphlist, localmeasureObj = locallist, SIMPLIFY = FALSE)
  assign(paste0("globmeas",name),globallist)
  assign(paste0("globalasslist",name),lapply(globallist, function(x){x$globalass}))
  assign(paste0("globalcluslist",name),lapply(globallist, function(x){x$globalclus}))
  assign(paste0("globalbarratcluslist",name),lapply(globallist, function(x){x$globalbarratclus}))
  assign(paste0("globalbarratweightcluslist",name),lapply(globallist, function(x){x$globalbarratweightclus}))
  assign(paste0("globaledenslist",name), lapply(globallist, function(x){x$edens}))
  assign(paste0("globaldiamLCClist",name), lapply(globallist,function(x){x$diameterLCC}))
  assign(paste0("globaldiamGlist",name), lapply(globallist,function(x){x$diameterG}))
  assign(paste0("meandistLCClist",name), lapply(globallist,function(x){x$meandistLCC}))
  assign(paste0("meandistGlist",name), lapply(globallist,function(x){x$meandistG}))
  assign(paste0("meanavpathlist",name), lapply(globallist,function(x){x$meanavpathlength}))
  globallist <- NULL
  
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



plotseqBN <- seq(6,length(knnkBN),length.out = 8)
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
numberofedgesLattices[plotseqLattices]

plotseqCN <- seq(16,76,length.out = 10)
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
# text(x= length(knnkRenyi[[lastRenyi]]), y=mean(knnkRenyi[[lastRenyi]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[lastRenyi])
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  # text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}
legend(x = "topright",legend = names(knnkRenyi)[plotseqRenyi], col = c(colrestRenyi,"black"), pch = 1, cex = 0.5 )



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
  # text(x= length(knnkLattices[[j]]), y=mean(knnkLattices[[j]],na.rm = TRUE), pos=4, labels= names(knnLattices)[j])
}
legend(x = "topright",legend = names(knnkLattices)[plotseqLattices], col = c(colrestLattices,"black"), pch = 1, cex = 0.5 )

plot(1:length(knnkBN[[lastBN]]),knnkBN[[lastBN]], 
     xlim = c(0,100),ylim = c(0,60), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
# text(x= length(knnkBN[[lastBN]]), y=mean(knnkBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(knnBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(knnkBN[[j]]),knnkBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  # text(x= length(knnkBN[[j]]), y=mean(knnkBN[[j]],na.rm = TRUE), pos=4, labels= names(knnBN)[j])
}
legend(x = "topright",legend = names(knnkBN)[plotseqBN], col = c(colrestBN, "black"), pch = 1, cex = 0.5 )
length(knnkCN)
plot(1:length(knnkCN[[lastCN]]),knnkCN[[lastCN]], 
      # xlim = c(0,100),ylim = c(20,120), 
     xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(knnkCN[[j]]),knnkCN[[j]],col = colrestCN[i])
  if(drawdegree == TRUE){lines(1:length(degreedistCN[[j]]),degreedistCN[[j]]*648,col = colrestCN[i])}
  # text(x= length(knnkCN[[j]])+20, y=mean(knnkCN[[j]],na.rm = TRUE), pos=4, labels= names(knnCN)[j])
}
legend(x = "bottomright",legend = names(knnkCN)[plotseqCN], col = c(colrestCN,"black"), pch = 1, cex = 0.5 )
names(knnCN)[plotseqCN]
dev.off()

#####################################################################################
# degree clustering (weighted)
# How do vertices of degree k cluster? 
#####################################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/degreeclustweight.pdf")
pdf(plotname)
par(mfrow = c(2,2))

drawdegree = FALSE


plot(1:length(degreeclusteringRenyi[[lastRenyi]]),degreeclusteringRenyi[[lastRenyi]], 
     # xlim = c(0,50),
     ylim = c(0,0.6), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Random Renyi")
if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[lastRenyi]]),degreedistRenyi[[lastRenyi]]*648)}
# text(x= length(degreeclusteringRenyi[[lastRenyi]]), y=mean(degreeclusteringRenyi[[lastRenyi]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[lastRenyi])
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(degreeclusteringRenyi[[j]]),degreeclusteringRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  # text(x= length(degreeclusteringRenyi[[j]]), y=mean(degreeclusteringRenyi[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[j])
}
legend(x = "topright",legend = names(degreeclusteringRenyi)[plotseqRenyi], col = c(colrestRenyi,"black"), pch = 1, cex = 0.5 )

plot(1:length(degreeclusteringLattices[[lastLattices]]),degreeclusteringLattices[[lastLattices]],
     # xlim = c(0,50),
     ylim = c(0,1), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Regular Lattices")
if(drawdegree == TRUE){lines(1:length(degreedistLattices[[lastLattices]]),degreedistLattices[[lastLattices]]*648)}
for (i in 1:length(restLattices)){
  j <- restLattices[i]
  points(1:length(degreeclusteringLattices[[j]]),degreeclusteringLattices[[j]],col = colrestLattices[i])
  if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
  # text(x= length(degreeclusteringLattices[[j]]), y=mean(degreeclusteringLattices[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringLattices)[j])
}
legend(x = "topright",legend = names(degreeclusteringLattices)[plotseqLattices], col = c(colrestLattices,"black"), pch = 1, cex = 0.5 )


plot(1:length(degreeclusteringBN[[lastBN]]),degreeclusteringBN[[lastBN]], 
     xlim = c(0,100),ylim = c(0,0.6), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Bayesian Networks")
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
# text(x= length(degreeclusteringBN[[lastBN]]), y=mean(degreeclusteringBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(GraphsBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(degreeclusteringBN[[j]]),degreeclusteringBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  # text(x= length(degreeclusteringBN[[j]]), y=mean(degreeclusteringBN[[j]],na.rm = TRUE), pos=4, labels= names(GraphsBN)[j])
}
legend(x = "topright",legend = names(knnBN)[plotseqBN], col = c(colrestBN,"black"), pch = 1, cex = 0.5 )


plot(1:length(degreeclusteringCN[[lastCN]]),degreeclusteringCN[[lastCN]], 
     xlim = c(0,120),ylim = c(0,1), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Correlation Network")
for (i in 1:length(restCN)){
  i <- 1
  
  j <- restCN[i]

  points(1:length(degreeclusteringCN[[j]]),degreeclusteringCN[[j]],col = colrestCN[i])
  
  # text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}
legend(x = "topright",legend = names(degreeclusteringCN)[plotseqCN], col = c(colrestCN,"black"), pch = 1, cex = 0.5 )
dev.off()



#####################################################################################
# strength per degree (weighted)
#####################################################################################
plot(1:length(meanstrengthBN[[lastBN]]),meanstrengthBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
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
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meanstrengthCN[[lastCN]]), y=mean(meanstrengthCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]]*648,col = colrestCN[i])}
  text(x= length(meanstrengthCN[[j]]), y=mean(meanstrengthCN[[j]],na.rm = TRUE), pos=4, labels= names(meanstrengthCN)[j])
}
###########################################################################
# mean distances (= mean strength per degree)
###########################################################################

plot(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Bayesian Networks")
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
      xlim = c(0,120),ylim = c(0,12000), 
     # xlab = "degree", ylab = "k_nn,k",
     main = "Correlation Networks")
points(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]],col = "red") 
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meandistanceCN[[lastCN]]), y=mean(meandistanceCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meandistanceCN[[j]]),meandistanceCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]]*648,col = colrestCN[i])}
  text(x= length(meandistanceCN[[j]]), y=mean(meandistanceCN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceCN)[j])
}




###########################################################################
# diameter
###########################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/diameter.pdf")
pdf(plotname)
par(mfrow = c(3,3))
plot(numberofedgesBN,globaldiamLCClistBN, col = "indianred", xlab = "edges", ylab = "diameter D")
points(numberofedgesBN,globaldiamGlistBN, col = "red")
points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,globaldiamGlistRenyi, col = "blue")
points(numberofedgesLattices,globaldiamLCClistLattices, col = "lightgreen")
points(numberofedgesLattices,globaldiamGlistLattices,col = "green")
text(numberofedgesLattices,globaldiamGlistLattices, pos = 4,labels = names(GraphsLattices))
m1 <- seq(0,16,4)
text(numberofedgesBN[m1],globaldiamLCClistBN[m1], pos = 4, labels = names(GraphsBN)[m1])

plot(numberofedgesCN,globaldiamLCClistCN, col = "grey", xlab = "edges", ylab = "diameter D")
points(numberofedgesCN,globaldiamGlistCN, col = "black")
points(numberofedgesBN,globaldiamGlistBN, col = "red")
# points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,globaldiamGlistRenyi, col = "blue")
points(numberofedgesLattices,globaldiamLCClistLattices, col = "lightgreen")
points(numberofedgesLattices,globaldiamGlistLattices,col = "green")
m2 <- seq(0,60,10)
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])

plot(numberofedgesCN,globaldiamLCClistCN, col = "grey", xlab = "edges", ylab = "diameter D",xlim = c(0,10000) )
points(numberofedgesCN,globaldiamGlistCN, col = "black")
points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,globaldiamGlistRenyi, col = "blue")
points(numberofedgesBN,globaldiamLCClistBN, col = "indianred", xlab = "edges", ylab = "diameter D")
points(numberofedgesBN,globaldiamGlistBN, col = "red")
points(numberofedgesLattices,globaldiamGlistLattices,col = "green")
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
dev.off()



#################################################################################################
# Mean av path length absolute  = mean distance (absolute) for unweighted. 
#################################################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/meandistance.pdf")
pdf(plotname)
par(mfrow = c(1,3))
plot(numberofedgesBN,meanavpathlistBN, col = "indianred", xlab = "edges", ylab = "meanavpathlength D", xlim = c(0,8000), ylim = c(0,20))
points(numberofedgesBN,meandistGlistBN, col = "red")
points(numberofedgesCN,meandistGlistCN, col = "black")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "lightblue")
points(numberofedgesRenyi,meandistGlistRenyi, col = "blue")
points(numberofedgesLattices,meanavpathlistLattices, col = "lightgreen")
points(numberofedgesLattices,meandistGlistLattices,col = "green")
text(numberofedgesLattices,meandistGlistLattices, pos = 4,labels = names(GraphsLattices))
m1 <- seq(0,16,4)
text(numberofedgesBN[m1],globaldiamLCClistBN[m1], pos = 4, labels = names(GraphsBN)[m1])

plot(numberofedgesCN,meanavpathlistCN, col = "grey", xlab = "edges", ylab = "meanavpathlength D", ylim = c(0,20))
points(numberofedgesCN,meandistGlistCN, col = "black")
points(numberofedgesBN,meanavpathlistBN, col = "red")
# points(numberofedgesRenyi,globaldiamLCClistRenyi, col = "lightblue")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "blue")
points(numberofedgesLattices,meanavpathlistLattices, col = "lightgreen")
points(numberofedgesLattices,meandistGlistLattices,col = "green")
m2 <- seq(0,60,10)
text(numberofedgesCN[m2],globaldiamLCClistCN[m2], pos = 4,labels = names(GraphsCN)[m2])

plot(numberofedgesCN,meanavpathlistCN, col = "grey", xlab = "edges", ylab = "meanavpathlength D",xlim = c(0,10000))
points(numberofedgesCN,meanavpathlistCN, col = "black")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "lightblue")
points(numberofedgesRenyi,meanavpathlistRenyi, col = "blue")
points(numberofedgesBN,meanavpathlistBN, col = "indianred", xlab = "edges")
points(numberofedgesBN,meanavpathlistBN, col = "red")
points(numberofedgesLattices,meanavpathlistLattices,col = "green")
text(numberofedgesCN[m2],meanavpathlistCN[m2], pos = 4,labels = names(GraphsCN)[m2])
dev.off()


#####################################################################################
# degree clustering 
# How do vertices of degree k cluster? 
#####################################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/degreeclustweight.pdf")
pdf(plotname)
par(mfrow = c(2,2))

drawdegree = FALSE


plot(1:length(degreeclusteringRenyi[[lastRenyi]]),degreeclusteringRenyi[[lastRenyi]], 
     # xlim = c(0,50),
     ylim = c(0,0.6), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Random Renyi")
if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[lastRenyi]]),degreedistRenyi[[lastRenyi]]*648)}
# text(x= length(degreeclusteringRenyi[[lastRenyi]]), y=mean(degreeclusteringRenyi[[lastRenyi]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[lastRenyi])
for (i in 1:length(restRenyi)){
  j <- restBN[i]
  points(1:length(degreeclusteringRenyi[[j]]),degreeclusteringRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistRenyi[[j]]),degreedistRenyi[[j]]*648,col = colrestRenyi[i])}
  # text(x= length(degreeclusteringRenyi[[j]]), y=mean(degreeclusteringRenyi[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringRenyi)[j])
}
legend(x = "topright",legend = names(degreeclusteringRenyi)[plotseqRenyi], col = colrestRenyi, pch = 1, cex = 0.5 )

plot(1:length(degreeclusteringLattices[[lastLattices]]),degreeclusteringLattices[[lastLattices]],
     # xlim = c(0,50),
     ylim = c(0,1), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Regular Lattices")
if(drawdegree == TRUE){lines(1:length(degreedistLattices[[lastLattices]]),degreedistLattices[[lastLattices]]*648)}
for (i in 1:length(restLattices)){
  j <- restLattices[i]
  points(1:length(degreeclusteringLattices[[j]]),degreeclusteringLattices[[j]],col = colrestLattices[i])
  if(drawdegree == TRUE){lines(1:length(degreedistLattices[[j]]),degreedistLattices[[j]]*648,col = colrestLattices[i])}
  # text(x= length(degreeclusteringLattices[[j]]), y=mean(degreeclusteringLattices[[j]],na.rm = TRUE), pos=4, labels= names(degreeclusteringLattices)[j])
}
legend(x = "topright",legend = names(degreeclusteringLattices)[plotseqLattices], col = colrestLattices, pch = 1, cex = 0.5 )


plot(1:length(degreeclusteringBN[[lastBN]]),degreeclusteringBN[[lastBN]], 
     xlim = c(0,100),ylim = c(0,0.6), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Bayesian Networks")
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
# text(x= length(degreeclusteringBN[[lastBN]]), y=mean(degreeclusteringBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(GraphsBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(degreeclusteringBN[[j]]),degreeclusteringBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(degreedistBN[[j]]),degreedistBN[[j]]*648,col = colrestBN[i])}
  # text(x= length(degreeclusteringBN[[j]]), y=mean(degreeclusteringBN[[j]],na.rm = TRUE), pos=4, labels= names(GraphsBN)[j])
}
legend(x = "topright",legend = names(degreeclusteringBN)[plotseqBN], col = colrestBN, pch = 1, cex = 0.5 )


plot(1:length(degreeclusteringCN[[lastCN]]),degreeclusteringCN[[lastCN]], 
     xlim = c(0,120),ylim = c(0,1), 
     xlab = "degree", ylab = "C(k) Barrat",
     main = "Correlation Network")
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(degreeclusteringCN[[j]]),degreeclusteringCN[[j]],col = colrestCN[i])
  
  # text(x= length(knnkRenyi[[j]]), y=mean(knnkRenyi[[j]],na.rm = TRUE), pos=4, labels= names(knnRenyi)[j])
}
legend(x = "topright",legend = names(degreeclusteringCN)[plotseqCN], col = colrestCN, pch = 1, cex = 0.5 )
dev.off()


#####################################################################################
# strength per degree
#####################################################################################
par(mfrow = c(1,2))
plot(1:length(meanstrengthBN[[lastBN]]),meanstrengthBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     xlab = "degree k", ylab = "strength s(k)",
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
legend(x = "topright",legend = names(meanstrengthBN)[plotseqBN], col = colrestBN, pch = 1, cex = 0.5 )





plot(1:length(meanstrengthCN[[lastCN]]),meanstrengthCN[[lastCN]], 
     # xlim = c(0,100),ylim = c(0,60), 
     xlab = "degree k", ylab = "strength s(k)",
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
legend(x = "bottomright",legend = names(meanstrengthCN)[plotseqCN], col = colrestCN, pch = 1, cex = 0.5 )

###########################################################################
# mean distances (= mean strength per degree)
###########################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen3/figures/degreevsstrength.pdf")
pdf(plotname)
par(mfrow = c(1,2))
plot(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]], 
     # xlim = c(0,100),ylim = c(0,60), 
      xlab = "degree", ylab = "Mean distance d(k)",
     main = "Bayesian Networks")
abline(0,1)
if(drawdegree == TRUE){lines(1:length(degreedistBN[[lastBN]]),degreedistBN[[lastBN]]*648)}
text(x= length(meandistanceBN[[lastBN]]), y=mean(meandistanceBN[[lastBN]],na.rm = TRUE), pos=4, labels= names(numberofedgesBN)[lastBN])
for (i in 1:length(restBN)){
  j <- restBN[i]
  points(1:length(meandistanceBN[[j]]),meandistanceBN[[j]],col = colrestBN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meandistanceBN[[j]]),meandistanceBN[[j]]*648,col = colrestBN[i])}
  # text(x= length(meandistanceBN[[j]]), y=mean(meandistanceBN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceBN)[j])
}
legend(x = "bottomright",legend = names(meandistanceBN)[plotseqBN], col = c(colrestBN,"black"), pch = 1, cex = 0.5 )

plot(1:length(meandistanceCN[[lastCN]]),meandistanceCN[[lastCN]], 
     # xlim = c(0,120),ylim = c(0,12000), 
     xlab = "degree", ylab = "Mean distance d(k)",
     main = "Correlation Networks")
# points(1:length(meandistanceBN[[lastBN]]),meandistanceBN[[lastBN]],col = "red") 
if(drawdegree == TRUE){lines(1:length(degreedistCN[[lastCN]]),degreedistCN[[lastCN]]*648)}
text(x= length(meandistanceCN[[lastCN]]), y=mean(meandistanceCN[[lastCN]],na.rm = TRUE), pos=4, labels= names(numberofedgesCN)[lastCN])
for (i in 1:length(restCN)){
  j <- restCN[i]
  points(1:length(meandistanceCN[[j]]),meandistanceCN[[j]],col = colrestCN[i])
  # points(1:length(knnkRenyi[[j]]),knnkRenyi[[j]],col = colrestRenyi[i])
  if(drawdegree == TRUE){lines(1:length(meanstrengthCN[[j]]),meanstrengthCN[[j]]*648,col = colrestCN[i])}
  # text(x= length(meandistanceCN[[j]]), y=mean(meandistanceCN[[j]],na.rm = TRUE), pos=4, labels= names(meandistanceCN)[j])
}
legend(x = "bottomright",legend = names(meandistanceCN)[plotseqCN], col = c(colrestCN,"black"), pch = 1, cex = 0.5 )
dev.off()

