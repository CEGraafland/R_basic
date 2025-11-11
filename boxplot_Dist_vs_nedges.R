#######################################################################################
# Complex measures BIC weighted BNS / cor weighted CNS
#######################################################################################
rm(list = ls())
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm1sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm2sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm3sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm4sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm5sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm1weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm2weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm3weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm4weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm5weights.rda")
library(visualizeR)
library(RColorBrewer)
library(bnlearn)
library(visualizeR)
library(grid)
library(gridExtra)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
colReds <- brewer.pal(9,"Reds")
colRedsamp <- colorRampPalette(colReds)(20)
colRainbow <- brewer.pal(9,"Set1")
colRainbowamp <- colorRampPalette(colReds)(20)
################################################
# Backpermutations BN
################################################
backpermutations <- list()
for (j in 1:length(datapermutations)){
  indback <- c()
  for (i in 1:ncol(datapermutations[[1]])){
    intback <- which(colnames(datapermutations[[j]]) == colnames(datapermutations[[1]][i]))
    indback[i] <- intback
  }
  backpermutations[[j]] <- indback
}
#######################################################################################################
# List as in Representations
#######################################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda") # <- no es gridGraph! list cn weights
length(listcnweights)
edgelistsCN <- lapply(listcnweights, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(listcnweights) <- as.character(numberofedgesCN)
#########################################################################################
# Plot CN 
#########################################################################################
distancesCN <- lapply(listcnweights, function(x) edge.attributes(x)$distances)
quantilesCN <- lapply(distancesCN, quantile)
minsCN <- lapply(quantilesCN, function(x) x[1])
firstsCN <- sapply(quantilesCN,function(x) x[2])
medsCN <- lapply(quantilesCN,function(x) x[3])
thirdsCN <- sapply(quantilesCN, function(x) x[4])
maxsCN <- sapply(quantilesCN, function(x) x[5])

par(mfrow = c(1,3))

graphics::plot(numberofedgesCN,maxsCN, type = "l", col = "lightblue")
graphics::lines(numberofedgesCN,thirdsCN,col = "blue")
graphics::lines(numberofedgesCN,medsCN,col = "black")
graphics::lines(numberofedgesCN,firstsCN,col = "blue")
graphics::lines(numberofedgesCN,minsCN,col = "lightblue")


graphics::plot(numberofedgesCN,maxsCN, xlim = c(0,30000), type = "l", col = "lightblue")
graphics::lines(numberofedgesCN,thirdsCN,col = "blue")
graphics::lines(numberofedgesCN,medsCN,col = "black")
graphics::lines(numberofedgesCN,firstsCN,col = "blue")
graphics::lines(numberofedgesCN,minsCN,col = "lightblue")
###########################################################################################
# Distances perm1weights
###########################################################################################
distancesBN <- lapply(perm1weights, function(x) edge.attributes(x)$distances)
numberofedgesBN <- sapply(perm1weights, function(x) length(E(x)))
quantilesBN <- lapply(distancesBN, quantile)
minsBN <- lapply(quantilesBN, function(x) x[1])
firstsBN <- sapply(quantilesBN,function(x) x[2])
medsBN <- lapply(quantilesBN,function(x) x[3])
thirdsBN <- sapply(quantilesBN, function(x) x[4])
maxsBN <- lapply(quantilesBN, function(x) x[5])

graphics::plot(numberofedgesBN,maxsBN, xlim = c(0,8000), ylim = c(0,20000), type = "l", col = "yellow")
graphics::lines(numberofedgesBN,thirdsBN,col = "orange")
graphics::lines(numberofedgesBN,medsBN,col = "red")
graphics::lines(numberofedgesBN,firstsBN,col = "orange")
graphics::lines(numberofedgesBN,minsBN,col = "yellow")

numberofedgesCN
join<-list(distancesBN[[6]],distancesCN[[85]],distancesBN[[14]],distancesBN[[18]], distancesBN[[32]],distancesCN[[50]],distancesCN[[35]],distancesCN[[30]],distancesCN[[3]])
boxplot.default(distancesBN[c(6,12,14,18)], y = numberofedgesBN)
boxplot.default(distancesCN[c(85,50,35,30,3)])
boxplot.default(join)

par(mfrow = c(1,2))
arrows

whoCN <- c(70,51,46,43,41,37,35,34,32,30,3)
whoBN <- c(5,10,15,18,20,25,30)

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/Distributionarclengths.pdf")
pdf(plotname)

plot(numberofedgesCN[whoCN], medsCN[whoCN], xlim = c(0,10000), ylim = c(0,20000), xlab = "Number of edges in network", ylab = "Geographical distance of edge", main = "Distribution of arc lengths in different networks")
points(numberofedgesBN[whoBN], medsBN[whoBN], col = "red")
arrows(numberofedgesCN[whoCN], thirdsCN[whoCN], numberofedgesCN[whoCN], firstsCN[whoCN], length=0.05, angle=90, code=3)
arrows(numberofedgesBN[whoBN], thirdsBN[whoBN], numberofedgesBN[whoBN], firstsBN[whoBN], length=0.05, angle=90, code=3, col = "red")
points(numberofedgesBN[whoBN],maxsBN[whoBN], xlim = c(0,8000), ylim = c(0,20000), pch = 0, col = "red")
points(numberofedgesCN[whoCN],maxsCN[whoCN], pch = 0, col = "grey")

dev.off()

plot(numberofedgesCN, medsCN, xlim = c(0,6000), ylim = c(0,5000))
points(numberofedgesBN, medsBN, col = "red")
arrows(numberofedgesCN, thirdsCN, numberofedgesCN, firstsCN, length=0.05, angle=90, code=3)
arrows(numberofedgesBN, thirdsBN, numberofedgesBN, firstsBN, length=0.05, angle=90, code=3, col = "red")

par(mfrow = c(2,1))
plot(numberofedgesBN, medsBN, xlim = c(8000,200000), ylim = c(0,10000), col = "red")
points(numberofedgesCN, medsCN, col = "black")
arrows(numberofedgesCN, thirdsCN, numberofedgesCN, firstsCN, length=0.05, angle=90, code=3)
arrows(numberofedgesBN, thirdsBN, numberofedgesBN, firstsBN, length=0.05, angle=90, code=3, col = "red")





#######################################################################################
# Same but different:
#######################################################################################
diffgraphs <- list()
diffgraphs[[1]] <- perm1weights$hc1_0_100i
for (i in 2:(length(perm1weights))){
  diffgraphs[[i]] <- difference(perm1weights[[i]], perm1weights[[i-1]])
}

xe <- seq(100,8700,100)
ye <- sapply(diffgraphs, function(x) mean(edge.attributes(x)$distances))
ye2 <- sapply(diffgraphs, function(x) max(edge.attributes(x)$distances))

plot(xe,ye2, 
     xlab = "Size of network", 
     ylab = "distance KM", 
     main = "Mean and max distance every new links",
     xlim = c(0,10000),
     ylim = c(0,20000), 
     col = "pink")
points(xe,ye, col = "red")


load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN1.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN2.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN3.rda")
load(file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphwCN4.rda")

gridGraphwCN <- c(rev(gridGraphwCN1),rev(gridGraphwCN2),rev(gridGraphwCN3),rev(gridGraphwCN4))
GraphsCN <- lapply(gridGraphwCN, function(x){x$graph}) # s = w?
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- lapply(edgelistsCN, length)
length(numberofedgesCN)
names(gridGraphwCN) <- as.character(numberofedgesCN)
edge.attributes(gridGraphwCN$`39`$graph)

diffgraphsCN <- list()
diffgraphsCN[[1]] <- GraphsCN$`39`
for (i in 2:(length(GraphsCN))){
  diffgraphsCN[[i]] <- difference(GraphsCN[[i]],GraphsCN[[i-1]])
}

xeCN <- numberofedgesCN
length(xeCN)
yeCN <- sapply(diffgraphsCN, function(x) mean(edge.attributes(x)$weight))
yeCN2 <- sapply(diffgraphsCN, function(x) max(edge.attributes(x)$weight))
length(yeCN)
points(xeCN,yeCN)
points(xeCN,yeCN2, col = "brown")
legend("bottomright",legend = c("BN mean","BN max","CN mean", "CN max"), pch = "0",col = c("red","pink","black", "brown"))


edgelistsdiffCN <- lapply(diffgraphsCN, E)
numberofedgesdiffCN <- sapply(edgelistsdiffCN, length)

###########################################################################################
# Betweenness vs networksize
###########################################################################################
# Without weights
betwBN <- lapply(perm1weights, edge.betweenness)
betwCN <- lapply(listcnweights, edge.betweenness)

par(mfrow= c(1,2))
dev.off()

plot(distancesCN$`1783` ,betwCN$`1783`, xlim = c(0,20000),ylim = c(0,60000), main = paste0("CN: ",length(betwCN$`1783`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_1700_1800i, betwBN$hc1_1700_1800i,xlim = c(0,20000),ylim = c(0,60000), main = paste0("BN: ",length(betwBN$hc1_1700_1800i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`2260` ,betwCN$`2260`, xlim = c(0,20000),ylim = c(0,60000), main = paste0("CN: ",length(betwCN$`2260`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_1700_1800i, betwBN$hc1_1700_1800i,xlim = c(0,20000),ylim = c(0,60000), main = paste0("BN: ",length(betwBN$hc1_2200_2300i)), xlab = "Edge distance",ylab =  "Edge betweenness")


plot(distancesCN$`8066`, betwCN$`8066`,xlim = c(0,20000),ylim = c(0,5000), main = paste0("CN: ",length(betwCN$`8066`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_8600_8700i, betwBN$hc1_8600_8700i,xlim = c(0,20000),ylim = c(0,5000), main = paste0("BN: ",length(betwBN$hc1_8600_8700i)), xlab = "Edge distance",ylab =  "Edge betweenness")



# With weights
bnweights <- lapply(perm1weights,function(x){edge.attributes(x)$weights})
cnweights <- lapply(listcnweights,function(x){edge.attributes(x)$strengths})
betwBN_weights <- mapply(edge.betweenness, graph = perm1weights,  weights = bnweights)
betwCN_weights <- mapply(edge.betweenness, graph = listcnweights, weights = cnweights)

par(mfrow= c(1,2))
plot(distancesCN$`1783` ,betwCN_weights$`1783`, xlim = c(0,20000), ylim = c(0,64000),main = paste0("CN: ",length(betwCN_weights$`1783`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_1700_1800i, betwBN_weights$hc1_1700_1800i,xlim = c(0,20000),ylim = c(0,64000), main = paste0("BN: ",length(betwBN_weights$hc1_1700_1800i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`2260` ,betwCN_weights$`2260`, xlim = c(0,20000),ylim = c(0,64000), main = paste0("CN: ",length(betwCN_weights$`2260`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_2200_2300i, betwBN_weights$hc1_2200_2300i,xlim = c(0,20000),ylim = c(0,64000), main = paste0("BN: ",length(betwBN_weights$hc1_2200_2300i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`2260` ,betwCN_weights$`2260`, xlim = c(0,20000),ylim = c(0,40000), main = paste0("CN: ",length(betwCN_weights$`2260`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_1700_1800i, betwBN_weights$hc1_1700_1800i,xlim = c(0,20000),ylim = c(0,40000), main = paste0("BN: ",length(betwBN_weights$hc1_1700_1800i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`8066`, betwCN_weights$`8066`,xlim = c(0,20000),ylim = c(0,5000), main = paste0("CN: ",length(betwCN_weights$`8066`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_8600_8700i, betwBN_weights$hc1_8600_8700i,xlim = c(0,20000),ylim = c(0,5000), main = paste0("BN: ",length(betwBN_weights$hc1_8600_8700i)), xlab = "Edge distance",ylab =  "Edge betweenness")

# With inverse weights
bnweightsinv <- lapply(perm1weights,function(x){1/edge.attributes(x)$weights})
cnweightsinv <- lapply(listcnweights,function(x){1/edge.attributes(x)$strengths})
betwBN_weights_inv <- mapply(edge.betweenness, graph = perm1weights,  weights = bnweightsinv)
betwCN_weights_inv <- mapply(edge.betweenness, graph = listcnweights, weights = cnweightsinv)

par(mfrow= c(1,2))
plot(distancesCN$`1783` ,betwCN_weights_inv$`1783`, xlim = c(0,20000), ylim = c(0,64000),main = paste0("CN: ",length(betwCN_weights_inv$`1783`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_1700_1800i, betwBN_weights_inv$hc1_1700_1800i,xlim = c(0,20000),ylim = c(0,64000), main = paste0("BN: ",length(betwBN_weights_inv$hc1_1700_1800i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`2260` ,betwCN_weights_inv$`2260`, xlim = c(0,20000),ylim = c(0,64000), main = paste0("CN: ",length(betwCN_weights_inv$`2260`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_2200_2300i, betwBN_weights_inv$hc1_2200_2300i,xlim = c(0,20000),ylim = c(0,64000), main = paste0("BN: ",length(betwBN_weights_inv$hc1_2200_2300i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`2260` ,betwCN_weights_inv$`2260`, xlim = c(0,20000),ylim = c(0,40000), main = paste0("CN: ",length(betwCN_weights_inv$`2260`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_1700_1800i, betwBN_weights_inv$hc1_1700_1800i,xlim = c(0,20000),ylim = c(0,40000), main = paste0("BN: ",length(betwBN_weights_inv$hc1_1700_1800i)), xlab = "Edge distance",ylab =  "Edge betweenness")

plot(distancesCN$`8066`, betwCN_weights_inv$`8066`,xlim = c(0,20000),ylim = c(0,5000), main = paste0("CN: ",length(betwCN_weights_inv$`8066`)), xlab = "Edge distance",ylab =  "Edge betweenness")
plot(distancesBN$hc1_8600_8700i, betwBN_weights_inv$hc1_8600_8700i,xlim = c(0,20000),ylim = c(0,5000), main = paste0("BN: ",length(betwBN_weights_inv$hc1_8600_8700i)), xlab = "Edge distance",ylab =  "Edge betweenness")

