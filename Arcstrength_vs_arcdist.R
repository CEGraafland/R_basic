##############################################################################
# Arcstrength by distances.
##############################################################################
rm(list = ls())
library(gridExtra)
library(gridGraphics)
library(grid)
library("bnlearn")
library(maps)
library(geosphere)
library(mapproj)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/world_pol.rda")
################################################
# For MAC:
# source("/Users/lisettegraafland/Desktop/mnt/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
################################################
# Load BNs
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm1sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm2sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm3sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm4sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10/perm5sort.rda")
# Load CNs
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
gridGraphsCN <- gridGraphs
rm(gridGraphs)
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(gridGraphsCN) <- as.character(numberofedgesCN)
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
# check backpermutations
proefDAG <- perm2sort$hc2_2000_2100i
datapermutations[[2]][,backpermutations[[2]]]
nodes(proefDAG)[backpermutations[[2]]]
arcs(proefDAG)
dag2 <-proefDAG
nodes(dag2) <- nodes(proefDAG)[backpermutations[[2]]]
arcs(dag2)
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
plot_long_strong_distances(dag2,data.dag = data, minimdist = 5000, which = "strongest", perm = permutations[[2]])


plot_long_distances(proefDAG, data, minimdist = 10000, smallcol = NA, perm = permutations[[2]])
plot_quantile_edges(proefDAG,data)
plot_long_strong_distances(proefDAG,data,minimdist = 10000, which = "normal", smallcol = NA, perm = permutations[[2]], k = NULL)
################################################################
# Strongest edges visualized BN All permutations
################################################################
size <- 30
sizename <- paste0(size,"00")
proefDAGs <- list(perm1sort[[size]], perm2sort[[size]],perm3sort[[size]],perm4sort[[size]],perm5sort[[size]])
dist <- 10000
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/tele_",dist,"_",sizename,".pdf")
pdf(plotname,height = 5, width = 10)
par(mfrow = c(2,3))
for (i in 1:5){
  dag <- proefDAGs[[i]]
  plot_long_strong_distances(dag,data,minimdist = dist, smallcol = NA, perm = permutations[[i]], k = NULL, title = "B", which = "strongest")
  }
dev.off()

########################################################################################
# Plots strengths/weights/zooms vs distances BN
########################################################################################
num1 <- 16
igraph.str <- bn_to_igraph.strengths(perm2sort[[num1]],perm = permutations[[2]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
igraph.str.dist <- igraph.distances(igraph = igraph.str, data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d),perm = permutations[[2]])

par(mfrow = c(2,2))
igraphmod <- igraph.weights(igraph = igraph.str.dist,type = "bn",fromdist = 0)
plotind <- which(E(igraphmod)$largeweights != 0)
plot(E(igraphmod)$distances[plotind],E(igraphmod)$largeweights[plotind],
     xlab = "distance", ylab = "weights",
     main = paste(length(E(igraphmod))," from  = ",min(E(igraphmod)$distances[plotind])),
     xlim = c(0,19000))

plot(E(igraphmod)$distances[plotind],E(igraphmod)$strengths[plotind],
     xlab = "distance", ylab = "strengths",
     main = paste(length(E(igraphmod))," from  = ",min(E(igraphmod)$distances[plotind])),
     xlim = c(0,19000), ylim = c(-1200,0))

igraphmod <- igraph.weights(igraph = igraph.str.dist, type = "bn",fromdist = 2000)
plotind <- which(E(igraphmod)$largeweights != 0)
plot(E(igraphmod)$distances[plotind],E(igraphmod)$largeweights[plotind],
     xlab = "distance", ylab = "weights",
     main = paste(length(E(igraphmod))," from  = ",min(E(igraphmod)$distances[plotind])),
     xlim = c(0,19000))

plot(E(igraphmod)$distances[plotind],E(igraphmod)$strengths[plotind],
     xlab = "distance", ylab = "strengths",
     main = paste(length(E(igraphmod))," from  = ",min(E(igraphmod)$distances[plotind])),
     xlim = c(0,19000), ylim = c(-1200,0))

plot_long_strong_distances(perm2sort[[16]], TimeCoordsAnom_from_Grid_rms(tas_ncep_10d), which = "normal", smallcol = NA,
                           minimdist = 10000,perm = permutations[[2]])

###########################################################################
# plot weights vs distances five permutations same size networks BN
###########################################################################
which <- 16
listselect <- list(perm1sort[[which]],perm2sort[[which]],perm3sort[[which]],perm4sort[[which]],perm5sort[[which]])
arcs <- lapply(listselect, arcs)
narcs <- lapply(arcs, nrow)
list <-mapply(str_vs_dis,dag = listselect,perm = permutations, MoreArgs = list(data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)), SIMPLIFY = FALSE)
par(mfrow = c(2,length(list)))
for (i in 1:length(list)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]))
}

############################################################################
# plot strengths in weights vs distances one permutation different sizes. BN
############################################################################

listselect <- perm4sort[6:25]
arcs <- lapply(listselect, arcs)
narcs <- sapply(arcs, nrow)
list <- lapply(listselect, str_vs_dis, perm = permutations[[4]], data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/strvsarc_bn_",min(narcs),"_",max(narcs),"_1.pdf")
pdf(plotname)
par(mfrow = c(2,length(list)/4))
for (i in 1:(length(list)/2)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]))
}
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/strvsarc_bn_",min(narcs),"_",max(narcs),"_2.pdf")
pdf(plotname)
par(mfrow = c(2,length(list)/4))
for (i in (length(list)/2+1):length(list)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]))
}
dev.off()
##########################################################################
# for Resumen 4 BN
############################################################################

listselectB1 <- perm1sort[c(5,10,15,20,25,30)]
arcsB1 <- lapply(listselectB1, arcs)
narcsB1 <- sapply(arcsB1, nrow)
listB1 <- lapply(listselectB1, str_vs_dis, perm = permutations[[1]], data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_bn_",min(narcsB1),"_",max(narcsB1),"_1.pdf")
pdf(plotname)
par(mfrow = c(2,length(listB1)/2))
for (i in 1:length(listB1)){
  graph <- listB1[[i]]
  plot(E(graph)$distances,(E(graph)$strengths)*-1,
       xlab = "distance",
       ylab = "+ BIC",
       main = as.character(narcsB1[[i]]),
       xlim = c(0,20000))
}
dev.off()





##########################################################################
# plot strengths in positive arcstrengths vs distances one permutation different sizes. BN
##########################################################################

listselect <- perm1sort[6:25]
arcs <- lapply(listselect, arcs)
narcs <- sapply(arcs, nrow)
list <- lapply(listselect, str_vs_dis, perm = permutations[[1]], data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/pos.arc.str.vsarc_bn_",min(narcs),"_",max(narcs),"_1.pdf")
pdf(plotname)
par(mfrow = c(2,length(list)/4))
for (i in 1:(length(list)/2)){
  graph <- list[[i]]
  plot(E(graph)$distances,(E(graph)$strengths)*-1,
       xlab = "distance",
       ylab = "+ BIC",
       main = as.character(narcs[[i]]),
       ylim = c(0,1200))
}
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/pos.arc.str.vsarc_bn_",min(narcs),"_",max(narcs),"_2.pdf")
pdf(plotname)
par(mfrow = c(2,length(list)/4))
for (i in (length(list)/2+1):length(list)){
  graph <- list[[i]]
  plot(E(graph)$distances,(E(graph)$strengths)*-1,
       xlab = "distance",
       ylab = "+ BIC",
       main = as.character(narcs[[i]]),
      ylim = c(0,1200))
}
dev.off()


###########################################################################
# plot five permutations distances vs weights same size networks BN
###########################################################################
which <- 16
listselect <- list(perm1sort[[which]],perm2sort[[which]],perm3sort[[which]],perm4sort[[which]],perm5sort[[which]])
arcs <- lapply(listselect, arcs)
narcs <- lapply(arcs, nrow)
list <-mapply(str_vs_dis,dag = listselect,perm = permutations, MoreArgs = list(data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)), SIMPLIFY = FALSE)
par(mfrow = c(2,length(list)))
for (i in 1:length(list)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]))
}

proef <- perm1sort[[which]]
valuesproef <- str_vs_dis(proef, data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
E(valuesproef)$strengths
E(valuesproef)$distances
larges <- which(E(valuesproef)$distances>2000)
plot(E(valuesproef)$distances[larges], E(valuesproef)$strengths[larges])

############################################################################
# Now the same voor correlation graphs.
############################################################################
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)

test <- str_vs_dis_2(GraphsCN[[50]],data)
plot(E(test)$distances,E(test)$strengths,
     xlab = "distance",
     ylab = "correlations",
     main = as.character(numberofedgesCN[50]))



listselect <- GraphsCN[41:70]
narcs <- numberofedgesCN[41:70]
list <- lapply(listselect, str_vs_dis_2, perm = permutations[[1]], data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/corvsarc_cn_",min(narcs),"_",max(narcs),"1.pdf")
pdf(plotname)
par(mfrow = c(2,length(list)/6))
for (i in 1:(length(list)/3)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$strengths,
       xlab = "distance",
       ylab = "correlation",
        main = as.character(narcs[[i]]),
       ylim = c(0,1))
}
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/corvsarc_cn_",min(narcs),"_",max(narcs),"2.pdf")
pdf(plotname)
par(mfrow = c(2,length(list)/6))
for (i in (length(list)/3+1):(2*length(list)/3)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$strengths,
       xlab = "distance",
       ylab = "correlation",
       main = as.character(narcs[[i]]),
       ylim = c(0,1))
}
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/corvsarc_cn_",min(narcs),"_",max(narcs),"3.pdf")
pdf(plotname)
for (i in (2*length(list)/3+1):length(list)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$strengths,
       xlab = "distance",
       ylab = "correlation",
       main = as.character(narcs[[i]]),
       ylim = c(0,1))
}
dev.off()


list <- lapply(list, igraph.weights, type = "cn", perm = permutations[[1]])
par(mfrow = c(2,length(list)/6))
for (i in 1:(length(list)/3)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]),
       ylim = c(0,1))
}
par(mfrow = c(2,length(list)/6))
for (i in (length(list)/3+1):(2*length(list)/3)){
  graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]),
       ylim = c(0,1))
}
for (i in (2*length(list)/3+1):length(list)){

    graph <- list[[i]]
  plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs[[i]]),
       ylim = c(0,1))
}

##############################################################################
# CN for resumen 4
##############################################################################
listselect2 <- GraphsCN[c(70,61,56,51,46,41)]
narcs2 <- numberofedgesCN[c(70,61,56,51,46,41)]
list2 <- lapply(listselect2, str_vs_dis_2, perm = permutations[[1]], data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/corvsarc_cn_",min(narcs2),"_",max(narcs2),".pdf")
pdf(plotname)
par(mfrow = c(2,length(list2)/2))
for (i in 1:(length(list2))){
  graph <- list2[[i]]
  plot(E(graph)$distances,E(graph)$strengths,
       xlab = "distance",
       ylab = "correlation",
       main = as.character(narcs2[[i]]),
       xlim = c(0,20000),
       ylim = c(0,1))
}
dev.off()



###############################################################################################
# To understand difference: correlation network and Bayesian network
# Plot network 
###############################################################################################
par(mfrow = c(1,2))
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
num <- 50
cornet <- GraphsCN[[num]]
cn.str.dist <- str_vs_dis_2(cornet,data) 
cn.str <- cn.strengths(cornet,data)
# cn.str.dist <- igraph.distances(cn.str, data.igraph = data, perm = NULL)
cn.str.dist.weight <- igraph.weights(cn.str.dist, type = "cn", fromdist = 3000, perm = permutations[[1]])
E(cn.str.dist.weight)$largeweights

plotgraph <- cn.str.dist.weight
plot(E(plotgraph)$distances,E(plotgraph)$strengths,
     xlab = "distance",
     ylab = "correlation",
     main = as.character(numberofedgesCN[num]),
     ylim = c(0,1))

plotind <- which(E(plotgraph)$largeweights != 0)
plot(E(plotgraph)$distances[plotind],E(plotgraph)$weights[plotind],
     xlab = "distance", ylab = "weights",
     main = paste(length(E(plotgraph))," from  = ",min(E(plotgraph)$distances[plotind])),
     xlim = c(0,19000))

par(mfrow = c(1,1))
perm = NULL
plot.spatial.igraph(network = cn.str.dist.weight, th = 0.3, by.dist = 5000, data.network = data, perm = NULL)
plotS.spatial.igraph(network = cn.str.dist.weight, data.network = data, type = "cn", th.type = "zoomweighted", th = 0.3, from.dist = 3000, by.dist = 5000, perm = NULL, remove = TRUE, shift = FALSE)
plotS.spatial.igraph(network = cn.str.dist.weight, data.network = data, type = "cn", th.type = "zoomweighted", th = 0.3, from.dist = 3000, by.dist = 5000, perm = NULL, remove = FALSE, shift = TRUE,curvature = 0.01)
# plot_long_distances(cn.str.dist.weight, data, 10000, smallcol = NA ,perm = NULL, title = "NA", remove = TRUE)
par(mfrow = c(1,2))

dag <- perm1sort$hc1_1500_1600i
perm <- permutations[[1]]
bn.str <- bn_to_igraph.strengths(dag,data.dag = data, perm = perm)
bn.str.dist <- igraph.distances(bn.str,data.igraph = data, perm = perm)
bn.str.dist.weight <- igraph.weights(bn.str.dist, type = "bn", fromdist = 3000, perm)
graph <- bn.str.dist.weight
plot(E(graph)$distances,E(graph)$weights,
       xlab = "distance",
       ylab = "weight",
       main = as.character(narcs(dag)),
      xlim = c(0,19000))

plotind <- which(E(graph)$largeweights != 0)
plot(E(graph)$distances[plotind],E(graph)$largeweights[plotind],
     xlab = "distance", ylab = "weights",
     main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
     xlim = c(3000,20000))
axis(1, c(3000), labels=c(3000))

par(mfrow = c(1,1))
plot.spatial.igraph(network = graph, th = 0.35, by.dist = 5000, data.network = data, perm = perm)
plotS.spatial.igraph(network = graph, data.network = data, type = "bn", th = 0.35, th.type = 'zoomweighted', by.dist = 5000, perm = perm, remove = FALSE, shift = TRUE)

title(paste0("th = ",th,"_>",by.dist))



###############################################################################################
# To understand difference: correlation network and Bayesian network
# Plot network ZOOM small Resumen 4
###############################################################################################
num <- 65
sizename <- as.character(numberofedgesCN[num])
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_zoom_cn",sizename,".pdf")
pdf(plotname)

par(mfrow = c(2,1))
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cornet <- GraphsCN[[num]]
cn.str.dist <- str_vs_dis_2(cornet,data) 
plotgraph <- cn.str.dist

plot(E(plotgraph)$distances,E(plotgraph)$strengths,
     xlab = "distance",
     ylab = "correlation",
     main = as.character(numberofedgesCN[num]),
     xlim = c(0,20000),
     ylim = c(0,1))

fromdist <- 3000
by.dist <- 5000
th <- 0.5
plotind <- which(E(plotgraph)$distances > fromdist)
col <- character(length(E(plotgraph)[plotind]))
indblack <- which(E(plotgraph)[plotind]$strengths <= th)
indpurple <- which(E(plotgraph)[plotind]$distances >  by.dist & E(plotgraph)[plotind]$strengths > th)
indyellow <- which(E(plotgraph)[plotind]$distances <= by.dist & E(plotgraph)[plotind]$strengths > th)
col[indblack] <- "black"
col[indpurple] <- "purple"
col[indyellow] <- "orange"
plot(E(plotgraph)$distances[plotind],E(plotgraph)$strengths[plotind],
     xlab = "distance", ylab = "correlation",
     main = paste(length(E(plotgraph))," from  = ",min(E(plotgraph)$distances[plotind])),
     xlim = c(fromdist,20000),
     col = col)

dev.off() 

par(mfrow = c(1,1))
from.distC = 3000
thC = 0.5
by.distC = 5000
# plotS.spatial.igraph(network = cn.str.dist.weight, data.network = data, type = "cn", th.type = "zoomweighted", th = 0.1, by.dist = 5000, perm = perm, remove = FALSE, shift = TRUE)
# plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = FALSE)
plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = TRUE)
title(paste0("|size| = ",sizename, " weight > ",thC," purple >",by.dist))
# plot_long_distances(cn.str.dist.weight, data, 10000, smallcol = NA ,perm = NULL, title = "NA", remove = TRUE)


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_zoom_bn",sizename,".pdf")
pdf(plotname)

par(mfrow = c(2,1))
dag <- perm1sort$hc1_1400_1500i
perm <- permutations[[1]]
nodes(dag)
bn.str <- bn_to_igraph.strengths(dag,data.dag = data, perm = perm)
bn.str.dist <- igraph.distances(bn.str,data.igraph = data, perm = perm)
bn.str.dist.weight <- igraph.weights(bn.str.dist, type = "bn", fromdist = 3000, perm)
graph <- bn.str.dist.weight
plot(E(graph)$distances,E(graph)$weights,
     xlab = "distance",
     ylab = "weight",
     main = as.character(narcs(dag)))



fromdist <- 3000
by.dist <- 5000
th <- 0.035
plotind <- which(E(graph)$distances > fromdist)
col <- character(length(E(graph)[plotind]))
indblack <- which(E(graph)[plotind]$weights <= th)
indpurple <- which(E(graph)[plotind]$distances >  by.dist & E(graph)[plotind]$weights > th)
indyellow <- which(E(graph)[plotind]$distances <= by.dist & E(graph)[plotind]$weights > th)
col[indblack] <- "black"
col[indpurple] <- "purple"
col[indyellow] <- "orange"
plot(E(graph)$distances[plotind],E(graph)$weights[plotind],
     xlab = "distance", ylab = "weight",
     main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
     xlim = c(fromdist,20000),
     col = col)

# plotind <- which(E(graph)$largeweights != 0)
# plot(E(graph)$distances[plotind],E(graph)$largeweights[plotind],
#      xlab = "distance", ylab = "weights",
#      main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
#      xlim = c(0,19000))

dev.off()

par(mfrow = c(1,1))
plotS.spatial.igraph(network = graph, th = 0.04, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = FALSE, shift = TRUE)
# plotS.spatial.igraph(network = graph, th = 0.035, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = TRUE, shift = FALSE)
title(paste0("|size| = ",length(E(graph)), " weight > ",th," purple >",by.dist))
dev.off()


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_plots_bn",sizename,".pdf")
pdf(plotname)

plotS.spatial.igraph(network = graph, th = 0.04, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = FALSE, shift = TRUE)
# title(paste0("|size| = ",length(E(graph)), " weight > ",th," purple >",by.dist))
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_plots_cn",sizename,".pdf")
pdf(plotname)
plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = TRUE)
# title(paste0("|size| = ",sizename, " weight > ",thC," purple >",by.dist))

dev.off()



###############################################################################################
# To understand difference: correlation network and Bayesian network
# Plot network ZOOM big
###############################################################################################
num <- 42
sizename <- as.character(numberofedgesCN[num])
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/strvsarc_zoom_cnbn",sizename,".pdf")
pdf(plotname)

par(mfrow = c(1,2))
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cornet <- GraphsCN[[num]]
cn.str.dist <- str_vs_dis_2(cornet,data) 
plotgraph <- cn.str.dist

plot(E(plotgraph)$distances,E(plotgraph)$strengths,
     xlab = "distance",
     ylab = "correlation",
     main = as.character(numberofedgesCN[num]),
     xlim = c(0,20000),
     ylim = c(0,1))

fromdist <- 3000
by.dist <- 5000
th <- 0.5
plotind <- which(E(plotgraph)$distances > fromdist)
col <- character(length(E(plotgraph)[plotind]))
indblack <- which(E(plotgraph)[plotind]$strengths <= th)
indpurple <- which(E(plotgraph)[plotind]$distances >  by.dist & E(plotgraph)[plotind]$strengths > th)
indyellow <- which(E(plotgraph)[plotind]$distances <= by.dist & E(plotgraph)[plotind]$strengths > th)
col[indblack] <- "black"
col[indpurple] <- "purple"
col[indyellow] <- "orange"
plot(E(plotgraph)$distances[plotind],E(plotgraph)$strengths[plotind],
     xlab = "distance", ylab = "correlation",
     main = paste(length(E(plotgraph))," from  = ",min(E(plotgraph)$distances[plotind])),
     xlim = c(fromdist,20000),
     col = col)

par(mfrow = c(1,1))
from.distC = 3000
thC = 0.5
by.distC = 5000
# plotS.spatial.igraph(network = cn.str.dist.weight, data.network = data, type = "cn", th.type = "zoomweighted", th = 0.1, by.dist = 5000, perm = perm, remove = FALSE, shift = TRUE)
# plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = FALSE)
plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = TRUE)
title(paste0("|size| = ",sizename, " weight > ",thC," purple >",by.dist))
# plot_long_distances(cn.str.dist.weight, data, 10000, smallcol = NA ,perm = NULL, title = "NA", remove = TRUE)
par(mfrow = c(1,2))

dag <- perm1sort$hc1_5000_5100i
perm <- permutations[[1]]
nodes(dag)
bn.str <- bn_to_igraph.strengths(dag,data.dag = data, perm = perm)
bn.str.dist <- igraph.distances(bn.str,data.igraph = data, perm = perm)
bn.str.dist.weight <- igraph.weights(bn.str.dist, type = "bn", fromdist = 3000, perm)
graph <- bn.str.dist.weight
plot(E(graph)$distances,E(graph)$weights,
     xlab = "distance",
     ylab = "weight",
     main = as.character(narcs(dag)))



fromdist <- 3000
by.dist <- 5000
th <- 0.035
plotind <- which(E(graph)$distances > fromdist)
col <- character(length(E(graph)[plotind]))
indblack <- which(E(graph)[plotind]$weights <= th)
indpurple <- which(E(graph)[plotind]$distances >  by.dist & E(graph)[plotind]$weights > th)
indyellow <- which(E(graph)[plotind]$distances <= by.dist & E(graph)[plotind]$weights > th)
col[indblack] <- "black"
col[indpurple] <- "purple"
col[indyellow] <- "orange"
plot(E(graph)$distances[plotind],E(graph)$weights[plotind],
     xlab = "distance", ylab = "weight",
     main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
     xlim = c(fromdist,20000),
     col = col)

# plotind <- which(E(graph)$largeweights != 0)
# plot(E(graph)$distances[plotind],E(graph)$largeweights[plotind],
#      xlab = "distance", ylab = "weights",
#      main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
#      xlim = c(0,19000))


par(mfrow = c(1,1))
plotS.spatial.igraph(network = graph, th = 0.035, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = FALSE, shift = TRUE)
# plotS.spatial.igraph(network = graph, th = 0.035, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = TRUE, shift = FALSE)
title(paste0("|size| = ",length(E(graph)), " weight > ",th," purple >",by.dist))
dev.off()


plotS.spatial.igraph(network = cn.str.dist.weight, data.network = data, type = "cn", th.type = "zoom", from.dist = 3000, th = 0.6, by.dist = 5000, perm = perm, remove = FALSE, shift = TRUE)
title(paste0("|Param| = ",length(E(cn.str.dist.weight)), " weight > ",th," purple >",by.dist))


###############################################################################################
# To understand difference: correlation network and Bayesian network
# Plot network ZOOM big Resumen 4
###############################################################################################
num <- 42
sizename <- as.character(numberofedgesCN[num])
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_zoom_cn",sizename,".pdf")
pdf(plotname)

par(mfrow = c(2,1))
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cornet <- GraphsCN[[num]]
cn.str.dist <- str_vs_dis_2(cornet,data) 
plotgraph <- cn.str.dist

plot(E(plotgraph)$distances,E(plotgraph)$strengths,
     xlab = "distance",
     ylab = "correlation",
     main = as.character(numberofedgesCN[num]),
     xlim = c(0,20000),
     ylim = c(0,1))

fromdist <- 3000
by.dist <- 5000
th <- 0.5
plotind <- which(E(plotgraph)$distances > fromdist)
col <- character(length(E(plotgraph)[plotind]))
indblack <- which(E(plotgraph)[plotind]$strengths <= th)
indpurple <- which(E(plotgraph)[plotind]$distances >  by.dist & E(plotgraph)[plotind]$strengths > th)
indyellow <- which(E(plotgraph)[plotind]$distances <= by.dist & E(plotgraph)[plotind]$strengths > th)
col[indblack] <- "black"
col[indpurple] <- "purple"
col[indyellow] <- "orange"
plot(E(plotgraph)$distances[plotind],E(plotgraph)$strengths[plotind],
     xlab = "distance", ylab = "correlation",
     main = paste(length(E(plotgraph))," from  = ",min(E(plotgraph)$distances[plotind])),
     xlim = c(fromdist,20000),
     col = col)

dev.off() 

par(mfrow = c(1,1))
from.distC = 3000
thC = 0.5
by.distC = 5000
# plotS.spatial.igraph(network = cn.str.dist.weight, data.network = data, type = "cn", th.type = "zoomweighted", th = 0.1, by.dist = 5000, perm = perm, remove = FALSE, shift = TRUE)
# plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = FALSE)
plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = TRUE)
title(paste0("|size| = ",sizename, " weight > ",thC," purple >",by.dist))
# plot_long_distances(cn.str.dist.weight, data, 10000, smallcol = NA ,perm = NULL, title = "NA", remove = TRUE)


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_zoom_bn",sizename,".pdf")
pdf(plotname)

par(mfrow = c(2,1))
dag <- perm1sort$hc1_5000_5100i
perm <- permutations[[1]]
nodes(dag)
bn.str <- bn_to_igraph.strengths(dag,data.dag = data, perm = perm)
bn.str.dist <- igraph.distances(bn.str,data.igraph = data, perm = perm)
bn.str.dist.weight <- igraph.weights(bn.str.dist, type = "bn", fromdist = 3000, perm)
graph <- bn.str.dist.weight
plot(E(graph)$distances,E(graph)$weights,
     xlab = "distance",
     ylab = "weight",
     main = as.character(narcs(dag)))



fromdist <- 3000
by.dist <- 5000
th <- 0.032
plotind <- which(E(graph)$distances > fromdist)
col <- character(length(E(graph)[plotind]))
indblack <- which(E(graph)[plotind]$weights <= th)
indpurple <- which(E(graph)[plotind]$distances >  by.dist & E(graph)[plotind]$weights > th)
indyellow <- which(E(graph)[plotind]$distances <= by.dist & E(graph)[plotind]$weights > th)
col[indblack] <- "black"
col[indpurple] <- "purple"
col[indyellow] <- "orange"
plot(E(graph)$distances[plotind],E(graph)$weights[plotind],
     xlab = "distance", ylab = "weight",
     main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
     xlim = c(fromdist,20000),
     col = col)

# plotind <- which(E(graph)$largeweights != 0)
# plot(E(graph)$distances[plotind],E(graph)$largeweights[plotind],
#      xlab = "distance", ylab = "weights",
#      main = paste(length(E(graph))," from  = ",min(E(graph)$distances[plotind])),
#      xlim = c(0,19000))

dev.off()

par(mfrow = c(1,1))
plotS.spatial.igraph(network = graph, th = .032, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = FALSE, shift = TRUE)
# plotS.spatial.igraph(network = graph, th = 0.035, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = TRUE, shift = FALSE)
title(paste0("|size| = ",length(E(graph)), " weight > ",th," purple >",by.dist))
dev.off()


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_plots_bn",sizename,".pdf")
pdf(plotname)

plotS.spatial.igraph(network = graph, th = 0.032, type = "bn", th.type ="zoom", from.dist = 3000, by.dist = 5000, data.network = data, perm = perm, remove = FALSE, shift = TRUE)
title(paste0("|size| = ",length(E(graph)), " weight > ",th," purple >",by.dist))
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/strvsarc_plots_cn",sizename,".pdf")
pdf(plotname)
plotS.spatial.igraph(network = plotgraph, data.network = data, type = "cn", th.type = "zoom", from.dist = from.distC, th = thC, by.dist = by.distC, perm = perm, remove = FALSE, shift = TRUE)
# title(paste0("|size| = ",sizename, " weight > ",thC," purple >",by.dist))

dev.off()


##################################################################
# 
###################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/strvsarc_cnbn",sizename,".pdf")
pdf(plotname)
# 1100
par(mfrow = c(2,2))
size <- 11
sizename <- paste0(size,"00")
cornet <- GraphsCN[[75]]
nedgescornet <- numberofedgesCN[75]
baynet <- perm2sort$hc2_1000_1100i
perm <- permutations[[2]]
# 1600
par(mfrow = c(2,2))
size <- 16
sizename <- paste0(size,"00")
cornet <- GraphsCN[[63]]
nedgescornet <- numberofedgesCN[63]
baynet <- perm2sort$hc2_1500_1600i
perm <- permutations[[2]]
# 3000
par(mfrow = c(2,2))
size <- 30
sizename <- paste0(size,"00")
cornet <- GraphsCN[[51]]
nedgescornet <- numberofedgesCN[51]
baynet <- perm2sort$hc2_2900_3000i
perm <- permutations[[2]]


# 8000
par(mfrow = c(2,2))
size <- 80
sizename <- paste0(size,"00")
cornet <- GraphsCN[[35]]
nedgescornet <- numberofedgesCN[35]
baynet <- perm2sort$hc2_8500_8600i
perm <- permutations[[2]]



cornetplus <- str_vs_dis_2(cornet,data)
baynetplus <- str_vs_dis(dag,perm = perm,data.dag = data)
# de sterkste:
E(baynetplus)[E(baynetplus)$strengths == min(E(baynetplus)$strengths)]
# zijn lengthe:
E(baynetplus)[E(baynetplus)$strengths == min(E(baynetplus)$strengths)]$distances

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Arcs_lengthsandstrengths/strvsarc_cnbn",sizename,".pdf")
pdf(plotname)
par(mfrow = c(2,2))
plot(E(cornetplus)$distances,E(cornetplus)$weights,
     xlab = "distance",
     ylab = "weight",
     main = as.character(nedgescornet),
    xlim = c(0,20000),
    ylim = c(0,1))



plot(E(baynetplus)$distances,E(baynetplus)$weights,
     xlab = "distance",
     ylab = "weight",
     main = as.character(narcs(baynet)))



plot(E(cornetplus)$distances,E(cornetplus)$strenghts,
     xlab = "distance",
     ylab = "correlation",
     main = as.character(nedgescornet),
     xlim = c(0,20000),
     ylim = c(0,1))

plot(E(baynetplus)$distances,E(baynetplus)$strengths,
     xlab = "distance",
     ylab = "arcstrength",
     main = as.character(narcs(baynet)))

dev.off()

#####################################################################
# high betweenness plot. 
#####################################################################
time.coords <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
meteodag <- listcnweights$`8066`
shift <- TRUE
weights <- TRUE
howmany <- 20

plot.edgebetw <- function(time.coords, meteodag, shift, howmany, weights = FALSE, strengths = FALSE){
  if (class(meteodag) == "igraph") {
    meteodag <- meteodag}
  
  if (class(meteodag) == "graphNEL"){ 
    meteodag <- igraph.from.graphNEL(meteodag)}
  
  if (class(meteodag) == "bn") {
    meteodag <- igraph.from.graphNEL(as.graphNEL(meteodag))}
  
  if(strengths == TRUE){
    E(meteodag)$betwe <- edge.betweenness(meteodag,weights = E(meteodag)$strengths)
  } else if (weights == TRUE){
    E(meteodag)$betwe <- edge.betweenness(meteodag,weights = E(meteodag)$weights)
    } else {E(meteodag)$betwe <- edge.betweenness(meteodag)}
  
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)
  if(shift == TRUE){x <- (x+360)%%360}
  
  points <- expand.grid(y, x)[2:1]
  
  betw_high2low <- sort(E(meteodag)$betwe, decreasing = TRUE)
  # howmany <- 20
  cutbetw <- betw_high2low[howmany]
  highbetw <- which(E(meteodag)$betwe >= cutbetw)
  # upp <- quantile(E(meteodag)$betwe,probs = seq(0,1,0.01))[99]
  E(meteodag)$color <- "NA"
  # E(meteodag)[E(meteodag)$betwe>upp]$color <- "red"
  E(meteodag)[highbetw]$color <- "red"

  if(shift == TRUE){map('world2',interior = FALSE,resolution = 0)}
  else {plot(wrld)}

  plot.igraph(meteodag,
            vertex.size = 100,
            vertex.color = NA,
            vertex.label = NA,
            # edge.curved = curvature,
            edge.arrow.size = 0.2,
            # edge.width = E(igraph)$width,
            lty = "dots",
            layout = as.matrix(points), add = TRUE, rescale = FALSE)
}


tc <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms =  TRUE)
meteodag <- listcnweights$`8066`
E(meteodag)$betwe <- edge.betweenness(meteodag,weights = E(meteodag)$strengths)
meteodag
meteodag <- perm1weights$hc1_1700_1800i
E(meteodag)$betwe <- edge.betweenness(meteodag,weights = E(meteodag)$weights)
lis <- FALSE
shift <- TRUE
TRUEedgecolor = "green"

plot.edgebetw(tc,listcnweights$`8066`,shift = TRUE,howmany = 50, weights = FALSE, strengths = TRUE )
plot.edgebetw(tc,listcnweights$`2260`,shift = TRUE,howmany = 50, weights = FALSE, strengths = TRUE )
plot.edgebetw(tc,perm1weights$hc1_1700_1800i,shift = TRUE,howmany = 50, weights = TRUE, strengths = FALSE )
