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
###################################################################################################
# (Same code is used in comunidades. )
# Use weighted igraphs: BN:
###################################################################################################
arcs1 <- lapply(perm1sort, arcs)
narcs1 <- sapply(arcs1, nrow)
nedgeslistbn <- narcs1[c(18,26,50,86)]

bnsel1 <- 18
bnsel2 <- 26
bnsel3 <- 50
bnsel4 <- 86


bn1 <- perm1weights[[bnsel1]]
bn2 <- perm1weights[[bnsel2]]
bn3 <- perm1weights[[bnsel3]]
bn4 <- perm1weights[[bnsel4]]
bnlist <- list(bn1,bn2,bn3,bn4)
#######################################################################################################
# List as in Representations
#######################################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda") # <- no es gridGraph! list cn weights
length(listcnweights)
edgelistsCN <- lapply(listcnweights, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(listcnweights) <- as.character(numberofedgesCN)
numberofedgesCN[[25]]

# cn:
nedgescnlist <- numberofedgesCN[c(61,56,35,32,1)]
numbersCN <- as.character(nedgescnlist)
cn1 <- listcnweights[[61]]
cn2 <- listcnweights[[56]]
cn3 <- listcnweights[[35]]
cn4 <- listcnweights[[32]]
cn00 <- listcnweights[[3]]
cnlist <- list(cn1,cn2,cn3,cn4,cn00)




nedgescnlist <- numberofedgesCN[c(69,42,35,32,3)]

cn1 <- listcnweights[[69]]
cn2 <- listcnweights[[42]]
cn3 <- listcnweights[[35]]
cn4 <- listcnweights[[32]]
cn00 <- listcnweights[[3]]

cnlist <- list(cn1,cn2,cn3,cn4,cn00)
# cnlist <- list(cn1,cn2,cn3)

graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(cnlist))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- cnlist[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(cnlist[[i]])
}
#####################################################################################
# Weighted normalized closeness
#####################################################################################
for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  bnlist[[i]]<-g
}
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  cnlist[[i]]<-g
}
weightscnlist <- lapply(cnlist, function(x)E(x)$weights)
weightsbnlist <- lapply(bnlist, function(x)E(x)$weights)
#####################################################################################
# Weighted normalized closeness
#####################################################################################
closenesscn <- mapply(closeness, graph = cnlist, weights = weightscnlist, normalized = TRUE, SIMPLIFY = FALSE)
closenessbn <-  mapply(closeness, graph = bnlist, weights = weightsbnlist, normalized = FALSE, SIMPLIFY = FALSE)

climclosenesscn <- lapply(closenesscn, quantity2clim, what = "closeness",ref.grid = tas_ncep_10d)
climclosenessbn <- lapply(closenessbn, quantity2clim, what = "closeness",ref.grid = tas_ncep_10d)

cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
attr(climclosenesscn[[1]]$Data, "climatology:fun")

plotscn <- list()
meas <- climclosenesscn
for (i in 1:length(meas)){
  plot <- spatialPlot(meas[[i]], main = list(paste0("CN:",nedgescnlist[[i]]), cex = 0.5),
                          backdrop.theme = "coastline", 
                          region = TRUE, col.regions= cb,
                          # set.max = 20000, 
                          colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  plotscn[[i]] <- plot
}
plotscn
do.call(grid.arrange, c(plotscn, nrow = 2, as.table = TRUE, top = paste0(attr(meas[[1]]$Data, "climatology:fun"))))


plotsbn <- list()
measbn <- climclosenessbn
for (i in 1:length(measbn)){
  plot <- spatialPlot(measbn[[i]], main = list(paste0("CN:",nedgeslistbn[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      region = TRUE, col.regions= cb,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)),
                      set.max = 0.000015)
  plotsbn[[i]] <- plot
}


do.call(grid.arrange, c(plotsbn, nrow = 2, as.table = TRUE, top = paste0(attr(meas[[1]]$Data, "climatology:fun"))))


#####################################################################################
# Weighted betweenness
#####################################################################################
betweennessbn <- mapply(betweenness,bnlist, weights = weightsbnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
betweennesscn <- mapply(betweenness,cnlist, weights = weightscnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
climbetwbn <- lapply(betweennessbn, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)
climbetwcn <- lapply(betweennesscn, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)





plotscn <- list()
meas <- climbetwcn
for (i in 1:length(climbetwcn)){
  plot <- spatialPlot(climbetwcn[[i]], main = list(paste0("CN:",nedgescnlist[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      region = TRUE, col.regions= cb,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)),
                      set.max = 20000
  )
  plotscn[[i]] <- plot
}
do.call(grid.arrange, c(plotscn, nrow = 3, as.table = TRUE, top = paste0(attr(meas[[1]]$Data, "climatology:fun"))))


plotsbn <- list()
measbn <- climbetwbn
for (i in 1:length(climbetwbn)){
  plot <- spatialPlot(climbetwbn[[i]], main = list(paste0("BN:",nedgeslistbn[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      region = TRUE, col.regions= cb,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)),
                      set.max = 20000
  )
  plotsbn[[i]] <- plot
}
do.call(grid.arrange, c(plotsbn, nrow = 2, as.table = TRUE, top = paste0(attr(measbn[[1]]$Data, "climatology:fun"))))

blank <- grid.rect(gp=gpar(col="white"))
plist <- list(plotsbn[[4]],plotscn[[2]],blank,plotscn[[3]],blank,plotscn[[4]], blank,plotscn[[5]])
do.call(grid.arrange, c(plist, ncol = 2))


##############################################################################
# Unweighted betweenness
#####################################################################################
weights1cnlist <- lapply(cnlist, function(x)E(x)$weight)
weights1bnlist <- lapply(bnlist, function(x)E(x)$weight)

betweennessbn <- mapply(betweenness,bnlist, weights = weights1bnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
betweennesscn <- mapply(betweenness,cnlist, weights = weights1cnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
logbetweennessbn <- lapply(betweennessbn, function(x)log(1+x))
logbetweennesscn <- lapply(betweennesscn, function(x)log(1+x))
# climbetwbn <- lapply(logbetweennessbn, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)
# climbetwcn <- lapply(logbetweennesscn, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)


movmed_betw_bn <- lapply(logbetweennessbn, movingmedias)
movmed_betw_bn_mean <- lapply(movmed_betw_bn, function(x) apply(x, MARGIN = 3, FUN = mean))

climbetwbn <- lapply(movmed_betw_bn_mean, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)

movmed_betw_cn <- lapply(logbetweennesscn, movingmedias)
movmed_betw_cn_mean <- lapply(movmed_betw_cn, function(x) apply(x, MARGIN = 3, FUN = mean))

climbetwcn <- lapply(movmed_betw_cn_mean, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)



plotscn <- list()
meas <- climbetwcn
for (i in 1:length(climbetwcn)){
  plot <- spatialPlot(climbetwcn[[i]], main = list(paste0("CN:",nedgescnlist[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      col = "Reds",
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
  )
  plotscn[[i]] <- plot
}
do.call(grid.arrange, c(plotscn, nrow = 3, as.table = TRUE, top = paste0("Moving Medias CN: Log(1 + ",attr(meas[[1]]$Data, "climatology:fun"),")")))


plotsbn <- list()
measbn <- climbetwbn
for (i in 1:length(climbetwbn)){
  plot <- spatialPlot(climbetwbn[[i]], main = list(paste0("BN:",nedgeslistbn[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      col = "Reds",
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
                      )
  plotsbn[[i]] <- plot
}
do.call(grid.arrange, c(plotsbn, nrow = 2, as.table = TRUE, top = paste0("Moving Medias 1st perm BN: Log(1 + ",attr(measbn[[1]]$Data, "climatology:fun"),")")))

blank <- grid.rect(gp=gpar(col="white"))
plist <- list(plotsbn[[1]],plotscn[[1]],plotsbn[[2]],plotscn[[2]],plotsbn[[4]],plotscn[[3]],blank,plotscn[[4]], blank,plotscn[[5]])
do.call(grid.arrange, c(plist, ncol = 2, top = paste0("Moving Medias perm 1 BN CN: Log(1 + ",attr(meas[[1]]$Data, "climatology:fun"),")")))

##############################################################################
# Inverse weighted betweenness
#####################################################################################
weights3cnlist <- lapply(cnlist, function(x)E(x)$weights)
invweights3cnlist <- lapply(weights3cnlist, function(x)1/x)
weights3bnlist <- lapply(bnlist, function(x)E(x)$weights)
invweights3bnlist <- lapply(weights3bnlist, function(x)1/x)

betweennessbn <- mapply(betweenness,bnlist, weights = invweights3bnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
betweennesscn <- mapply(betweenness,cnlist, weights = invweights3cnlist, SIMPLIFY = FALSE) # HAS TO BE POSITIVE
logbetweennessbn <- lapply(betweennessbn, function(x)log(1+x))
logbetweennesscn <- lapply(betweennesscn, function(x)log(1+x))
# climbetwbn <- lapply(logbetweennessbn, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)
# climbetwcn <- lapply(logbetweennesscn, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)


movmed_betw_bn <- lapply(logbetweennessbn, movingmedias)
movmed_betw_bn_mean <- lapply(movmed_betw_bn, function(x) apply(x, MARGIN = 3, FUN = mean))

climbetwbn <- lapply(movmed_betw_bn_mean, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)

movmed_betw_cn <- lapply(logbetweennesscn, movingmedias)
movmed_betw_cn_mean <- lapply(movmed_betw_cn, function(x) apply(x, MARGIN = 3, FUN = mean))

climbetwcn <- lapply(movmed_betw_cn_mean, quantity2clim, what = "betweenness",ref.grid = tas_ncep_10d)



plotscn <- list()
meas <- climbetwcn
for (i in 1:length(climbetwcn)){
  plot <- spatialPlot(climbetwcn[[i]], main = list(paste0("CN:",nedgescnlist[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      col = "Reds",
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
  )
  plotscn[[i]] <- plot
}
do.call(grid.arrange, c(plotscn, nrow = 3, as.table = TRUE, top = paste0("Moving Medias CN: Log(1 + ",attr(meas[[1]]$Data, "climatology:fun"),")\ninverse weights: Strong correlation = short distance")))


plotsbn <- list()
measbn <- climbetwbn
for (i in 1:length(climbetwbn)){
  plot <- spatialPlot(climbetwbn[[i]], main = list(paste0("BN:",nedgeslistbn[[i]]), cex = 0.5),
                      backdrop.theme = "coastline", 
                      lonCenter = 180,
                      col = "Reds",
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
  )
  plotsbn[[i]] <- plot
}
do.call(grid.arrange, c(plotsbn, nrow = 2, as.table = TRUE, top = paste0("Moving Medias 1st perm BN: Log(1 + ",attr(measbn[[1]]$Data, "climatology:fun"),")\ninverse weights: Strong BIC = short distance")))

blank <- grid.rect(gp=gpar(col="white"))
plist <- list(plotsbn[[1]],plotscn[[1]],plotsbn[[2]],plotscn[[2]],plotsbn[[4]],plotscn[[3]],blank,plotscn[[4]], blank,plotscn[[5]])
do.call(grid.arrange, c(plist, ncol = 2, top = paste0("Moving Medias perm 1 BN CN: Log(1 + ",attr(meas[[1]]$Data, "climatology:fun"),")\ninverse weights: Strong BIC and COR = short distance")))

