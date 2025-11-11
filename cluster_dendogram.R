#######################################################################################
# Cluster dendogram BN and CN
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

################################################################################
# BN load and order
################################################################################
arcs1 <- lapply(perm1sort, arcs)
narcs1 <- sapply(arcs1, nrow)

nedgeslistbn <- narcs1[c(7,10,18)]
numbersBN <- as.character(nedgeslistbn)
bnsel1 <- 7
bnsel2 <- 10
bnsel3 <- 18

nedgeslistbn <- narcs1[c(18,26,87)]
numbersBN <- as.character(nedgeslistbn)
bnsel1 <- 18
bnsel2 <- 26
bnsel3 <- 87

bn1 <- perm1weights[[bnsel1]]
bn2 <- perm1weights[[bnsel2]]
bn3 <- perm1weights[[bnsel3]]

bnlist <- list(bn1,bn2,bn3)
commu_betsbn <- list()

##############################################################
# Weighted (Do not execute) ------------- 1
##############################################################
for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weights
  )
  commu_betsbn[[i]]<- commu_bet
}


name <- paste0("commu_BN_",numbersBN[[1]],"_",numbersBN[[2]],"_",numbersBN[[3]])
save(list = "commu_betsbn",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))
##############################################################
# Unweighted (Do not execute) -------------- 2
##############################################################
for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weight
  )
  commu_betsbn[[i]]<- commu_bet
}

name <- paste0("commu_BN_unw_",numbersBN[[1]],"_",numbersBN[[2]],"_",numbersBN[[3]])
commu_betsbn_unw <- commu_betsbn
save(list = "commu_betsbn_unw",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))
##############################################################
# Inverse weighted (Do not execute) ------------ 3
##############################################################
# inverse weighted
for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = 1/E(g)$weights
  )
  commu_betsbn[[i]]<- commu_bet
}

name <- paste0("commu_BN_invw_",numbersBN[[1]],"_",numbersBN[[2]],"_",numbersBN[[3]])
commu_betsbn_invw <- commu_betsbn
save(list = "commu_betsbn_invw",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))

###############################################################
# 1 weighted
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_BN_1795_2592_7862.rda")
# 2 unweighted
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_BN_unw_1795_2592_7862.rda")
# 3 inverse weighted 
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_BN_invw_1795_2592_7862.rda")

plotsbn <- list()
for(i in 1:length(bnlist)){
  ncom <- length(levels(factor(commu_betsbn[[i]]$membership)))
  mem <- membership(commu_betsbn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d, backperm = backpermutations[[1]])
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                              set.min = NULL, set.max = NULL, 
                              lonCenter = 180, 
                              regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                              main = paste0("Comunities BN ",nedgeslistbn[i]),
                              colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                              lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, c(plotsbn, ncol = 1))
#######################################################################################################
# CN load and order:  List as in Representations
#######################################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda") # <- no es gridGraph! list cn weights

edgelistsCN <- lapply(listcnweights, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(listcnweights) <- as.character(numberofedgesCN)


# This one is used in borrador. 999 of betweenness not yet. 
nedgescnlist <- numberofedgesCN[c(61,56,35,32,1)]
numbersCN <- as.character(nedgescnlist)
cn1 <- listcnweights[[61]]
cn2 <- listcnweights[[56]]
cn3 <- listcnweights[[35]]
cn4 <- listcnweights[[32]]

cnlist <- list(cn1,cn2,cn3,cn4)
commu_betscn <- list()
plotscn <- list()


##################################################################
# weighted  (Do not execute) 4
##################################################################
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weights
  )
  commu_betscn[[i]]<- commu_bet}

# For the two extra
name <- paste0("commu_CN_",numbersCN[[1]],"_",numbersCN[[2]])
save(list = "commu_betscn",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))
# Which we already used in borrador
name <- paste0("commu_CN_",numbersCN[[1]],"_",numbersCN[[2]],"_",numbersCN[[3]],"_",numbersCN[[4]])
save(list = "commu_betscn",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))

################################################################
# unweighted (Do not execute) 5
################################################################
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weight
  )
  commu_betscn[[i]]<- commu_bet}

commu_betscn_unw <- commu_betscn
name <- paste0("commu_CN_unw_",numbersCN[[1]],"_",numbersCN[[2]],"_",numbersCN[[3]],"_",numbersCN[[4]])
save(list = "commu_betscn_unw",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))

################################################################
# inverse weighted (Do not execute) 6
################################################################
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = 1/E(g)$weights
  )
  commu_betscn[[i]]<- commu_bet}

commu_betscn_invw <- commu_betscn
name <- paste0("commu_CN_invw_",numbersCN[[1]],"_",numbersCN[[2]],"_",numbersCN[[3]],"_",numbersCN[[4]])
save(list = "commu_betscn_invw",file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))

################################################################################
# 4
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_CN_1783_2260_8066_10017.rda")
# 5
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_CN_unw_1783_2260_8066_10017.rda")
# 6
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_CN_invw_1783_2260_8066_10017.rda")

commu_betscn <- commu_betscn_unw
for (i in 1:length(cnlist)){
  ncom <- length(levels(factor(commu_betscn[[i]]$membership)))
  mem <- membership(commu_betscn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotscn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                              set.min = NULL, set.max = NULL, 
                              lonCenter = 180, 
                              regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                              main = paste0("Comunities CN ",nedgescnlist[i]),
                              colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                              lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, plotscn)  
########################################################################################
# Investigate dendograms.
########################################################################################
cn_dends <- lapply(commu_betscn, as.dendrogram)
plot(cn_dends[[2]],leaflab = "none")
cn_dends_cut20 <- lapply(cn_dends, cut, h = 648-14)
cn_dends_cut20
as.dendrogram(cn_dends_cut20[[3]]$upper)
cn_clust_cut20_3 <- cn_dends_cut20[[3]]$upper
plot(cn_clust_cut20_3)
plot(cn_dends_cut20[[3]]$upper,leaflab = "none", ylim = c(0,648))
cn_clusts <- lapply(commu_betscn, as.hclust)

plot(cn_clusts[[2]],labels = FALSE)
membership(c1)

#commu_betsbn <- commu_betsbn_unw
bn_dends <- lapply(commu_betsbn, as.dendrogram)
plot(bn_dends[[1]])
bn_dends_cut20 <- lapply(bn_dends, cut, h = 648-14)
as.dendrogram(bn_dends_cut20[[1]]$upper)
bn_clust_cut20_3 <- bn_dends_cut20[[1]]$upper
plot(bn_clust_cut20_3)
plot(bn_dends_cut20[[1]]$upper,leaflab = "none", ylim = c(0,648))
plot(bn_dends_cut20[[3]]$upper)
bn_clusts <- lapply(commu_betsbn, as.hclust)
plot(bn_clusts[[2]])
membership(c1)

nobs(bn_dends[[1]])
bn_dends[[1]][c()]
nobs(bn_dends_cut20[[3]]$upper)
bn_dends_cut20[[3]]$lower
length(bn_dends[[1]][[3]])
nobs(bn_dends[[1]][[1]])
nobs(bn_dends[[1]][[2]])
bn_dends[[1]][[2]]
commu_betsbn[[1]]

############################################################################
# Plot dendograms as belong to communities
############################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/dends_BN",nedgeslistbn[1],"_CN",nedgescnlist[2],"_",nedgescnlist[3],".pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/dends_BN",nedgeslistbn[1],"_CN",nedgescnlist[2],"_",nedgescnlist[3],"_unw.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/dends_BN",nedgeslistbn[1],"_CN",nedgescnlist[2],"_",nedgescnlist[3],"_invw.pdf")
pdf(plotname)
par(mfrow = c(3,1))
yas <- 20

ncom1 <- 8
ncom2 <- 12
ncom3 <- 20
plot(bn_dends[[1]],leaflab = "none", ylim = c(649-yas-0.5,649), main = paste0("BN: ",nedgeslistbn[1]),yaxt = "n")
axis(2,at = 649 - c(1,ncom1,ncom2,yas)-0.5, labels = as.character(c(1,ncom1,ncom2,yas)))
abline(h = 649-1-0.5, lty = 2)
abline(h = 649-ncom1-0.5, lty = 2)
abline(h = 649-ncom2-0.5, lty = 2)
abline(h = 649-ncom3-0.5, lty = 2)
plot(cn_dends[[2]],leaflab = "none", ylim = c(649-yas-0.5,649), main = paste0("CN: ",nedgescnlist[2]), yaxt = "n")
axis(2,at = 649 - c(1,ncom1,ncom2,yas)-0.5, labels = as.character(c(1,ncom1,ncom2,yas)))
abline(h = 649-1-0.5, lty = 2)
abline(h = 649-ncom1-0.5, lty = 2)
abline(h = 649-ncom2-0.5, lty = 2)
abline(h = 649-ncom3-0.5, lty = 2)
plot(cn_dends[[3]],leaflab = "none", ylim = c(649-yas-0.5,649), main = paste0("CN: ",nedgescnlist[3]), yaxt = "n")
axis(2,at = 649 - c(1,ncom1,ncom2,yas)-0.5, labels = as.character(c(1,ncom1,ncom2,yas)))
abline(h = 649-1-0.5, lty = 2)
abline(h = 649-ncom1-0.5, lty = 2)
abline(h = 649-ncom2-0.5, lty = 2)
abline(h = 649-ncom3-0.5, lty = 2)
dev.off()
#############################################################################
# p log p dendograms
############################################################################

cn_dends <- lapply(commu_betscn, as.dendrogram)
bn_dends <- lapply(commu_betsbn, as.dendrogram)


commu.hist <- function(commuobject, numb.commu){
  dend <- as.dendrogram(commuobject)
  heightdend <- cut(dend, h = 648 - numb.commu)
  branchsize <- sapply(heightdend$lower, nobs)
  branchsizef <- as.factor(branchsize)
  
  
  logn <- log(length(levels(branchsizef)))
  probs <- hist(branchsize, 
                breaks = seq(0,max(branchsize),1), 
                probability = TRUE,
                plot = FALSE)
  
  return(probs)
}

commu.plogp(cn_dends[[1]],647)

commuobject <- bn_dends[[1]]
c5 <- cut_at(commu_betsbn[[1]],no =3)
commuobject <- commu_betscn[[3]]
tcs <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
lattitudes <- attr(tcs,"VertexCoords")[,"y"]
which(c5 ==1)
totarea <- sum(cos(y/180*pi))
no <- 3
cno <- cut_at(commu_betsbn[[1]],no =no)
cnoarreas <- c()
i <- 1
for(i in 1:no){
  commuarrea <- sum(cos(y/180*pi)[which(cno==i)])
  cnoarreas[i] <- commuarrea
}
sum(cnoarreas)
length(heightdend$lower)
members(heightdend$lower)
numb.commu <- 2

commu.plogp.aw <- function(commuobject, numb.commu, ycoords){
  y <- ycoords
  cno <- cut_at(commuobject, no = numb.commu)
  cnoarreas <- c()
  for(i in 1:numb.commu){
    commuarrea <- sum(cos(y/180*pi)[which(cno==i)])
    cnoarreas[i] <- commuarrea
  }
  totarea <- sum(cos(y/180*pi))
  
  logn <- log2(numb.commu)
  probs <- cnoarreas/totarea
  plogps <- probs*(log2(probs))
  
  plogp <- -sum(plogps)
  entr <- plogp/logn
  return(list(plogp,entr))
}
commu.plogp.aw(commu_betsbn[[1]],2, y)
commu.plogp.aw(commuobject = commu_betscn[[3]],ycoords = lattitudes, numb.commu = 2)

commu.plogp <- function(commuobject, numb.commu){
  dend <- as.dendrogram(commuobject)
  heightdend <- cut(dend, h = 648  - numb.commu)
  branchsize <- sapply(heightdend$lower, nobs)
  branchsizef <- as.factor(branchsize)
  
  logn <- log2(numb.commu)
  probs <- branchsize/648
  plogps <- probs*(log2(probs))
  
  plogp <- -sum(plogps)
  entr <- plogp/logn
  return(list(plogp,entr))
}


commu.plogp.2 <- function(commuobject, numb.commu){
  dend <- as.dendrogram(commuobject)
  heightdend <- cut(dend, h = 648 -1 - numb.commu)
  branchsize <- sapply(heightdend$lower, nobs)
  branchsizef <- as.factor(branchsize)
  
  
  logn <- log(length(levels(branchsizef)))
  probs <- hist(branchsize, 
                breaks = seq(0,max(branchsize),1), 
                probability = TRUE,
                plot = FALSE)
  plogps <- c()
  for(i in 1:max(branchsize)){
    pi <- probs$density[i]
    binwidthi <- probs$breaks[2]-probs$breaks[1]
    plogp <- pi*log(pi/binwidthi)
    plogps[i] <- plogp
  }
  
  plogp <- sum(plogps,na.rm = TRUE)
  entr <- -plogp/logn
  return(list(plogp,entr))
}


nedgescnlistchar <- as.character(nedgescnlist[1:4])
dim(nedgescnlistchar) <- (c(1,length(nedgescnlistchar))) 
dfcnnames <- apply(X = nedgescnlistchar, MARGIN = 2,FUN = function(x) paste0("CN",x))
ncom <- 648

dfcn <- array(data = NA, dim = c(ncom-1,length(dfcnnames),2), dimnames = list(as.character(2:ncom),dfcnnames,c("plogp","entropy"))) 
for(i in 1:length(dfcnnames)){
  output<- apply(X = matrix(2:ncom,nrow = 1, ncol = ncom-1) ,MARGIN = 2, FUN = commu.plogp, commuobject = cn_dends[[i]])
  dfcn[,i,1] <- sapply(output,function(x){x[[1]]})
  dfcn[,i,2] <- sapply(output,function(x){x[[2]]})
  }

dfcn

nedgesbnlistchar <- as.character(nedgeslistbn[1:3])
dim(nedgesbnlistchar) <- (c(1,length(nedgesbnlistchar))) 
dfbnnames <- apply(X = nedgesbnlistchar, MARGIN = 2,FUN = function(x) paste0("BN",x))
ncom <- 648

dfbn <- array(data = NA, dim = c(ncom-1,length(dfbnnames),2), dimnames = list(as.character(2:ncom),dfbnnames,c("plogp","entropy"))) 
for(i in 1:length(dfbnnames)){
  output<- apply(X = matrix(2:ncom,nrow = 1, ncol = ncom-1) ,MARGIN = 2, FUN = commu.plogp, commuobject = bn_dends[[i]])
  dfbn[,i,1] <- sapply(output,function(x){x[[1]]})
  dfbn[,i,2] <- sapply(output,function(x){x[[2]]})
}





# plogp plot 
bluebn <- colorRampPalette(c("light blue","navy"), alpha = TRUE)(3)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/plogp.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/plogp_smallBNs.pdf")
pdf(plotname)
dev.off()

plot(log2(x=2:ncom),y = dfcn[,"CN1783","plogp"], col = blackcn[[1]], xlab = "log2(number of communities)", ylab = "I = - sum pi·log2(pi)", xlim = c(0,10), main = "average Information gain communities") 
axis(1,seq(1,9,1),round(2^seq(1,9,1)))
points(log2(x=2:ncom),y = dfbn[,"BN1795","plogp"],col = bluebn[[1]]) 
points(log2(x=2:ncom),y = dfbn[,"BN700","plogp"],col = bluebn[[1]]) 
points(log2(x=2:ncom),y = dfcn[,"CN2260","plogp"], ylim = c(-3,0),col = blackcn[[2]])   
points(log2(x=2:ncom),y = dfbn[,"BN999","plogp"], ylim = c(-3,0),col = bluebn[[2]]) 
points(log2(x=2:ncom),y = dfbn[,"BN2592","plogp"], ylim = c(-3,0),col = bluebn[[2]]) 
points(log2(x=2:ncom),y = dfcn[,"CN8066","plogp"], ylim = c(-3,0),col = blackcn[[3]])   
points(log2(x=2:ncom),y = dfbn[,"BN1795","plogp"], ylim = c(-3,0),col = bluebn[[3]]) 
points(log2(x=2:ncom),y = dfbn[,"BN7862","plogp"], ylim = c(-3,0),col = bluebn[[3]]) 
points(log2(x=2:ncom),y = dfcn[,"CN10017","plogp"], ylim = c(-3,0),col = blackcn[[4]]) 
lines(log2(x=1:ncom),log2(x=1:ncom),lty = 2)
legend("bottomright",legend = c(dfbnnames,dfcnnames),fill = c(bluebn,blackcn))
legend("bottomright",legend = c(dfbnnames[1],dfbnnames[3],dfcnnames[2],dfcnnames[3]),fill = c(bluebn,blackcn))
# entropy plot
plot(x=1:nrow(dfcn),y = dfcn[,"CN1783","entropy"], ylim = c(0,1), col = blackcn[[1]], xlab = "number of communities", ylab = "entropy") 
points(x=1:nrow(dfbn),y = dfbn[,"BN1795","entropy"], col = bluebn[[1]]) 
points(x=1:nrow(dfcn),y = dfcn[,"CN2260","entropy"], col = blackcn[[2]])   
points(x=1:nrow(dfbn),y = dfbn[,"BN2592","entropy"], col = bluebn[[2]]) 
points(x=1:nrow(dfcn),y = dfcn[,"CN8066","entropy"], col = blackcn[[3]])   
points(x=1:nrow(dfbn),y = dfbn[,"BN7862","entropy"], col = bluebn[[3]]) 
points(x=1:nrow(dfcn),y = dfcn[,"CN10017","entropy"], col = blackcn[[4]]) 
legend("topright",legend = c(dfbnnames,dfcnnames),fill = c(bluebn,blackcn))



# Entropy plot
par(mfrow = c(3,2))
plot(x=1:nrow(dfcn),y = dfcn$CN1783, ylim = c(0,1)) 
points(x=1:nrow(dfbn),y = dfbn$BN1795, ylim = c(0,1),col = "blue") 

plot(x=1:nrow(dfcn),y = dfcn$CN2260, ylim = c(0,1))   
points(x=1:nrow(dfbn),y = dfbn$BN2592,ylim = c(0,1),col = "blue") 

plot(x=1:nrow(dfcn),y = dfcn$CN8066, ylim = c(0,1)) 
points(x=1:nrow(dfbn),y = dfbn$BN7862, ylim = c(0,1), col = "blue")  
dev.off()

# Entropyplot 8066 ft 1795
plot(x=1:nrow(dfcn),y = dfcn$CN8066, ylim = c(0,1)) 
points(x=1:nrow(dfbn),y = dfbn$BN1795, ylim = c(0,1),col = "blue") 



#######################################################################################
# area weighted plogp AFMAKEN!!!
#######################################################################################

nedgescnlistchar <- as.character(nedgescnlist[1:4])
dim(nedgescnlistchar) <- (c(1,length(nedgescnlistchar))) 
dfcnnames <- apply(X = nedgescnlistchar, MARGIN = 2,FUN = function(x) paste0("CN",x))
ncom <- 648
timecoords <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
lattitudes <- attr(timecoords,"VertexCoords")$y

dfcn <- array(data = NA, dim = c(ncom-1,length(dfcnnames),2), dimnames = list(as.character(2:ncom),dfcnnames,c("plogp","entropy"))) 
i <- 3
for(i in 1:length(dfcnnames)){
  output<- apply(X = matrix(2:ncom,nrow = 1, ncol = ncom-1) ,MARGIN = 2, FUN = commu.plogp.aw, commuobject = commu_betscn[[i]], ycoords = lattitudes)
  dfcn[,i,1] <- sapply(output,function(x){x[[1]]})
  dfcn[,i,2] <- sapply(output,function(x){x[[2]]})
}

dfcn

nedgesbnlistchar <- as.character(nedgeslistbn[1:3])
dim(nedgesbnlistchar) <- (c(1,length(nedgesbnlistchar))) 
dfbnnames <- apply(X = nedgesbnlistchar, MARGIN = 2,FUN = function(x) paste0("BN",x))
ncom <- 648

dfbn <- array(data = NA, dim = c(ncom-1,length(dfbnnames),2), dimnames = list(as.character(2:ncom),dfbnnames,c("plogp","entropy"))) 
for(i in 1:length(dfbnnames)){
  output<- apply(X = matrix(2:ncom,nrow = 1, ncol = ncom-1) ,MARGIN = 2, FUN = commu.plogp.aw, commuobject = commu_betsbn[[i]], ycoords = lattitudes)
  dfbn[,i,1] <- sapply(output,function(x){x[[1]]})
  dfbn[,i,2] <- sapply(output,function(x){x[[2]]})
}





# plogp plot 
bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(3)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

plot(log2(x=2:ncom),y = dfcn[,"CN1783","plogp"], col = blackcn[[1]], xlab = "log2(number of communities)", ylab = "I = - sum pi·log2(pi)", xlim = c(0,10),main = "average Information gain communities AREA WEIGHTED", xaxt= "n") 
axis(1,seq(1,9,1),round(2^seq(1,9,1)))
points(log2(x=2:ncom),y = dfbn[,"BN1795","plogp"],col = bluebn[[1]]) 
points(log2(x=2:ncom),y = dfcn[,"CN2260","plogp"], ylim = c(-3,0),col = blackcn[[2]])   
points(log2(x=2:ncom),y = dfbn[,"BN2592","plogp"], ylim = c(-3,0),col = bluebn[[2]]) 
points(log2(x=2:ncom),y = dfcn[,"CN8066","plogp"], ylim = c(-3,0),col = blackcn[[3]])   
points(log2(x=2:ncom),y = dfbn[,"BN7862","plogp"], ylim = c(-3,0),col = bluebn[[3]]) 
points(log2(x=2:ncom),y = dfcn[,"CN10017","plogp"], ylim = c(-3,0),col = blackcn[[4]]) 
lines(log2(x=1:ncom),log2(x=1:ncom),lty = 2)
legend("bottomright",legend = c(dfbnnames,dfcnnames),fill = c(bluebn,blackcn))
legend("bottomright",legend = c(dfbnnames[1],dfbnnames[3],dfcnnames[2],dfcnnames[3]),fill = c(bluebn,blackcn))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/plogpWEIGHTED.pdf")
pdf(plotname)
dev.off()

# entropy plot
plot(x=1:nrow(dfcn),y = dfcn[,"CN1783","entropy"], ylim = c(0,1), col = blackcn[[1]], xlab = "number of communities", ylab = "entropy") 
points(x=1:nrow(dfbn),y = dfbn[,"BN1795","entropy"], col = bluebn[[1]]) 
points(x=1:nrow(dfcn),y = dfcn[,"CN2260","entropy"], col = blackcn[[2]])   
points(x=1:nrow(dfbn),y = dfbn[,"BN2592","entropy"], col = bluebn[[2]]) 
points(x=1:nrow(dfcn),y = dfcn[,"CN8066","entropy"], col = blackcn[[3]])   
points(x=1:nrow(dfbn),y = dfbn[,"BN7862","entropy"], col = bluebn[[3]]) 
points(x=1:nrow(dfcn),y = dfcn[,"CN10017","entropy"], col = blackcn[[4]]) 
legend("topright",legend = c(dfbnnames,dfcnnames),fill = c(bluebn,blackcn))



###################################################################################
# histograms
###################################################################################
nedgescnlistchar <- as.character(nedgescnlist[1:4])
dim(nedgescnlistchar) <- (c(1,length(nedgescnlistchar))) 
dfcnnames <- apply(X = nedgescnlistchar, MARGIN = 2,FUN = function(x) paste0("CN",x))
ncom <- 648

histcn <- list()
for(i in 1:length(dfcnnames)){
  output<- apply(X = matrix(1:ncom,nrow = 1, ncol = ncom) ,MARGIN = 2, FUN = commu.hist, commuobject = cn_dends[[i]])
  histcn[[i]] <- output
}

plot(histcn[[1]][[100]])

nedgesbnlistchar <- as.character(nedgeslistbn[1:3])
dim(nedgesbnlistchar) <- (c(1,length(nedgesbnlistchar))) 
dfbnnames <- apply(X = nedgesbnlistchar, MARGIN = 2,FUN = function(x) paste0("BN",x))
ncom <- 648

histbn <- list()
for(i in 1:length(dfbnnames)){
  output<- apply(X = matrix(1:ncom,nrow = 1, ncol = ncom) ,MARGIN = 2, FUN = commu.hist, commuobject = bn_dends[[i]])
  histbn[[i]] <- output
}

bluebn<- colorRampPalette(c("light blue","navy"), alpha = TRUE)(3)
blackcn <- colorRampPalette(c("grey","black"), alpha = TRUE)(4)

########################################################################
# Histogram Cn BN own scale
########################################################################
which <- 100
plot(histcn[[4]][[which]], col = blackcn[[4]], main= paste0("Histogram of branchsize ",which," communities"))
plot(histcn[[3]][[which]], col = blackcn[[3]],add = TRUE)
plot(histcn[[2]][[which]], col = blackcn[[2]], add = TRUE)
plot(histcn[[1]][[which]], col = blackcn[[1]], add = TRUE)

legend("topright",legend = c(dfcnnames),fill = c(blackcn))

plot(histbn[[3]][[which]], col = bluebn[[3]], main= paste0("Histogram of branchsize ",which," communities"))
plot(histbn[[2]][[which]], col = bluebn[[2]],add = TRUE)
plot(histbn[[1]][[which]], col = bluebn[[1]], add = TRUE)

legend("topright",legend = c(dfbnnames),fill = c(bluebn))
########################################################################
# histogram same small scale
########################################################################
which <- 14
plot(histcn[[4]][[which]], col = blackcn[[4]], main= paste0("Histogram of branchsize ",which," communities"), xlim = c(0,50))
plot(histcn[[3]][[which]], col = blackcn[[3]], add = TRUE)
plot(histcn[[2]][[which]], col = blackcn[[2]], add = TRUE)
plot(histcn[[1]][[which]], col = blackcn[[1]], add = TRUE)

legend("topright",legend = c(dfcnnames),fill = c(blackcn))

plot(histbn[[3]][[which]], col = bluebn[[3]], main= paste0("Histogram of branchsize ",which," communities"), xlim = c(0,50))
plot(histbn[[2]][[which]], col = bluebn[[2]], add = TRUE)
plot(histbn[[1]][[which]], col = bluebn[[1]], add = TRUE)

legend("topright",legend = c(dfbnnames),fill = c(bluebn))

#######################################################################
# Plot dendograms
#######################################################################

i <- 1
for(i in 1:3){
# plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/dends_BN",nedgeslistbn[i],"_CN",nedgescnlist[i],"_unw.pdf")
 plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/dends_BN",nedgeslistbn[i],"_CN",nedgescnlist[i],"_invw.pdf")
  # plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/dends_BN",nedgeslistbn[i],"_CN",nedgescnlist[i],".pdf")
  pdf(plotname)
par(mfrow = c(2,1))
plot(bn_dends[[i]], main = paste0("BN: ",nedgeslistbn[i]),leaflab = c("none"))
plot(cn_dends[[i]], main = paste0("CN: ",nedgescnlist[i]),leaflab = c("none"))
dev.off()
}

#######################################################################
# plot clusters
#######################################################################
par(mfrow = c(2,1))
for(i in 1:3){
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/dendograms/clusts_BN",nedgeslistbn[i],"_CN",nedgescnlist[i],".pdf")
  pdf(plotname)
  par(mfrow = c(2,1))
  plot(bn_clusts[[i]], main = paste0("BN: ",nedgeslistbn[i]), labels = FALSE)
  plot(cn_clusts[[i]], main = paste0("CN: ",nedgescnlist[i]), labels = FALSE)
  dev.off()
}



##############################################################################
#
##############################################################################
rect.hclust
## S3 method for class 'communities'
print(c1)

## S3 method for class 'communities'
modularity(x, ...)

## S3 method for class 'communities'
length(x)

sizes(c1)

algorithm(communities)

merges(c1)

crossing(communities, graph)

code_len(communities)

is_hierarchical(communities)

## S3 method for class 'communities'
dend <- as.dendrogram(c1, hang = 0.1,
                      use.modularity = FALSE)
str(dend)
plot(dend, "triangle")
dendh <- cut(dend,610)

plot(dendh$upper,edge.root = TRUE)
str(dend)
## S3 method for class 'communities'
clus <- as.hclust(c1)
plot(clus)

as_phylo(x, ...)

## S3 method for class 'communities'
as_phylo(x, use.modularity = FALSE, ...)

cut_at(communities, no, steps)

show_trace(communities)

## S3 method for class 'communities'
plot(x, y, col = membership(x),
     mark.groups = communities(x), edge.color = c("black", "red")[crossing(x,
                                                                           y) + 1], ...)