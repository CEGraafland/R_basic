library(igraph)
library(RColorBrewer)
library(MASS)
library(RColorBrewer)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
############################################################################
# Betweenness representation CN BN equal loglikelihoods. 
############################################################################
############################################################################
# Calculate betweenness CN networks
############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda")
############################################################################
# Calculate weighted betweenness CN
############################################################################
# select weights: In case of CN: correlation
CNweights <- lapply(listcnweights, function(x){E(x)$weights})
# Betweenness can not be calculated for edges with weight zero. Lets give them minimal weight. 
CNbetind <- sapply(listcnweights, function(x){which(E(x)$weights==0)})
listcnweightsplus <- list()
for(i in 1:length(listcnweights)){
  graph <- listcnweights[[i]]
  E(graph)$weights[CNbetind[i]] <- 0.000000001
  listcnweightsplus[[i]] <- graph
}
CNweightsplus <- lapply(listcnweightsplus, function(x){E(x)$weights})
CNbetw <- mapply(betweenness, graph = listcnweightsplus, weights = CNweightsplus, MoreArgs= list(normalized = TRUE))
############################################################################
# Calculate betweenness perm1sort (perm2 / perm3 / perm4 / perm 5)
############################################################################
permweight1_10 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10", full.names = T)
lapply(permweight1_10, load, envir=.GlobalEnv)
# or manual
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm1weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm2weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm3weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm4weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm5weights.rda")
# Create one list with all permutations.
BNpermweightslist <- list(perm1weights, perm2weights, perm3weights, perm4weights, perm5weights)
# select weights: In case of BN: eBIC
# Betweenness can not be calculated for edges with weight zero. Lets give them minimal weight. 
BNbetwperms <- list()
BNsperms <- list()
names <- c()

for (p in 1:length(BNpermweightslist)){
  names[p] <- paste0("BNbetw",p)
  listused <- BNpermweightslist[[p]]
  weights <- lapply(listused, function(x){E(x)$weights})
  betind <- sapply(listused, function(x){which(E(x)$weights==0)})
  listusedweightsplus <- list()
  for(i in 1:length(listused)){
    graph <- listused[[i]]
    E(graph)$weights[betind[i]] <- 0.000000001
    listusedweightsplus[[i]] <- graph
  }
  weightsplus <- lapply(listusedweightsplus, function(x){E(x)$weights})
  betwWeightedlist <- mapply(betweenness, graph = listusedweightsplus, weights = weightsplus, MoreArgs= list(normalized = TRUE))
  BNbetwperms[[p]] <- betwWeightedlist
  BNsperms[[p]] <- listusedweightsplus
} 
names(BNbetwperms) <- names
names(BNsperms) <- names
##################################################################################
# We have two lists: BNbetwperms (list of 5 BNbetw), CNbetw (one CNbetw)
# Choose from which you want a presentation
###################################################################################
# Bayesian list:
perm <- 1
permu <- permutations[[perm]]
i <-18
igraph <- BNsperms[[perm]][[i]]
betw <- BNbetwperms[[perm]][,i] 
# betw2 <- betw[order(betw, decreasing = TRUE)]
igraph <- set_vertex_attr(igraph, "bet", value = betw)
hist(V(igraph)$bet,freq = TRUE, xlim = c(0,0.15), breaks = 3)
# Correlation list
i <- 10
perm <- 1
permu <- permutations[[perm]]
igraph <- listcnweights[[i]]
length(E(igraph))
betw <- CNbetw[,i] 
igraph <- set_vertex_attr(igraph, "bet", value = betw)
hist(V(igraph)$bet, freq = TRUE, breaks = 3)
################################################################################
# Code to visualize.
# Choose how many clusters you want. 
# Vertices with equal betweenness are clustered. 
################################################################################
probs <- seq(0,1,(0.25))
lprobs <- length(probs)
probints <- lprobs -1 
quantiles<- quantile(betw,probs = probs)
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}
# Get sufficient colors for all different vertices in all clusters. 
if(probints <= 9){
  colReds <- brewer.pal(probints,"Reds")
  colReds <- rev(colReds)
} else {
  colReds <- brewer.pal(9,"Reds")
  # colReds <- colorRampPalette(colReds, alpha = TRUE)(probints)
  colReds <- colorRampAlpha(colReds, n=probints, alpha=0.5)
  colReds <- rev(colReds)}

# Coloring / membermaking / vertices. 
# First cluster
V(igraph)[V(igraph)$bet <= quantiles[2]]$color <- colReds[1]
V(igraph)[V(igraph)$bet <= quantiles[2]]$memb <- 1
# middle clusters and last cluster
for(i in 2:(probints)){
  V(igraph)[V(igraph)$bet> quantiles[i] & V(igraph)$bet <=quantiles[i+1]]$color <- colReds[i]
  V(igraph)[V(igraph)$bet> quantiles[i] & V(igraph)$bet <=quantiles[i+1]]$memb <- i
}

ordall <- order(V(igraph)$bet,decreasing = TRUE)
permordall <- c()
for(j in 1:length(ordall)){
  permordall[ordall[j]] <- j
}
igraphord <- permute(igraph,permordall)



# Make igraph clusters by members. 
# And plot all clusters. 
clust <- make_clusters(igraphord, membership = V(igraphord)$memb, algorithm = "betw quantile", merges = NULL,
                       modularity = TRUE)

# Subgraph
# Choose how many clusters you want to be presented. 
# We present only strongest. 
nclust <- 1
idsvector <- c()
idslist <- list()
for(i in 1:nclust){
  idsvector <- c(idsvector, V(igraphord)[V(igraphord)$memb== (probints-(i-1))])
  idslist[[i]] <- V(igraphord)[V(igraphord)$memb== (probints-(i-1))]
}

ng <- induced_subgraph(igraphord,idsvector)
clustng <- make_clusters(ng, membership = V(ng)$memb, algorithm = "betw quantile", merges = NULL,
                         modularity = TRUE)
plot(clustng,ng,vertex.size = 3, edge.arrow.size = 0.2, vertex.label = V(ng)$name, layout = layout_in_circle)
plot.igraph(ng,
            vertex.size = 1,
            vertex.label = V(ng)$name,
            layout = layout_in_circle,
            mark.groups = idslist,
            mark.col = colorRampAlpha(colReds, n=length(idslist), alpha=0.5),
            #mark.col = rev(brewer.pal(length(idslist),"Reds")),
            mark.border = NA,
            main = paste0(length(E(igraphord))))
legend("left",legend = V(ng)$name, pch = NA)
plot.igraph(ng,
            vertex.size = 1,
            vertex.label = V(ng)$name,
            layout = layout_with_dh,
            mark.groups = idslist,
            mark.col = colorRampAlpha(colReds, n=length(idslist), alpha=0.5),
            mark.border = NA,
            main = paste0(length(E(igraphord))))
legend("left",legend = V(ng)$name, pch = NA)

nclust <- 1
idsvector <- c()
idslist <- list()
for(i in 1:nclust){
  idsvector <- c(idsvector, V(igraph)[V(igraph)$memb== (probints-(i-1))])
  idslist[[i]] <- V(igraph)[V(igraph)$memb== (probints-(i-1))]
}

ng <- induced_subgraph(igraph,idsvector)
V(ng)

# Get sufficient colors for all different vertices in all clusters. 
if(nclust <= 9){
  colReds <- brewer.pal(nclust,"Reds")
  colReds <- rev(colReds)
} else {
  colReds <- brewer.pal(9,"Reds")
  # colReds <- colorRampPalette(colReds, alpha = TRUE)(probints)
  colReds <- colorRampAlpha(colReds, n=nclust, alpha=0.5)
  colReds <- rev(colReds)}
# Coloring
for(i in 1:(nclust)){
  V(ng)[V(ng)$memb == (probints-(i-1))]$color <- colReds[i]
}


data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
shift = TRUE

x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)

if(shift == TRUE){x <- (x+360)%%360}

points <- expand.grid(y, x)[2:1]

if (!is.null(permu)){
  points <- points[permu,]
}

# load("/home/catharina/Documents/lisworld.rda")
# load("/Volumes/ubuntu/Documents/lisworld.rda")
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
if(shift == TRUE){raster::plot(crop(l1[[2]],m), main = paste0(length(E(igraphord))))} else {raster::plot(lisword, main = paste0(length(E(igraphord))))}

plot.igraph(ng,
            vertex.size = 400,
            # vertex.color = "blue",
            vertex.label = NA,
            edge.arrow.size = 0.2,
            # edge.color = NA,
            # mark.groups = idslist,
            # edge.width = E(igraph)$width,
            lty = "dots",
            main = paste0(length(E(igraphord))),
            layout = as.matrix(points[idsvector,]), add = TRUE, rescale = FALSE)


plot(clustng,ng,vertex.size = 3, edge.arrow.size = 0.2, vertex.label = V(ng)$name, layout = layout_in_circle)

##############################################################################
# OTher extra
##############################################################################
plot(clust,igraphord,vertex.size = 3,edge.arrow.size = 0.2,layout = layout_on_grid)
plot(clust, delete_edges(igraphord, edges = E(igraphord)),vertex.size = 3,edge.arrow.size = 0.2, layout = layout_in_circle)
# WHy Does plot(clust,igraph,vertex.size = 3,edge.arrow.size = 0.2, layout = layout_on_sphere) works??

