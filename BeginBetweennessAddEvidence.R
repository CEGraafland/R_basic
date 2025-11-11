library(igraph)
library(RColorBrewer)
############################################################################
# Betweenness representation CN BN equal loglikelihoods. 
############################################################################
############################################################################
# Calculate betweenness CN networks
############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda")
############################################################################
# put weights on igraphs cn
############################################################################

weights <- lapply(listcnweights, function(x){E(x)$weights})
betwlist <- sapply(listcnweights, betweenness, weights = NULL, normalized = TRUE)

weights <- lapply(listcnweights, function(x){E(x)$weights})
betind <- sapply(listcnweights, function(x){which(E(x)$weights==0)})
listcnweightsplus <- list()
for(i in 1:length(listcnweights)){
  graph <- listcnweights[[i]]
  E(graph)$weights[betind[i]] <- 0.000000001
  listcnweightsplus[[i]] <- graph
}

weightsplus <- lapply(listcnweightsplus, function(x){E(x)$weights})

betwWeightedlist <- mapply(betweenness, graph = listcnweightsplus, weights = weightsplus, MoreArgs= list(normalized = TRUE))

betwWeightedlist2 <- list()
for(i in 1:length(listcnweightsplus)){
  graph <- listcnweightsplus[[i]]
  betw <- betweenness(graph = graph, weights = weightsplus[[i]], normalized = TRUE)
  betwWeightedlist2[[i]] <- betw
}

all.equal(betwWeightedlist2[[70]],betwWeightedlist[,70])  # TRUE!

############################################################################
# Color vertexes: Quantile Betweenesss
############################################################################
i <-35
igraph <- listcnweightsplus[[i]]
length(E(igraph))
betw <- betwWeightedlist[,i] 
igraph <- set_vertex_attr(igraph, "bet", value = betw)
V(igraph)$bet

probs <- seq(0,1,(0.01))
lprobs <- length(probs)
probints <- lprobs -1 

quantiles<- quantile(betw,probs = probs)


if(probints <= 9){
  colReds <- brewer.pal(probints,"Reds")
} else {
  colReds <- brewer.pal(9,"Reds")
  colReds <- colorRampPalette(colReds)(probints)}

V(igraph)[betw <= quantiles[2]]$color <- colReds[1]
V(igraph)[betw <= quantiles[2]]$memb <- 1
V(igraph)[betw<= quantiles[2]]$label <- NA
for(i in 2:(probints)){
  V(igraph)[betw> quantiles[i] & betw <=quantiles[i+1]]$color <- colReds[i]
  V(igraph)[betw> quantiles[i] & betw <=quantiles[i+1]]$memb <- i
}
for(i in 2:(probints-1)){
  
  V(igraph)[betw > quantiles[2] & betw <=quantiles[i+1]]$label <- NA
}
V(igraph)[betw > quantiles[probints] & betw <=quantiles[probints+1]]
V(igraph)[memb==probints]
  # V(igraph)[betw > quantiles[3] & betw <=quantiles[4]]$color <- "red"
  # V(igraph)[betw > quantiles[4] & betw <=quantiles[5]]$color <- "blue"
  # 
  # 
  # V(igraph)[betwWeightedlist[,i] <= quantiles[2]]$memb <- 1
  # V(igraph)[betwWeightedlist[,i] > quantiles[2] & betwWeightedlist[,i] <=quantiles[3]]$memb <- 2
  # V(igraph)[betwWeightedlist[,i] > quantiles[3] & betwWeightedlist[,i] <=quantiles[4]]$memb <- 3
  # V(igraph)[betwWeightedlist[,i] > quantiles[4] & betwWeightedlist[,i] <=quantiles[5]]$memb <- 4
  # 
  # V(igraph)[betwWeightedlist[,i] <= quantiles[2]]$label <- NA
  # V(igraph)[betwWeightedlist[,i] > quantiles[2] & betwWeightedlist[,i] <=quantiles[3]]$label <- NA
  # V(igraph)[betwWeightedlist[,i] > quantiles[3] & betwWeightedlist[,i] <=quantiles[4]]$label <- NA
  # V(igraph)[betwWeightedlist[,i] > quantiles[4] & betwWeightedlist[,i] <=quantiles[5]]$label 
  
clust <- make_clusters(igraph, membership = V(igraph)$memb, algorithm = "betw quantile", merges = NULL,
                modularity = TRUE)
plot(clust,igraph,vertex.size = 3)
  
  
# weights <- ifelse(crossing(cl, graph), 1, 100)
  # layout <- layout_with_kk(igraph, weights=weights)
  
  grouplist <- list(which(V(igraph)$color == "blue"),
                    which(V(igraph)$color == "red"),
                    which(V(igraph)$color == "yellow"),
                    which(V(igraph)$color == "green"))
  
  clp <- cluster_label_prop(listcnweightsplus[[1]])
  class(clp)
  
  layout_with_kk(igraph, weights = betw)
  
  plot.igraph(igraph,
              vertex.size = 2,
              edge.arrow.size = 0.2, 
              mark.groups = list(which(V(igraph)$color == "blue"),
                               which(V(igraph)$color == "red"),
                               which(V(igraph)$color == "yellow"),
                               which(V(igraph)$color == "green")),
              layout = layout_with_drl)
  
############################################################################
# Save only strongest betweenness edges. 
############################################################################
ids <- V(igraph)[memb==probints]
ids <- c(V(igraph)[memb==probints-2],V(igraph)[memb==probints-1],V(igraph)[memb==probints])

ng <- induced_subgraph(igraph,ids)

plot(ng,
     vertex.size = 2,
     vertex.label = V(ng)$name)

ord <- order(V(ng)$bet,decreasing = TRUE)
permord <- c()
for(j in 1:length(ord)){
  permord[ord[j]] <- j
}
ngord <- permute(ng,permord)
V(ngord)
############################################################################
# Calculate layout (plot i.graph)
############################################################################

clustngord <- make_clusters(ngord, membership = V(ngord)$memb, algorithm = "betw quantile", merges = NULL,
                            modularity = TRUE)
ids2 <- list(V(ngord)[memb==probints-2],V(ngord)[memb==probints-1],V(ngord)[memb==probints])
plot.igraph(ngord,vertex.size = 1, edge.arrow.size = 0.2, vertex.label = V(ngord)$name, mark.groups = ids2,layout = layout_in_circle(ngord))
V(ngord)$name
V(ng)$name
layout_on_sphere(ngord)
id21 <- V(ngord)[memb==probints-2]
plot.igraph(ngord,
            vertex.size = 1,
            vertex.label = V(ngord)$name,
            # layout = layout_on_sphere(ngord),
            mark.groups = ids2,
            mark.col = brewer.pal(length(ids2),"Reds"))
plot.igraph(ngord,vertex.size = 1,mark.groups = ids2,layout = layout_in_circle(ngord))

plot.igraph(ngord,
            vertex.size = 1,
            vertex.label = V(ngord)$name,
            layout = layout_with_dh,
            mark.groups = ids2,
            mark.col = brewer.pal(length(ids2),"Reds"))


############################################################################
# Calculate betweenness perm1sort (perm2 / perm3 / perm4 / perm 5)
############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm2weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm3weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm4weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm5weights.rda")

listused <- perm2weights
weights <- lapply(listused, function(x){E(x)$weights})
betwlist <- sapply(listused, betweenness, weights = NULL, normalized = TRUE)

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

betwWeightedlist2 <- list()

for(i in 1:length(listusedweightsplus)){
  graph <- listusedweightsplus[[i]]
  betw <- betweenness(graph = graph, weights = weightsplus[[i]], normalized = TRUE)
  betwWeightedlist2[[i]] <- betw
}

all.equal(betwWeightedlist2[[70]],betwWeightedlist[,70])  # TRUE!



######################################################################################
i <-13
igraph <- listusedweightsplus[[i]]
betw <- betwWeightedlist[,i] 
betw2 <- betw[order(betw, decreasing = TRUE)]
igraph <- set_vertex_attr(igraph, "bet", value = betw)
V(igraph)$bet

probs <- seq(0,1,(0.01))
lprobs <- length(probs)
probints <- lprobs -1 

quantiles<- quantile(betw,probs = probs)
quantiles2 <- quantile(betw2,probs = probs)
quantiles
quantiles2
if(probints <= 9){
  colReds <- brewer.pal(probints,"Reds")
} else {
  colReds <- brewer.pal(9,"Reds")
  colReds <- colorRampPalette(colReds)(probints)}

V(igraph)[V(igraph)$bet <= quantiles[2]]$color <- colReds[1]
V(igraph)[V(igraph)$bet <= quantiles[2]]$memb <- 1
V(igraph)[V(igraph)$bet<= quantiles[2]]$label <- NA

for(i in 2:(probints)){
  V(igraph)[V(igraph)$bet> quantiles[i] & V(igraph)$bet <=quantiles[i+1]]$color <- colReds[i]
  V(igraph)[V(igraph)$bet> quantiles[i] & V(igraph)$bet <=quantiles[i+1]]$memb <- i
}
for(i in 2:(probints-1)){
  
  V(igraph)[betw > quantiles[2] & betw <=quantiles[i+1]]$label <- NA
}
length(unique(V(igraph)$memb))

V(igraph)[V(igraph)$bet > quantiles[100]]
V(igraph)[V(igraph)$memb==probints]
# V(igraph)[betw > quantiles[3] & betw <=quantiles[4]]$color <- "red"
# V(igraph)[betw > quantiles[4] & betw <=quantiles[5]]$color <- "blue"
# 
# 
# V(igraph)[betwWeightedlist[,i] <= quantiles[2]]$memb <- 1
# V(igraph)[betwWeightedlist[,i] > quantiles[2] & betwWeightedlist[,i] <=quantiles[3]]$memb <- 2
# V(igraph)[betwWeightedlist[,i] > quantiles[3] & betwWeightedlist[,i] <=quantiles[4]]$memb <- 3
# V(igraph)[betwWeightedlist[,i] > quantiles[4] & betwWeightedlist[,i] <=quantiles[5]]$memb <- 4
# 
# V(igraph)[betwWeightedlist[,i] <= quantiles[2]]$label <- NA
# V(igraph)[betwWeightedlist[,i] > quantiles[2] & betwWeightedlist[,i] <=quantiles[3]]$label <- NA
# V(igraph)[betwWeightedlist[,i] > quantiles[3] & betwWeightedlist[,i] <=quantiles[4]]$label <- NA
# V(igraph)[betwWeightedlist[,i] > quantiles[4] & betwWeightedlist[,i] <=quantiles[5]]$label 

clust <- make_clusters(igraph, membership = V(igraph)$memb, algorithm = "betw quantile", merges = NULL,
                       modularity = TRUE)
clust$membership
clust$membership[205]
plot(clust,igraph,vertex.size = 3,edge.arrow.size = 0.2)

# Subgraph
# Choose how many clusters you want to be presented. 
# We present only strongest. 
nclust <- 4
idsvector <- c()
idslist <- list()
for(i in 1:nclust){
  idsvector <- c(idsvector, V(igraph)[V(igraph)$memb== (probints-(i-1))])
  idslist[[i]] <- V(igraph)[V(igraph)$memb== (probints-(i-1))]
}

ng <- induced_subgraph(igraph,idsvector)
V(ng)


ids <- V(igraph)[V(igraph)$memb==probints]
ids <- c(V(igraph)[V(igraph)$memb==probints-2],V(igraph)[V(igraph)$memb==probints-1],V(igraph)[V(igraph)$memb==probints])
ng <- induced_subgraph(igraph,ids)
clustng <- make_clusters(ng, membership = V(ng)$memb, algorithm = "betw quantile", merges = NULL,
                         modularity = TRUE)
plot(clustng,ng,vertex.size = 3, edge.arrow.size = 0.2, vertex.label = V(ng)$name, layout = layout_with_dh)
V(ng)$memb
V(ng)$bet
V(igraph)[ids]$bet
ord <- order(V(ng)$bet,decreasing = TRUE)
permord <- c()
for(j in 1:length(ord)){
  permord[ord[j]] <- j
}
ngord <- permute(ng,permord)
clustngord <- make_clusters(ngord, membership = V(ngord)$memb, algorithm = "betw quantile", merges = NULL,
                            modularity = TRUE)
plot(clustngord,ngord,vertex.size = 3, edge.arrow.size = 0.2, vertex.label = V(ngord)$name, layout = layout_in_circle(ngord))
plot(clustngord,ngord,vertex.size = 3, edge.arrow.size = 0.2, vertex.label = V(ngord)$name, layout = layout_on_grid(ngord))
plot(clustng,ng,vertex.size = 3, edge.arrow.size = 0.2, vertex.label = V(ng)$name, layout = layout_with_dh)

V(ng)$memb
ids2 <- list(V(ngord)[memb==probints-2],V(ngord)[memb==probints-1],V(ngord)[memb==probints])

plot.igraph(ngord,
     vertex.size = 1,
     vertex.label = V(ngord)$name,
     layout = layout_in_circle,
     mark.groups = ids2,
     mark.col = brewer.pal(length(ids2),"Reds"),
     mark.border = NA)

plot.igraph(ngord,
            vertex.size = 1,
            vertex.label = V(ngord)$name,
            layout = layout_on_grid,
            mark.groups = ids2,
            mark.col = brewer.pal(length(ids2),"Reds"),
            mark.border = NA)

plot.igraph(ngord,
            vertex.size = 1,
            vertex.label = V(ngord)$name,
            layout = layout_with_dh,
            mark.groups = ids2,
            mark.col = brewer.pal(length(ids2),"Reds"))

plot.igraph(ngord,
            vertex.size = 1,
            vertex.label = V(ngord)$name,
            layout = layout_with_dh,
            mark.groups = ids2,
            mark.col = brewer.pal(length(ids2),"Reds"))


data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
shift = TRUE

x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)

if(shift == TRUE){x <- (x+360)%%360}

points <- expand.grid(y, x)[2:1]

if (!is.null(perm)){
  points <- points[perm,]
}

# load("/home/catharina/Documents/lisworld.rda")
# load("/Volumes/ubuntu/Documents/lisworld.rda")
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
if(shift == TRUE){raster::plot(crop(l1[[2]],m))} else {raster::plot(lisworld)}

plot.igraph(igraph,
            vertex.size = 300,
            vertex.color = "blue",
            vertex.label = NA,
            edge.arrow.size = 0.2,
            mark.groups = idslist,
            # edge.width = E(igraph)$width,
            lty = "dots",
            layout = as.matrix(points[idsvector,]), add = TRUE, rescale = FALSE)

