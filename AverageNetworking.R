#################################################################################
# Average networking
#################################################################################
library(bnlearn)
library(igraph)
library(RColorBrewer)
library(transformeR)
library(visualizeR)
#################################################################################
# Load sorted BNs permutations 1-25
#################################################################################
perm1_10 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort1_10", full.names = T)
perm1_10sort <- perm1_10
perm1_10sort[11] <- perm1_10[1] 
perm1_10sort <- perm1_10sort[2:11]
perm11_25 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permsort11_25", full.names = T)
permsort <- list()
namespermsort <- c()

le1 <- length(perm1_10sort)
le2 <- length(perm11_25)
for (i in 1:(le1+le2)){
  temp.space <- new.env()
  if(i <= le1){
    variable <- get(load(perm1_10sort[i], temp.space), envir = temp.space)
    namespermsort[i] <- ls(envir = temp.space)
    permsort[[i]] <- variable
    } else { 
    variable <- get(load(perm11_25[i-le1], temp.space), envir = temp.space)
    namespermsort[i] <- ls(envir = temp.space)
    permsort[[i]] <- variable
    }
  rm(temp.space)
}


names(permsort) <- namespermsort
#################################################################################
# Load BNS converted in igraphs with edge weight igual to BIC
#################################################################################
perm1_10weight <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10", full.names = T)
perm1_10weightsort <- perm1_10weight
perm1_10weightsort[11] <- perm1_10weight[1] 
perm1_10weightsort <- perm1_10weightsort[2:11]


perm11_25 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights11_25", full.names = T)
permweight <- list()
namespermweight<- c()

le1 <- length(perm1_10weightsort)
le2 <- length(perm11_25)
for (i in 1:(le1+le2)){
  temp.space <- new.env()
  if(i <= le1){
    variable <- get(load(perm1_10weightsort[i], temp.space), envir = temp.space)
    namespermweight[i] <- ls(envir = temp.space)
    permweight[[i]] <- variable
  } else { 
    variable <- get(load(perm11_25[i-le1], temp.space), envir = temp.space)
    namespermweight[i] <- ls(envir = temp.space)
    permweight[[i]] <- variable
  }
  rm(temp.space)
}


names(permweight) <- namespermweight

#############################################################################
# Average of highest betweenness nodes in all permutations
#############################################################################
nodos <- V(bnsw1800[[2]])
# Select all 1800 networks
bnsw1800_1 <- lapply(permweight[1:10], function(x) x[[18]])
bnsw1800_2 <- lapply(permweight[11:25], function(x) x[[7]])
bnsw1800 <- c(bnsw1800_1,bnsw1800_2)
# Add to the graphs a betweenness attribute
weights <- lapply(bnsw1800, function(x){E(x)$weights})
betind <- lapply(bnsw1800, function(x) which(E(x)$weights ==0))
listusedweightsplus <- list()
for(i in 1:length(bnsw1800)){
  x <- bnsw1800[[i]]
  y <- betind[[i]]
  E(bnsw1800[[i]])$weights[y] <- 0.000000001

}
weightsplus <- lapply(bnsw1800, function(x) E(x)$weights)
invweightsplus <- lapply(weightsplus, function(x) 1/x)

betwlist <- list()
# large weight = large distance
for (i in 1:length(bnsw1800)){
  betwlist[[i]] <- log(1+betweenness(graph = bnsw1800[[i]], weights = weightsplus[[i]], normalized = FALSE))}
# unweighted
for (i in 1:length(bnsw1800)){
  betwlist[[i]] <- log(1+betweenness(graph = bnsw1800[[i]], weights = NULL, normalized = FALSE))}
# large weight = short distance
for (i in 1:length(bnsw1800)){
  betwlist[[i]] <- log(1+betweenness(graph = bnsw1800[[i]], weights = invweightsplus[[i]], normalized = FALSE))}


for (i in 1:length(bnsw1800)){
bnsw1800[[i]] <- set_vertex_attr(bnsw1800[[i]], "bet", value = betwlist[[i]])}
# select upper 83 procent of every graph (Every graphs has other values, therefore select upper strongest)
probs <- seq(0,1,(1/4))
lprobs <- length(probs)
probints <- lprobs -1 
quantiles <- lapply(betwlist, quantile, probs = probs)

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

# Choose how strong you want to be the network
listhigh <- list()
i <- 1
for (i in 1:length(bnsw1800)){
  igraph <- bnsw1800[[i]]
  listhigh[[i]] <- V(igraph)[V(igraph)$bet >= quantiles[[i]][4]]}

listhighvs <- lapply(listhigh,names)
merges <- listhighvs[[1]]
for(i in 2:length(listhighvs)){
  merges <- c(merges,listhighvs[[i]])
}

# Choose how frequent you want to be the nodes 
mergesf <- factor(merges)
mergesf
tablemergesf <- table(mergesf)
dfmergesf <- as.data.frame(tablemergesf)
ind <- which(dfmergesf[,2]>15)
ind
highv <- as.character(dfmergesf[ind,1])
length(highv)
data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
shift = TRUE

x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)

if(shift == TRUE){x <- (x+360)%%360}

points <- expand.grid(y, x)[2:1]


# load("/home/catharina/Documents/lisworld.rda")
# load("/Volumes/ubuntu/Documents/lisworld.rda")
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
if(shift == TRUE){raster::plot(crop(l1[[2]],m), main = paste0("Average: ",length(E(bnsw1800[[1]]))))} else {raster::plot(lisword, main = paste0("Average: ",length(E(bnsw1800[[1]]))))}
igraph <- bnsw1800[[2]]
V(igraph)$color <- NA
V(igraph)[highv]$color <- colReds[1]
V(igraph)$label <- NA
V(igraph)[highv]$label <- as.character(dfmergesf[ind,2])


plot.igraph(igraph,
            vertex.size = 20,
            # vertex.color = "blue",
            # vertex.label = NA,
            vertex.label.color = "red",
            edge.arrow.size = 0.2,
             edge.color = NA,
            # mark.groups = idslist,
            # edge.width = E(igraph)$width,
            # lty = "points",
            main = paste0(length(E(igraph))),
            layout = as.matrix(points), add = TRUE, rescale = FALSE)


##################################################################################
# METHOD 2 (LOAD perm 1 to 5)
##################################################################################
# Bayesian list:
perm <- 1
permu <- permutations[[perm]]
i <-16
igraph <- BNsperms[[perm]][[i]]
betw <- BNbetwperms[[perm]][,i] 
# betw2 <- betw[order(betw, decreasing = TRUE)]
igraph <- set_vertex_attr(igraph, "bet", value = betw)


load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm2sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm3sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm4sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm5sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm2weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm3weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm4weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm5weights.rda")

bns1600 <- list(perm1sort[[16]],perm2sort[[16]],perm3sort[[16]],perm4sort[[16]],perm5sort[[16]])
bns1300 <- list(perm1sort[[13]],perm2sort[[13]],perm3sort[[13]],perm4sort[[13]],perm5sort[[13]])
bns1800 <- list(perm1sort[[18]],perm2sort[[18]],perm3sort[[18]],perm4sort[[18]],perm5sort[[18]])
nodos <- nodes(bns1800[[2]])
perm1weights

########################################################################################
# METHOD 2 (LOAD BN 1-25)
########################################################################################
# Select all 1800 networks
bns1800_1 <- lapply(permsort[1:10], function(x) x[[18]])
bns1800_2 <- lapply(permsort[11:25], function(x) x[[7]])
bns1800 <- c(bns1800_1,bns1800_2)

# Custom.strength measures arc.strength as frequency in networks. Does not measure it as real.arc.strength.
# Make average network
nodos <- nodes(bns1800[[2]])
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d,rms = TRUE)
arcstrengths <- lapply(bns1800, arc.strength,data = as.data.frame(dataRMS))
avstrength <- custom.strength(networks = bns1800, nodes = nodos, weights = NULL, cpdag = FALSE, debug = FALSE)
max(avstrength$strength)
media <- averaged.network(avstrength, nodes = nodos, threshold = 0.88)
media <- cextend(media)

dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
mediafit <- bn.fit(media, dataRMS)

igraphdir <- igraph.from.graphNEL(as.graphNEL(mediafit))
igraphske<- as.undirected(igraphdir)
edgelist <- E(igraphske)
nedgeslist <- length(edgelist)

graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)

# Calculate now weights (EBIC)
mediastrenghts <- bn_to_igraph.strengths(media, perm = permutations[[i]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
mediadists <- igraph.distances(mediastrenghts, perm = permutations[[i]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
mediaweights <- igraph.weights(mediadists, type = "bn", fromdist = 2000)
  
# save(perm1sort, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1_5media.rda")
# save(gridGraphsBN, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBN.rda")

# Get features of the network (only posible for METHOD 2)
# comunidades
bn1 <- mediaweights
betind1 <- which(E(bn1)$weights==0)
E(bn1)$weights[betind1] <- 0.000000001
bnbetcom1 <- cluster_edge_betweenness(as.undirected(bn1),
                                      weights = E(bn1)$weights
)

ncom <- length(levels(factor(bnbetcom1$membership)))
mem <- membership(bnbetcom1) # posibility 1 (ALL COMUNITIES)
# mem <- cutat(bnbetcom1,16) # posibility 2 (less comunities)

memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
if(ncom <= 9){
  colRainbow <- brewer.pal(ncom,"Set1")
} else {
  colRainbow <- brewer.pal(9,"Set1")
  colRainbow<- colorRampPalette(colRainbow)(ncom)
}
memClim$Data
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter = 180, regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                     main = paste0("Comunities BN ",length(E(bn1))),
                     colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                     lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
)
plot1

igraph <- bn1
betw<- betweenness(graph = igraph, weights = E(igraph)$weights, normalized = TRUE)
igraph <- set_vertex_attr(igraph, "bet", value = betw)
hist(V(igraph)$bet,freq = TRUE,  breaks = 3)



# Plot high edge betweenness nodes 
probs <- seq(0,1,(0.01))
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
nclust <- 10
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

nclust <- 10
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
} else {
  colReds <- brewer.pal(9,"Reds")
  # colReds <- colorRampPalette(colReds, alpha = TRUE)(probints)
  colReds <- colorRampAlpha(colReds, n=nclust, alpha=0.5)}

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
