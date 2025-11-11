#######################################################################################
# Investigating topological structure BN
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
########################################################################################
# number of reversible arcs
# number of relative reversible arcs
########################################################################################
list <- perm1sort
# oef <- perm1sort$hc1_1400_1500i
# oef <- perm2sort$hc2_1700_1800i
# oef <- perm3sort$hc3_1700_1800i
# oef <- perm4sort$hc4_5000_5100i
oef <- perm5sort$hc5_1700_1800i
revarcs <- lapply(list, reversible.arcs)

length(reversible.arcs(oef))
revarcs <- lapply(list, reversible.arcs)
nrevarcs <- sapply(revarcs,nrow)
arcs <- lapply(list,arcs)
narcs <- sapply(arcs,nrow)
plot(narcs,nrevarcs, main = "number of reversible arcs")
which.max(nrevarcs)
nrelrevarcs <- nrevarcs/narcs
plot(narcs,nrelrevarcs, main = "reversible arcs / number of arcs")
###########################################################################################
# number of v structures vs arcs
# number of relative v structures vs arcs
#########################################################################################
vstructmoral <- lapply(list, vstructs, moral = TRUE)
vstructimmoral <- lapply(list,vstructs,moral = FALSE)
nvstructmoral <- sapply(vstructmoral, nrow)
nvstructimmoral <- sapply(vstructimmoral, nrow)
plot(narcs,nvstructmoral, main = "number of v-structures")
plot(narcs,nvstructimmoral, main = "number of v-structures")
nvstructmoral - nvstructimmoral
# relative v structures
degmat <- matrix(nrow = length(nodes(list[[1]])), ncol = length(list))
twopathsmat <- matrix(nrow = length(nodes(list[[1]])), ncol = length(list))
for(j in 1:length(list)){
  g <- list[[j]]
  deg<- c()
  twopaths <- c()
  for(i in 1:length(nodes(g))){
    node <- nodes(g)[i]
    deg[i] <- bnlearn::degree(g,node)
    if(deg[i]>1){
    twopaths[i] <- ncol(combn(deg[i],2))
    } else{twopaths[i] <- 0}
  }
  degmat[,j] <- deg
  twopathsmat[,j] <- twopaths
}

possible2pathslist <- colSums(twopathsmat)
relvstructlist <- nvstructimmoral/possible2pathslist
plot(narcs,relvstructlist, main = "number of v-structures / possible paths of length 2")

par(mfrow = c(2,2))
plot(narcs,nrevarcs, main = "number of reversible arcs")
plot(narcs,nrelrevarcs, main = "reversible arcs / number of arcs")
plot(narcs,nvstructimmoral, main = "number of v-structures")
plot(narcs,relvstructlist, main = "number of v-structures / possible paths of length 2")
##########################################################
# shielders local: ONE CASE ONLY
##########################################################
oef <- perm1sort$hc1_1400_1500i
names <- c()
names
for(i in 1:648) names <- append(names,paste0("V",i))
# obtain vstructs immoral
voefimmor <- vstructs(oef, moral = FALSE)
# Which of those vstructs leads to dseperation?
dsepvimmor <- logical()
for(i in 1:nrow(voefimmor)){
  answer <- dsep(oef,voefimmor[i,1],voefimmor[i,3])
  dsepvimmor[i] <- answer
}
blockedim <- which(dsepvimmor)
blockedim <- blocked
voefimmor[blockedim,]
dsep(oef,"V34","V205")

# obtain vstructs moral
voefmor <- vstructs(oef, moral = TRUE)
# Which of those vstructs leads to dseperation?
# Given shielder all are not dseperated: GOOD
# Without evidence, moral and immoral give same seperation: GOOD
dsepvmor <- logical()
for(i in 1:nrow(voefmor)){
  answer <- dsep(oef,voefmor[i,1],voefmor[i,3])
  dsepvmor[i] <- answer
}
blockedm <- which(dsepvmor)
voefmor[blockedm,]

nrow(voefmor)
nrow(voefimmor)

# Neither of the rows are equal; THis is normal: can not be moral ánd immoral. 
all.equal(voefmor,voefimmor)
j <- 1
cc <- c()
for(j in 1:nrow(voefimmor)){
  ans <- which(apply(voefimmor, 1, function(x) all.equal(voefimmor, voefmor[j,])) == "TRUE")
  cc[j] <- ans
}

# Indicate shielders/ blockers / middle nodes of v structures
shielders <- voef[,2]
voeff <- factor(shielders, levels = names)
voeff
plot(voeff)
table(voeff)
which.max(table(voeff))
or <- order(table(voeff))
table(voeff)[or]


ins <- c()
outs <- c()
allcomb <- c()
vcomb <- c()
betwcomb <- c()
for (i in 1:length(names)){
  node <- names[i]
  ins[i] <- in.degree(oef,node)
  outs[i] <- out.degree(oef,node)
  if(ins[i]+outs[i] >1){
  allcomb[i] <- ncol(combn(ins[i]+outs[i],2))
  } else allcomb[i] <- 0
  if(ins[i]>1){
  vcomb[i] <- ncol(combn(ins[i],2))
  } else vcomb[i] <- 0
  
  betwcomb[i] <- allcomb[i] -vcomb[i]
}

shields <- table(voeff)
ins
outs
allcomb
vcomb
betwcomb

# twopaths <- c()
# for (i in 1:length(names)){
#   indegree <- ins[i]
#   outdegree <- outs[i]
#   pos1 <- ncol(combn(outs,2))
#   
# }
quantile(betwcomb)
all.equal(as.vector(shields),as.vector(vcomb))


shieldsClim<- quantity2clim(quantity = shields, what = "shields", ref.grid = tas_ncep_10d) 
relshieldsClim <- quantity2clim(quantity = shields/ncol(combn(ins+outs,2)), what = "relshields", ref.grid = tas_ncep_10d) 
betwdirClim <- quantity2clim(quantity = betwcomb, what = "betwDir", ref.grid = tas_ncep_10d) 
relbetwdirClim <- quantity2clim(quantity = betwcomb/ncol(combn(ins+outs,2)), what = "relbetwDir", ref.grid = tas_ncep_10d) 
vdirClim <- quantity2clim(quantity = vcomb, what = "vDir", ref.grid = tas_ncep_10d) 
insClim <- quantity2clim(quantity = ins, what = "indegree", ref.grid = tas_ncep_10d) 
outsClim <- quantity2clim(quantity = outs, what = "outdegree", ref.grid = tas_ncep_10d) 


par(mfrow =c(3,2))
a <- spatialPlot(grid = shieldsClim, backdrop.theme = "coastline", set.min = NULL,
  set.max = NULL, lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE, main = "Conditioners")
b <- spatialPlot(grid = relshieldsClim, backdrop.theme = "coastline", set.min = NULL,
                 set.max = NULL, lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE, main = "rel Conditioners")
c <- spatialPlot(grid = betwdirClim, backdrop.theme = "coastline", set.min = NULL,
             set.max = 50,
            lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE, main = "Passers")
d <- spatialPlot(grid = relbetwdirClim, backdrop.theme = "coastline", set.min = NULL,
                 set.max = 50,
                 lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE, main = "rel Passers")
# d <- spatialPlot(grid = vdirClim, backdrop.theme = "coastline", set.min = NULL,
            # set.max = NULL, lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE)
e <- spatialPlot(grid = insClim, backdrop.theme = "coastline", set.min = NULL,
                 set.max = NULL, lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE,  main = "indegree")
f <- spatialPlot(grid = outsClim, backdrop.theme = "coastline", set.min = NULL,
                 set.max = NULL, lonCenter = 180, col.regions = colRedsamp,  rev.colors = FALSE, main = "outdegree")

do.call(grid.arrange,c(list(a,c,b,d,e,f), ncol = 2))  

max(ins)
which.max(outs)
in.degree(oef, node)
out.degree(x, node)
########################################################################
# Shielders difference in probabilities. Theoreticaly not possible. 
########################################################################
voefmor
datafit <- data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d,rms = TRUE))
names(datafit) <- names
oeffit <- bn.fit(oef,data = datafit)
oeffit$V1$fitted.values
###########################################################################################
# number of v structures vs arcs all permutations. 
# number of relative v structures vs arcs
#########################################################################################
permsorts <- list(perm1sort,perm2sort,perm3sort,perm4sort,perm5sort)
degarray <- array(dim = c(length(permsorts),length(nodes(permsorts[[1]][[1]])),length(permsorts[[1]])))
twopathsarray <- array(dim = c(length(permsorts),length(nodes(permsorts[[1]][[1]])),length(permsorts[[1]])))


for (k in 1:length(permsorts)){
 
  for(j in 1:length(permsorts[[1]])){
    g <- permsorts[[k]][[j]]
    deg<- c()
    twopaths <- c()
    for(i in 1:length(nodes(g))){
      node <- nodes(g)[i]
      deg[i] <- bnlearn::degree(g,node)
      if(deg[i]>1){
        twopaths[i] <- ncol(combn(deg[i],2))
      } else{twopaths[i] <- 0}
    }
    degarray[k,,j] <- deg
    twopathsarray[k,,j] <- twopaths
  }
}

possible2pathsmatrix<-  apply(degarray, c(1,3), sum)
possible2pathsmatrixperm <- possible2pathsmatrix[,c(8,9,10,13,16,26,50,80)]


##########################################################
# shielders local: list case in TERMINAL
##########################################################

# for(k in c(8,9,10,13,16,26,50,80)){
# 
#   oefs <- lapply(list(perm1sort,perm2sort,perm3sort,perm4sort,perm5sort), function(x) x[[k]])
#   names <- c()
#   for(m in 1:648) names <- append(names,paste0("V",m))
#   # obtain vstructs immoral
#   voefsimmor <- lapply(oefs,vstructs,  moral = FALSE)
#   # Which of those vstructs leads to dseperation?
#   dsepsvimmor <- list()
#   blockeds <- list()
#   blockedsimmor <- list()
#   for(j in 1:length(oefs)){
#     voefimmor <- voefsimmor[[j]]
#     dsepvimmor <- logical()
#     for(i in 1:nrow(voefimmor)){
#       answer <- dsep(oefs[[j]],voefimmor[i,1],voefimmor[i,3])
#       dsepvimmor[i] <- answer
#     }
#     dsepsvimmor[[j]] <- dsepvimmor
#     blockedsim <- which(dsepvimmor)
#     blockedsimmor[[j]] <- voefimmor[blockedsim,]
#   }
#   assign(paste0("dseppingvstructs",k),blockedsimmor)
#   save(list = paste0("dseppingvstructs",k),file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dseppingvstructs/dseppingvstructs",k,".rda"))
# }
# for(k in c(80)){
#   
#   oefs <- lapply(list(perm1sort,perm2sort,perm3sort,perm4sort,perm5sort), function(x) x[[k]])
#   names <- c()
#   for(m in 1:648) names <- append(names,paste0("V",m))
#   # obtain vstructs immoral
#   voefsimmor <- lapply(oefs,vstructs,  moral = FALSE)
#   # Which of those vstructs leads to dseperation?
#   dsepsvimmor <- list()
#   blockeds <- list()
#   blockedsimmor <- list()
#   for(j in 1:length(oefs)){
#     voefimmor <- voefsimmor[[j]]
#     dsepvimmor <- logical()
#     for(i in 1:nrow(voefimmor)){
#       answer <- dsep(oefs[[j]],voefimmor[i,1],voefimmor[i,3])
#       dsepvimmor[i] <- answer
#     }
#     dsepsvimmor[[j]] <- dsepvimmor
#     blockedsim <- which(dsepvimmor)
#     blockedsimmor[[j]] <- voefimmor[blockedsim,]
#   }
#   assign(paste0("dseppingvstructs",k),blockedsimmor)
#   save(list = paste0("dseppingvstructs",k),file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dseppingvstructs/dseppingvstructs",k,".rda"))
# }

#####################################################################################################
# Visualizing vstructrures that d-separate the end points. 
#####################################################################################################
  files<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dseppingvstructs", full.names = T)
  filesnames <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/dseppingvstructs")
  filesnames <- gsub(".rda", "", filesnames)
  filesnames
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))
  dseplist <- list()

  for (i in 1:length(files)){
    load(files[i])
    variablepos <- get(load(files[i]))
    dseplist[[i]]<- variablepos
  }
  
  names(dseplist) <- filesnames

dseplist$dseppingvstructs8
dseppingvstructs8
dseppingvstructs9
dseppingvstructs10
dseppingvstructs13
dseppingvstructs16
dseppingvstructs26
dseppingvstructs50
dseppingvstructs80

par(mfrow = c(2,3))
for (i in 1:length(dseppingvstructs80)){
  plot.vstruct.bn(dseppingvstructs80[[i]],TimeCoordsAnom_from_Grid_rms(tas_ncep_10d), shift = TRUE)
}
dev.off()

par(mfrow = c(2,3))
for (i in 1:length(dseppingvstructs16)){
  plot.vstruct.bn(dseppingvstructs50[[i]],TimeCoordsAnom_from_Grid_rms(tas_ncep_10d), shift = TRUE)
}
dev.off()
#########################################################################################
# relative d seperating vstructs
#########################################################################################
permsorts <- list(perm1sort,perm2sort,perm3sort,perm4sort,perm5sort)
tamanos <- list(dseppingvstructs8, dseppingvstructs9,dseppingvstructs10, dseppingvstructs13, dseppingvstructs16, dseppingvstructs26, dseppingvstructs50,dseppingvstructs80)
  matx <- matrix(nrow = length(permsorts), ncol = length(tamanos))
  maty <- matrix(nrow = length(permsorts), ncol = length(tamanos))
  maty2 <-  matrix(nrow = length(permsorts), ncol = length(tamanos))
  matfout <- matrix(nrow = length(permsorts), ncol = length(tamanos))

for(i in 1:length(permsorts)){
  permv <- lapply(tamanos, function(x) x[[i]])
  maty[i,]<- sapply(permv, nrow)
  matx[i,]<- sapply(permsorts[[i]][c(8,9,10,13,16,26,50,80)],narcs)
  matfout[i,] <- maty[i,]/matx[i,]
  maty2[i,] <- maty[i,]/possible2pathsmatrixperm[i,]
}
plot(matx,maty, col = rainbow(5), xlab = "|E|",  ylab = "dsepvstructs")
dev.off()
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen6/figures/numberdsepvstructs.pdf"
pdf(plotname)
plot(matx,maty2, col = rainbow(5), xlab = "|E|", ylab = "dsepvstructs/All 2 paths", main = "relative number of dseparating vstructures")
legend(x = "topright",legend = c("perm 1","perm 2","perm 3","perm 4","perm 5"), fill =  rainbow(5),col = rainbow(5))
dev.off()
plot(matx,matfout, col = rainbow(5), xlab = "|E|",  ylab = "dsepvstructs/|E|")


#####################################################################
# overlpping algorithm based on cliques
#####################################################################
clique.community <- function(graph, m) {
  graph <- cn0
  k <- NULL
  clq <- cliques(graph, min=k, max=k)
  k <- m
  edges <- c()
  for (i in seq_along(clq)) {
    for (j in seq_along(clq)) {
      if ( length(unique(c(clq[[i]], clq[[j]]))) == k+1 ) {
        edges <- c(edges, c(i,j))
      }
    }
  }
  clq.graph <- simplify(graph(edges))
  V(clq.graph)$name <- seq_len(vcount(clq.graph))
  comps <- decompose.graph(clq.graph)
  
  lapply(comps, function(x) {
    unique(unlist(clq[ V(x)$name ]))
  })
}
cliques(bn4)
clique.community(bn1,1)
clique.commun
###################################################################################################
# comunidades. 
# Use weighted igraphs: BN:
###################################################################################################
arcs1 <- lapply(perm1sort, arcs)
narcs1 <- sapply(arcs1, nrow)

# train tamanos
nedgeslistbn <- narcs1[c(10,27)]
nedgeslistbn

bnsel1 <- 10
bnsel2 <- 27


# used in borrador
nedgeslistbn <- narcs1[c(18,26,87)]
nedgeslistbn

bnsel1 <- 18
bnsel2 <- 26
bnsel3 <- 87


nedgeslistbn <- narcs1[c(13,50,86,18)]

bnsel1 <- 13
bnsel2 <- 50
bnsel3 <- 86
bnsel4 <- 18


bn1 <- perm1weights[[bnsel1]]
bn2 <- perm1weights[[bnsel2]]
bn3 <- perm1weights[[bnsel3]]
bn4 <- perm1weights[[bnsel4]]

betind1 <- which(E(bn1)$weights==0)
betind2 <- which(E(bn2)$weights==0)
betind3 <- which(E(bn3)$weights==0)
betind4 <- which(E(bn4)$weights==0)

E(bn1)$weights[betind1] <- 0.000000001
E(bn2)$weights[betind2] <- 0.000000001
E(bn3)$weights[betind3] <- 0.000000001
E(bn4)$weights[betind4] <- 0.000000001

bnbetcom1 <- cluster_edge_betweenness(as.undirected(bn1),
                        weights = E(bn1)$weights
                       )

bnbetcom2 <- cluster_edge_betweenness(as.undirected(bn2),
                                      weights = E(bn2)$weights
                        )
bnbetcom3 <- cluster_edge_betweenness(as.undirected(bn3),
                                      weights = E(bn3)$weights
)

ncom <- length(levels(factor(bnbetcom1$membership)))
mem <- membership(bnbetcom1)
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





bnlist <- list(bn1,bn2,bn3,bn4)
# THis one used in borrador
bnlist <- list(bn1,bn2,bn3)
bnlist
# Train tamanos 
bnlist <- list(bn1,bn2)
commu_betsbn <- list()
plotsbn <- list()

# weighted
for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weights
  )
  commu_betsbn[[i]]<- commu_bet
}

# Unweighted
for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weight
  )
  commu_betsbn[[i]]<- commu_bet
}

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
#############################################################################################
# comunidades.
# corelation network.
#############################################################################################
# load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
# gridGraphsCN <- gridGraphs
# rm(gridGraphs)
# str(gridGraphsCN[[1]])
# GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
# edgelistsCN <- lapply(GraphsCN, E)
# numberofedgesCN <- sapply(edgelistsCN, length)
# names(gridGraphsCN) <- as.character(numberofedgesCN)
# listcn <- GraphsCN[30:45]
# narcs <- numberofedgesCN[30:45]
# 
# listcnstrenghts <- lapply(listcn, str_vs_dis_2, perm = NULL, data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
# # save(listcnstrengths, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
# listcnweights <- lapply(listcnstrenghts, igraph.weights, type = "cn", fromdist = 0, perm = permutations[[1]])
# # save(listcnweights, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")


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
nedgescnlist <- numberofedgesCN[c(69,42,35,32,3)]
cn00 <- listcnweights[[3]]
cn1 <- listcnweights[[69]]
cn2 <- listcnweights[[42]]
cn3 <- listcnweights[[35]]
cn4 <- listcnweights[[32]]

# Train tamanos
nedgescnlist <- numberofedgesCN[c(78,53)]
cn1 <- listcnweights[[78]]
cn2 <- listcnweights[[53]]

# This one is used in borrador. 999 of betweenness not yet. +78 + 53
nedgescnlist <- numberofedgesCN[c(61,56,35,32,1)]
cn1 <- listcnweights[[61]]
cn2 <- listcnweights[[56]]
cn3 <- listcnweights[[35]]
cn4 <- listcnweights[[32]]

cnlist <- list(cn1,cn2,cn3,cn4,cn00)
cnlist <- list(cn1,cn2,cn3,cn4)
cnlist<- list(cn1,cn2)
commu_betscn <- list()
plotscn <- list()

for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                           weights = E(g)$weights
  )
  commu_betscn[[i]]<- commu_bet}

# unweighted
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weight
  )
  commu_betscn[[i]]<- commu_bet}

# inverse weighted
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  ind <- which(E(g)$weights==0)
  E(g)$weights[ind] <- 0.000000001
  commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = 1/E(g)$weights
  )
  commu_betscn[[i]]<- commu_bet}



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
##############################################################
#  join;
##############################################################
blank <- grid.rect(gp=gpar(col="white"))
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_betw_all.pdf"
pdf(plotname)
plotsall <- list(plotsbn[[1]],plotscn[[1]],plotsbn[[2]],plotscn[[2]],plotsbn[[3]],plotscn[[3]],plotsbn[[4]],plotscn[[4]])
plotsall <- list(plotsbn[[1]],plotscn[[1]],plotsbn[[2]],plotscn[[2]],plotsbn[[3]],plotscn[[3]])
do.call(grid.arrange,c(plotsall,nrow = 4, top = "edge betweenness"))
dev.off()
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_betw.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_betw_1800.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen6/figures/commu.betw.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_betw_CN_optimums.pdf"
pdf(plotname)
pdf(plotname, height = 4, width = 11)
do.call(grid.arrange,c(list(plotsbn[[4]],plotscn[[3]]),nrow = 1, top = "edge betweenness"))
dev.off()
##############################################################
# infomap Not a hierarchical community structure!!
##############################################################
# bn
###############################################################
commu_infosbn <- list()
plotsinfosbn <- list()

# cluster_infomap(graph, e.weights = NULL, v.weights = NULL, nb.trials = 10,
                # modularity = TRUE)

cutoff = 15

for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  commu_info <- cluster_infomap(as.undirected(g), e.weights = E(g)$weights)
  commu_infosbn[[i]]<- commu_info
}
for(i in 1:length(bnlist)){
  ncom <- length(levels(factor(commu_infosbn[[i]]$membership)))
  mem <- membership(commu_infosbn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)

  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsinfosbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                              set.min = NULL, set.max = NULL, 
                              lonCenter = 180, 
                              regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                              main = paste0("Comunities BN ",nedgeslistbn[i]),
                              colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                              lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}
do.call(grid.arrange, plotsinfosbn)
##############################################################
# cn
###############################################################
commu_infoscn <- list()
plotsinfoscn <- list()

for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  commu_info <- cluster_infomap(as.undirected(g), e.weights = E(g)$weights)
  commu_infoscn[[i]]<- commu_info}
for (i in 1:length(cnlist)){
  ncom <- length(levels(factor(commu_infoscn[[i]]$membership)))
  mem <- membership(commu_infoscn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsinfoscn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                              set.min = NULL, set.max = NULL, 
                              lonCenter = 180, 
                              regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                              main = paste0("Comunities CN ",nedgescnlist[i]),
                              colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                              lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, plotsinfoscn)  
##############################################################
#  join infomap
##############################################################
plotsinfosall <- list(plotsinfosbn[[1]],plotsinfoscn[[1]],plotsinfosbn[[2]],plotsinfoscn[[2]],plotsinfosbn[[3]],plotsinfoscn[[3]],plotsinfosbn[[4]],plotsinfoscn[[4]])
do.call(grid.arrange,c(plotsinfosall, nrow = 4, top = "infomap"))


##############################################################
# random walktrap
##############################################################
# bn
###############################################################
commu_fgsbn <- list()
plotsfgsbn <- list()


for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  commu_fg <- cluster_walktrap(as.undirected(g), weights = E(g)$weights, steps =3)
  commu_fgsbn[[i]]<- commu_fg
}
for(i in 1:length(bnlist)){
  ncom <- length(levels(factor(commu_fgsbn[[i]]$membership)))
  mem <- membership(commu_fgsbn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsfgsbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                   set.min = NULL, set.max = NULL, 
                                   lonCenter = 180, 
                                   regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                   main = paste0("Comunities BN ",nedgeslistbn[i]),
                                   colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                   lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}
do.call(grid.arrange, plotsfgsbn)
##############################################################
# cn
###############################################################
commu_fgscn <- list()
plotsfgscn <- list()

for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  commu_fg <- cluster_walktrap(as.undirected(g), weights = E(g)$weights, steps = 4)
  commu_fgscn[[i]]<- commu_fg}
for (i in 1:length(cnlist)){
  ncom <- length(levels(factor(commu_fgscn[[i]]$membership)))
  mem <- membership(commu_fgscn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsfgscn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                   set.min = NULL, set.max = NULL, 
                                   lonCenter = 180, 
                                   regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                   main = paste0("Comunities CN ",nedgescnlist[i]),
                                   colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                   lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, c(plotsfgscn, top = "CN random walk with walktrap"))  
##############################################################
#  join fgs
##############################################################
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_walktrap_all.pdf"
pdf(plotname)
plotsfgsall <- list(plotsfgsbn[[1]],plotsfgscn[[1]],plotsfgsbn[[2]],plotsfgscn[[2]],plotsfgsbn[[3]],plotsfgscn[[3]],plotsfgsbn[[4]],plotsfgscn[[4]])
do.call(grid.arrange,c(plotsfgsall,nrow = 4, top = "random walk with walktrap"))
dev.off()
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_walktrap.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/Resumen6/figures/commu.walktrap.pdf"
pdf(plotname, height = 4, width = 11)
do.call(grid.arrange,c(list(plotsfgsbn[[4]],plotsfgscn[[2]]),nrow = 1, top = "random walk with walktrap"))
dev.off()


##############################################################
# Eigen vector clustering
##############################################################
# bn Failed to converge
###############################################################
commu_eigensbn <- list()
plotseigensbn <- list()
cutoff = NULL

for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  commu_eigenbn <- cluster_leading_eigen(as.undirected(g), weights = E(g)$weights)
  commu_eigensbn[[i]]<- commu_eigenbn
}
for(i in 1:length(bnlist)){
  ncom <- length(levels(factor(commu_eigensbn[[i]]$membership)))
  mem <- membership(commu_eigensbn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)

  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotseigensbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                 set.min = NULL, set.max = NULL, 
                                 lonCenter = 180, 
                                 regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                 main = paste0("Comunities BN ",nedgeslistbn[i]),
                                 colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                 lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}
do.call(grid.arrange, c(plotseigensbn, nrow = 4, top = "Eigen vector"))
##############################################################
# cn
###############################################################
commu_eigenscn <- list()
plotseigenscn <- list()
i <- 4
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  commu_eigen <- cluster_leading_eigen(as.undirected(g), weights = E(g)$weights)
  commu_eigenscn[[i]]<- commu_eigen}
for (i in 1:length(cnlist)){
  ncom <- length(levels(factor(commu_eigenscn[[i]]$membership)))
  mem <- membership(commu_eigenscn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotsfgscn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                 set.min = NULL, set.max = NULL, 
                                 lonCenter = 180, 
                                 regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                 main = paste0("Comunities CN ",nedgescnlist[i]),
                                 colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                 lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, plotseigenscn)  



##############################################################
# Label propagation clustering
##############################################################
# bn is not a hierarchical community structure!
###############################################################
commu_labelsbn <- list()
plotslabelsbn <- list()

for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  commu_labelbn <- cluster_label_prop(as.undirected(g), weights = E(g)$weights)
  commu_labelsbn[[i]]<- commu_labelbn
}

for(i in 1:length(bnlist)){
  ncom <- length(levels(factor(commu_labelsbn[[i]]$membership)))
  mem <- membership(commu_labelsbn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)

  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotslabelsbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                    set.min = NULL, set.max = NULL, 
                                    lonCenter = 180, 
                                    regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                    main = paste0("Comunities BN ",nedgeslistbn[i]),
                                    colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                    lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}
do.call(grid.arrange, c(plotslabelsbn, nrow = 4, top = "label propagation"))
##############################################################
# cn
###############################################################
commu_labelscn <- list()
plotslabelscn <- list()
i <- 4
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  commu_label <- cluster_label_prop(as.undirected(g), weights = E(g)$weights)
  commu_labelscn[[i]]<- commu_label}
for (i in 1:length(cnlist)){
  ncom <- length(levels(factor(commu_labelscn[[i]]$membership)))
  mem <- membership(commu_labelscn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotslabelscn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                 set.min = NULL, set.max = NULL, 
                                 lonCenter = 180, 
                                 regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                 main = paste0("Comunities CN ",nedgescnlist[i]),
                                 colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                 lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, c(plotslabelscn, top = "CN label propagation"))
##########################################################################
# Join label propagations
########################################################################## 
blank <- grid.rect(gp=gpar(col="white"))
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_label_prop.pdf"
pdf(plotname)

plotslabelsall <- list(plotslabelsbn[[4]],plotslabelscn[[3]],blank,plotslabelscn[[4]],blank,plotslabelscn[[5]])
do.call(grid.arrange,c(plotslabelsall,ncol = 2))
dev.off()




##############################################################
# Louvain clustering not hierarchical. 
##############################################################
# bn 
###############################################################
commu_louvainsbn <- list()
plotslouvainsbn <- list()

for (i in 1:length(bnlist)){
  g <- bnlist[[i]]
  commu_louvainbn <- cluster_louvain(as.undirected(g), weights = E(g)$weights)
  commu_louvainsbn[[i]]<- commu_louvainbn
}

for(i in 1:length(bnlist)){
  ncom <- length(levels(factor(commu_louvainsbn[[i]]$membership)))
  mem <- membership(commu_louvainsbn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)

  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotslouvainsbn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                    set.min = NULL, set.max = NULL, 
                                    lonCenter = 180, 
                                    regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                    main = paste0("Comunities BN ",nedgeslistbn[i]),
                                    colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                    lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}
do.call(grid.arrange, c(plotslouvainsbn, nrow = 4, top = "BN louvain"))
##############################################################
# cn
###############################################################
commu_louvainscn <- list()
plotslouvainscn <- list()
i <- 4
for (i in 1:length(cnlist)){
  g <- cnlist[[i]]
  commu_louvain <- cluster_louvain(as.undirected(g), weights = E(g)$weights)
  commu_louvainscn[[i]]<- commu_louvain}
for (i in 1:length(cnlist)){
  ncom <- length(levels(factor(commu_louvainscn[[i]]$membership)))
  mem <- membership(commu_louvainscn[[i]])
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plotslouvainscn[[i]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                                    set.min = NULL, set.max = NULL, 
                                    lonCenter = 180, 
                                    regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                                    main = paste0("Comunities CN ",nedgescnlist[i]),
                                    colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                                    lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )
}


do.call(grid.arrange, c(plotslouvainscn, top = "CN louvain"))



##############################################################
#  join multilevels (louvian)
##############################################################
blank <- grid.rect(gp=gpar(col="white"))
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/communities/commu_louvain.pdf"
pdf(plotname)

plotslouvainsall <- list(plotslouvainsbn[[4]],plotslouvainscn[[3]],blank,plotslouvainscn[[4]],blank,plotslouvainscn[[5]])
do.call(grid.arrange,c(plotslouvainsall,ncol = 2))
dev.off()


#################################################################
# large CN communities
#################################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")


load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda") # <- no es gridGraph! list cn weights
edgelistsCN <- lapply(listcnweights, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(listcnweights) <- as.character(numberofedgesCN)

nedgescnlist <- numberofedgesCN[c(69,42,35,32)]
cn00 <- listcnweights[[10]]
g <- cn00
length(E(g))

ind <- which(E(g)$weights==0)
E(g)$weights[ind] <- 0.000000001
commu_bet <- cluster_edge_betweenness(as.undirected(g),
                                        weights = E(g)$weights
)
assign(paste0("commu_bet_",length(E(g))),g)
name <- paste0("commu_bet_",length(E(g)))

save(list = name,file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/",name,".rda"))

commu_bet_est <- estimate_edge_betweenness(as.undirected(g), weights = E(g)$weights, cutoff = -1)
bet <- edge_betweenness(as.undirected(g), weights = E(g)$weights)
#######################################################################
# Less merges / cut communities: Used to generate evolution plots
#######################################################################
commu_bet$merges
commu_betscn <- list()
i <- 1
commu_betscn[[i]]<- commu_bet
commu_bet12 <- cut_at(commu_betsbn[[1]],9)


par(mfrow = c(1,3))
k <-1
plots <- list()
cuts <- c(3,5,9,56,65,80)
cuts <- c(3,5,7,9,12,14,16)
cuts <- 1:16
cuts <- c(44,50,56,65,80)
cuts <- c(56)
cuts <- c(5,6,9,20,30,60)
cuts <- c(4,8,12,16,20,30)
for(j in 1:length(cuts)){
i <- cuts[[j]]
commu_bet12 <- cut_at(commu_betsbn[[k]],i)
commu_bet12
g <- bnlist[[k]]


  ncom <- length(levels(factor(commu_bet12)))
  mem <- commu_bet12
  memClim <- quantity2clim(mem, "membership", tas_ncep_10d)
  if(ncom <= 9){
    colRainbow <- brewer.pal(ncom,"Set1")
  } else {
    colRainbow <- brewer.pal(9,"Set1")
    colRainbow<- colorRampPalette(colRainbow)(ncom)
  }
  if(ncom > 15){colRainbow <- sample(colRainbow,length(colRainbow))}
  plots[[j]] <- spatialPlot(grid = memClim, backdrop.theme = "coastline", 
                              set.min = NULL, set.max = NULL, 
                              lonCenter = 180, 
                              regions = TRUE,col.regions = colRainbow, at = 0:ncom, rev.colors = FALSE, 
                              main = paste0("Comunities BN ",length(E(g))),
                              colorkey = list(col = colRainbow, width = 0.6, at = 0:ncom,
                                              lables = list(cex = 0.5, labels =as.character(0:ncom),at = 0:ncom))
  )

}
plots
do.call(grid.arrange,c(plots))


