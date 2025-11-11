# install.packages(pkgs = "/home/catharina/Documents/Installations/modified_bnlearn.tar.gz", repos = NULL,type ='source')
# install.packages(pkgs = "/home/catharina/Documents/Installations/transformeR-1.2.0.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(igraph)
library(transformeR)
library(gridExtra)
library(visualizeR)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
library(mopa)
data(wrld)

# devtools::install_github('SantanderMetGroup/visualizeR')
######################################################
# Investigate alphas PC
#######################################################
View(loglik10d_19_17)
View(loglik10d_17_15)
View(loglik10d_21_13)
library(BiocGenerics)
x <- BiocGenerics::unique(loglik10d_21_13)
ncol(loglik10d_21_13)
ncol(x)
head(loglik10d_21_13)
########################################################
# make range alphas: 603 below 1 e -13, 
#######################################################
myalphasrest <- c(2:10 %o% 10^(-14:0))

length(myalphasrest) # 135

myalphas603 <- c(2:10 %o% 10^(-80:-14))
myalphas603
length(myalphas603) # 603

myalphaslow <- c(2:10 %o% 10^(-109:-81))
myalphaslow
length(myalphaslow) # 261

myalphaslowlow <- c(2:10 %o% 10^(-130:-110))
length(myalphaslowlow) # 189

myalphaslowlowlow <- c(2:10 %o% 10^(-150:-131)) # NOT YET -200 still give dags. (15 edges)
length(myalphaslowlowlow)

c(myalphaslowlow, myalphaslow)
######################################################################################
# calculate values PC algorithm for different alphas range.
######################################################################################
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
corrspear10d <- cor(data10d, method = "spearman")
loglik10d_rangehigh <- sapply(myalphasrest, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_range603 <- sapply(myalphas603, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_rangelow <- sapply(myalphaslow, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_rangelowlow <- sapply(myalphaslowlow, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_all <- sapply(c(myalphaslowlow,myalphaslow,myalphas603), alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_all[1:5,1:5]
######################################################################################
# Plot values PC algorithm vs HC for different alphas range
#######################################################################################
plot.nedge.vs.loglik_combi(loglik10d_range603,hc_edges_loglik_10d_1000i)
plot.nedge.vs.loglik_combi(loglik10d_rangelow,hc_edges_loglik_10d_1000i)
dev.off()
# loglik10d_range603[,seq(1,600,3)]
plot.nedge.vs.loglik_combi(cbind(loglik10d_rangelowlow,loglik10d_rangelow,loglik10d_range603),hc_edges_loglik_10d_1000i)
plot.nedge.vs.loglik_combi(t(arcset3),hc_edges_loglik_10d_1000i)
#######################################################################################################
# Convert pc-edge-loglik-datasets to dataset with correct number of undirected arcs.
# (As the amount of undirected arcs was calculated by length(undirected.arcs(cpdag)))
# And the correct formula is nrow(undirected.arcs(cpdag))/2 (Each undirected edge appears twice)
# We have to divide the undirected arcs by 4. (2 to correct "length", 2 to correct "dubbel appearancy")
#######################################################################################################
arcset <- cbind(loglik10d_rangelowlow,loglik10d_rangelow,loglik10d_range603)
arcset[,1:5]
arcset2 <- t(arcset)
arcset2

nrow(arcset2)
ncol(arcset2)

arcset2[,1]
unlist(arcset2[,1])
z <- matrix(data = NA, nrow = nrow(arcset2),ncol = ncol(arcset2))
for(i in 1:ncol(arcset2)) {
  z[,i] <- unlist(arcset2[,i])
}
z
z[,2]
z[,2]/4
arcset3 <- z
arcset3[,2] <- z[,2]/4
arcset3
colnames(arcset3) <- c("alpha",
                    "nundeqclass",
                    "nedgeseqclass",
                    "nedgesdag",
                    "likgbn")
pcstructure <- arcset3
save(pcstructure, file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/PCstructuredata.rda")
rm(arcset3,pcstructure)
load(file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/PCstructuredata.rda")

# for introducing in plot Combi use:
t(arcset3)
###############################################################################
# check arcset3
###############################################################################
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
class(data10d)
nrow(z)
z[1000:1053,]
# 8e-111  164   67   67 -323486.4
# 4e-86  300  140  140 -303983.6
# 6e-64  280  244  244 -281205.7
# 8e-42  332  399  399 -251108.5
# 1e-19  612  551  551 -227449.0
# 4e-18  672  564   NA        NA
# 3e-16  672  582   NA        NA
# 1e-13  772  603  603 -222533.7

alpha_loglik_pc3_withdag(data = data10d, alpha = 8e-111 ) #inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 2e-130) #inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 4e-86) #inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 6e-64) #inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 8e-42) #inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 1e-19) #inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 4e-18) #inderdaad doet het niet
alpha_loglik_pc3_simple(data = data10d, alpha = 4e-18, method = "spearman") # inderdaad delen door 4
alpha_loglik_pc3_simple(data = data10d, alpha = 3e-16, method = "spearman") # inderdaad delen door 4
alpha_loglik_pc3_withdag(data = data10d, alpha = 1e-13, method = "spearman") # inderdaad delen door 4
alpha_loglik_pc3_simple(data = data10d, alpha = 2e-13, method = "spearman") # inderdaad delen door 4
alpha_loglik_pc3_simple(data = data10d, alpha = 2e-10, method = "spearman") # inderdaad delen door 4

# check the function
data <- data10d
data <- data10d
alpha <- 2e-10

if (is.null(correlation)) {correlation <- cor(data, method = method)} 
degreestat <- list(C = correlation, n = nrow(data))
eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
amata <- as(eqclass,"amat")
bn.eqclass <- as.bn(eqclass@graph)
bn.eqclass
nedgeseqclass <- narcs(bn.eqclass)
undirected.arcs(bn.eqclass)
nrow(undirected.arcs(bn.eqclass))
nundeqclass <- nrow(undirected.arcs(bn.eqclass))/2
nundeqclass

if(isValidGraph(amata, type = "pdag")){
  dag <- pcalg::pdag2dag(eqclass@graph)
  bn.dag <- as.bn(dag$graph)
  nedgesdag <- narcs(bn.dag)
  
  newdata <- as.data.frame(data)
  names(newdata) <- nodes(bn.dag)
  gbn <- bn.fit(bn.dag, newdata)
  likgbn <- logLik(gbn, newdata)
}
else {likgbn <- NA
nedgesdag <- NA}


return(list(networkdata = data.frame(alpha,
                                     nundeqclass,
                                     nedgeseqclass,
                                     nedgesdag,
                                     likgbn), networks = list(eqclass = bn.eqclass, dag = bn.dag)))




#######################################################################################################
# New loglik - plots: extended version of those in resumen script
#######################################################################################################
plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Plotloglikpc1000hc.pdf"
pdf(plotname)
plot.nedge.vs.loglik_combi(cbind(loglik10d_rangelowlow,loglik10d_rangelow,loglik10d_range603),hc_edges_loglik_10d_1000i)
abline(v = 648, col = "red")
dev.off()

plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Plotloglikpc8695hc.pdf"
pdf(plotname)
plot.nedge.vs.loglik_combi(cbind(loglik10d_rangelowlow,loglik10d_rangelow,loglik10d_range603),hc_edges_loglik_10d_8695i)
abline(v = 648, col = "red")
dev.off()


save(loglik10d_15_13,loglik10d_17_15,loglik10d_19_17,loglik10d_21_19,loglik10d_rangelowlow,loglik10d_rangelow,loglik10d_range603,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_pcnetworks10d.rda")

##############################################################################################
# Compare HC and PC 
# With long distances.
###############################################################################################
plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HCvsPCminimdistance2000.pdf"
pdf(plotname)
10 * pi*6000/180
par(mfrow = c(1,2), cex = 0.5)
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plot_long_distances(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 2000)
plot_long_distances(pc_10d_1e_13$networks$dag, data10d, minimdist = 2000)
dev.off()



plot_long_distances(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 3000)
plot_long_distances(pc_10d_1e_13$networks$dag, data10d, minimdist = 6741)

##########################################################################################
# Find maximum distance in particular graph.
##########################################################################################
# dag, data.dag, minimdist
dag <- hc_edges_loglik_10d_400_600i$networks
dag <- pc_10d_1e_13$networks$dag
data.dag <- data10d

igraph <- igraph.from.graphNEL(as.graphNEL(dag))
adjmat <- as.matrix(as_adjacency_matrix(igraph))

nodenames <- nodes(dag) 
length(nodenames)

edgeindices <- which(adjmat == 1, arr.ind = TRUE)
fromV <- nodenames[edgeindices[,1]]
toV <- nodenames[edgeindices[,2]]
fromtoV <- cbind(fromV,toV)

longitude <- attributes(data.dag)$VertexCoords$x
lattitude <- attributes(data.dag)$VertexCoords$y

distances <- c()
for (i in 1:nrow(edgeindices)){
  x1Lat <- lattitude[edgeindices[i,1]]
  x1Lon <- longitude[edgeindices[i,1]]
  x2Lat <- lattitude[edgeindices[i,2]]
  x2Lon <- longitude[edgeindices[i,2]]
  
  disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
  distances[i] <- disti}

arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
max(as.numeric(arcdistances[,3]))

#################################################################################################
# Visualize teleconnections HC
#################################################################################################
hc_list <- list(#hc_edges_loglik_10d_600_800i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1400_1600i,
                hc_edges_loglik_10d_1800_2000i)

# select networks
networks <- lapply(hc_list, function(m) m[["networks"]])
nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]
minimdist = 10000
#calculate area weighted dag degree for all hillclimbing dags
plotname <- paste0("/home/catharina/Documents/PlotsResumen/HClongdists",firstel,"_",lastel,"_",minimdist,"km.pdf")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HClongdists",firstel,"_",lastel,"_",minimdist,"km.pdf")
plotname
pdf(plotname)
par(mfcol = c(3,1),cex = 0.5)
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
lapply(networks, plot_long_distances_3, data.dag = data10d, minimdist = minimdist, type = "HC")
dev.off()

###################################################################################################
# Draw graphs PC and HC 200/400/600 and draw longest edges in 600. (*)
###################################################################################################
load(file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Data_graphsPC.rda")
# In the following also include the graph of 400 edges of PC
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plotname <- "/home/catharina/Documents/PlotsResumen/HCvsPCnormalgraphs.pdf"
plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HCvsPCnormalgraphswithlarge.pdf"
pdf(plotname)
par(mfcol = c(3,2), cex = 0.5)
plot.Meteodag(data10d, hc_edges_loglik_10d_200i$networks)
title(main = paste0("dag from HC \n number of edges = ",narcs(hc_edges_loglik_10d_200i$networks)))
plot.Meteodag(data10d, hc_edges_loglik_10d_200_400i$networks)
title(main = paste0("dag from HC \n number of edges = ",narcs(hc_edges_loglik_10d_200_400i$networks)))
plot_long_distances_2(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 3000, type = "HC")


plot.Meteodag(data10d, pc_10d_4e_70$networks$dag)
title(main = paste0("dag from PC \n number of edges = ",narcs(pc_10d_4e_70$networks$dag)))
plot.Meteodag(data10d, pc_10d_12e_42$networks$dag)
title(main = paste0("dag from PC \n number of edges = ",narcs(pc_10d_12e_42$networks$dag)))
plot_long_distances_2(pc_10d_1e_13$networks$dag, data10d, minimdist = 3000, type = "PC")
dev.off()

plot_long_distances_2(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 3000, type = "HC")
plot_long_distances_2(pc_10d_1e_13$networks$dag, data10d, minimdist = 3000, type = "PC")

#################################################################################################
# Visualize one graph
#################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda")
# visualize 1 graph: with plot.Meteodag (time.coords only uses X and Y coords attributes)
dag <- hc_edges_loglik_10d_1400_1600i$networks
numberofnodes <- length(nodes(dag))
numberofedges <- narcs(dag)

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Bayesiangraph",narcs(dag),".pdf")
pdf(plotname)

edgedensity <- round(numberofedges/(choose(numberofnodes,2)), digits = 4)
TimeCoordsAnom_from_Grid(tas_ncep_10d)
plot.Meteodag(dag, time.coords = TimeCoordsAnom_from_Grid(tas_ncep_10d))
title(main = paste0("dag from HC \n number of edges = ", numberofedges))
dev.off()


###########################################################################################
# New long distance function like plotmeteodag for coloring largest link. (*)
#########################################################################################
plot_long_distances_2 <- function(dag, data.dag, minimdist,type){
  
  
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  
  nodenames <- nodes(dag) 
  length(nodenames)
  
  edgeindices <- which(adjmat == 1, arr.ind = TRUE)
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  longitude <- attributes(data.dag)$VertexCoords$x
  lattitude <- attributes(data.dag)$VertexCoords$y
  
  distances <- c()
  for (i in 1:nrow(edgeindices)){
    x1Lat <- lattitude[edgeindices[i,1]]
    x1Lon <- longitude[edgeindices[i,1]]
    x2Lat <- lattitude[edgeindices[i,2]]
    x2Lon <- longitude[edgeindices[i,2]]
    
    disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
    distances[i] <- disti}
  
  arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
  arcdistances
  
  as_edgelist(igraph)[,1]
  arcdistances[,1]
  
  #identify which indices igraph correspond to indices in bn
  indsame <- c()
  
  for (i in 1:nrow(arcdistances)){
    int <- intersect(which(as_edgelist(igraph)[,1] == arcdistances[i,1]),
                     which(as_edgelist(igraph)[,2] == arcdistances[i,2]))
    indsame[i] <- int
  }
  indsame
  
  #create distance vector in igraph object
  newdistances <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    newdistances[ind] <- distances[i]
  }
  newdistances
  
  E(igraph)$distances <- newdistances
  E(igraph)[E(igraph)$distances<= minimdist]$color <- "green"
  E(igraph)[E(igraph)$distances > minimdist]$color <- "red"
  
  #maak het plaatje:
  x <- attr(data.dag, "Xcoords", exact = FALSE)
  y <- attr(data.dag, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  data(wrld)
  plot(wrld)
  plot.igraph(igraph, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              #edge.lty = 1+(E(igraph)$normweights)*2,
              main = paste0("dag from ",type,"\n number of edges = ",narcs(dag)),
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  title(main = paste0("dag from ",type,"\n number of edges = ", narcs(dag)))
}
plot_long_distances_2(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 3000, type = "HC")
plot_long_distances_2(pc_10d_1e_13$networks$dag, data10d, minimdist = 3000, type = "PC")



plot_long_distances_3 <- function(dag, data.dag, minimdist,type){
  
  
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  
  nodenames <- nodes(dag) 
  length(nodenames)
  
  edgeindices <- which(adjmat == 1, arr.ind = TRUE)
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  longitude <- attributes(data.dag)$VertexCoords$x
  lattitude <- attributes(data.dag)$VertexCoords$y
  
  distances <- c()
  for (i in 1:nrow(edgeindices)){
    x1Lat <- lattitude[edgeindices[i,1]]
    x1Lon <- longitude[edgeindices[i,1]]
    x2Lat <- lattitude[edgeindices[i,2]]
    x2Lon <- longitude[edgeindices[i,2]]
    
    disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
    distances[i] <- disti}
  
  arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
  arcdistances
  
  as_edgelist(igraph)[,1]
  arcdistances[,1]
  
  #identify which indices igraph correspond to indices in bn
  indsame <- c()
  
  for (i in 1:nrow(arcdistances)){
    int <- intersect(which(as_edgelist(igraph)[,1] == arcdistances[i,1]),
                     which(as_edgelist(igraph)[,2] == arcdistances[i,2]))
    indsame[i] <- int
  }
  indsame
  
  #create distance vector in igraph object
  newdistances <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    newdistances[ind] <- distances[i]
  }
  newdistances
  
  E(igraph)$distances <- newdistances
  E(igraph)[E(igraph)$distances<= minimdist]$color <- "NA"
  E(igraph)[E(igraph)$distances > minimdist]$color <- "red"
  
  #maak het plaatje:
  x <- attr(data.dag, "Xcoords", exact = FALSE)
  y <- attr(data.dag, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  data(wrld)
  plot(wrld)
  plot.igraph(igraph, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              #edge.lty = 1+(E(igraph)$normweights)*2,
              main = paste0("dag from ",type,"\n number of edges = ",narcs(dag)),
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  title(main = paste0("dag from ",type,"\n number of edges = ", narcs(dag),"\n minimum distance = ",minimdist))
}
plot_long_distances_3(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 3000, type = "HC")




##################################################################################################
# Probeersel commandos met nieuwe bic gt test marco
##################################################################################################
library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")



load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_1.rda")
str(tabu_1$networkdata)
360/1801
tabu_1$networkdata[1154,]
log10(sum(tabu_1$networkdata[1:1154,"tests"]))
log10(test$networkdata$tests[3])
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_1.rda")
alphas <- c(2:10 %o% 10^(-112:-110))
alphas

data <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
correlation <- NULL
alphas <- c(1*10^-150, 1*10^-50, 1*10^-13, 1*10^-10, 1*10^-5, 1*10^-1)
alphas <- c(2:10 %o% 10^(-112:-2))
alpha <- alphas


test <- iteration_pc(data = data, correlation = NULL, alpha = alphas)
test$networks$pc_10d_0.2
test$networkdata


try <- mmpc(gaussian.test, test = "bic-gt")
arcs(try)
skeleton(gaussian.test)
bnlearn:::elist2arcs

step1 <- mmhc(gaussian.test, maximize.args = list(max.iter = 1))
step2a <- rsmax2(gaussian.test, restrict = "mmpc", restrict.args = list(test ="bic-gt"), maximize = "hc", maximize.args = list(max.iter = 1))
step2 <- hc(gaussian.test, start = step1, max.iter = 1)
step3a <-  mmhc(gaussian.test, maximize.args = list(max.iter = 3))
step3 <- hc(gaussian.test, start = step2, max.iter = 1)
step1
step2a
step2
step3a
step3

mmhc

if (bnlearn:::data.type(data) == "factor") {test = "bic-t"}
if (bnlearn:::data.type(data) == "continuous") {test = "bic-gt"}

net = pc.stable(gaussian.test, test = test)
net$learning$from = attr(data, "network")
net$learning$label = attr(data, "label")
net$learning$np = attr(data, "np")
pc.stable
bnlearn



data <- gaussian.test
iterations <- 5
netmmhc <- iteration_mmhc(data, 3)

save(netmmhc, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/netmmhc.rda")


data <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))

try <- mmpc(data, test = "bic-gt")
a <- make_full_graph(ncol(gaussian.test))
full <- as.bn(igraph.to.graphNEL(a))
nodes(full) <- names(gaussian.test)

tryigraph1 <- igraph.from.graphNEL(as.graphNEL(full))
tryigraph2 <- igraph.from.graphNEL(as.graphNEL(try))

dif <- difference(tryigraph1, tryigraph2)
difbn <- as.bn(igraph.to.graphNEL(dif))

check <- hc(gaussian.test, blacklist = arcs(difbn), max.iter = 3)
checkb <- mmhc(gaussian.test, maximize.args = list(max.iter = 3))

all.equal(check,checkb)

acyclic(net)

##################################################################################################
# Concept 
##################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_2.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_1.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_1.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g10.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g9.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g8.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g7.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g6.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g25.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g0.5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2_eBIC_g0.25.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_1.rda")


load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable2_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstable_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstablem_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pcstablem_4.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple5.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_align2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_align3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_align4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_2.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gsm_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gsm_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gsm_1.rda")

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gslow_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gslow_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gslow_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gslow_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gslow_1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gs_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gs_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gs_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gs_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs/gs_1.rda")
######################################################################################################
# Add g0 for mmhc. 
######################################################################################################
mmhcg0s <- list(mmhc_1,mmhc_2,mmhc_3,mmhc_4,mmhc_5)
i <- 1
for (i in 1:5){
networkdata <- mmhcg0s[[i]]$networkdata[nrow(mmhcg0s[[i]]$networkdata),]
networkdata$tests <- sum(mmhcg0s[[i]]$networkdata$tests)+ mmhcg0s[[i]]$begin$learning$ntests
networkdata$gamma <- 0
networks <- list(mmhcg0s[[i]]$networks[[2]])
names(networks) <- paste0("mmhc_",i,"_eBIC_g0")
assign(paste0("mmhc_",i,"_eBIC_g0"), list(networks = networks, networkdata = networkdata, begin = mmhcg0s[[i]]$begin))
save(list = paste0("mmhc_",i,"_eBIC_g0"), file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_",i,"_eBIC_g0.rda"))
}
######################################################################################################
# hcg0s <- list(hc_1,hc_2,hc_3,hc_4,hc_5)
# hc_4$networkdata[7801:8000,]
# for (i in 1:5){
#   networkdata <- hcg0s[[i]]$networkdata[nrow(hcg0s[[i]]$networkdata),]
#   networkdata$tests <- sum(hcg0s[[i]]$networkdata$tests)
#   assign(paste0("hc_",i,"_eBIC_g0"), list(networks = hcg0s[[i]]$networks[[2]], networkdata = networkdata))
#   save(list = paste0("hc_",i,"_eBIC_g0"), file = paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_",i,"_eBIC_g0.rda"))
# }

save(mmhc_1_eBIC_g0,mmhc_2_eBIC_g0,mmhc_3_eBIC_g0,mmhc_4_eBIC_g0,mmhc_5_eBIC_g0, file = )
mmhc_2_eBIC_g0
mmhc_1$
lmmhc1 <- length(mmhc_1$networkdata$logliks)
mmhc_1$networkdata[(lmmhc1-10):lmmhc1,]
log10(mmhc_1$begin$learning$ntests)
tabu_1$networkdata
sum(mmhc_1$networkdata$tests)+ mmhc_1$begin$learning$ntests
log10(sum(tabu_1$networkdata))

sum(c(1,1,1))
# pcstable2_1$networkdata[nrow(pcstable2_1$networkdata),]
# tabu_2_eBIC_g5$networkdata[nrow(tabu_2_eBIC_g5$networkdata),]
# tabu_2_eBIC_g5$networks$tabu_10d_494i
# pcstable2_2
# pcstable2_3
# pcstable2_4
# tabu_2_eBIC_g0.5$networks$tabu_10d_1554i
# tabu_2_eBIC_g0.5$networkdata$logliks[700:750] 
# tabu_2$networkdata$logliks[700:750]
# networkdat
# 
# 
# gs_2$networkdata
# gsm_1$networkdata
# tabu_4$networkdata[888,]
# fit4 <- bn.fit(tabu_mmhc_simple4,as.data.frame(datapermutations[[4]]))
# c(narcs(fit4),logLik(fit4,as.data.frame(datapermutations[[4]])),nparams(fit4))
# 
# tabu_2$networkdata[600:700,]
# tabu_mmhc_2$networkdata[190:205,]
# narcs(tabu_2$networks$tabu_10d_505i)
# narcs(tabu_mmhc_align2$networks$tabu_10d_302i)
# tabu_mmhc_align2$networkdata
# narcs(mmhc_2$networks[[3]])
# tabu_mmhc_align3$networks$tabu_10d_290i
# tabu_mmhc_align4$networks$tabu_10d_383i
##############################################################################################
# get graph
##############################################################################################
attributes(tabu_1$networks$tabu_10d_8000i)
360/(648*648)
360/1301
nodes
bn.net(pc_2$networks$pc_10d_0.2)
pcstable_1$networkdata
pcstable_2
pcstable_3
pcstable_4
pcstable_5
tabu_2$networkdata[7800:8000,]
360/8500

cextend
par1 <- nparams(bn.fit(hc_edges_loglik_10d_200_400i$networks, as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))))
arcs1 <- narcs(bn.fit(hc_edges_loglik_10d_200_400i$networks, as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))))
(tabu_4$networks$tabu_10d_6112i)
mmhc_1$networkdata
par2 <- nparams(mmhc_1$networks[[1]],as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)))
arcs2 <- narcs(mmhc_1$networks[[1]])
(tabu_4$networks$tabu_10d_6112i)
nparams(mmhc_2$networks[[2]],as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)))
mmhc_2$networkdata[nrow(mmhc_2$networkdata),]
narcs(mmhc_2$networks[[2]])
narcs(mmhc_4$networks[[2]])
narcs(mmhc_5$networks[[2]])
648*2
360/(648*2)
648*2+898
narcs(tabu_1$networks$tabu_10d_553i)
narcs(tabu_2$networks$tabu_10d_505i)
tabu_1$networkdata[553,]
tabu_2$networkdata[505,]
nparams(tabu_1$networks$tabu_10d_553i,as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, TRUE))) -553
nparams(tabu_2$networks$tabu_10d_505i,as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, TRUE)))
##############################################################################################
# Plot parametersratio vs speed 
# Old method: with iterations in HC en Tabu en MMHC. Without gammas.
#############################################################################################
testssum <- function(x){
  f1 <- function(x, n){ sum(head(x,n))}
  sums <- c()
  for (i in 1:length(x)){
    sums[i] <- f1(x,i)
  }
  return(sums)
}

# nets <- mmhc_2
# attributes(nets$networks$pc_10d_0.2)
# mmhc_2$networkdata$testsplus <- mmhc_2$networkdata$tests +mmhc_2$begin$learning$ntests
# check <- function(x){networkdata <- x$networkdata; return(networkdata)}
# check(pc_2)

ejesoud <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- 360/networkdata$params
  miss <- !is.na(xeje)
  if (length(nets) == 3){
    sumtests <- testssum(networkdata$tests)
    sumtests <- sumtests + nets$begin$learning$ntests
    yeje <- log10(sumtests)
  } else if (length(networkdata) == 5){
    sumtests <- testssum(networkdata$tests)
    yeje <- log10(sumtests)
  } else {yeje <- log10(networkdata$tests)}
  return(list(xeje,yeje))
}

nets <- tabu_1
tabu_1$networkdata$logliks[7800:8000]
mmhc_2_eBIC_g0.1$networks$mmhc_10d_g0.1

ejesedges <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- networkdata$edges
  miss <- !is.na(xeje)
  if (length(nets) == 3){
    sumtests <- testssum(networkdata$tests)
    sumtests <- sumtests + nets$begin$learning$ntests
    yeje <- log10(sumtests)
  } else if (length(networkdata) == 5){
    sumtests <- testssum(networkdata$tests)
    yeje <- log10(sumtests)
  } else {yeje <- log10(networkdata$tests)}
  return(list(xeje,yeje))
}


ejes <- function(nets){
networkdata <- nets$networkdata
network <- nets$networks
xeje <- 360/networkdata$params
miss <- !is.na(xeje)
if (length(nets) == 3){
  sumtests <- testsum(networkdata$testsplus)
  sumtests <- sumtests + nets$begin$learning$ntests
  yeje <- log10(sumtests)
} else if (length(networkdata) == 5){
  sumtests <- testsum(networkdata$tests)
  yeje <- log10(sumtests)
} else {yeje <- log10(networkdata$tests)}
return(list(xeje,yeje))
}
#####################################################################################
# Iterations plots: WATCH OUT: ORGANIZE ejesfunction and plot region. And check nparams dataset.
#######################################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/testsandloglik.pdf")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/testsvsedgesTabuMmhc.pdf")
pdf(plotname, width = 10, height = 5) 
par(mfrow = c(1,2))
########## choose edges funtion.
ejesmmhc <- lapply(list(mmhc_1, mmhc_2, mmhc_3, mmhc_4, mmhc_5),ejesedges)
ejestabu <- lapply(list(tabu_1,tabu_2, tabu_3, tabu_4, tabu_5),ejesedges)
ejestabuBIC <- lapply(list(tabu_2_eBIC_g0.25,tabu_2_eBIC_g0.5,tabu_2_eBIC_g5,tabu_2_eBIC_g10,tabu_2_eBIC_g25),ejes)
ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4),ejes)
# ejespc <- lapply(list(pc_2, pc_3),ejes)
ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejes)
# ejeslist <- lapply(list(mmhc_4,tabu_4, pcstablem_4), ejes)
# ejeslist <- lapply(list(tabu_2,tabu_3,tabu_4),ejes)

plot(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], type = "l",
     lwd = 1,
     col = rgb(0,255,0, alpha = 255, maxColorValue = 255),
     main = "tests vs number of params",
     xlab = expression(paste("n/",theta)), 
     ylab = "log10(calls to the statistical criterion)",
     xlim = c(0,0.3),
     ylim = c(2,10))


col2rgb("green")
plot(ejesmmhc[[1]][[1]], ejesmmhc[[1]][[2]], type = "l",
     lwd = 1,
     col = rgb(0,255,0, alpha = 255, maxColorValue = 255),
     main = "iteration process: tests vs number of edges, gamma = 0",
     xlab = "narcs", 
     ylab = "log10(calls to the statistical criterion)",
     xlim = c(0,7500),
     ylim = c(2,10))
points(ejesmmhc[[2]][[1]], ejesmmhc[[2]][[2]], lwd = 1,col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
points(ejesmmhc[[3]][[1]], ejesmmhc[[3]][[2]], col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
points(ejesmmhc[[4]][[1]], ejesmmhc[[4]][[2]], col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
points(ejesmmhc[[5]][[1]], ejesmmhc[[5]][[2]], col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
col2rgb("red")

lines(ejestabu[[1]][[1]], ejestabu[[1]][[2]], col = rgb(255,0,0, alpha = 255, maxColorValue = 255))
points(ejestabu[[2]][[1]], ejestabu[[2]][[2]], col = rgb(255,0,0, alpha = 20, maxColorValue = 255))
points(ejestabu[[3]][[1]], ejestabu[[3]][[2]], col = rgb(255,0,0, alpha = 20, maxColorValue = 255))
points(ejestabu[[4]][[1]], ejestabu[[4]][[2]], col = rgb(255,0,0, alpha = 20, maxColorValue = 255))
points(ejestabu[[5]][[1]], ejestabu[[5]][[2]], col = rgb(255,0,0, alpha = 20, maxColorValue = 255))
dev.off()
col2rgb("navyblue")
col2rgb("black")
lines(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,0, alpha = 255, maxColorValue = 255))
points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,0, alpha = 125, maxColorValue = 255))
points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,0, alpha = 125, maxColorValue = 255))
points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,0, alpha = 125, maxColorValue = 255))
# col2rgb("skyblue")
# lines(ejespc[[1]][[1]], ejespc[[1]][[2]], col = rgb(135,206,235, alpha = 255, maxColorValue = 255))
# points(ejespc[[2]][[1]], ejespc[[2]][[2]], col = rgb(135,206,235, alpha = 1, maxColorValue = 255))
col2rgb("blue")
lines(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 100, maxColorValue = 255))
points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
col2rgb("darkgrey")
lines(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], col = rgb(169,169,169, alpha = 255, maxColorValue = 255))
lines(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(0,0,0, alpha = 255, maxColorValue = 255))
legend(x="left", c("mmhc",expression(paste("tabu/hc ",gamma,"= 0,0.5,5")),expression(paste("pcstable ",gamma,"\u22655")),expression(paste("gs ",gamma,"\u22655"))), fill = c("green","red","navyblue","blue"),  col = c("green","red","navyblue","skyblue","blue"))
pcstable2_2$networkdata
##############################################################################################
# Plot parametersratio vs loglik
#############################################################################################
ejeslog <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- 360/networkdata$params
  miss <- !is.na(xeje)
  if (length(networkdata) == 5){
    logliks <- networkdata$logliks
    yeje <- logliks
  } else if (length(nets) == 3){
    logliks <- networkdata$logliks
    yeje <- logliks
  } else {yeje <- networkdata$logliks}
  return(list(xeje,yeje))
}

ejesmmhc <- lapply(list(mmhc_1, mmhc_2, mmhc_3, mmhc_4, mmhc_5),ejeslog)
ejestabu <- lapply(list(tabu_2, tabu_3, tabu_4, tabu_5),ejeslog)
ejestabuBIC <- lapply(list(tabu_2_eBIC_g0.25,tabu_2_eBIC_g0.5,tabu_2_eBIC_g5,tabu_2_eBIC_g6,tabu_2_eBIC_g7,tabu_2_eBIC_g8,tabu_2_eBIC_g9,tabu_2_eBIC_g10,tabu_2_eBIC_g25),ejeslog)
ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4),ejeslog)
# ejespc <- lapply(list(pc_2, pc_3),ejeslog)
ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejeslog)

col2rgb("green")



plot(ejesmmhc[[1]][[1]], ejesmmhc[[1]][[2]], type = "l",
     lwd = 1,
     col = rgb(0,255,0, alpha = 255, maxColorValue = 255),
     main = "logliklihood vs sample-parameter ratio",
     xlab = expression(paste("n/",theta)), 
     ylab = expression(paste("P(G|(x,",theta,")")),
     xlim = c(0,0.3),
     ylim = c(-380000,-100000)
     )
points(ejesmmhc[[2]][[1]], ejesmmhc[[2]][[2]], lwd = 1,col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
points(ejesmmhc[[3]][[1]], ejesmmhc[[3]][[2]], col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
points(ejesmmhc[[4]][[1]], ejesmmhc[[4]][[2]], col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
points(ejesmmhc[[5]][[1]], ejesmmhc[[5]][[2]], col = rgb(0,255,0, alpha = 1, maxColorValue = 255))
col2rgb("red")
lines(ejestabu[[1]][[1]], ejestabu[[1]][[2]], col = rgb(255,0,0, alpha = 255, maxColorValue = 255))
points(ejestabu[[2]][[1]], ejestabu[[2]][[2]], col = rgb(255,0,0, alpha = 1, maxColorValue = 255))
points(ejestabu[[3]][[1]], ejestabu[[3]][[2]], col = rgb(255,0,0, alpha = 1, maxColorValue = 255))
points(ejestabu[[4]][[1]], ejestabu[[4]][[2]], col = rgb(255,0,0, alpha = 1, maxColorValue = 255))
col2rgb("navyblue")
col2rgb("black")
lines(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,0, alpha = 255, maxColorValue = 255))
points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,0, alpha = 125, maxColorValue = 255))
points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,0, alpha = 1, maxColorValue = 255))
points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,0, alpha = 1, maxColorValue = 255))
# col2rgb("skyblue")
# lines(ejespc[[1]][[1]], ejespc[[1]][[2]], col = rgb(135,206,235, alpha = 255, maxColorValue = 255))
# points(ejespc[[2]][[1]], ejespc[[2]][[2]], col = rgb(135,206,235, alpha = 1, maxColorValue = 255))
col2rgb("blue")
lines(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 100, maxColorValue = 255))
points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))

lines(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], col = rgb(169,169,169, alpha = 255, maxColorValue = 255))
lines(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(0,0,0, alpha = 255, maxColorValue = 255))
legend(x="left", c("mmhc","tabu/hc","pcstable","gs"), fill = c("green","red","navyblue","blue"),  col = c("green","red","navyblue","skyblue","blue"))


plot(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], type = "l",
     lwd = 1,
     col = rgb(0,255,0, alpha = 255, maxColorValue = 255),
     main = "logliklihood vs sample-parameter ratio",
     xlab = expression(paste("n/",theta)), 
     ylab = expression(paste("P(G|(x,",theta,")")),
     xlim = c(0.202,0.203),
     ylim = c(-228000,-226000)
)
lines(ejestabu[[1]][[1]], ejestabu[[1]][[2]], col = rgb(255,0,0, alpha = 255, maxColorValue = 255))
lines(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
lines(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 255, maxColorValue = 255))
lines(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 255, maxColorValue = 255))

all.equal(tabu_2_eBIC_g0.5$networks$tabu_10d_901i, tabu_2$networks$tabu_10d_2317i)
all.equal(tabu_2_eBIC_g0.5$networkdata[600:900,],tabu_2$networkdata[600:900,])

ejestabuBIC[[4]][[1]]
###########################################################################################
# Permutations (done in Terminalscript)
###########################################################################################
permutations <- list()
for (j in 1:length(datapermutations)){
  indsame <- c()
  for (i in 1:ncol(datapermutations[[1]])){
    int <- which(colnames(datapermutations[[1]]) == colnames(datapermutations[[j]][i]))
    indsame[i] <- int
  }
permutations[[j]] <- indsame
}

permutations
save(permutations, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
##############################################################################################
# Compare longdistances different 
###############################################################################################


nodes(mmhc_3$networks[[3]])
all.equal(datapermutations[[1]],TimeCoordsAnom_from_Grid_std(tas_ncep_10d))
plot_long_distances(mmhc_3$networks[[1]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[3]])
plot_long_distances(mmhc_3$networks[[3]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[3]])
plot_long_distances(mmhc_2$networks[[1]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[2]])
plot_long_distances(hc_3$networks[[3]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist =5000, perm = permutations[[3]])
plot_long_distances(tabu_1$networks[[4]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[1]])
plot_long_distances(tabu_2_eBIC_g1$networks$tabu_10d_1123i, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[2]])
plot_long_distances(tabu_2_eBIC_g5$networks$tabu_10d_1123i, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[2]])
plot_long_distances(tabu_2_eBIC_g5$networks$tabu_10d_493i, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[2]])
plot_long_distances(tabu_2_eBIC_g0.5$networks$tabu_10d_1554i, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[2]])
plot_long_distances(tabu_mmhc_simple2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 1000, perm = permutations[[2]])

plot_long_distances(hc_edges_loglik_10d_600_800i$networks, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[1]])
plot_long_distances(pc_10d_1e_13$networks$dag, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[1]])

plot.Meteodag(time.coords = TimeCoordsAnom_from_Grid_std(tas_ncep_10d), hc_edges_loglik_10d_600_800i$networks)

mmhc_3$networks[[3]]
all.equal(mmhc_3$networks[[2]],mmhc_3$networks[[3]])
attributes(mmhc_3$networks[[1]])
View(mmhc_3$networks)

##################################################################################################
# Compare propagation mmhc. 
##################################################################################################
# load mmhc propagation
pattern <- "mmhc.*V205"
filesmmhc <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", full.names = T, pattern = pattern)
filesmmhcnames <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", pattern = pattern)
filesmmhcnames <- gsub(".rda", "", filesmmhcnames)
filesmmhcnames

variablelistmmhc <- list()
for (i in 1:length(filesmmhc)){
  variablepos <- get(load(filesmmhc[i]))
  variablelistmmhc[[i]] <- variablepos
}
names(variablelistmmhc) <- filesmmhcnames
str(variablelistmmhc)
# load Tabu propagation
pattern <- "tabu.*V205"
filestabu <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", full.names = T, pattern = pattern)
filestabunames <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", pattern = pattern)
filestabunames <- gsub(".rda", "", filestabunames)

variablelisttabu <- list()
for (i in 1:length(filestabu)){
  variablepos <- get(load(filestabu[i]))
  variablelisttabu[[i]] <- variablepos
}
names(variablelisttabu) <- filestabunames
str(variablelisttabu)

mmhcclims <- list()
tabuclims <- list()
for(i in 1:length(variablelisttabu)){
mmhcclims[[i]] <- quantity2clim(variablelistmmhc[[i]]$with, attr(variablelistmmhc[[i]]$with, "probability"), tas_ncep_10d)
tabuclims[[i]] <- quantity2clim(variablelisttabu[[i]]$with, attr(variablelisttabu[[i]]$with, "probability"), tas_ncep_10d)
}

plotsclim <- list()
plotstabu <- list()
for( i in 1:length(variablelisttabu)){
  plotsclim[[i]] <- plotClimatology(mmhcclims[[i]],backdrop.theme = "coastline", main = list(attr(mmhcclims[[i]]$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  plotstabu[[i]] <- plotClimatology(tabuclims[[i]],backdrop.theme = "coastline", main = list(attr(tabuclims[[i]]$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
}

totallist <- c(plotsclim,plotstabu)
lay = rbind(c(1,6),
            c(2,7),
            c(3,8),
            c(4,9),
            c(5,10))
grid.arrange(grobs = totallist, layout_matrix = lay)
##################################################################################################
pattern <- "tabu_1600.*V81"
filestabu <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", full.names = T, pattern = pattern)
filestabunames <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", pattern = pattern)
filestabunames <- gsub(".rda", "", filestabunames)

variablelisttabu <- list()
for (i in 1:length(filestabu)){
  variablepos <- get(load(filestabu[i]))
  variablelisttabu[[i]] <- variablepos
}
names(variablelisttabu) <- filestabunames
str(variablelisttabu)

# mmhcclims <- list()
tabuclims <- list()
for(i in 1:length(variablelisttabu)){
  # mmhcclims[[i]] <- quantity2clim(variablelistmmhc[[i]]$with, attr(variablelistmmhc[[i]]$with, "probability"), tas_ncep_10d)
  tabuclims[[i]] <- quantity2clim(variablelisttabu[[i]]$with, attr(variablelisttabu[[i]]$with, "probability"), tas_ncep_10d)
}

# plotsclim <- list()
plotstabu <- list()
for( i in 1:length(variablelisttabu)){
  # plotsclim[[i]] <- plotClimatology(mmhcclims[[i]],backdrop.theme = "coastline", main = list(attr(mmhcclims[[i]]$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  plotstabu[[i]] <- plotClimatology(tabuclims[[i]],backdrop.theme = "coastline", main = list(attr(tabuclims[[i]]$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
}

lay = rbind(c(1),
            c(2),
            c(3),
            c(4),
            c(5))
grid.arrange(grobs = plotstabu, layout_matrix = lay)

plotClimatology(test1quantity,backdrop.theme = "coastline", main = list(attr(test1quantity$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
test2quantity <- quantity2clim(prop_tabu1_844_V81_equal2$with, attr(prop_tabu1_844_V81_equal2$with, "probability"), tas_ncep_10d)
plotClimatology(test2quantity,backdrop.theme = "coastline", main = list(attr(test2quantity$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
bn_mmhc_3_1 <- bn.fit(mmhc_3$networks[[1]],as.data.frame(initialdata[permutations[[3]]]))
bn_mmhc_3_3 <- bn.fit(mmhc_3$networks[[3]],as.data.frame(initialdata[permutations[[3]]]))
test1 <- PropagationExactGeneralPerm(baysnet = bn_mmhc_3_1,
                            nodesEvents = 1:648,
                            valueEvent = ">=1",
                            nodesEvidence = c(81,280),
                            valueEvidence = c(2,2),
                            perm = permutations[[3]])
test3 <- PropagationExactGeneralPerm(baysnet = bn_mmhc_3_3,
                                     nodesEvents = 1:648,
                                     valueEvent = ">=1",
                                     nodesEvidence = c(81,280),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[3]])

test1quantity <- quantity2clim(test1$with, attr(test1$with, "probability"), tas_ncep_10d)
plotClimatology(test1quantity,backdrop.theme = "coastline", main = list(attr(test1quantity$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
test3quantity <- quantity2clim(test3$with, attr(test3$with, "probability"), tas_ncep_10d)
plotClimatology(test3quantity,backdrop.theme = "coastline", main = list(attr(test3quantity$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))


bn_hc_1_1 <- bn.fit(hc_edges_loglik_10d_600_800i$networks, initialdata)
testhc1 <- PropagationExactGeneralPerm(baysnet = bn_hc_1_1,
                                       nodesEvents = 1:648,
                                       valueEvent = ">=1",
                                       nodesEvidence = c(81,280),
                                       valueEvidence = c(2,2),
                                       perm = permutations[[1]])





##################################################################################################
# Orden and compare Plots LONG DISTANCES
##################################################################################################
# mmhc
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_1.rda")
# Tabu simples
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_mmhc_simple5.rda")
# Tabus large
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizegraphs/tabu_simple_1600_5.rda")
# gs size graphs
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/tabu_gs_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/tabu_gs_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/tabu_gs_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/tabu_gs_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/tabu_gs_simple5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/mmhc_gs_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/mmhc_gs_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/mmhc_gs_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/mmhc_gs_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/mmhc_gs_simple5.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/gs_simple1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/gs_simple2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/gs_simple3.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/gs_simple4.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gssizegraphs/gs_simple5.rda")
#permutations
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")

mmhcs <- list(mmhc_1$networks[[3]], 
              mmhc_2$networks[[3]], 
              mmhc_3$networks[[3]], 
              mmhc_4$networks[[3]], 
              mmhc_5$networks[[3]])

tabus <- list(tabu_mmhc_simple1,
              tabu_mmhc_simple2, 
              tabu_mmhc_simple3,
              tabu_mmhc_simple4,
              tabu_mmhc_simple5)

tabuslarge <- list(tabu_simple_1600_1,
                  tabu_simple_1600_2,
                  tabu_simple_1600_3,
                  tabu_simple_1600_4,
                  tabu_simple_1600_5)



load(file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
library(gridExtra)
library(gridGraphics)
load(file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/utilnetworks_hcnetworks10d.rda")
load(file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Data_graphsPC.rda")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/comparegraphsdistances.pdf")
pdf(plotname, width = 14, height = 8)
op <- par(no.readonly = TRUE)
lay <- rbind(c(1,0,0),
             c(2,3,0),
             c(4,5,6))
graphics::layout(as.matrix(lay))
par(mar = c(0.0, 0.0, 0.0, 0.0))
par(mgp = c(0, 0, 0))
par(oma = c(0, 0, 0, 0))


# best of pcstable
plot_long_distances(cextend(gs_simple2), data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(0,0,255,alpha = 100, maxColorValue = 255), perm = permutations[[2]],title = TRUE)

# pcstablbest and best of mmhc
plot_long_distances(mmhc_gs_simple2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(0,255,0,alpha = 100, maxColorValue = 255) ,perm = permutations[[2]], title = TRUE)
plot_long_distances(mmhc_2$networks[[2]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(0,255,0,alpha = 100, maxColorValue = 255), perm = permutations[[2]],title = TRUE)

# pcstablebest, mmhcbest, hcbien. 
plot_long_distances(tabu_gs_simple2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(200,100,100, alpha = 100, maxColorValue = 255), perm = permutations[[2]],title = TRUE, remove = TRUE)
plot_long_distances(tabu_mmhc_simple2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(200,100,100, alpha = 50,maxColorValue = 255), perm = permutations[[2]], title = TRUE, remove = TRUE)
plot_long_distances(tabu_simple_1600_2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = NA ,perm = permutations[[2]], title = FALSE, remove = TRUE)
par(op)
dev.off()

#################################################################### Compare
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/longmmhctabu.pdf")
pdf(plotname, width = , height = 13)
lay <- rbind(c(1,2),
             c(3,4),
             c(5,6),
             c(7,8),
             c(9,10))
graphics::layout(as.matrix(lay))


for(i in 1:length(mmhcs)){
#plot_long_distances(mmhcs[[i]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[i]])
plot_long_distances(tabus[[i]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[i]])
plot_long_distances(tabuslarge[[i]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, perm = permutations[[i]])
}
dev.off()
##################################################################################################
# Orden and compare Plots PROPAGATION
##################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_1_tabu_10d_g0.1_V81_equal2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_1_tabu_10d_g0.5_V81_equal2.rda")

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/comparegraphs.pdf")
png(plotname,width = 400,height = 500)
lay = rbind(c(1,NA,NA),
            c(2,3,NA),
            c(4,5,6))
grid.arrange(grobs = plots, layout_matrix = lay)
dev.off()
setwd

iteration_mmhc_split()
iteration_tabu(as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)), 18)



# load mmhc propagation
pattern <- "mmhc5.*V81"
filesmmhc <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", full.names = T, pattern = pattern)
filesmmhcnames <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", pattern = pattern)
filesmmhcnames <- gsub(".rda", "", filesmmhcnames)
filesmmhcnames

variablelistmmhc <- list()
for (i in 1:length(filesmmhc)){
  variablepos <- get(load(filesmmhc[i]))
  variablelistmmhc[[i]] <- variablepos
}
names(variablelistmmhc) <- filesmmhcnames
str(variablelistmmhc)
# load Tabu propagation
pattern <- "tabu5.*V81|tabu_1600_5"
filestabu <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", full.names = T, pattern = pattern)
filestabunames <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", pattern = pattern)
filestabunames <- gsub(".rda", "", filestabunames)

variablelisttabu <- list()
for (i in 1:length(filestabu)){
  variablepos <- get(load(filestabu[i]))
  variablelisttabu[[i]] <- variablepos
}
names(variablelisttabu) <- filestabunames
str(variablelisttabu)

tabu_1_eBIC_g0.5
variablelisttabu <- list(prop_2_tabu_10d_g0.25_V81_equal2)
mmhcclims <- list()
tabuclims <- list()
for(i in 1:length(variablelisttabu)){
 tabuclims[[i]] <- quantity2clim(variablelisttabu[[i]]$with-variablelisttabu[[i]]$without, paste0(attr(variablelisttabu[[i]]$with, "probability"),"-",attr(variablelisttabu[[i]]$without, "probability")), tas_ncep_10d)
}
for(i in 1:length(variablelistmmhc)){
  mmhcclims[[i]] <- quantity2clim(variablelistmmhc[[i]]$with-variablelistmmhc[[i]]$without,paste0(attr(variablelisttabu[[i]]$with, "probability"),"-",attr(variablelisttabu[[i]]$without, "probability")), tas_ncep_10d)
}
plotsclim <- list()
plotstabu <- list()
for( i in 1:length(variablelisttabu)){
  plotstabu[[i]] <- plotClimatology(tabuclims[[i]],backdrop.theme = "coastline", main = list(attr(tabuclims[[i]]$Data,"climatology:fun"),cex = 0.7),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0.05,1,0.05))
}
for( i in 1:length(variablelistmmhc)){
plotsclim[[i]] <- plotClimatology(mmhcclims[[i]],backdrop.theme = "coastline", main = list(attr(mmhcclims[[i]]$Data,"climatology:fun"),cex = 0.7),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0.05,1,0.05))
}
totallist <- c(plotsclim,plotstabu)

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/comparepropgraphs.pdf")
pdf(plotname, width = 11, height = 6)
lay = rbind(c(1,NA),
            c(3,2))
grid::gpar(mar = c(0.0, 0.0, 0.0, 0.0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0))
grid.arrange(grobs = totallist, layout_matrix = lay, gp =  gpar(mar = c(0.0, 0.0, 0.0, 0.0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0)))

dev.off()

##################################################################################################
# Orden and compare Plots PROPAGATION Extendent version with H2pc. 
##################################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_1_tabu_10d_g0.1_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_1_tabu_10d_g0.25_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_1_tabu_10d_g0.5_V81_equal2.rda")
cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/comparegraphs.pdf")
png(plotname,width = 400,height = 500)
lay = rbind(c(1,NA,NA),
            c(2,3,NA),
            c(4,5,6))
grid.arrange(grobs = plots, layout_matrix = lay)
dev.off()
setwd

# load mmhc propagation
pattern <- "mmhc5.*V81"
filesmmhc <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", full.names = T, pattern = pattern)
filesmmhcnames <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhcsizeprop", pattern = pattern)
filesmmhcnames <- gsub(".rda", "", filesmmhcnames)
filesmmhcnames

variablelistmmhc <- list()
for (i in 1:length(filesmmhc)){
  variablepos <- get(load(filesmmhc[i]))
  variablelistmmhc[[i]] <- variablepos
}
names(variablelistmmhc) <- filesmmhcnames
str(variablelistmmhc)

# load h2pc propagation
pattern <- "h2pc_g0.1_2"
filesh2pc <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC", full.names = T, pattern = pattern)
filesh2pcnames <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC", pattern = pattern)
filesh2pccnames <- gsub(".rda", "", filesh2pcnames)
filesh2pccnames

variablelisth2pc <- list()
for (i in 1:length(filesh2pc)){
  variablepos <- get(load(filesh2pc[i]))
  variablelisth2pc[[i]] <- variablepos
}
names(variablelisth2pc) <- filesh2pcnames
str(variablelisth2pc)

# load Tabu propagation



pattern <- "tabu_10d_g0.25_V81_"
filestabu <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/", full.names = T, pattern = pattern)
filestabunames <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/", pattern = pattern)
filestabunames <- gsub(".rda", "", filestabunames)

variablelisttabu <- list()
for (i in 1:length(filestabu)){
  variablepos <- get(load(filestabu[i]))
  variablelisttabu[[i]] <- variablepos
}
names(variablelisttabu) <- filestabunames
str(variablelisttabu)

h2pcclims <- list()
tabuclims <- list()
mmhcclims <- list()
for(i in 1:length(variablelisttabu)){
  tabuclims[[i]] <- quantity2clim(variablelisttabu[[i]]$with-variablelisttabu[[i]]$without, paste0(attr(variablelisttabu[[i]]$with, "probability"),"-",attr(variablelisttabu[[i]]$without, "probability")), tas_ncep_10d)
}
for(i in 1:length(variablelisth2pc)){
  h2pcclims[[i]] <- quantity2clim(variablelisth2pc[[i]]$with-variablelisth2pc[[i]]$without, paste0(attr(variablelisth2pc[[i]]$with, "probability"),"-",attr(variablelisth2pc[[i]]$without, "probability")), tas_ncep_10d)
}
plotsh2pc <- list()
plotstabu <- list()
for( i in 1:length(variablelisttabu)){
  plotstabu[[i]] <- spatialPlot(tabuclims[[i]],backdrop.theme = "coastline", main = list(paste0("Tabu: |E| = ",narcs(tabu_1_eBIC_g0.25$networks$tabu_10d_g0.25),"\n",attr(tabuclims[[i]]$Data,"climatology:fun")),cex = 0.7), colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0.05,1,0.05),lonCenter = 180, region = TRUE, col.regions= cb)
}

for( i in 1:length(variablelistmmhc)){
  plotsh2pc[[i]] <- spatialPlot(h2pcclims[[i]],backdrop.theme = "coastline", main = list(paste0("H2PC: |E| = ",narcs(h2pc_2_eBIC_g0.1$networks$h2pc_10d_g0.1),"\n",attr(h2pcclims[[i]]$Data,"climatology:fun")),cex = 0.7),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0.05,1,0.05), lonCenter = 180, region = TRUE, col.regions= cb)
}
totallist <- c(plotsh2pc,plotstabu)
totallist[[1]]

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/MarcoNewPlots/propH2pcTabu.pdf")
pdf(plotname)
lay = rbind(c(1,NA),
            c(2,NA))

grid::gpar(mar = c(0.0, 0.0, 0.0, 0.0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0))
grid.arrange(grobs = totallist[c(1,2)], layout_matrix = lay, gp =  gpar(mar = c(0.0, 0.0, 0.0, 0.0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0)))

dev.off()




#############################################################################################
# h2pc
#############################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC/prop_h2pc_g0.1_1_1130_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC/prop_h2pc_g0.1_2_1121_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC/prop_h2pc_g0.1_3_1110_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC/prop_h2pc_g0.1_4_1120_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/h2pc_eBIC/prop_h2pc_g0.1_5_1122_V81_equal2.rda")

variablelisth2pc <- list(prop_h2pc_g0.1_1_1130_V81_equal2,
                         prop_h2pc_g0.1_2_1121_V81_equal2,
                         prop_h2pc_g0.1_3_1110_V81_equal2,
                         prop_h2pc_g0.1_4_1120_V81_equal2,
                         prop_h2pc_g0.1_5_1122_V81_equal2)

h2pcclims <- list()
for(i in 1:length(variablelisth2pc)){
  h2pcclims[[i]] <- quantity2clim(variablelisth2pc[[i]]$with-variablelisth2pc[[i]]$without, paste0(attr(variablelisth2pc[[i]]$with, "probability"),"-",attr(variablelisth2pc[[i]]$without, "probability")), tas_ncep_10d)
}
plotsh2pc<- list()
for( i in 1:length(variablelisth2pc)){
  plotsh2pc[[i]] <- spatialPlot(h2pcclims[[i]],backdrop.theme = "coastline", main = list(attr(h2pcclims[[i]]$Data,"climatology:fun"),cex = 0.7),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0.05,1,0.05),lonCenter = 180, rev.colors = TRUE)
}
do.call(grid.arrange, c(plotsh2pc, nrow = 2))

        