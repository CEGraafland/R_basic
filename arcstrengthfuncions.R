library(rJava)
library(loadeR)
library(transformeR)
library(magrittr)
library(abind)
library(igraph)
library(devtools)
library(mopa)
data(wrld)
str(wrld)
library(sp)
library(graph)
library(rbmn)

library("pcalg")
library(bnlearn)
library(gRain)
library(Rgraphviz)


dag <- hc_edges_loglik_10d_200_400i$networks #gemaakt met tas_ncep_10d
dag <- pc_10d_1e_13$networks$dag
time.coords <- TimeCoordsAnom_from_Grid(tas_ncep_10d) # behoudt grid atrributes

logLik(dag,as.data.frame(time.coords))

#dag learned from bnlearn package
plot_quantile_edges <- function(dag,time.coords, cex = 1){

  dag.data <- as.data.frame(time.coords)
  arcsbn <- arcs(dag)
  linkstrength <- arc.strength(dag, dag.data)
  
  #create igraph from bn graph for plotting purpose
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))

  #identify which indices igraph correspond to indices in bn
  indsame <- c()
  
  for (i in 1:nrow(linkstrength)){
    int <- intersect(which(as_edgelist(igraph)[,1] == arcsbn[i,1]),
                     which(as_edgelist(igraph)[,2] == arcsbn[i,2]))
    indsame[i] <- int
  }
  indsame
  
  #create strengths vector in igraph object
  #low strength corresponds to high weight
  strengths <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    strengths[ind] <- linkstrength[i,3]
  }


  #stel de strengths en weights in:
  E(igraph)$strengths <- strengths
  max(strengths)
  min(strengths)
  normalize <- max(strengths)-min(strengths)
  E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  E(igraph)$weights
  
  #Stel colors in
  c_scale <- rev(topo.colors(5))

  #for HC manual
  #E(igraph)[weights <0.1 ]$color <- c_scale[1]
  #E(igraph)[(weights >= 0.1) & (weights < 0.2) ]$color <- c_scale[2]
  #E(igraph)[(weights >= 0.2) & (weights < 0.3) ]$color <- c_scale[3]
  #E(igraph)[(weights >= 0.3) & (weights < 0.4) ]$color <- c_scale[4]
  #E(igraph)[weights >=0.4 ]$color <- c_scale[5]
  #E(igraph)[weights >= 0.5]$color <- "red"

  #quantiles set colors to quantiles
  # E(igraph)[weights <= quantile(E(igraph)$weights)[2]]$color <- c_scale[2]
  # E(igraph)[weights > quantile(E(igraph)$weights)[2] & weights <=quantile(E(igraph)$weights)[3]]$color <- NA
  # E(igraph)[weights > quantile(E(igraph)$weights)[3] & weights <=quantile(E(igraph)$weights)[4]]$color <- NA
  # E(igraph)[weights > quantile(E(igraph)$weights)[4] & weights <=quantile(E(igraph)$weights)[5]]$color <- NA

  #quantiles set colors to quantiles
  E(igraph)[weights <= quantile(E(igraph)$weights)[2]]$color <- c_scale[2]
  E(igraph)[weights > quantile(E(igraph)$weights)[2] & weights <=quantile(E(igraph)$weights)[3]]$color <- c_scale[3]
  E(igraph)[weights > quantile(E(igraph)$weights)[3] & weights <=quantile(E(igraph)$weights)[4]]$color <- c_scale[4]
  E(igraph)[weights > quantile(E(igraph)$weights)[4] & weights <=quantile(E(igraph)$weights)[5]]$color <- c_scale[5]
  
  #maak het plaatje:
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)

  points <- expand.grid(y, x)[2:1]
  data(wrld)
  plot(wrld)
  plot.igraph(igraph, 
            vertex.size = 100,
            vertex.color = "blue",
            vertex.label = NA,
            edge.arrow.size = 0.2,
            edge.lty = 1+(E(igraph)$normweights)*2,
            layout = as.matrix(points), add = TRUE, rescale = FALSE)


  legend(-260,90,legend=c(paste0("<",round(quantile(E(igraph)$weights)[2],digits = 2)),
                         paste0("<",round(quantile(E(igraph)$weights)[3],digits = 2)),
                         paste0("<",round(quantile(E(igraph)$weights)[4],digits = 2)),
                         paste0("<",round(quantile(E(igraph)$weights)[5],digits = 2))),
        fill=c_scale[1:4],
        cex=cex)
}

dev.off()
#hcdag <- hc_edges_loglik_10d_200_400i$networks
#data <- TimeCoordsAnom_from_Grid(tas_ncep_10d)

# quantile plots HC
par(mfrow = c(2,2), cex =0.3)
plot_quantile_edges(hc_edges_loglik_10d_200i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_200_400i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_400_600i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_600_800i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_800_1000i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_1000_1200i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_1200_1400i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_1400_1600i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
plot_quantile_edges(hc_edges_loglik_10d_2800_3000i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
dev.off()
# quantile plots PC (NOT YET)
pc_10d_1e_13 <- alpha_loglik_pc3_withdag(graph_from_Grid(tas_ncep_10d),1e-13)
nodes(pc_10d_1e_13$networks$dag)<-nodes(hc_edges_loglik_10d_200_400i$networks)

plot_quantile_edges(pc_10d_1e_13$networks$dag, TimeCoordsAnom_from_Grid(tas_ncep_10d))



### ---- Function Stronweak egdes ------------------------------------------------------------
plot_strongweak_edges <- function(dag,time.coords,treshold,which = c("both","strong"),cex =1, criterion = NULL, k = NULL ){

  dag.data <- as.data.frame(time.coords)
  names(dag.data) <- nodes(dag)
  arcsbn <- arcs(dag)
  linkstrength <- arc.strength(dag, dag.data, criterion = criterion ,k  = k)
  
  #groter dan treshold:
  strongarcs <- as.matrix(linkstrength[linkstrength["strength"] <= treshold,,])[,1:2]
  strongarcs
  
  #create igraph from bn graph for plotting purpose
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))

  ##identify which indices igraph correspond to indices in bn
  #indsame <- c()

  #for (i in 1:nrow(linkstrength)){
  #  int <- intersect(which(as_edgelist(igraph)[,1] == arcsbn[i,1]),
  #                  which(as_edgelist(igraph)[,2] == arcsbn[i,2]))
  #  indsame[i] <- int
  #}
  #indsame
  
  ##create strengths vector in igraph object
  ##low strength corresponds to high weight
  #strengths <- numeric(nrow(as_edgelist(igraph)))
  #for (i in 1:nrow(as_edgelist(igraph))){
  #  ind <- indsame[i]
  #  as_edgelist(igraph)[ind,]
  #  strengths[ind] <- linkstrength[i,3]
  #}
  
  #identify indices sequence with 'strong arcs' in igraph600
  indstrong <- c()
  nrow(strongarcs)
  strongarcs
  for (i in 1:nrow(strongarcs)){
    int <- intersect(which(as_edgelist(igraph)[,1] == unname(strongarcs[i,1])),
                     which(as_edgelist(igraph)[,2] == unname(strongarcs[i,2])))
    indstrong[i] <- int
  }
  indstrong
  sort(indstrong)
  as_edgelist(igraph)[sort(indstrong),]
  
  #Zwakke links:
  no.indstrong <- as.vector(1:nrow(as_edgelist(igraph)))
  no.indstrong <- no.indstrong[-sort(indstrong)]
  no.indstrong
  
  #STel de kleuren in voor sterk en zwak:
  if(which == "both"){
    weak <- "blue"
  }
  if(which == "strong"){
    weak <- "NA"
  }
  E(igraph)[ sort(no.indstrong) ]$color <- weak
  E(igraph)[ sort(indstrong) ]$color <- "red"
  #maak het plaatje:
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  data(wrld)
  plot(wrld)
  plot.igraph(igraph, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              edge.lty = 1+(E(igraph)$normweights)*2,
              main = paste(deparse(substitute(dag))),
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  
  legend(-250,90,legend=c(paste0("<",treshold),
                          paste0(">",treshold)),
         fill=c("red",weak),
         x.intersp = 0.2,
         xjust = 0.5,
         cex = cex)
}






# TEST HC
plot_strongweak_edges(hc_edges_loglik_10d_200_400i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d),
                      treshold = -230.0000, which = "strong")
plot_strongweak_edges(hc_edges_loglik_10d_400_600i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d),
                      treshold = -230.0000, which = "strong")
plot_strongweak_edges(hc_edges_loglik_10d_1400_1600i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d),treshold = -230.0000,
                      which = "strong")
plot_strongweak_edges(hc_edges_loglik_10d_2800_3000i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d),treshold = -230.0000)

# TEST PC
data.dag <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
names(data.dag) <- nodes(pc_10d_1e_13$networks$dag)
nedges <- narcs(pc_10d_1e_13$networks$dag)
nedges
quantile1pc <-signif(quantile(arc.strength(pc_10d_1e_13$networks$dag, data.dag)[,3]), digits =3)
quantile1pc[2]
plot_strongweak_edges(pc_10d_1e_13$networks$dag, TimeCoordsAnom_from_Grid(tas_ncep_10d),treshold = quantile1pc[3], which = "strong")
title(main = paste0("dag from PC 25 procent \n number of edges = ",nedges))

# Test tabu
plot_strongweak_edges(tabu_1$networks$tabu_10d_1154i, TimeCoordsAnom_from_Grid(tas_ncep_10d),
                      treshold = -10, which = "strong", criterion = "bic-g", k = 1)
dag <- tabu_1$networks$tabu_10d_6848i
dag.data <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
names(dag.data) <- nodes(dag)
arcsbn <- arcs(dag)
linkstrength <- arc.strength(dag, dag.data, criterion = "bic-g", k = 5)
linkstrength




#####################################################################################################################
# make different plots with 50 procent quantile first plot.
#####################################################################################################################
load("sftp://ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
# set variable names in data
nodes(hc_edges_loglik_10d_200i$networks)
data.dag <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
names(data.dag) <- nodes(hc_edges_loglik_10d_200i$networks)
# calculate quantiles first 200 nodes
quantiles <-quantile(arc.strength(hc_edges_loglik_10d_3600_3800i$networks,data.dag)[,3])
quantiles
par(mfrow = c(3, 2))
data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
for (i in list(hc_edges_loglik_10d_200i,
               hc_edges_loglik_10d_200_400i,
               hc_edges_loglik_10d_400_600i,
               hc_edges_loglik_10d_600_800i,
               hc_edges_loglik_10d_800_1000i,
               hc_edges_loglik_10d_1000_1200i,
               hc_edges_loglik_10d_1200_1400i,
               hc_edges_loglik_10d_1400_1600i,
               hc_edges_loglik_10d_1600_1800i,
               hc_edges_loglik_10d_1800_2000i,
               hc_edges_loglik_10d_2000_2200i,
               hc_edges_loglik_10d_2200_2400i,
               hc_edges_loglik_10d_2400_2600i,
               hc_edges_loglik_10d_2600_2800i,
               hc_edges_loglik_10d_2800_3000i,
               hc_edges_loglik_10d_3000_3200i,
               hc_edges_loglik_10d_3200_3400i,
               hc_edges_loglik_10d_3400_3600i,
               hc_edges_loglik_10d_3600_3800i)){
  #quantiles <-stats::quantile(arc.strength(i$networks,as.data.frame(data.dag))[,3],probs = seq(0,1,0.05))
  plot_strongweak_edges(i$networks,data.dag, 
                        treshold = round(quantiles[2]),
                        which ="strong",
                        cex =0.2)
  title(paste0("number of edges = ",length(i$networks$arcs)/2))
}
dev.off()
quantiles
#####################################################################################################################
# make different quantile plots
#####################################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo Doctorado/R_practice/Data/tas_ncep_10d.rda")
# set variable names in data
nodes(hc_edges_loglik_10d_200i$networks)
data.dag <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
names(data.dag) <- nodes(hc_edges_loglik_10d_200i$networks)

par(mfrow = c(3, 2))
data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
for (i in list(hc_edges_loglik_10d_200i,
               hc_edges_loglik_10d_200_400i,
               hc_edges_loglik_10d_400_600i,
               hc_edges_loglik_10d_600_800i,
               hc_edges_loglik_10d_800_1000i,
               hc_edges_loglik_10d_1000_1200i,
               hc_edges_loglik_10d_1200_1400i,
               hc_edges_loglik_10d_1400_1600i,
               hc_edges_loglik_10d_1600_1800i,
               hc_edges_loglik_10d_1800_2000i,
               hc_edges_loglik_10d_2000_2200i,
               hc_edges_loglik_10d_2200_2400i,
               hc_edges_loglik_10d_2400_2600i,
               hc_edges_loglik_10d_2600_2800i,
               hc_edges_loglik_10d_2800_3000i)){
  
  plot_quantile_edges(i$networks, data.dag, cex = 0.5 )
  title(paste0("number of edges = ",length(i$networks$arcs)/2))
}
dev.off()


######################################################################################################
# degree area weighted 
# degree large distance.
# degree important links
######################################################################################################

# Calculate awconnectivity for DAGs  and return as Climatology object
areaweighted_dag_degree <- function(dag,grid,data.dag = NULL, what = c("both","in","out")){
  # distinguish between grid and time.coords.matrix:
  if (is.null(data.dag)) {data.dag <- TimeCoordsAnom_from_Grid(grid)} 

  lattitude <- attributes(data.dag)$VertexCoords$y
  adjmat <- amat(dag)
  sumArea <- sum(cos(lattitude/180*pi))
  
  # indegree
  if (what == "in") {awconnectivity <- as.vector(cos(lattitude/180*pi)%*%adjmat) / sumArea}
  # outdegree
  if (what == "out") {awconnectivity <- as.vector(adjmat%*%cos(lattitude/180*pi)) / sumArea}
  # in and out degree 
  if (what == "both") {awconnectivity <- (as.vector(adjmat%*%cos(lattitude/180*pi)) + as.vector(adjmat%*%cos(lattitude/180*pi))) / sumArea}
  
               
  mat <- matrix(awconnectivity, nrow = 1)

  grid$Data <- mat2Dto3Darray(mat, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))

  attr(grid$Data, "climatology:fun") <- paste0("awconnectivity of ",what," degree")
  return(grid)
}

################################################# Apply on test.
test <- hc_edges_loglik_10d_2800_3000i
 
testin <- areaweighted_dag_degree(test$networks, tas_ncep_10d, data10d, "in")
testout <-areaweighted_dag_degree(test$networks, tas_ncep_10d, data10d, "out")
testboth <- areaweighted_dag_degree(test$networks, tas_ncep_10d, data10d, "both")

Multi <- makeMultiGrid(testin,testout,testboth)
plotClimatology(Multi, backdrop.theme = "coastline", 
                names.attr = c("testin","testout","testboth"), 
                main = "area weighted degree")
plotClimatology(testout, backdrop.theme = "coastline", main = paste0(narcs(test$networks)))

############## apply area weighted connectivity on all iterations
hc_list <- list(hc_edges_loglik_10d_200i,
                 hc_edges_loglik_10d_200_400i,
                 hc_edges_loglik_10d_400_600i,
                 hc_edges_loglik_10d_600_800i,
                 hc_edges_loglik_10d_800_1000i,
                 hc_edges_loglik_10d_1000_1200i,
                 hc_edges_loglik_10d_1200_1400i,
                 hc_edges_loglik_10d_1400_1600i,
                 hc_edges_loglik_10d_1600_1800i,
                 hc_edges_loglik_10d_1800_2000i,
                 hc_edges_loglik_10d_2000_2200i,
                 hc_edges_loglik_10d_2200_2400i,
                 hc_edges_loglik_10d_2400_2600i,
                 hc_edges_loglik_10d_2600_2800i,
                 hc_edges_loglik_10d_2800_3000i,
                hc_edges_loglik_10d_3000_3200i,
                hc_edges_loglik_10d_3200_3400i,
                hc_edges_loglik_10d_3400_3600i,
                hc_edges_loglik_10d_3600_3800i)
hc_list <- list(hc_edges_loglik_10d_200i,
                hc_edges_loglik_10d_400_600i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1200_1400i,
                hc_edges_loglik_10d_1600_1800i,
                hc_edges_loglik_10d_2000_2200i,
                hc_edges_loglik_10d_2400_2600i,
                hc_edges_loglik_10d_2800_3000i)
#select networks
networks <- lapply(hc_list, function(m) m[["networks"]])
#calculate area weighted dag degree for all hillclimbing dags
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
awdds <-lapply(networks, areaweighted_dag_degree, grid =tas_ncep_10d , data.dag = data10d, what = "both")


# make multi grids, give names, and plot with plotClimatology
awddsclim <-makeMultiGrid(awdds[])
arcscounts <- sapply(networks, narcs)
what <- "both"
plotClimatology(awddsclim, backdrop.theme = "coastline", 
                names.attr = arcscounts, 
                main = paste0("area weighted",what,"degree")
                #,at = seq(0,0.08,0.005)
                )

#other methodç
plotClimatology(awdds[[1]])
for (i in 1:length(awdds)){
  assign(paste0("plot",i),plotClimatology(awdds[[i]],backdrop.theme = "coastline",main = paste0(arcscounts[i])))
}

grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6)
grid.arrange(plot6,plot7,plot8,plot9,plot10,plot11)
grid.arrange(plot11,plot12,plot13,plot14,plot15,plot16)

################### compare with others xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
clim.awcon.test10d <- measure2clim(measures10d, what = "awconnectivity", ref.grid = tas_ncep_10d)
plotClimatology(clim.awcon.test10d, backdrop.theme = "coastline")

################## long distance degree.


############################# HAVERSINE : introduce longitudes and latitudes
haversine <- function(x1Lat,x1Lon,x2Lat,x2Lon){
    earthR <- 6371 #using mean radius
    mLat <- as.double(x1Lat)
    bLat <- as.double(x2Lat)
    mLong <- as.double(x1Lon)
    bLong <- as.double(x2Lon)
    changeLat <- (mLat - bLat)/180*pi
    changeLong <- (mLong - bLong)/180*pi
    a <- sin(changeLat/2) * sin(changeLat/2) + cos((mLat)/180*pi) * 
      cos(bLat/180*pi) * sin(changeLong/2) * sin(changeLong/2)
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    distKm <- earthR * c
    distMi <- as.double(distKm * 0.621371192)
    output <- c(x1Lat,x1Lon,x2Lat,x2Lon,distKm,distMi)
    return(output[5])}             

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/utilnetworks_hcnetworks10d.rda")

plot_long_distances <- function(dag, data.dag, minimdist, smallcol = rgb(0,0,255, alpha = 125, maxColorValue = 255),perm = NULL, title = TRUE, remove = TRUE){
  dag <- gridGraphsBN$hc1_1800_1900i$graph
  # dag <- hc_2$networks$tabu_10d_505i
    # dag <- tabu_mmhc_simple2
    data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
    minimdist <- 10000
    perm <- permutations[[2]]
  require(data.table)  
  # Check if class is igraph or bn, add characteristics for plotting purpose  
  if (class(dag) == "igraph") {
    names <- c()
    names
    for(i in 1:648) names <- append(names,paste0("V",i))
    V(dag)$name <- names
    nodenames <- V(dag)$name
    igraph <- dag
    as.directed(igraph, mode = "arbitrary")
    }
  if (class(dag) == 'bn') {
    igraph <- igraph.from.graphNEL(as.graphNEL(dag))
    nodenames <- nodes(dag) }
 
  
  # Make edgelist in igraph class
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  edgeindices <- which(adjmat == 1, arr.ind = TRUE)
  edgeindices
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  # Make coordinates for every variable
  longitude <- attributes(data.dag)$VertexCoords$x
  lattitude <- attributes(data.dag)$VertexCoords$y
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }

  # estimate distance of all edges in igraph-edgelist.
  distances <- c()
  i <- 17
  for (i in 1:nrow(edgeindices)){
    x1Lat <- lattitude[edgeindices[i,1]]
    x1Lon <- longitude[edgeindices[i,1]]
    x2Lat <- lattitude[edgeindices[i,2]]
    x2Lon <- longitude[edgeindices[i,2]]

    disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
    distances[i] <- disti}

  # Make dataframe with departing variable, end variable and distance
  arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
  arcdistances

  #identify which indices in igraph edgelist correspond to indices in bn
  indsame <- c()
  nrow(as_edgelist(igraph))
  for (i in 1:nrow(arcdistances)){
    int <- intersect(which(as_edgelist(igraph)[,1] == arcdistances[i,1]),
                  which(as_edgelist(igraph)[,2] == arcdistances[i,2]))
    indsame[i] <- int
  }

    #permutate distance vector to find belong distance vector for igraph object
    newdistances <- numeric(nrow(as_edgelist(igraph)))
    for (i in 1:nrow(as_edgelist(igraph))){
      ind <- indsame[i]
      as_edgelist(igraph)[ind,]
      newdistances[ind] <- distances[i]
    }
    newdistances
    
    # Remove limit neighbours
    names <- c()
    for(i in 1:ncol(data.dag)) {names <- append(names,paste0("V",i))}
    nvertpts <- sqrt(ncol(data.dag)/2)
    limitsl <- names[1:nvertpts]
    limitsl <- t(limitsl)
    limitsr <- names[(ncol(data.dag)-nvertpts+1):ncol(data.dag)]
    limitsr <- t(limitsr)
    
    black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
    black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
    blacks <- rbind(black,black2)
    blacks
    
    m1 <- as_edgelist(igraph)
    m2 <- blacks
    colnames(m2) <- paste0("V", seq(len=ncol(m2)))
    DT1 <- data.table(m1)
    DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
    setnames(DT2, c(head(names(DT2), -1L), "found"))
    a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
    a <-as.matrix(a)
    
    # Color edges:
    # Large edges, small edges, and border edges
    E(igraph)$distances <- newdistances
    E(igraph)[E(igraph)$distances<= minimdist]$color <- smallcol
    E(igraph)[E(igraph)$distances > minimdist]$color <- "black"
    if(remove == TRUE) {
      E(igraph)[which(a[,1] == 1, arr.ind = TRUE)]$color <- NA}

    
    
    # draw the graph
    x <- attr(data.dag, "Xcoords", exact = FALSE)
    y <- attr(data.dag, "Ycoords", exact = FALSE)
    
    points <- expand.grid(y, x)[2:1]
    if (!is.null(perm)){
      points <- points[perm,]
    }
   
    load("/home/catharina/Documents/lisworld.rda")
    plot(lisworld)
    plot.igraph(igraph, 
            vertex.size = 100,
            vertex.color = "blue",
            vertex.label = NA,
            edge.arrow.size = 0.2,
            layout = as.matrix(points), add = TRUE, rescale = FALSE)
    if(title == TRUE){
      title(mar = c(0,0,0,0), main = paste0("minimdist = ",minimdist," |E| = ",nrow(arcs(dag))))
      }
}
######################################################
# Now in basic Network Functions 
######################################################
plot_long_distances <- function(dag, data.dag, minimdist, smallcol = rgb(0,0,255, alpha = 125, maxColorValue = 255),perm = NULL, title = "NA", remove = TRUE){
  # dag <- gridGraphsLattices$`1242`$graph
  # dag <- gridGraphsBN$hc1_1600_1700i$graph
  # dag <- hc_2$networks$tabu_10d_505i
  # dag <- tabu_mmhc_simple2
  # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
  # minimdist <- 10000
  # perm <- permutations[[2]]
  require(data.table)  
  # Check if class is igraph or bn, add characteristics for plotting purpose  
  if (class(dag) == "igraph") {
    names <- c()
    names
    for(i in 1:648) names <- append(names,paste0("V",i))
    V(dag)$name <- names
    nodenames <- V(dag)$name
    igraph <- dag
    as.directed(igraph, mode = "arbitrary")
  }
  if (class(dag) == 'bn') {
    igraph <- igraph.from.graphNEL(as.graphNEL(dag))
    nodenames <- nodes(dag) }
  
  
  # Make edgelist in igraph class

  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  adjmattri[lower.tri(adjmattri)] <- 0
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  edgeindices
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  # Make coordinates for every variable
  longitude <- attributes(data.dag)$VertexCoords$x
  lattitude <- attributes(data.dag)$VertexCoords$y
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }
  
  # estimate distance of all edges in igraph-edgelist.
  distances <- c()
  for (i in 1:nrow(edgeindices)){
    x1Lat <- lattitude[edgeindices[i,1]]
    x1Lon <- longitude[edgeindices[i,1]]
    x2Lat <- lattitude[edgeindices[i,2]]
    x2Lon <- longitude[edgeindices[i,2]]
    
    disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
    distances[i] <- disti}
  
  # Make dataframe with departing variable, end variable and distance
  arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
  arcdistances
  
  #identify which indices in igraph edgelist correspond to indices in bn
  indsame <- c()
  nrow(as_edgelist(igraph))
  for (i in 1:nrow(arcdistances)){
    int <- intersect(which(as_edgelist(igraph)[,1] == arcdistances[i,1]),
                     which(as_edgelist(igraph)[,2] == arcdistances[i,2]))
    indsame[i] <- int
  }
  
  #permutate distance vector to find belong distance vector for igraph object
  newdistances <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    newdistances[ind] <- distances[i]
  }
  newdistances
  
  # Remove limit neighbours
  names <- c()
  for(i in 1:ncol(data.dag)) {names <- append(names,paste0("V",i))}
  nvertpts <- sqrt(ncol(data.dag)/2)
  limitsl <- names[1:nvertpts]
  limitsl <- t(limitsl)
  limitsr <- names[(ncol(data.dag)-nvertpts+1):ncol(data.dag)]
  limitsr <- t(limitsr)
  
  black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
  black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
  blacks <- rbind(black,black2)
  blacks
  
  m1 <- as_edgelist(igraph)
  m2 <- blacks
  colnames(m2) <- paste0("V", seq(len=ncol(m2)))
  DT1 <- data.table(m1)
  DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
  setnames(DT2, c(head(names(DT2), -1L), "found"))
  a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
  a <-as.matrix(a)
  
  # Color edges:
  # Large edges, small edges, and border edges
  E(igraph)$distances <- newdistances
  E(igraph)[E(igraph)$distances<= minimdist]$color <- smallcol
  E(igraph)[E(igraph)$distances > minimdist]$color <- "black"
  if(remove == TRUE) {
    E(igraph)[which(a[,1] == 1, arr.ind = TRUE)]$color <- NA}
  
  
  
  # draw the graph
  x <- attr(data.dag, "Xcoords", exact = FALSE)
  y <- attr(data.dag, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  if (!is.null(perm)){
    points <- points[perm,]
  }
  
  # load("/home/catharina/Documents/lisworld.rda")
  # load("/Volumes/ubuntu/Documents/lisworld.rda")
  load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  raster::plot(lisworld)
  plot.igraph(igraph, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  if(title == "A"){
    title(mar = c(0,0,0,0), main = paste0("minimdist = ",minimdist," |E| = ",nrow(as_edgelist(igraph))))
  }
  if(title == "B"){
    title(mar = c(0,0,0,0), main = paste0(" |E| = ",nrow(as_edgelist(igraph))))
  }
}
################################################################
# Function set_distances: in Basic Network Functions
################################################################

# set_distances <- function(graphObj){
#   # graphObj <- gridGraphsBN[[3]]
#   # dag <- gridGraphsLattices$`1242`$graph
#   # dag <- gridGraphsBN$hc1_1600_1700i$graph
#   # dag <- hc_2$networks$tabu_10d_505i
#   # dag <- tabu_mmhc_simple2
#   # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
#   # minimdist <- 10000
#   # perm <- permutations[[2]]
#   require(data.table)  
#   # Check if class is igraph or bn, add characteristics for plotting purpose  
#   # if (class(dag) == "igraph") {
#   #   names <- c()
#   #   names
#   #   for(i in 1:648) names <- append(names,paste0("V",i))
#   #   V(dag)$name <- names
#   #   nodenames <- V(dag)$name
#   #   igraph <- dag
#   #   as.directed(igraph, mode = "arbitrary")
#   # }
#   # if (class(dag) == 'bn') {
#   #   igraph <- igraph.from.graphNEL(as.graphNEL(dag))
#   #   nodenames <- nodes(dag) }
#   # 
#   
#   # Make edgelist in igraph class
#   
#   adjmat <- as.matrix(as_adjacency_matrix(graphObj$graph))
#   adjmattri <- adjmat
#   adjmattri[lower.tri(adjmattri)] <- 0
#   edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
#   edgeindices
#   fromV <- nodenames[edgeindices[,1]]
#   toV <- nodenames[edgeindices[,2]]
#   fromtoV <- cbind(fromV,toV)
#   
#   # Make coordinates for every variable
#   longitude <- graphObj$VertexCoords$x
#   lattitude <- graphObj$VertexCoords$y
#   # if (!is.null(perm)){
#   #   longitude <- longitude[perm]
#   #   lattitude <- lattitude[perm]
#   # }
#   
#   # estimate distance of all edges in igraph-edgelist.
#   distances <- c()
#   for (i in 1:nrow(edgeindices)){
#     x1Lat <- lattitude[edgeindices[i,1]]
#     x1Lon <- longitude[edgeindices[i,1]]
#     x2Lat <- lattitude[edgeindices[i,2]]
#     x2Lon <- longitude[edgeindices[i,2]]
#     
#     disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
#     distances[i] <- disti}
#   
#   # Make dataframe with departing variable, end variable and distance
#   arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
#   arcdistances
#   
#   #identify which indices in igraph edgelist correspond to indices in bn
#   indsame <- c()
#   nrow(as_edgelist(graphObj$graph))
#   for (i in 1:nrow(arcdistances)){
#     int <- intersect(which(as_edgelist(graphObj$graph)[,1] == arcdistances[i,1]),
#                      which(as_edgelist(graphObj$graph)[,2] == arcdistances[i,2]))
#     indsame[i] <- int
#   }
#   
#   #permutate distance vector to find belong distance vector for igraph object
#   newdistances <- numeric(nrow(as_edgelist(graphObj$graph)))
#   for (i in 1:nrow(as_edgelist(graphObj$graph))){
#     ind <- indsame[i]
#     as_edgelist(graphObj$graph)[ind,]
#     newdistances[ind] <- distances[i]
#   }
#   newdistances
#   
#   # # Remove limit neighbours
#   # names <- c()
#   # for(i in 1:ncol(graphObj$time.coords.matrix)) {names <- append(names,paste0("V",i))}
#   # nvertpts <- sqrt(ncol(graphObj$time.coords.matrix)/2)
#   # limitsl <- names[1:nvertpts]
#   # limitsl <- t(limitsl)
#   # limitsr <- names[(ncol(graphObj$time.coords.matrix)-nvertpts+1):ncol(graphObj$time.coords.matrix)]
#   # limitsr <- t(limitsr)
#   # 
#   # black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
#   # black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
#   # blacks <- rbind(black,black2)
#   # blacks
#   # 
#   # m1 <- as_edgelist(graphObj$graph)
#   # m2 <- blacks
#   # colnames(m2) <- paste0("V", seq(len=ncol(m2)))
#   # DT1 <- data.table(m1)
#   # DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
#   # setnames(DT2, c(head(names(DT2), -1L), "found"))
#   # a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
#   # a <-as.matrix(a)
#   # 
#   # Color edges:
#   # Large edges, small edges, and border edges
#   E(graphObj$graph)$weight <- newdistances
#   return(graphObj)
# }

withweight3 <- set_distances(gridGraphsBN[[3]])
E(withweight3$graph)$weight
E(gridGraphsBN[[3]]$graph)$weight
diameter(withweight3$graph)
diameter(gridGraphsBN[[3]]$graph)

dev.off()
10 * pi*6000/180
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plot_long_distances(gridGraphsBN$hc1_900_1000i$graph, data10d, minimdist = 1000)
plot_long_distances(pc_10d_1e_13$networks$dag, data10d, minimdist = 1000)

hc_list <- list(hc_edges_loglik_10d_200i,
                hc_edges_loglik_10d_400_600i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1200_1400i,
                hc_edges_loglik_10d_1600_1800i,
                hc_edges_loglik_10d_2000_2200i,
                hc_edges_loglik_10d_2400_2600i,
                hc_edges_loglik_10d_2800_3000i)
hc_list <- list(hc_edges_loglik_10d_1400_1600i,
                hc_edges_loglik_10d_1600_1800i,
                hc_edges_loglik_10d_1800_2000i,
                hc_edges_loglik_10d_2000_2200i,
                hc_edges_loglik_10d_2200_2400i,
                hc_edges_loglik_10d_2400_2600i,
                hc_edges_loglik_10d_2600_2800i,
                hc_edges_loglik_10d_2800_3000i)
hc_list <- list(hc_edges_loglik_10d_200i,
                hc_edges_loglik_10d_200_400i,
                hc_edges_loglik_10d_400_600i,
                hc_edges_loglik_10d_600_800i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1000_1200i,
                hc_edges_loglik_10d_1200_1400i,
                hc_edges_loglik_10d_1400_1600i)
#select networks
networks <- lapply(hc_list, function(m) m[["networks"]])
#calculate area weighted dag degree for all hillclimbing dags
par(mfrow = c(4,2),cex = 0.4)
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
lapply(networks, plot_long_distances, data.dag = data10d, minimdist = 12000)

dev.off()
############################################################################ Voor Complex graph
igraph <- graph10d$graph
fromV <- V(igraph)[edgeindices[,1]]
toV <- V(igraph)[edgeindices[,2]]

################################# degree long distances.


mean_distance_edges <- function(dag,grid,data.dag, what = c("out","in","both")){
  
  
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
  which(E(igraph)$distances < 2000)
  E(igraph)[E(igraph)$distances < 2000]
  
  if (what == "out"){
  meandist <- unlist(lapply(V(igraph), function(x) mean(E(igraph)[from(x)]$distances)))}
  if (what == "in"){
  meandist <- unlist(lapply(V(igraph), function(x) mean(E(igraph)[to(x)]$distances)))}
  if (what == "both"){
  meandist <- unlist(lapply(V(igraph), function(x) mean(E(igraph)[from(x) | to(x)]$distances)))}
  
  mat <- matrix(meandist, nrow = 1)
  mat
  
  
  grid$Data <- mat2Dto3Darray(mat, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
  
  attr(grid$Data, "climatology:fun") <- "meandist"
  return(grid)}
  

data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
meandist <- mean_distance_edges(hc_edges_loglik_10d_400_600i$networks,tas_ncep_10d,data.dag, what = "out")
plotClimatology(meandist, backdrop.theme = "coastline")

hc_list <- list(hc_edges_loglik_10d_200i,
                hc_edges_loglik_10d_200_400i,
                hc_edges_loglik_10d_400_600i,
                hc_edges_loglik_10d_600_800i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1000_1200i,
                hc_edges_loglik_10d_1200_1400i,
                hc_edges_loglik_10d_1400_1600i)
hc_list <- list(hc_edges_loglik_10d_1400_1600i,
                hc_edges_loglik_10d_1600_1800i,
                hc_edges_loglik_10d_1800_2000i,
                hc_edges_loglik_10d_2000_2200i,
                hc_edges_loglik_10d_2200_2400i,
                hc_edges_loglik_10d_2400_2600i,
                hc_edges_loglik_10d_2600_2800i,
                hc_edges_loglik_10d_2800_3000i)
#select networks
networks <- lapply(hc_list, function(m) m[["networks"]])
nedgesnetworks <- as.character(sapply(networks, narcs))


grids <- lapply(networks, mean_distance_edges, grid = tas_ncep_10d, data.dag = data10d)
grids[[1]]
plotClimatology(grids[[1]])

multis = list()
for (i in 1:length(grids)){
plot <- plotClimatology(grids[[i]], main = paste0("number of edges = ",nedgesnetworks[i]) ,backdrop.theme = "coastline")
multis[[i]] <- plot
}
multis[[2]]
n <- length(multis)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(multis, ncol=nCol))

Multigrid <- makeMultiGrid(grids)
plotClimatology(Multigrid, backdrop.theme = "coastline")
# Make MultiGrid as in areaweighted dag_degree

# For degree long distanc:
#newgraph <- delete.edges(igraph, E(igraph)[E(igraph)$distances < 2000])
#newgraph


#as.vector(igraph::degree(newgraph))
#degree.distribution(newgraph)
