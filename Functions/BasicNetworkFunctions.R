###############################################################################
# Functions: 
# plot.Meteodag
# Time_Coords_Anom_from_Grid_RMS
# Time_Coords_Anom_from_Grid_std
# Time_Coords_Anom_from_Grid_RMS
# Time_Coords_Anom_from_Grid_aslist
#
if(grep("Untitled",getwd())==1){
  load("/data/Untitled/Trabajo/R_practice/Data/lisworld.rda")
  load("/data/Untitled/Trabajo/R_practice/Data/lisworldm.rda")
  load("/data/Untitled/Trabajo/R_practice/Data/coastline.rda")
} else {
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworldm.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/coastline.rda")
  }

library(data.table)
library(sp)
library(rgeos)
library(transformeR)
###############################################################################
plot.Meteodag <- function(time.coords, meteodag, lis = FALSE, edgecolor = "green"){
  if (class(meteodag) == "igraph") igraphDegree <- meteodag
  
  if (class(meteodag) == "graphNEL") igraphDegree <- igraph.from.graphNEL(meteodag)
  
  if (class(meteodag) == "bn") igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))
  
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  
  if (lis == TRUE){plot(lisworld)}
  else {plot(wrld)}

  plot.igraph(igraphDegree, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.color= edgecolor,
              edge.arrow.size = 0.3,
              edge.lty = 2,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  
}

TimeCoordsAnom_from_Grid <- function(grid, subind = NULL) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]}
  out <- time.coords.matrix
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}

TimeCoordsAnom_from_Grid_std <- function(grid, subind = NULL) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    #WATCH OUT SCALEGRID(scale = TRUE) changed in scaleGRID
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE) })
  #localScaling(grid,time.frame = "monthly")
  
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]}
  
  out <- time.coords.matrix
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}


TimeCoordsAnom_from_Grid_rms <- function(grid, subind = NULL, rms = FALSE,mask =NULL) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE) })
  #localScaling(grid,time.frame = "monthly")
  
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  if (rms == TRUE) {
    time.coords.matrix <- scale(time.coords.matrix, center = FALSE, scale = TRUE)}
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]
    row.names(time.coords.matrix) <- as.character(subind)}
  if (!is.null(mask)) {
    landid <- which(!is.na(time.coords.matrix),arr.ind=TRUE)
    landnodes<- unique(landid[,2])
    time.coords.matrix <- time.coords.matrix[,landnodes]
    }
  
  out <- time.coords.matrix
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}

TimeCoordsAnom_from_Grid_stand1_2times <- function(grid, subind = NULL,twotimes = FALSE) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid(type = "standardize") %>% redim(drop = TRUE) })
  #localScaling(grid,time.frame = "monthly")
  
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  if(twotimes == TRUE){
    time.coords.matrix <- scale(time.coords.matrix, center = TRUE, scale = TRUE)}
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]
    row.names(time.coords.matrix) <- as.character(subind)}
  
  out <- time.coords.matrix
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}

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



TimeCoordsAnom_from_Grid_aslist <- function(grid) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% scaleGrid() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <-  bindGrid(seas.list, dimension = "time")
  aux <-  redim(aux, drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  
  out <- list("data_coords" = time.coords.matrix,
              "VertexCoords" = ref.coords)
  
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}

set_distances <- function(graphObj){
  # Make edgelist in igraph class
  if (!is.null(names(graphObj$graph))){
  names <- c()
  for(i in 1:648) names <- append(names,paste0("V",i))
  nodenames <- names
  }
  
  adjmat <- as.matrix(as_adjacency_matrix(graphObj$graph))
  adjmattri <- adjmat
  adjmattri[lower.tri(adjmattri)] <- 0
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  edgeindices
  if (!is.null(names(graphObj$graph))){
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  } else {
    fromV <- edgeindices[,1]
    toV <- edgeindices[,2]
    fromtoV <- cbind(fromV,toV)
  }
  
  # Make coordinates for every variable
  longitude <- graphObj$VertexCoords$x
  lattitude <- graphObj$VertexCoords$y
  # if (!is.null(perm)){
  #   longitude <- longitude[perm]
  #   lattitude <- lattitude[perm]
  # }
  
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
  nrow(as_edgelist(graphObj$graph))

  for (i in 1:nrow(arcdistances)){
    int <- intersect(which(as_edgelist(graphObj$graph)[,1] == arcdistances[i,1]),
                     which(as_edgelist(graphObj$graph)[,2] == arcdistances[i,2]))
    indsame[i] <- int
  }
  
  #permutate distance vector to find belong distance vector for igraph object
  newdistances <- numeric(nrow(as_edgelist(graphObj$graph)))
  for (i in 1:nrow(as_edgelist(graphObj$graph))){
    ind <- indsame[i]
    as_edgelist(graphObj$graph)[ind,]
    newdistances[ind] <- distances[i]
  }
  newdistances
  
  E(graphObj$graph)$weight <- newdistances
  return(graphObj)
}



movingmedias <- function(x){
  vect <- x
  names <- character(length = 648)
  for (i in 1:648){names[i] <- paste0("V",i)}
  
  dimnames <-list("row" = c("r1","r2","r3"),"column"=c("c1",'c2',"c3"),"coord" = names)
  matrixlist <- array(data =NA, dim = c(3,3,648), dimnames = dimnames)
  
  for (i in 1:648){
    if (i<=18|i%%18 == 0) {
      if (i%%18 == 0){ # j1 in matrix
        matrixlist[1,1,i] <- NA
      } else {matrixlist[1,1,i] <- vect[i-17 + 648]} 
    } else {matrixlist[1,1,i] <- vect[i-17]} 
    
    if (i%%18 == 0){ # j2 in matrix
      matrixlist[1,2,i] <- NA
    } else {matrixlist[1,2,i] <-vect[i+1]}
    
    if (i>= 631|i%%18 == 0){ # j3 in matrix
      if (i%%18 == 0){
        matrixlist[1,3,i] <- NA 
      } else  {matrixlist[1,3,i] <- vect[i+19 - 648]}
    } else {matrixlist[1,3,i] <- vect[i+19]}
    
    if (i<=18){ # j4 in matrix
      matrixlist[2,1,i] <- vect[i-18 + 648]
    } else {matrixlist[2,1,i] <- vect[i-18]}
    
    matrixlist[2,2,i] <- vect[i] # j5 in matrix
    
    if (i>= 631){ # j6 in matrix
      matrixlist[2,3,i] <- vect[i+18 -648]
    } else matrixlist[2,3,i] <- vect[i+18]
    
    if (i<=18|i%%18 == 1){ # j7 in matrix
      if (i%%18 == 1) {matrixlist[3,1,i] <- NA
      } else {matrixlist[3,1,i] <- vect[i-19 + 648]}
    } else {matrixlist[3,1,i] <- vect[i-19]}
    
    if (i%%18 == 1){ # j8 in matrix
      matrixlist[3,2,i] <- NA
    } else {matrixlist[3,2,i] <-vect[i-1]}
    
    if (i%%18 == 1 |i>= 631){ # j9 in matrix
      if (i%%18 == 1){
        matrixlist[3,3,i] <- NA
      } else {matrixlist[3,3,i] <- vect[i+17-648]}
    } else {matrixlist[3,3,i] <- vect[i+17]}
  }
  return(matrixlist)
}



#movingmedias <- function(x){
  #vect <- x
  #names <- character(length = 648)
  #for (i in 1:648){names[i] <- paste0("V",i)}
  
  #dimnames <-list("row" = c("r1","r2","r3"),"column"=c("c1",'c2',"c3"),"coord" = names)
  #matrixlist <- array(data =NA, dim = c(3,3,648), dimnames = dimnames)
  
#   for (i in 1:648){
#     if (i<=18|i%%18 == 0) {
#       if (i%%18 == 0){ # j1 in matrix
#         matrixlist[1,1,i] <- 0
#       } else {matrixlist[1,1,i] <- vect[i-17 + 648]} 
#     } else {matrixlist[1,1,i] <- vect[i-17]} 
#     
#     if (i%%18 == 0){ # j2 in matrix
#       matrixlist[1,2,i] <- 0
#     } else {matrixlist[1,2,i] <-vect[i+1]}
#     
#     if (i>= 631|i%%18 == 0){ # j3 in matrix
#       if (i%%18 == 0){
#         matrixlist[1,3,i] <- 0 
#       } else  {matrixlist[1,3,i] <- vect[i+19 - 648]}
#     } else {matrixlist[1,3,i] <- vect[i+19]}
#     
#     if (i<=18){ # j4 in matrix
#       matrixlist[2,1,i] <- vect[i-18 + 648]
#     } else {matrixlist[2,1,i] <- vect[i-18]}
#     
#     matrixlist[2,2,i] <- vect[i] # j5 in matrix
#     
#     if (i>= 631){ # j6 in matrix
#       matrixlist[2,3,i] <- vect[i+18 -648]
#     } else matrixlist[2,3,i] <- vect[i+18]
#     
#     if (i<=18|i%%18 == 1){ # j7 in matrix
#       if (i%%18 == 1) {matrixlist[3,1,i] <- 0
#       } else {matrixlist[3,1,i] <- vect[i-19 + 648]}
#     } else {matrixlist[3,1,i] <- vect[i-19]}
#     
#     if (i%%18 == 1){ # j8 in matrix
#       matrixlist[3,2,i] <- 0
#     } else {matrixlist[3,2,i] <-vect[i-1]}
#     
#     if (i%%18 == 1 |i>= 631){ # j9 in matrix
#       if (i%%18 == 1){
#         matrixlist[3,3,i] <- 0
#       } else {matrixlist[3,3,i] <- vect[i+17-648]}
#     } else {matrixlist[3,3,i] <- vect[i+17]}
#   }
#   return(matrixlist)
# }


# movingmedias <- function(x){
#   vect <- x
#   nnodes <- length(vect)
#   nlat <- sqrt(nnodes/2)
#   nlong <- 2*nlat
#   
#   names <- character(length = nnodes)
#   for (i in 1:nnodes){names[i] <- paste0("V",i)}
#   
#   dimnames <-list("row" = c("r1","r2","r3"),"column"=c("c1",'c2',"c3"),"coord" = names)
#   matrixlist <- array(data =NA, dim = c(3,3,nnodes), dimnames = dimnames)
#   
#   for (i in 1:nnodes){
#     if (i<=nlat|i%%nlat == 0) {
#       if (i%%nlat == 0){ # j1 in matrix
#         matrixlist[1,1,i] <- 0
#       } else {matrixlist[1,1,i] <- vect[i-(nlat-1) + nnodes]} 
#     } else {matrixlist[1,1,i] <- vect[i-(nlat-1)]} 
#     
#     if (i%%nlat == 0){ # j2 in matrix
#       matrixlist[1,2,i] <- 0
#     } else {matrixlist[1,2,i] <-vect[i+1]}
#     
#     if (i>= (nnodes-(nlat-1))|i%%nlat == 0){ # j3 in matrix
#       if (i%%nlat == 0){
#         matrixlist[1,3,i] <- 0 
#       } else  {matrixlist[1,3,i] <- vect[i+(nlat+1) - nnodes]}
#     } else {matrixlist[1,3,i] <- vect[i+(nlat+1)]}
#     
#     if (i<=nlat){ # j4 in matrix
#       matrixlist[2,1,i] <- vect[i-nlat + nnodes]
#     } else {matrixlist[2,1,i] <- vect[i-nlat]}
#     
#     matrixlist[2,2,i] <- vect[i] # j5 in matrix
#     
#     if (i>= (nnodes-(nlat-1))){ # j6 in matrix
#       matrixlist[2,3,i] <- vect[i+nlat -nnodes]
#     } else matrixlist[2,3,i] <- vect[i+nlat]
#     
#     if (i<=nlat|i%%nlat == 1){ # j7 in matrix
#       if (i%%nlat == 1) {matrixlist[3,1,i] <- 0
#       } else {matrixlist[3,1,i] <- vect[i-(nlat+1) + nnodes]}
#     } else {matrixlist[3,1,i] <- vect[i-(nlat+1)]}
#     
#     if (i%%nlat == 1){ # j8 in matrix
#       matrixlist[3,2,i] <- 0
#     } else {matrixlist[3,2,i] <-vect[i-1]}
#     
#     if (i%%nlat == 1 |i>= (nnodes-(nlat-1))){ # j9 in matrix
#       if (i%%nlat == 1){
#         matrixlist[3,3,i] <- 0
#       } else {matrixlist[3,3,i] <- vect[i+(nlat-1)-nnodes]}
#     } else {matrixlist[3,3,i] <- vect[i+(nlat-1)]}
#   }
#   return(matrixlist)
# }




#########################################################################
# plot_long_distances
# original from arcstrengthfunctions
# uses haversine. 
#########################################################################
plot_long_distances <- function(dag, data.dag, minimdist, smallcol = rgb(0,0,255, alpha = 125, maxColorValue = 255),bigcol = "black",perm = NULL, title = "NA", remove = TRUE){
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
  E(igraph)[E(igraph)$distances > minimdist]$color <- bigcol
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
  # load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  load("../Data/lisworld.rda")
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

plot_undirected_distances <- function(dag, data.dag, minimdist, smallcol = rgb(0,0,255, alpha = 125, maxColorValue = 255),perm = NULL, title = "NA", remove = TRUE, arrowmode = FALSE){
  # dag <- gridGraphsLattices$`1242`$graph
  # dag <- gridGraphsBN$hc1_1600_1700i$graph
  # dag <- pcstable2_withnet_2$networks$pcstable_10d_g4
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
  
  # # estimate distance of all edges in igraph-edgelist.
  # distances <- c()
  # for (i in 1:nrow(edgeindices)){
  #   x1Lat <- lattitude[edgeindices[i,1]]
  #   x1Lon <- longitude[edgeindices[i,1]]
  #   x2Lat <- lattitude[edgeindices[i,2]]
  #   x2Lon <- longitude[edgeindices[i,2]]
  #   
  #   disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
  #   distances[i] <- disti}
  # 
  # # Make dataframe with departing variable, end variable and distance
  # arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
  # arcdistances
  # 
  # #identify which indices in igraph edgelist correspond to indices in bn
  # indsame <- c()
  # nrow(as_edgelist(igraph))
  # for (i in 1:nrow(arcdistances)){
  #   int <- intersect(which(as_edgelist(igraph)[,1] == arcdistances[i,1]),
  #                    which(as_edgelist(igraph)[,2] == arcdistances[i,2]))
  #   indsame[i] <- int
  # }
  # 
  # #permutate distance vector to find belong distance vector for igraph object
  # newdistances <- numeric(nrow(as_edgelist(igraph)))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   ind <- indsame[i]
  #   as_edgelist(igraph)[ind,]
  #   newdistances[ind] <- distances[i]
  # }
  # newdistances
  
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

  E(igraph)[which_mutual(igraph)]$color <- smallcol
  E(igraph)[!which_mutual(igraph)]$color <- "darkgrey"
  
  if(arrowmode == TRUE){
  arrowmode <- rep(2,length(E(igraph)))
  arrowmode[which_mutual(igraph)] <- 0} else {arrowmode <- NULL}
  
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
  # load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  raster::plot(lisworld)
  plot.igraph(igraph, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              edge.arrow.mode = arrowmode,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  if(title == "A"){
    title(mar = c(0,0,0,0), main = paste0("minimdist = ",minimdist," |E| = ",nrow(as_edgelist(igraph))))
  }
  if(title == "B"){
    title(mar = c(0,0,0,0), main = paste0(" |E| = ",nrow(as_edgelist(igraph))))
  }
}

plot_long_strong_distances <- function(dag, data.dag, minimdist, smallcol = rgb(0,0,255, alpha = 125, maxColorValue = 255),k = NULL, criterion = NULL, perm = NULL, which = c("normal","strongest","half"),title = "NA", remove = TRUE){
    # dag <- proefDAG
    # criterion <- NULL
    # k <- NULL
   # minimweight <- 0.037
   # smallcol <- NA
   # remove <- TRUE
  # dag <- gridGraphsLattices$`1242`$graph
  # dag <- gridGraphsBN$hc1_1600_1700i$graph
  # dag <- hc_2$networks$tabu_10d_505i
  # dag <- tabu_mmhc_simple2
    # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
    # minimdist <- 10000
     # perm <- permutations[[5]]
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
    arcsbn <- arcs(dag)
    linkstrength <- arc.strength(dag, as.data.frame(data.dag), criterion = criterion ,k  = k)
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

  # Make edgelist in igraph class method 2
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

  # permutate distance vector to find belong distance vector for igraph object
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
  c_scale <- rev(topo.colors(5))
  E(igraph)$distances <- newdistances
  quantiles<- quantile(E(igraph)[E(igraph)$distances > minimdist]$weights)

  if(which == "strongest"){
    E(igraph)[weights <= quantiles[2]]$color <- "NA"
    E(igraph)[weights > quantiles[2] & weights <=quantiles[3]]$color <- "NA"
    E(igraph)[weights > quantiles[3] & weights <=quantiles[4]]$color <- "NA"
    E(igraph)[weights > quantiles[4] & weights <=quantiles[5]]$color <- c_scale[5]
  } else if(which == "half"){
    E(igraph)[weights <= quantiles[2]]$color <- "NA"
    E(igraph)[weights > quantiles[2] & weights <=quantiles[3]]$color <- "NA"
    E(igraph)[weights > quantiles[3] & weights <=quantiles[4]]$color <- c_scale[4]
    E(igraph)[weights > quantiles[4] & weights <=quantiles[5]]$color <- c_scale[5]
  } else if (which == "normal"){
    E(igraph)[weights <= quantiles[2]]$color <- c_scale[2]
    E(igraph)[weights > quantiles[2] & weights <=quantiles[3]]$color <- c_scale[3]
    E(igraph)[weights > quantiles[3] & weights <=quantiles[4]]$color <- c_scale[4]
    E(igraph)[weights > quantiles[4] & weights <=quantiles[5]]$color <- c_scale[5]
  }

  # E(igraph)[E(igraph)$weights <= minimweight]$color <- "orange"
  # E(igraph)[E(igraph)$weights > minimweight]$color <- "black"
  E(igraph)[E(igraph)$distances<= minimdist]$color <- smallcol

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
  # load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
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

plot_long_strong_v_distances <- function(dag, data.dag, minimdist, smallcol = rgb(0,0,255, alpha = 125, maxColorValue = 255),k = NULL, criterion = NULL, perm = NULL, title = "NA", remove = TRUE){
   dag <- proefDAG
  # criterion <- NULL
  # k <- NULL
  # minimweight <- 0.037
  # smallcol <- NA
  # remove <- TRUE
  # dag <- gridGraphsLattices$`1242`$graph
  # dag <- gridGraphsBN$hc1_1600_1700i$graph
  # dag <- hc_2$networks$tabu_10d_505i
  # dag <- tabu_mmhc_simple2
  # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
  # minimdist <- 10000
  # perm <- permutations[[1]]
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
    arcsbn <- arcs(dag)
    linkstrength <- arc.strength(dag, as.data.frame(data.dag), criterion = criterion ,k  = k)
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
  
  # Make edgelist in igraph class method 2
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
  
  # permutate distance vector to find belong distance vector for igraph object
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
  c_scale <- rev(topo.colors(5))
  E(igraph)$distances <- newdistances
  quantiles<- quantile(E(igraph)[E(igraph)$distances > minimdist]$weights)
  
  E(igraph)[weights <= quantiles[2]]$color <- c_scale[2]
  E(igraph)[weights > quantiles[2] & weights <=quantiles[3]]$color <- c_scale[3]
  E(igraph)[weights > quantiles[3] & weights <=quantiles[4]]$color <- c_scale[4]
  E(igraph)[weights > quantiles[4] & weights <=quantiles[5]]$color <- c_scale[5]
  
  # E(igraph)[E(igraph)$weights <= minimweight]$color <- "orange"
  # E(igraph)[E(igraph)$weights > minimweight]$color <- "black"
  E(igraph)[E(igraph)$distances<= minimdist]$color <- smallcol
  
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
  # load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
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




str_vs_dis <- function(dag, data.dag, k = NULL, criterion = NULL, perm = NULL){
  # dag <- proefDAG
  # criterion <- NULL
  # k <- NULL
  # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
  # perm <- permutations[[1]]
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
    arcsbn <- arcs(dag)
    linkstrength <- arc.strength(dag, as.data.frame(data.dag), criterion = criterion ,k  = k)
    igraph <- igraph.from.graphNEL(as.graphNEL(dag))
    nodenames <- nodes(dag) }
  
  
  
  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  # adjmattri[lower.tri(adjmattri)] <- 0
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  edgeindices
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  # Make edgelist in igraph class method 2
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
  normalize <- max(strengths)-min(strengths)
  E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  
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
  
  # identify which indices in igraph edgelist correspond to indices in bn
  indsame <- c()
  nrow(as_edgelist(igraph))
  for (i in 1:nrow(arcdistances)){
    int <- intersect(which(as_edgelist(igraph)[,1] == arcdistances[i,1]),
                     which(as_edgelist(igraph)[,2] == arcdistances[i,2]))
    indsame[i] <- int
  }
  
  # permutate distance vector to find belong distance vector for igraph object
  newdistances <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    newdistances[ind] <- distances[i]
  }
  newdistances
  E(igraph)$distances <- newdistances
  
  return(igraph)
}

str_vs_dis_2 <- function(dag, data.dag, k = NULL, criterion = NULL, perm = NULL, othermatrix = NULL){
     # dag <- GraphsCN[[50]]
  # # criterion <- NULL
  # # k <- NULL
   # data.dag <- data
  # perm <- permutations[[1]]
  # # require(data.table)
  
  # Check if class is igraph or bn, add characteristics for plotting purpose  
  if (class(dag) == "igraph") {
    names <- c()
    names
    for(i in 1:ncol(data.dag)) names <- append(names,paste0("V",i))
    V(dag)$name <- names
    nodenames <- V(dag)$name
    igraph <- dag
  }
  if (class(dag) == 'bn') {
    arcsbn <- arcs(dag)
    linkstrength <- arc.strength(dag, as.data.frame(data.dag), criterion = criterion ,k  = k)
    igraph <- igraph.from.graphNEL(as.graphNEL(dag))
    nodenames <- nodes(dag) }
  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  adjmattri[lower.tri(adjmattri)] <- 0
  cormatrix <- abs(cor(data.dag))
  rownames(cormatrix) <- colnames(cormatrix) <-  names
  cormatrix[lower.tri(cormatrix)] <- 0
  
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  noedgeindices <-  which(adjmattri == 0, arr.ind = TRUE)
  cormatrix[noedgeindices] <- 0
  
  edgeindices2 <- which(cormatrix != 0, arr.ind = TRUE)
  all.equal(edgeindices, edgeindices2)
  
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  # edgeindices[1,]
  # cormatrix[edgeindices[,]]
  # 
  
  # identify which indices in igraph edgelist correspond to indices in adjmatrix
  # in order to add edgeattributes to igraph
  indsame <- c()
  for (i in 1:nrow(fromtoV)){
    int <- intersect(which(fromtoV[i,1] == as_edgelist(igraph)[,1] ),
                     which(fromtoV[i,2] == as_edgelist(igraph)[,2] ))
    indsame[i] <- int
  }
  
  if (is.null(othermatrix)) {
    strengths <- cormatrix[edgeindices[,]]
    } else {strengths <- othermatrix[edgeindices[,]]}

  
  # strengthsedgelist <- strengths[indsame]
  
  # i <- 1
  # fromtoV[59,]
  # as_edgelist(igraph)[1,]
  # newstrengths <- numeric(nrow(as_edgelist(igraph)))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   ind <- indsame[i]
  #   as_edgelist(igraph)[ind,]
  #   newstrengths[i] <- strengths[ind]
  # }
  
  # E(igraph)$strengths <- strengthsedgelist
  # normalize <- max(strengths)-min(strengths)
  # E(igraph)$weights <- (strengthsedgelist-min(strengthsedgelist))/(max(strengthsedgelist) - min(strengthsedgelist))
  
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
  
  
  # test <- arcdistances[,3]
  # newdistances <- round(distances, digits = 3)[indsame]
  as_edgelist(igraph)[1,]
    # permutate distance vector to find belong distance vector for igraph object
    newdistances <- numeric(nrow(as_edgelist(igraph)))
    for (i in 1:nrow(as_edgelist(igraph))){
      ind <- indsame[i]
      as_edgelist(igraph)[ind,]
      newdistances[ind] <- distances[i]
    }
    newdistances
  E(igraph)$distances <- newdistances
  
  strengthsedgelist <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    strengthsedgelist[ind] <- strengths[i]
  }

  E(igraph)$strengths <- strengthsedgelist
  # Make vector with mean strengths vs distances 
  
  # values <- sort(unique(distances))
  # weightvsdis <- numeric(length(values))
  
  # weightvsdis <- numeric(length(distances))
  # 
  # for( i in 1:length(values)){
  #   k <- distances[i]
  #   kaas <- E(igraph)[E(igraph)$distances == k]
  #   c_k <- sum(E(igraph)$weights[kaas])/length(kaas)
  #   weightvsdis[i] <- c_k
  # } 
  # 
  # return(list(weights = weightvsdis, distances = values ))
  
  return(igraph)
}

cn.strengths <- function(igraph, data.igraph, k = NULL, criterion = NULL, perm = NULL){
  # igraph <- graphs_int_CN[[100]]
  # criterion <- NULL
  # k <- NULL
  # data.igraph<- TimeCoordsAnom_from_Grid_rms(tas_interim_10d, rms = TRUE)
  # perm <- permutations[[1]]
  require(data.table)  
  
  # Check if class is igraph or bn, add characteristics for plotting purpose  

    names <- c()
    names
    for(i in 1:648) names <- append(names,paste0("V",i))
    V(igraph)$name <- names
    nodenames <- V(igraph)$name


  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  adjmattri[lower.tri(adjmattri)] <- 0
  cormatrix <- abs(cor(data.igraph))
  cormatrix
  rownames(cormatrix) <- colnames(cormatrix) <-  names
  cormatrix[lower.tri(cormatrix)] <- 0
  
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  noedgeindices <-  which(adjmattri == 0, arr.ind = TRUE)
  cormatrix[noedgeindices] <- 0
  
  edgeindices2 <- which(cormatrix != 0, arr.ind = TRUE)
  all.equal(edgeindices, edgeindices2)
  
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  edgeindices[1,]
  cormatrix[edgeindices[,]]
  
  
  # identify which indices in igraph edgelist correspond to indices in adjmatrix
  # in order to add edgeattributes to igraph
  indsame <- c()
  nrow(as_edgelist(igraph))

  fromtoV[4,]
  as_edgelist(igraph)[7,]
  for (i in 1:nrow(fromtoV)){
    int <- intersect(which(fromtoV[i,1] == as_edgelist(igraph)[,1] ),
                     which(fromtoV[i,2] == as_edgelist(igraph)[,2] ))
    indsame[i] <- int
  }

  
  
  strengths <- cormatrix[edgeindices]
  strengths
  strengthsedgelist <- strengths[indsame]
  
  # i <- 1
  # fromtoV[59,]
  # as_edgelist(igraph)[1,]
  # newstrengths <- numeric(nrow(as_edgelist(igraph)))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   ind <- indsame[i]
  #   as_edgelist(igraph)[ind,]
  #   newstrengths[i] <- strengths[ind]
  # }
  
  E(igraph)$strenghts <- strengthsedgelist
  
  return(igraph)
}

cn.strengths_2 <- function(igraph, data.igraph, k = NULL, criterion = NULL, perm = NULL){
  # igraph <- graphs_int_CN[[100]]
  # criterion <- NULL
  # k <- NULL
  # data.igraph<- TimeCoordsAnom_from_Grid_rms(tas_interim_10d, rms = TRUE)
  # perm <- permutations[[1]]
  require(data.table)  
  
  # Check if class is igraph or bn, add characteristics for plotting purpose  
  
  names <- c()
  names
  for(i in 1:648) names <- append(names,paste0("V",i))
  V(igraph)$name <- names
  nodenames <- V(igraph)$name
  
  
  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  adjmattri[lower.tri(adjmattri)] <- 0
  cormatrix <- abs(cor(data.igraph))
  cormatrix
  rownames(cormatrix) <- colnames(cormatrix) <-  names
  cormatrix[lower.tri(cormatrix)] <- 0
  
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  noedgeindices <-  which(adjmattri == 0, arr.ind = TRUE)
  cormatrix[noedgeindices] <- 0
  
  edgeindices2 <- which(cormatrix != 0, arr.ind = TRUE)
  all.equal(edgeindices, edgeindices2)
  
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  edgeindices[1,]
  cormatrix[edgeindices[,]]
  
  
  # identify which indices in igraph edgelist correspond to indices in adjmatrix
  # in order to add edgeattributes to igraph
  indsame <- c()
  nrow(as_edgelist(igraph))
  i <- 3
  fromtoV[3,]
  as_edgelist(igraph)[7,]
  for (i in 1:nrow(fromtoV)){
    int <- intersect(which(as_edgelist(igraph)[i,1] == fromtoV[,1] ),
                     which(as_edgelist(igraph)[i,2] == fromtoV[,2] ))
    indsame[i] <- int
  }
  
  fromtoV[indsame,]
  all.equal(fromtoV[indsame,],as_edgelist(igraph), check.attributes = FALSE)
  
  strengths <- cormatrix[edgeindices]
  strengths
  strengthsedgelist <- strengths[indsame]
  
  # i <- 1
  # fromtoV[59,]
  # as_edgelist(igraph)[1,]
  # newstrengths <- numeric(nrow(as_edgelist(igraph)))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   ind <- indsame[i]
  #   as_edgelist(igraph)[ind,]
  #   newstrengths[i] <- strengths[ind]
  # }
  
  E(igraph)$strenghts <- strengthsedgelist
  
  return(igraph)
}

bn_to_igraph.strengths <- function(dag, data.dag, k = NULL, criterion = NULL, perm = NULL){
  # dag <- proefDAG
  # criterion <- NULL
  # k <- NULL
  # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
  # perm <- permutations[[1]]
  require(data.table)  

  arcsbn <- arcs(dag)
  linkstrength <- arc.strength(dag, as.data.frame(data.dag), criterion = criterion ,k  = k)
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))
  nodenames <- nodes(dag) 
  
  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  # adjmattri[lower.tri(adjmattri)] <- 0
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  edgeindices
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  # Make edgelist in igraph class method 2
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
  
  return(igraph)
}



igraph.distances <- function(igraph, data.igraph, k = NULL, criterion = NULL, perm = NULL){

  # criterion <- NULL
  # k <- NULL
  # data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
  # perm <- permutations[[1]]
  require(data.table)  

  nodenames <- V(igraph)$name
  # Check if class is igraph or bn, add characteristics for plotting purpose  
  if (is.null(V(igraph)$name)) {
    names <- c()
    names
    for(i in 1:length(V(igraph))) names <- append(names,paste0("V",i))
    V(igraph)$name <- names
    nodenames <- V(igraph)$name

    # as.directed(igraph, mode = "arbitrary")
  }
  # dir <- is.directed(igraph)
  # if(!dir){as.directed(igraph, mode = "arbitrary")}


  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  # if(!dir){adjmattri[upper.tri(adjmattri)] <- 0}
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)
  edgeindices
  as_edgelist(igraph)
  fromV <- nodenames[edgeindices[,1]]
  toV <- nodenames[edgeindices[,2]]
  fromtoV <- cbind(fromV,toV)
  
  # # Make edgelist in igraph class method 2
  # indsame <- c()
  # for (i in 1:nrow(linkstrength)){
  #   int <- intersect(which(as_edgelist(igraph)[,1] == arcsbn[i,1]),
  #                    which(as_edgelist(igraph)[,2] == arcsbn[i,2]))
  #   indsame[i] <- int
  # }
  # indsame
  
  # Make coordinates for every variable
  longitude <- attributes(data.igraph)$VertexCoords$x
  lattitude <- attributes(data.igraph)$VertexCoords$y
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
  nrow(arcdistances)
  
  # identify which indices in igraph edgelist correspond to indices in bn
  indsame <- c()
  nrow(as_edgelist(igraph))
  for (i in 1:nrow(as_edgelist(igraph))){
    int <- intersect(which(arcdistances[,1] == as_edgelist(igraph)[i,1]),
                     which(arcdistances[,2] == as_edgelist(igraph)[i,2]))
    indsame[i] <- int
  }
  
  # permutate distance vector to find belong distance vector for igraph object
  newdistances <- distances[indsame]
  
  # newdistances <- numeric(nrow(as_edgelist(igraph)))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   ind <- indsame[i]
  #   as_edgelist(igraph)[ind,]
  #   newdistances[ind] <- distances[i]
  # }
  # newdistances
  E(igraph)$distances <- newdistances
  # if(!dir){as.undirected(igraph)}
  return(igraph)
}
igraph.distances.2 <- function(igraph, data.igraph, k = NULL, criterion = NULL, perm = NULL){
  #gfG <- graph_from_Grid(grid = gridused, method = "pearson", th = 0.35)
  #igraph <- gfG$graph

  require(data.table)  
  
  nodenames <- V(igraph)$name
  # Check if class is igraph or bn, add characteristics for plotting purpose  
  if (is.null(V(igraph)$name)) {
    names <- c()
    names
    for(i in 1:length(V(igraph))) names <- append(names,paste0("V",i))
    V(igraph)$name <- names
    nodenames <- V(igraph)$name
    
    # as.directed(igraph, mode = "arbitrary")
  }
  # dir <- is.directed(igraph)
  # if(!dir){as.directed(igraph, mode = "arbitrary")}
  
  
  # Make edgelist in igraph class by edgeindices adjacancy matrix
  # for selecting edgeindices when applying Haversine
  adjmat <- as.matrix(as_adjacency_matrix(igraph))
  adjmattri <- adjmat
  adjmattri[upper.tri(adjmattri)] <- 0
  edgeindices <- which(adjmattri == 1, arr.ind = TRUE)

  fromV <- nodenames[edgeindices[,2]]
  toV <- nodenames[edgeindices[,1]]
  fromtoV <- cbind(fromV,toV)

  # Make coordinates for every variable
  longitude <- attributes(data.igraph)$VertexCoords$x
  lattitude <- attributes(data.igraph)$VertexCoords$y
  
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }
  
  # estimate distance of all edges in igraph-edgelist.
  distances <- c()
  for (i in 1:nrow(edgeindices)){
    x1Lat <- lattitude[edgeindices[i,2]]
    x1Lon <- longitude[edgeindices[i,2]]
    x2Lat <- lattitude[edgeindices[i,1]]
    x2Lon <- longitude[edgeindices[i,1]]
    
    disti <- haversine(x1Lat,x1Lon,x2Lat,x2Lon)
    distances[i] <- disti}
  
  # # Make dataframe with departing variable, end variable and distance
  # arcdistances <- cbind(unname(fromtoV), round(distances, digits = 3))
  # nrow(arcdistances)
  # 
  # # identify which indices in igraph edgelist correspond to indices in bn
  # indsame <- c()
  # nrow(as_edgelist(igraph))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   int <- intersect(which(arcdistances[,1] == as_edgelist(igraph)[i,1]),
  #                    which(arcdistances[,2] == as_edgelist(igraph)[i,2]))
  #   indsame[i] <- int
  # }
  # 
  # # permutate distance vector to find belong distance vector for igraph object
  # newdistances <- distances[indsame]
  
  # newdistances <- numeric(nrow(as_edgelist(igraph)))
  # for (i in 1:nrow(as_edgelist(igraph))){
  #   ind <- indsame[i]
  #   as_edgelist(igraph)[ind,]
  #   newdistances[ind] <- distances[i]
  # }
  # newdistances
  E(igraph)$distances <- distances
  # if(!dir){as.undirected(igraph)}
  return(igraph)
}

igraph.weights <- function(igraph, type = c("bn","cn"), fromdist = 2000, k = NULL, criterion = NULL, perm = NULL){
  # type = "bn"
  # igraph <- perm1strenghts[[1]]
  # fromdist <- 0
  strengths <- E(igraph)$strengths
  distances <- E(igraph)$distances
  
  if (type == "bn"){
  E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  } else if (type == "cn"){
  E(igraph)$weights <- (strengths - min(strengths))/(max(strengths) - min(strengths))  
  }
  
  larges <- which(E(igraph)$distances > fromdist)
  largedistances <- E(igraph)$distances[larges]
  largestrengths <- E(igraph)$strengths[larges]
  
  normalize <- max(largestrengths)-min(largestrengths)
  if (type == "bn"){
  largeweights <- (max(largestrengths) - largestrengths)/normalize
  } else if (type == "cn"){
  largeweights <- (largestrengths - min(largestrengths))/normalize
  }
  
  largeweights.edge <- numeric(length = length(E(igraph)))
  largeweights.edge[larges] <- largeweights
  
  E(igraph)$largeweights <- largeweights.edge
  
  return(igraph)
}
  
igraph.zoom <- function(igraph, type = c("bn","cn"), fromdist = 2000, perm = NULL){
  #  type = "cn"
  # igraph <- cn.str.dist
  strengths <- E(igraph)$strengths
  distances <- E(igraph)$distances
  
  if (type == "bn"){
    E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  } else if (type == "cn"){
    E(igraph)$weights <- (strengths - min(strengths))/(max(strengths) - min(strengths))  
  }
  
  larges <- which(E(igraph)$distances > fromdist)
  largedistances <- E(igraph)$distances[larges]
  largestrengths <- E(igraph)$strengths[larges]
  
  normalize <- max(largestrengths)-min(largestrengths)
  if (type == "bn"){
    largeweights <- (max(largestrengths) - largestrengths)/normalize
  } else if (type == "cn"){
    largeweights <- (largestrengths - min(largestrengths))/normalize
  }
  
  largeweights.edge <- numeric(length = length(E(igraph)))
  largeweights.edge[larges] <- largeweights
  
  E(igraph)$largeweights <- largeweights.edge
  
  return(igraph)
}

plot.spatial.igraph <- function(network, data.network, th = 0.2, by.dist = NULL, perm = NULL, title = "NA", remove = TRUE){
  # data.network <- data
  # remove <- TRUE
  # perm <- permutations[[1]]
  # require(data.table)
  # # Check if class is igraph or bn, add characteristics for plotting purpose
  # th <- 0.2
  # by.dist <- 3000
  # network <- cn.str.dist.weight
  # distances <- E(network)$distances
  # largeweights <- E(network)$largeweights
 
  # Make coordinates for every variable
  longitude <- attributes(data.network)$VertexCoords$x
  lattitude <- attributes(data.network)$VertexCoords$y
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }
  
  # Remove limit neighbours
  names <- c()
  for(i in 1:ncol(data.network)) {names <- append(names,paste0("V",i))}
  nvertpts <- sqrt(ncol(data.network)/2)
  limitsl <- names[1:nvertpts]
  limitsl <- t(limitsl)
  limitsr <- names[(ncol(data.network)-nvertpts+1):ncol(data.network)]
  limitsr <- t(limitsr)
  
  black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
  black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
  blacks <- rbind(black,black2)
  blacks
  
  m1 <- as_edgelist(network)
  m2 <- blacks
  colnames(m2) <- paste0("V", seq(len=ncol(m2)))
  DT1 <- data.table(m1)
  DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
  setnames(DT2, c(head(names(DT2), -1L), "found"))
  a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
  a <-as.matrix(a)
  
  # Color edges:
  # Large edges, small edges, and border edges
  c_scale <- rev(topo.colors(5))
  # quantiles<- quantile(E(igraph)[E(igraph)$distances > minimdist]$weights)
  
 
  if (is.null(by.dist)){
  E(network)[largeweights <= th]$color <- NA
  E(network)[largeweights > th]$color <- c_scale[5]
  # E(network)[largeweights <= th]$width <- 0
  # E(network)[largeweights > th]$width<- 1
  } else { 
  E(network)[E(network)$largeweights <= th]$color <- NA
  E(network)[which(E(network)$distances >  by.dist & E(network)$largeweights > th)]$color <- c_scale[5]
  E(network)[which(E(network)$distances <= by.dist & E(network)$largeweights > th)]$color <- c_scale[1]
  #E(igraph)[largeweights <= th]$width <- 0
  #E(igraph)[distances <= by.dist & largeweights > th]$width <- 0.5
  #E(igraph)[distances > by.dist & largeweights > th]$width <- 1
  }

  if(remove == TRUE) {
    E(network)[which(a[,1] == 1, arr.ind = TRUE)]$color <- NA}


  
  # draw the graph
  x <- attr(data.network, "Xcoords", exact = FALSE)
  y <- attr(data.network, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  if (!is.null(perm)){
    points <- points[perm,]
  }
  
  # load("/home/catharina/Documents/lisworld.rda")
  # load("/Volumes/ubuntu/Documents/lisworld.rda")
  # load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
  raster::plot(lisworld)
  plot.igraph(network,
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              # edge.width = E(igraph)$width,
              lty = "dots",
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
}

plotS.spatial.igraph <- function(network, data.network, type = c("cn","bn"), smallcol = "purple",th = 0.2, th.type = c("zoom","zoomweighted"), from.dist = 0, by.dist = NULL, perm = NULL, title = "NA", remove = TRUE, shift = FALSE, curvature = FALSE){
  # data.network <- data
  # remove <- TRUE
  # perm <- permutations[[1]]
  # require(data.table)
  # # Check if class is igraph or bn, add characteristics for plotting purpose
  # th <- 0.2
  # by.dist <- 3000
  # network <- cn.str.dist.weight
  # distances <- E(network)$distances
  # largeweights <- E(network)$largeweights
  
  # Make coordinates for every variable
  longitude <- attributes(data.network)$VertexCoords$x
  lattitude <- attributes(data.network)$VertexCoords$y
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }
  
  # Remove limit neighbours
  names <- c()
  for(i in 1:ncol(data.network)) {names <- append(names,paste0("V",i))}
  nvertpts <- sqrt(ncol(data.network)/2)
  limitsl <- names[1:nvertpts]
  limitsl <- t(limitsl)
  limitsr <- names[(ncol(data.network)-nvertpts+1):ncol(data.network)]
  limitsr <- t(limitsr)
  
  black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
  black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
  blacks <- rbind(black,black2)
  blacks
  
  m1 <- as_edgelist(network)
  m2 <- blacks
  colnames(m2) <- paste0("V", seq(len=ncol(m2)))
  DT1 <- data.table(m1)
  DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
  setnames(DT2, c(head(names(DT2), -1L), "found"))
  a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
  a <-as.matrix(a)
  
  c_scale <- rev(topo.colors(5))
  # Color edges:
  # Large edges, small edges, and border edges
  
  if(th.type == "zoomweighted"){
    if (is.null(by.dist)){
      E(network)[largeweights <= th]$color <- NA
      E(network)[largeweights > th]$color <- "purple"
      # E(network)[largeweights <= th]$width <- 0
      # E(network)[largeweights > th]$width<- 1
    } else { 
      E(network)[E(network)$largeweights <= th]$color <- NA
      E(network)[which(E(network)$distances >  by.dist & E(network)$largeweights > th)]$color <- "orange"
      E(network)[which(E(network)$distances <= by.dist & E(network)$largeweights > th)]$color <- smallcol
      #E(igraph)[largeweights <= th]$width <- 0
      #E(igraph)[distances <= by.dist & largeweights > th]$width <- 0.5
      #E(igraph)[distances > by.dist & largeweights > th]$width <- 1
    }
  }
  
  if(th.type == "zoom" & type == "cn"){
    if (is.null(by.dist)){
      E(network)[strengths <= th]$color <- NA
      E(network)[strengths > th]$color <- "purple"
      E(network)[distances <= from.dist]$color <- NA
      # E(network)[largeweights <= th]$width <- 0
      # E(network)[largeweights > th]$width<- 1
    } else { 
      E(network)[E(network)$strengths <= th]$color <- NA
      E(network)[which(E(network)$distances >  by.dist & E(network)$strengths > th)]$color <- "purple"
      E(network)[which(E(network)$distances <= by.dist & E(network)$strengths > th)]$color <- "orange"
      E(network)[E(network)$distances <= from.dist]$color <- NA
      #E(igraph)[largeweights <= th]$width <- 0
      #E(igraph)[distances <= by.dist & largeweights > th]$width <- 0.5
      #E(igraph)[distances > by.dist & largeweights > th]$width <- 1
    }
  }
  
  if(th.type == "zoom" & type == "bn"){
    if (is.null(by.dist)){
      E(network)[weights <= th]$color <- NA
      E(network)[weights > th]$color <- c_scale[5]
      E(network)[distances <= from.dist]$color <- NA
      # E(network)[largeweights <= th]$width <- 0
      # E(network)[largeweights > th]$width<- 1
    } else { 
      E(network)[E(network)$weights <= th]$color <- NA
      E(network)[which(E(network)$distances >  by.dist & E(network)$weights > th)]$color <- "purple"
      E(network)[which(E(network)$distances <= by.dist & E(network)$weights > th)]$color <- "orange"
      E(network)[E(network)$distances <= from.dist]$color <- NA
      #E(igraph)[largeweights <= th]$width <- 0
      #E(igraph)[distances <= by.dist & largeweights > th]$width <- 0.5
      #E(igraph)[distances > by.dist & largeweights > th]$width <- 1
    }
  }
  
  
  if(remove == TRUE) {
    E(network)[which(a[,1] == 1, arr.ind = TRUE)]$color <- NA}
  
  # draw the graph
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
  # if(shift == TRUE){raster::plot(m)}
  if(shift == TRUE){map('world2',interior = FALSE,resolution = 0)}
  else {raster::plot(lisworld)}
  
  plot.igraph(network,
              vertex.size = 100,
              vertex.color = NA,
              vertex.label = NA,
              edge.curved = curvature,
              edge.arrow.size = 0.2,
              # edge.width = E(igraph)$width,
              lty = "dots",
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  title(main = title)
}


plotS.cluster.igraph <- function(network, data.network, type = c("cn","bn"), th = 0.2, th.type = c("zoom","zoomweighted"), from.dist = 0, by.dist = NULL, perm = NULL, title = "NA", remove = TRUE, shift = FALSE){
  # data.network <- data
  # remove <- TRUE
  # perm <- permutations[[1]]
  # require(data.table)
  # # Check if class is igraph or bn, add characteristics for plotting purpose
  # th <- 0.2
  # by.dist <- 3000
  # network <- cn.str.dist.weight
  # distances <- E(network)$distances
  # largeweights <- E(network)$largeweights
  
  # Make coordinates for every variable
  longitude <- attributes(data.network)$VertexCoords$x
  lattitude <- attributes(data.network)$VertexCoords$y
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }
  
  # Remove limit neighbours
  names <- c()
  for(i in 1:ncol(data.network)) {names <- append(names,paste0("V",i))}
  nvertpts <- sqrt(ncol(data.network)/2)
  limitsl <- names[1:nvertpts]
  limitsl <- t(limitsl)
  limitsr <- names[(ncol(data.network)-nvertpts+1):ncol(data.network)]
  limitsr <- t(limitsr)
  
  black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
  black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
  blacks <- rbind(black,black2)
  blacks
  
  m1 <- as_edgelist(network)
  m2 <- blacks
  colnames(m2) <- paste0("V", seq(len=ncol(m2)))
  DT1 <- data.table(m1)
  DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
  setnames(DT2, c(head(names(DT2), -1L), "found"))
  a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
  a <-as.matrix(a)
  
  # Color edges:
  # Large edges, small edges, and border edges
  
  
  
  
  if(remove == TRUE) {
    E(network)[which(a[,1] == 1, arr.ind = TRUE)]$color <- NA}
  
  # draw the graph
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
  if(shift == TRUE){raster::plot(m)}
  else {raster::plot(lisworld)}
  
  plot.igraph(network,
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              # edge.width = E(igraph)$width,
              lty = "dots",
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
}


plot.vstruct.bn <- function(dseplist, data.network, perm = NULL, title = "NA", remove = FALSE, shift = FALSE){
   # data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
   # bn <- perm3sort$hc3_500_600i
  # remove <- TRUE
   # perm <- permutations[[1]]
   # shift <- TRUE
  # require(data.table)
  # # Check if class is igraph or bn, add characteristics for plotting purpose
  # th <- 0.2
  # by.dist <- 3000
  # network <- cn.str.dist.weight
  # distances <- E(network)$distances
  # largeweights <- E(network)$largeweights
  # dseplist <- dseppingvstructs13[[5]]
  names <- c()
  for(i in 1:ncol(data.network)) {names <- append(names,paste0("V",i))}
  bn <- empty.graph(names)
  dsep1 <- dseplist[,1:2]
  dsep2 <- dseplist[,c(3,2)]
  dsep <- rbind(dsep1,dsep2)
  colnames(dsep) <- c("from","to")
  arcs(bn) <- dsep
  network <- igraph.from.graphNEL(as.graphNEL(bn))
  
  # Make coordinates for every variable
  longitude <- attributes(data.network)$VertexCoords$x
  lattitude <- attributes(data.network)$VertexCoords$y
  if (!is.null(perm)){
    longitude <- longitude[perm]
    lattitude <- lattitude[perm]
  }
  
  # Remove limit neighbours
  nvertpts <- sqrt(ncol(data.network)/2)
  limitsl <- names[1:nvertpts]
  limitsl <- t(limitsl)
  limitsr <- names[(ncol(data.network)-nvertpts+1):ncol(data.network)]
  limitsr <- t(limitsr)
  
  black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
  black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
  blacks <- rbind(black,black2)
  blacks
  
  m1 <- as_edgelist(network)
  m2 <- blacks
  colnames(m2) <- paste0("V", seq(len=ncol(m2)))
  DT1 <- data.table(m1)
  DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
  setnames(DT2, c(head(names(DT2), -1L), "found"))
  a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
  a <-as.matrix(a)
  

  
  
  if(remove == TRUE) {
    E(network)[which(a[,1] == 1, arr.ind = TRUE)]$color <- NA}
  
  # draw the graph
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
  if(shift == TRUE){plot(crop(l1[[2]],m))} else {raster::plot(lisworld)}
  plot.igraph(network,
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              # edge.width = E(igraph)$width,
              lty = "dots",
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
}#########################################################################
# CMIP utils
#########################################################################
IPweight <- function(bdists,ref.perf,sD,sS){
  refID <- which(rownames(bdists) == ref.perf)
  D <- bdists[ref.perf,-refID]
  S <- bdists[-refID,-refID]
  W <- numeric(length=length(D))
  for(i in 1:length(names(D))){
    name.i <-names(D)[i]
    nameID.S <- which(rownames(S) == name.i)
    nameID.D <- which(names(D) == name.i)
    Di <- D[nameID.D]
    Sij <- S[nameID.S,][-nameID.S]
    sumS <- sum(exp(-Sij^2/sS^2))
    W[i] <- exp(-Di^2/sD^2)/(1+sumS)
  }
  W <- W/sum(W)
  names(W) <- names(D)
  return(W)
}