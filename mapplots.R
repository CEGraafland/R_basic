library(maps)
library(maptools)
library(geosphere)
library(mapproj)
library(ggplot2)
library(broom)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/world_pol.rda")

par(mfrow = c(2,2), mar=c(0,0,0,0))

map("usa", col="tomato",  border="gray10", fill=TRUE, bg="gray30")
map("state", col="orange",  border="gray10", fill=TRUE, bg="gray30")
map("county", col="palegreen",  border="gray10", fill=TRUE, bg="gray30")
map("world", col="skyblue",  border="gray10", fill=TRUE, bg="gray30")
map("world")

plot(lisworld)
m2 <- crop(l1[[2]],m)
dim(l1[[2]]@lines)
class(m2)
plot(m2)
SpatialLinesDataFrame(lisworld,data)
SpatialLines2map(l1[[2]])
l1[[2]]@lines
dim(l1[[2]]) <- c(128,1)
l1[[2]][3][1]
plot(l1[[1]][1])
SpatialLines2map


map2SpatialPolygons()
mapworld2 <- map("world2", fill = TRUE)
mapworld2$range
IDs <- mapworld2$names
length(mapworld2)
class(mapworld2)
str(mapworld2)
map2SpatialPolygons(mapworld2, IDs = IDs)
IDs <- sapply(strsplit(mapworld2$names, ":"), function(x) x[1])
mapworld <- SpatialPolygons2map(world_pol)
class(mapworld)

database <- crop(l1[[2]],lisworld)
database <- crop(l1[[2]],m)
namefield <- NULL
# function (database, namefield = NULL) 
# {
if (!inherits(database, "SpatialLines")) 
  stop("database must be a SpatialLines[DataFrame] object.")
line.names <- NULL
if (inherits(database, "SpatialLinesDataFrame") & !is.null(namefield)) {
  namcol <- lapply(namefield, function(x) which(tolower(names(database)) == 
                                                  tolower(x)))
  if (any(lapply(namcol, length) != 1)) {
    zz <- which(lapply(namcol, length) != 1)
    warning(paste0("database does not (uniquely) contain the field '", 
                   namefield[zz], "'."))
  }
  else {
    zz <- as.data.frame(lapply(database@data[unlist(namcol)], 
                               as.character), stringsAsFactors = FALSE)
    line.names <- vapply(1:dim(zz)[1], function(x) paste(zz[x, 
                                                            ], collapse = ":"), FUN.VALUE = "a")
  }
}
if (is.null(line.names)) 
  line.names <- unlist(lapply(database@lines, function(x) x@ID))
nlines <- length(line.names)
nseg <- vapply(1:nlines, FUN = function(i) length(database@lines[[i]]@Lines), 
               FUN.VALUE = 1)
line.names <- unlist(lapply(1:dim(database)[1], function(i) {
  if (nseg[i] == 1) 
    line.names[i]
  else paste(line.names[i], 1:nseg[i], sep = ":")
}))
allpoly <- lapply(database@lines, function(x) lapply(x@Lines, 
                                                     function(y) y@coords))
mymap <- do.call(rbind, lapply(do.call(c, allpoly), function(x) rbind(c(NA, 
                                                                        NA), x)))[-1, ]
result <- list(x = mymap[, 1], y = mymap[, 2], names = line.names, 
               range = c(range(mymap[, 1], na.rm = TRUE), range(mymap[, 
                                                                      2], na.rm = TRUE)))
class(result) <- "map"
result
plot(result)

col.1 <- adjustcolor("orange red", alpha=0.4)
col.2 <- adjustcolor("orange", alpha=0.4)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)


map(mapworld,col="skyblue",  border="gray10", fill=TRUE, bg="gray30")
map(result, col="skyblue",  border="gray10", fill=TRUE, bg="gray30")
data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)
x <- (x+360)%%360
p <- expand.grid(y, x)[2:1]
rm(points)
x[x < 0] = x[x<0] + 360


if (!is.null(perm)){
  points <- points[perm,]
}

mapply(points, x=p$Var2, y=p$Var1, pch=19, col="orange")
points(x=points$Var2, y=points$Var1, pch=19, col="orange")

graph <- listcnweights$`2260`
graph <- perm1weights$hc1_1700_1800i
edgelist <- as_edgelist(graph) 
row.names(p) <- names(V(graph))
gridpoints <- names(V(graph))

dists <- edge.attributes(graph)$distances

largegraph <- subgraph.edges(graph,E(graph)[edge.attributes(graph)$distances>10000], delete.vertices = FALSE)
edgelist <- as_edgelist(largegraph)
i <- 1
arcslist <- list()
i <- 1
for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]

    node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=10, addStartEnd=TRUE, breakAtDateLine = FALSE)
  # edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
  
  negs <- which(arc[,'lon']<0)

  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
   # lines(arc, col=edge.col[edge.ind], lwd=edge.ind)
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}
naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes
dev.off()
warnings()


install.packages("rgdal")
install.packages("mapproj")
library(rgdal)

coastline_pol <- readOGR("path", layer = "filename without extension")

map('world', fill = TRUE, col = 1:10, wrap=c(-180,180) )
map('world', fill = TRUE, col = 1:10, wrap=c(-180,180) )

map('coastline')

world2MapEnv
map(mapworld, projection="rectangular", parameter=0, 
    orientation=c(90,0,0), wrap = TRUE, fill=TRUE, resolution=0,col=0)
map('world2', interior = FALSE,projection="rectangular", parameter=0, 
    orientation=c(90,0,0))
plot(world2MapEnv)
plot(m)
map.grid(c(0,360,-90,90))

plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

plot.map("world", center=180, col="white",bg="gray",
         fill=TRUE,ylim=c(-90,90),mar=c(0,0,0,0))



meteodag <- perm1weights$hc1_1500_1600i

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



library(ggmap)
library(broom)

# 
# names(arcslist[[1]])
# lapply(arcslist, function(x) {
#   cbind.data.frame(route=x, lon = x[,'lon'], lat = x[,'lat'],stringsAsFactors=FALSE)})
# lapply(arcslist, function(x) {
#   cbind.data.frame(route=x,stringsAsFactors=FALSE)})
# 
# names(arcslist) <- as.character(1:length(arcslist))
# lapply(arcslist, function(x) {})
# 
# dflist <- list()
# for(i in 1:length(arcslist)){
# eerste <- cbind.data.frame( x = arcslist[[i]])
# tweede <- eerste
# tweede$naam <- attributes(arcslist)$names[i]
# dflist[[i]] <- tweede
# }
    # listpaths <- lapply(arcslist, function(x) {cbind.data.frame(route=x,stringsAsFactors=FALSE)})




##################################################################
# Plot maps with ggplot and lines as got by gcIntermediate
# Adjust to CN and BN in name of plot
##################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda") # <- no es gridGraph! list cn weights
edgelistsCN <- lapply(listcnweights, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(listcnweights) <- as.character(numberofedgesCN)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/permweights1_10/perm1weights.rda")
mapworld2 <- map("world2", fill = TRUE)
dev.off()
data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
x <- attr(data.network, "Xcoords", exact = FALSE)
y <- attr(data.network, "Ycoords", exact = FALSE)
# x <- (x+360)%%360
p <- expand.grid(y, x)[2:1]


graphCN <- listcnweights$`2260`
graphCN <- listcnweights$`8066`
graphBN <- perm1weights$hc1_1700_1800i
graph <- graphCN
edgelist <- as_edgelist(graph) 
nrow(edgelist)
row.names(p) <- names(V(graph))
gridpoints <- names(V(graph))

dists <- edge.attributes(graphBN)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

arcslist <- list()
for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  
  edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)

  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  arc <- cbind(arc,edge.ind)
  # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}

naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes


world_pol_df <- tidy(mapworld2, IDs = "region")
    
do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
      cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
    })) -> all_paths
    


bla <- ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  ggtitle(paste0("BN: ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

bla

#####################################################################
# high betweenness maps plot. 
#####################################################################
mapworld2 <- map("world2", fill = TRUE)
dev.off()

graphCN <- listcnweights$`1783`
graphCN <- listcnweights$`2260`
graphCN <- listcnweights$`8066`
graphBN <- perm1weights$hc1_1700_1800i
graph <- graphBN
edgelist <- as_edgelist(graph) 
nrow(edgelist)
row.names(p) <- names(V(graph))
gridpoints <- names(V(graph))

invweighted <- TRUE
weighted <- FALSE
strengths <- FALSE
weights <- TRUE


if(invweighted == TRUE){
  if(strengths == TRUE){
    E(graph)$betwe <- edge.betweenness(graphCN,weights = 1/E(graphCN)$strengths)
  } else if (weights == TRUE){
    E(graph)$betwe <- edge.betweenness(graphBN,weights = 1/E(graphBN)$weights)}
} else if (weighted == TRUE){
  if(strengths == TRUE){
    E(graph)$betwe <- edge.betweenness(graphCN,weights = E(graphCN)$strengths)
  } else if (weights == TRUE){
    E(graph)$betwe <- edge.betweenness(graphBN,weights = E(graphBN)$weights)}
} else {E(graph)$betwe <- edge.betweenness(graph)}




dists <- edge.attributes(graphBN)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

col.1 <- adjustcolor("light blue", alpha=0.2)
col.2 <- adjustcolor("red", alpha=1)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

arcslist <- list()
for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  
  edge.ind <- ceiling(edge_attr(graph,"betwe",E(graph)[i])*100 / max(edge_attr(graph,"betwe",E(graph))))
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)
  
  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  arc <- cbind(arc,edge.ind)
  # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}
naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes


world_pol_df <- tidy(mapworld2, IDs = "region")

do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
  cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
})) -> all_paths



bla <- ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
  scale_color_gradient(low = col.1, high = col.2, guide = guide_colorbar(title = "Edge Betweenness")) +
  ggtitle(paste0("BN: ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

bla

##################################################################
# Plot maps with ggplot and lines as got by gcIntermediate
# Adjust to CN and BN in name of plot
##################################################################
mapworld2 <- map("world2", fill = TRUE)
dev.off()

graphCN <- listcnweights$`2260`
graphCN <- listcnweights$`8066`
graphBN <- perm1weights$hc1_1700_1800i
graph <- graphCN
edgelist <- as_edgelist(graph) 
nrow(edgelist)
row.names(p) <- names(V(graph))
gridpoints <- names(V(graph))

dists <- edge.attributes(graphCN)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

col.1 <- adjustcolor("light blue", alpha=0.8)
col.2 <- adjustcolor("navy blue", alpha=0.8)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

arcslist <- list()
for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  
  edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)
  
  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  arc <- cbind(arc,edge.ind)
  # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}
naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes


world_pol_df <- tidy(mapworld2, IDs = "region")

do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
  cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
})) -> all_paths



bla <- ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
  scale_color_gradient(low = col.1, high = col.2,limits = c(0,100), breaks = c(0,25,50,75,100), guide = guide_colorbar(title = "Distance"), labels = distsc) +
  ggtitle(paste0("BN: ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

bla

#####################################################################
# removed edges maps plot. 
#####################################################################
###############################################################
# 1 weighted
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_BN_1795_2592_7862.rda")
# 2 unweighted
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_BN_unw_1795_2592_7862.rda")
# 3 inverse weighted 
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_BN_invw_1795_2592_7862.rda")
# 4
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_CN_1783_2260_8066_10017.rda")
# 5
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_CN_unw_1783_2260_8066_10017.rda")
# 6
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/large_communities/commu_CN_invw_1783_2260_8066_10017.rda")


mapworld2 <- map("world2", fill = TRUE)
dev.off()

removed <- commu_betsbn_unw[[1]]$removed.edges
removed <- commu_betscn_unw[[2]]$removed.edges
removed <- commu_betscn_unw[[3]]$removed.edges

graphCN <- listcnweights$`1783`
graphCN <- listcnweights$`2260`
graphCN <- listcnweights$`8066`
graphBN <- perm1weights$hc1_1700_1800i
graph <- graphBN
edgelist <- as_edgelist(graph) 
nrow(edgelist)
row.names(p) <- names(V(graph))
gridpoints <- names(V(graph))
invweighted <- FALSE
weighted <- FALSE
strengths <- FALSE
weights <- TRUE


if(invweighted == TRUE){
  if(strengths == TRUE){
    E(graph)$betwe <- edge.betweenness(graphCN,weights = 1/E(graphCN)$strengths)
  } else if (weights == TRUE){
    E(graph)$betwe <- edge.betweenness(graphBN,weights = 1/E(graphBN)$weights)}
} else if (weighted == TRUE){
  if(strengths == TRUE){
    E(graph)$betwe <- edge.betweenness(graphCN,weights = E(graphCN)$strengths)
  } else if (weights == TRUE){
    E(graph)$betwe <- edge.betweenness(graphBN,weights = E(graphBN)$weights)}
} else {E(graph)$betwe <- edge.betweenness(graph)}


dists <- edge.attributes(graphBN)$distances
distsv <- round(c(0,max(dists)/4,max(dists)/2,max(dists)/4*3, max(dists)))
distsc <- as.character(distsv)

col.1 <- adjustcolor("yellow", alpha=0.2)
col.2 <- adjustcolor("red", alpha = 0.8)
col.3 <- adjustcolor("red", alpha=1)
edge.pal <- colorRampPalette(c(col.1, col.2,col.3), alpha = TRUE)
edge.col <- edge.pal(100)

arcslist <- list()
i <- 1
for(i in 1:nrow(edgelist))  {
  node1 <- gridpoints[gridpoints == edgelist[i,1]]
  
  node2 <- gridpoints[gridpoints == edgelist[i,2]]
  
  arc <- gcIntermediate(c(p[node1,]$Var2, p[node1,]$Var1), 
                        c(p[node2,]$Var2, p[node2,]$Var1),
                        n=100, addStartEnd=TRUE, breakAtDateLine = FALSE)
  
  edge.ind <- which(edgelist[i,1] == edgelist[removed,][,1] & edgelist[i,2] == edgelist[removed,][,2])
  if(length(edge.ind)==0){edge.ind < 1} else {edge.ind <- edge.ind*100 / length(removed)}
  col <- as.character(edge.col[edge.ind])
  negs <- which(arc[,'lon']<0)
  
  if(!is.null(negs)){arc[negs] <- (arc[negs]+360)%%360}
  arc <- cbind(arc,edge.ind)
  # lines(arc, col=edge.col[edge.ind], lwd=edge.ind) interesting for edgewidth!
  # points(x = arc[1,1], y = arc[1,2],col = "pink", pch = 19)
  arcslist[[i]]<- arc
  
}
naampjes <- character()
for(i in 1:length(arcslist)){naampjes[i] <- paste0("arc",i)}
names(arcslist) <- naampjes


world_pol_df <- tidy(mapworld2, IDs = "region")

do.call(rbind.data.frame, lapply(names(arcslist), function(x) {
  cbind.data.frame(route=x, arcslist[[x]], stringsAsFactors=FALSE)
})) -> all_paths



bla <- ggplot() +
  geom_polygon(data = world_pol_df, mapping = aes(x = long, y = lat, group = group)) +
  geom_line(data=all_paths, aes(x=lon, y=lat, group=route, col = edge.ind)) +
  scale_color_gradient(low = col.1, high = col.2, guide = guide_colorbar(title = "Removing order")) +
  ggtitle(paste0("BN: ",nrow(edgelist))) +
  theme(plot.title = element_text(hjust = 0.5))

bla





