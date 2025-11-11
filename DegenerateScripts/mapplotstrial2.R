library(maps)
library(geosphere)
par(mfrow = c(2,2), mar=c(0,0,0,0))

map("usa", col="tomato",  border="gray10", fill=TRUE, bg="gray30")
map("state", col="orange",  border="gray10", fill=TRUE, bg="gray30")
map("county", col="palegreen",  border="gray10", fill=TRUE, bg="gray30")
map("world", col="skyblue",  border="gray10", fill=TRUE, bg="gray30")
map("world")
SpatialLines2map(m$)

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
]

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

  
  col.1 <- adjustcolor("orange red", alpha=0.4)
  col.2 <- adjustcolor("orange", alpha=0.4)
  edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
  edge.col <- edge.pal(100)
  
  
  
  map(result, col="skyblue",  border="gray10", fill=TRUE, bg="gray30")
  data.network <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
  x <- attr(data.network, "Xcoords", exact = FALSE)
  y <- attr(data.network, "Ycoords", exact = FALSE)
  x <- (x+360)%%360
  points <- expand.grid(y, x)[2:1]
  
  
  if (!is.null(perm)){
    points <- points[perm,]
  }
 
 
  points(x=points$Var2, y=points$Var1, pch=19, col="orange")

  graph <- listcnweights$`2260`
  graph <- perm1weights$hc1_1700_1800i
edgelist <- as_edgelist(graph) 
row.names(points) <- names(V(graph))
gridpoints <- names(V(graph))

dists <- edge.attributes(graph)$distances

largegraph <- subgraph.edges(graph,E(graph)[edge.attributes(graph)$distances>10000], delete.vertices = FALSE)
edgelist <- as_edgelist(largegraph)
i <- 1
  for(i in 1:nrow(edgelist))  {
    node1 <- gridpoints[gridpoints == edgelist[i,1]]
    node2 <- gridpoints[gridpoints == edgelist[i,2]]
    
    arc <- gcIntermediate( c(points[node1,]$Var2, points[node1,]$Var1), 
                           c(points[node2,]$Var2, points[node2,]$Var1),
                           n=1000, addStartEnd=TRUE )
    edge.ind <- ceiling(edge_attr(graph,"distances",E(graph)[i])*100 / max(dists))
    
    lines(arc, col=edge.col[edge.ind], lwd=edge.ind/30)

  }
dev.off()
warnings()








dummy_graph <-make_full_graph(5) * 2 +edge(1, 6)
merge_matrix <-rbind(c(1, 2),c(3, 4),c(6, 7),c(9, 10),c(12, 5),c(8, 14),c(11, 15),c(13, 16),c(17, 18))
dendrogram <-make_clusters(dummy_graph, merges = merge_matrix)
dendrogram
