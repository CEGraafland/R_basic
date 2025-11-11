library(igraph)
library(transformeR)
library(devtools)
install_github('SantanderMetGroup/mopa')
library(mopa)
data(wrld)
library(sp)

load(file = "/Users/lisettegraafland/Documents/R_practice/Data/ncep/tas_ncep.rda")
str(tas_ncep)
tas_ncep.clim <- climatology(tas_ncep)
plotClimatology(tas_ncep.clim, backdrop.theme = "countries")


#interpolation 5d met afgesneden polen
lon <- tas_ncep$xyCoords$x
lat <- tas_ncep$xyCoords$y

#After interpGrid in Console: (see ncepTASanalyse)
load(file = "/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_5d2.rda")

tas_ncep_5d2.clim <- climatology(tas_ncep_5d2)
plotClimatology(tas_ncep_5d2.clim, backdrop.theme = "countries")

#adjacency matrix function
adj.matrix.function5 <- function(eengrid){
  months <- list()
  for(i in 1:12){
    sea <- subsetGrid(eengrid, season = i)
    months[[i]] <- climatology(sea)
  }
  
  grid2 <- redim(bindGrid.time(months), drop = T)
  
  matref <- array3Dto2Dmat(grid2$Data)
  mat <- array3Dto2Dmat(eengrid$Data)
  ind <- rep(1:12,30)
  ind
  
  outmat <- mat 
  
  for (i in 1:12){
    month <- which(ind == i)
    outmat[month,] <- mat[month,] - matref[i,]
  }
  
  
  cor.matrix <- cor(outmat,method="spearman")
  
  adj.matrix5 <- cor.matrix
  diag(adj.matrix5) <- 0
  adj.matrix5[adj.matrix5 <= .99 ] <- 0
  adj.matrix5[adj.matrix5 > .99] <- 1

  return(adj.matrix5)}

bewaar <- adje
adje <- adj.matrix.function5(tas_ncep_5d2)
identical(adj.matrix.function2(tas_ncep_5d2),adj.matrix.function5(tas_ncep_5d2))
identical(adje,adj.matrix.function2(tas_ncep_5d2)) #FALSE!! 


tas5d2.adj.matrix2 <- adj.matrix.function2(tas_ncep_5d2)

#Find NA values
# in adjacency matrix:
which(is.na(tas5d2.adj.matrix2), arr.ind=TRUE) #hence none

#Graph
tas5d2.graph2<- graph_from_adjacency_matrix(tas5d2.adj.matrix2, mode = "undirected")
adje.graph <- graph_from_adjacency_matrix(adje, mode = "undirected")
# edge density
edge_density(tas5d2.graph2, loops = TRUE)
edge_density(adje.graph, loops = TRUE)

# degree distribution
degree_distribution(tas5d2.graph2)
plot(degree_distribution(tas5d2.graph2))
     
#betweenness
betw.tas5d2.graph2 <- betweenness(tas5d2.graph2)
betw.adje.graph <- betweenness(adje.graph)

log.betw <- log(betw.tas5d2.graph2 + 1)
log.betw.adje <- log(betw.adje.graph)
#betweenness edges
edge.betw.tas5d2.graph2 <- edge.betweenness(tas5d2.graph2)


betw.sp <- SpatialPixelsDataFrame(points,data.frame(betw.matrix[,3]))
spplot(betw.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

?SpatialPoints
mat2 <- array3Dto2Dmat(tas_ncep_5d2$Data)
mat2b <- grid3Dto2Dmatrix(tas_ncep_5d2)
str(mat2b)

mat2[1:360,1:3]
?array3Dto2Dmat
points2 <- expand.grid(tas_ncep_5d2$xyCoords$y, tas_ncep_5d2$xyCoords$x)[2:1]
betw2.sp <- SpatialPixelsDataFrame(points2, data.frame(betw.tas5d2.graph2))
betw.adje.sp <- SpatialPixelsDataFrame(points2, data.frame(betw.adje.graph))
spplot(betw2.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))
spplot(betw.adje.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

log.betw.sp <- SpatialPixelsDataFrame(points2, data.frame(log.betw))
spplot(log.betw.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

log.betw.adje.sp <- SpatialPixelsDataFrame(points2, data.frame(log.betw.adje))
spplot(log.betw.adje.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

betw.adje.log.plot <- spplot(log.betw.adje.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))
save(betw.adje.log.plot, file = "/Users/lisettegraafland/Documents/R_practice/Data/betw_adje_log_plot.rda")
betw.adje.log.plot

adj.matrix.function6 <- function(eengrid){
  months <- list()
  for(i in 1:12){
    sea <- subsetGrid(eengrid, season = i)
    months[[i]] <- climatology(sea)
  }
  
  grid2 <- redim(bindGrid.time(months), drop = T)
  
  matref <- array3Dto2Dmat(grid2$Data)
  mat <- array3Dto2Dmat(eengrid$Data)
  ind <- rep(1:12,30)
  ind
  
  outmat <- mat 
  
  for (i in 1:12){
    month <- which(ind == i)
    outmat[month,] <- mat[month,] - matref[i,]
  }
   
  
  cor.matrix <- abs(cor(outmat,method="spearman"))
  
  adj.matrix5 <- cor.matrix
  diag(adj.matrix5) <- 0
  adj.matrix5[adj.matrix5 <= .99 ] <- 0
  adj.matrix5[adj.matrix5 > .99] <- 1
  
  return(adj.matrix5)}

tas5d2.adj.matrix6 <- adj.matrix.function6(tas_ncep_5d2)
tas5d2.graph6<- graph_from_adjacency_matrix(tas5d2.adj.matrix6, mode = "undirected")

betw.tas5d2.graph6 <- betweenness(tas5d2.graph6)
log.betw6 <- log(betw.tas5d2.graph6 + 1)

close.tas5d2.graph6 <- closeness(tas5d2.graph6, mode ="all")

points2 <- expand.grid(tas_ncep_5d2$xyCoords$y, tas_ncep_5d2$xyCoords$x)[2:1]
grid3D <- mat2Dto3Darray(betw.tas5d2.graph6,
                         x = tas_ncep_5d2$xyCoords$x,
                         y=tas_ncep_5d2$xyCoords$y)

str(betw.tas5d2.graph6)

betw6.sp <- SpatialPixelsDataFrame(points2, data.frame(betw.tas5d2.graph6))
spplot(betw6.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

log.betw6.sp <- SpatialPixelsDataFrame(points2, data.frame(log.betw6))
spplot(log.betw6.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))

close6.sp <- SpatialPixelsDataFrame(points2, data.frame(close.tas5d2.graph6))
spplot(close6.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))







#### powerlaw attempt -----------------

load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_10d2.rda")


plplot <- function(grid, tau){
  graphObj <- graph_from_Grid(grid, th = tau)
  list_graphObj <- graph2measure(graphObj)
  
  par(mfrow = c(1, 2))
  graph1 <- plot(log(1:length(list_graphObj$ddist)),
     log(list_graphObj$ddist), main = paste0("degree_distribution\nloglog alpha = ", tau), 
     type = "l", 
     xlab = "log (degree K)" , 
     ylab = "log (p(K))")
  graph2 <- plot(1:length(list_graphObj$ddist),
       list_graphObj$ddist, 
     main = paste0("degree_distribution\n alpha = ", tau), 
     type = "l", 
     xlab = "degree K" , 
     ylab = "p(K)")
}

grid = tas_ncep_10d
tau = c(0.1,0.15,0.2,0.3,0.4,0.6,0.7,0.8)
tau[1]

plplotlist <- function(grid, tau){
  par(mar=c(4,4,2,2))
  par(mfrow = c(length(tau)/2, 4))
  
  for (i in 1:length(tau)){
    graphObj <- graph_from_Grid(grid, th = tau[i])
    list_graphObj <- graph2measure(graphObj)
  
    plot(log(1:length(list_graphObj$ddist)),
                 log(list_graphObj$ddist), main = paste0("degree_distribution\nloglog alpha = ", tau[i]), 
                 type = "l", 
                 xlab = "log (degree K)" , 
                 ylab = "log (p(K))")
    plot(1:length(list_graphObj$ddist),
                 list_graphObj$ddist, 
                 main = paste0("degree_distribution\n alpha = ", tau[i]), 
                 type = "l", 
                 xlab = "degree K" , 
                 ylab = "p(K)")
  }
}


plplotlist(tas_ncep_10d,c(0.1,0.15,0.2,0.3,0.4,0.6,0.7,0.8))
plplotlist(tas_ncep_5d2,c(0.1,0.15,0.2,0.3,0.4,0.6,0.7,0.8))
dev.off()




# Show powerlaw plots quickly in RSTudio 
for( a in c(0.2,0.4,0.6,0.7,0.8)){
  par(mfrow = c(5, 1))
  plplot(tas_ncep_10d2, a)
}
for( a in c(0.2,0.4,0.6,0.7,0.8)){
  plplot(tas_ncep_5d2, a)
}
# create powerlaw plots
for(i in c(0.25,0.3,0.35,0.45,0.5,0.55,0.65)){
  tau<-i
  plplotname <- paste0("/Users/lisettegraafland/Documents/R_practice/plots/powerlaw/powlaw10da",tau,".pdf")
  pdf(plplotname)
  plplot(tas_ncep_10d2,tau)
  dev.off()
}


#Pogingen tot Fit exact powerlaw

#Attempt power low fit from highest peak
#Daarna powerlaws fitten.(Data is fout, maar idee plaatje is goed)
plot(3:length(list_cn_0.6_10d2$ddist),
     list_cn_0.6_10d2$ddist[3:length(list_cn_0.6_10d2$ddist)], main = paste0("degree_distribution\nloglog alpha = ", tau), 
     type = "l", 
     xlab = "log (degree K)" , 
     ylab = "log (p(K))")

list_cn_0.6_10d2$ddist
which.max(list_cn_0.6_10d2$ddist)
fit_10d_tau0.6 <- fit_power_law(list_cn_0.6_10d2$ddist[3:length(list_cn_0.6_10d2$ddist)])
fit_10d_tau0.6

plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d2), cn_0.6_10d2$graph)

cn_0.6_10d2$graph

Mode(igraph::degree(f$graph))
f$ddist
str(f$ddist)

which.max(f$ddist)
f$ddist[49]
max(f$ddist*18048)

plot(log(1:length(f$ddist)),log(f$ddist), main = "degree_distribution", type = "l", xlab = "degree K" , ylab = "p(K)")
plot(1:length(f$ddist),f$ddist, main = "degree_distribution", type = "l", xlab = "degree K" , ylab = "p(K)")
fit <- fit_power_law(f$ddist, xmin =f$ddist[49], implementation = "plfit")
str(fit)
fit2 <- fit_power_law(igraph::degree(f$graph))
fit3 <- fit_power_law(f$ddist)
str(fit3)
str(fit2)

x <- as.vector(48:500)
x
fit2$alpha
pl <- function(x,tau){
  y <- x^(-tau)
  return(y)}
pl(log(50))
2^2
2^(-(fit2$alpha)*1000) 

lapply(log(x),function(x){y <- log((x)^(-(fit2$alpha)*18048)) return(y)})
plot(x,lapply(x,pl))

#Oude script voor functie plplot
cn_0.6_10d2 <- graph_from_Grid(tas_ncep_10d2, th = 0.6)
list_cn_0.6_10d2 <- graph2measure(cn_0.6_10d2)
list_cn_0.6_10d2$ddist
length(list_cn_0.6_10d2$ddist)
par(mfrow = c(1, 1))
plot(log(1:length(list_cn_0.6_10d2$ddist)),
     log(list_cn_0.6_10d2$ddist), main = "degree_distribution", 
     type = "l", 
     xlab = "log (degree K)" , 
     ylab = "log (p(K))")
plot(1:length(list_cn_0.6_10d2$ddist),
     list_cn_0.6_10d2$ddist, 
     main = "degree_distribution", 
     type = "l", 
     xlab = "degree K" , 
     ylab = "p(K)")

fit_10d_tau0.6 <- fit_power_law(list_cn_0.6_10d2$ddist)
fit_10d_tau0.6



## tas_Ncep helemaal ------------------

load(file = "/Users/lisettegraafland/Documents/R_practice/Data/ncep/tas_ncep.rda")

tas_ncep.clim <- climatology(tas_ncep)
plotClimatology(tas_ncep.clim, backdrop.theme = "countries")

b <- betweennessGrid(tas_ncep, th = .65) 
b$edens


climb <- betweenness2clim(b, ref.grid = tas_ncep)

plotClimatology(climb, backdrop.theme = "coastline", at = seq(0,17,1.5), col.regions = colores(69))

colores <- colorRampPalette(c("white","lightblue","blue","blue","olivedrab","olivedrab","green", "gold", "orange", "red","red"))
colores2 <- colorRampPalette(c("grey80","blue","olivedrab","gold","red1","red2"))
c <- betweennessGrid(tas_ncep) 
c$edens

climc <- betweenness2clim(c, ref.grid = tas_ncep)

plotClimatology(climc, backdrop.theme = "coastline", at = seq(0,20,1.5), col.regions = colores(69))

d <- betweennessGrid(tas_ncep, th = .73) 
d$edens

climd <- betweenness2clim(d, ref.grid = tas_ncep)

plotClimatology(climd, backdrop.theme = "coastline", at = seq(0,17,1), col.regions = colores2(69))


e<- betweennessGrid(tas_ncep, th = .75) 
e$edens

f<- betweennessGrid(tas_ncep, th = .71) 
str(f)
f$edens
str(f$graph)

igraph::degree(f$graph)

climf <- betweenness2clim(f, what = "betweenness", ref.grid = tas_ncep) 
climf.close <- betweenness2clim(f, what = "closeness", ref.grid = tas_ncep) 
climf.AWcen <- betweenness2clim(f, what = "AWcentrality", ref.grid = tas_ncep)
climf.localclus <- betweenness2clim(f, what = "localclustering", ref.grid = tas_ncep) 




plotClimatology(climf, backdrop.theme = "coastline", at = seq(0,20,.5), col.regions = colores2(69))

plotClimatology(climf.close, backdrop.theme = "coastline",at = seq(0.0000021,0.0000029,0.00000005)) #NA overrekenen goede schaal uitzoeken

plotClimatology(climf.AWcen,backdrop.theme = "coastline")

plotClimatology(climf.localclus,backdrop.theme = "coastline")
## NUEVA -------------------------------

library(transformeR)
library(magrittr)
library(abind)
library(igraph)

# Use of CorFuncNueva functies 20d (voorloper op GridComplexNetworks functies)
degree20 <- betweennessGrid(tas_ncep_20d2, th = .62) 
degree20$edens
clim20d <- betweenness2clim(degree20, what = "betweenness", ref.grid = tas_ncep_20d2) 
clim20d.close <- betweenness2clim(degree20, what = "closeness", ref.grid = tas_ncep_20d2) 
clim20d.AWcen <- betweenness2clim(degree20, what = "AWcentrality", ref.grid = tas_ncep_20d2)
clim20d.localclus <- betweenness2clim(degree20, what = "localclustering", ref.grid = tas_ncep_20d2) 

plotClimatology(clim20d, backdrop.theme = "coastline")
plotClimatology(clim20d.close, backdrop.theme = "coastline") 
plotClimatology(clim20d.AWcen,backdrop.theme = "coastline")
plotClimatology(clim20d.localclus,backdrop.theme = "coastline")

#Use of CorFuncNueva functies 5d (voorloper op GridComplexNetworks functies)
load(file = "/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_5d2.rda")

a <- betweennessGrid(tas_ncep_5d2, th = .68) 
a$edens
str(a$closeness)
str(a$adjacency)
str(a$VertexCoords$y)
str(a$AWcentrality)

#Area weighted connectivity
a$VertexCoords$y # labdas
str(a$adjacency[1,]) # A_vi voor v = 1 en alle i
cos

cos((a$VertexCoords$y)/(180)*pi) # cos(labda_i) voor alle i
sum(cos((a$VertexCoords$y)/(180)*pi)) # sum cos(labda_i) = total area
a$adjacency[3,]%*%cos((a$VertexCoords$y)/(180)*pi) # AWC_1
as.vector(cos((a$VertexCoords$y)/(180)*pi)%*%a$adjacency) #AWC_v voor alle v

clim.AWcen <- betweenness2clim(a, what = "AWcentrality", ref.grid = tas_ncep_5d2)

clim <- betweenness2clim(a, ref.grid = tas_ncep_5d2)
clima.close <- betweenness2clim(a, what = "closeness", ref.grid = tas_ncep_5d2)
clima.locclus <- betweenness2clim(a, what = "localclustering", ref.grid = tas_ncep_5d2)

plotClimatology(clim, backdrop.theme = "coastline", at = seq(0,13,.5), col.regions = topo.colors(69))
plotClimatology(clim, backdrop.theme = "coastline", at = seq(0,15,0.5), col.regions = colores2(69))

plotClimatology(clima.close, backdrop.theme = "coastline",at = seq(0.0000025,0.000003,0.00000005))

plotClimatology(clim.AWcen,backdrop.theme = "coastline")
plotClimatology(clima.locclus,backdrop.theme = "coastline")
str(out)


betweennessGrid <- function(grid, th = 0.8) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% localScaling() %>% redim(drop = TRUE)
  })
  grid <- NULL
  aux <- do.call("bindGrid.time", seas.list) %>% redim(drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  # Correlation matrix
  cor.matrix <- cor(time.coords.matrix, method = "spearman")
  abs.cor.matrix <- abs(cor.matrix)
  adj.matrix <- abs.cor.matrix
  # Adjacency matrix
  diag(adj.matrix) <- 0
  adj.matrix[adj.matrix <= th ] <- 0
  adj.matrix[adj.matrix > th ] <- 1
  # Graph
  graph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  # edge density
  edens <- edge_density(graph, loops = TRUE)
  # degree distrib
  ddist <- degree_distribution(graph)
  # betweenness
  betw <- log1p(betweenness(graph))
  # betweenness edges
  edge.betw <- edge.betweenness(graph)
  # closeness edges
  close <- closeness(graph, mode ="all")
  # AWcentrality
  sumArea <- sum(cos(ref.coords$y/(180)*pi))
  AWcentrality <- as.vector(cos(ref.coords$y/(180)*pi)%*%adj.matrix) / sumArea
  # transitivity (local clustering)
  localclus <- transitivity(graph, type = "local")
  
  out <- list("graph" = graph,
              "data_coords" = time.coords.matrix,
              "corelation" = cor.matrix,
              "edens" = edens, 
              "VertexCoords" = ref.coords ,
              "adjacency" = adj.matrix ,
              "ddist" = ddist, 
              "betweenness" = betw, 
              "closeness" = close, 
              "AWcentrality" = AWcentrality,
              "localclustering" = localclus)
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  return(out)
}
  
betweenness2clim <- function(betObj, what = c("betweenness","closeness","AWcentrality","localclustering"), ref.grid) {
    if (what == "betweenness") {
      mat <- matrix(betObj$betweenness, nrow = 1)
    }
    if (what == "closeness") {
    mat <- matrix(betObj$closeness, nrow = 1)
    }
    if (what == "AWcentrality") {
    mat <- matrix(betObj$AWcentrality, nrow = 1)
    }
    if (what == "localclustering") {
    mat <- matrix(betObj$localclustering, nrow = 1)  
    }
    ref.grid$Data <- mat2Dto3Darray(mat, x = attr(betObj, "Xcoords"), y = attr(betObj, "Ycoords"))
    attr(ref.grid$Data, "climatology:fun") <- what
    return(ref.grid)
}
  
    # points2 <- expand.grid(tas_ncep_5d2$xyCoords$y, tas_ncep_5d2$xyCoords$x)[2:1]
    # betw2.sp <- SpatialPixelsDataFrame(points2, data.frame(betw.tas5d2.graph2))
    # betw.adje.sp <- SpatialPixelsDataFrame(points2, data.frame(betw.adje.graph))
    # spplot(betw2.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))
    # spplot(betw.adje.sp, sp.layout=list(wrld, first =F), col.regions = rev(heat.colors(16)))



plotClimatology(climatology(anom1))


# adjacency matrix function with localscale:
adj.matrix.function4 <- function(grid){
  Rescale2 <-localScaling(grid, time.frame = "monthly")
  Rescale <-redim(Rescale2, drop = T)
  time.coords.matrix <- array3Dto2Dmat(Rescale$Data)
  
  cor.matrix4 <- cor(time.coords.matrix,method="spearman")
  
  
  
  adj.matrix4 <- cor.matrix4
  diag(adj.matrix4) <- 0
  adj.matrix4[adj.matrix4 <= .8 ] <- 0
  adj.matrix4[adj.matrix4 > .8 ] <- 1
  
  return(adj.matrix4)
}

