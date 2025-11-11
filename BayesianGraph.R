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


## [graphs creations] --------------
## 10 d graph
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")

plotClimatology(climatology(tas_ncep_10d), backdrop.theme = "coastline")

x <- getCoordinates(tas_ncep_10d)$x
y <- getCoordinates(tas_ncep_10d)$y

indx <- which(findInterval(x, vec = c(120,172), left.open = TRUE)==1)
indy <- which(findInterval(y, vec = c(-10,12), left.open = TRUE)==1)

# 10degree without el nino basin
tas_ncep_10dnonino <- tas_ncep_10d
for(i in 1:length(indx)){
  for(l in 1:length(indy)){
    tas_ncep_10dnonino$Data[,indy[l],indx[i]] <- NA
  }
}

plotClimatology(climatology(tas_ncep_10dnonino), backdrop.theme = "coastline")

### Time.coords.matrix from graph_from_Grid---------
## 10 d graph
degree10d <- graph_from_Grid(tas_ncep_10d, th = .68) #Check of nu consistent DAG
CNdegree10d<- graph2measure(degree10d)
CNdegree10d$edens

## 20 d graph 
str(tas_ncep_20d2)
degree20d <- graph_from_Grid(tas_ncep_20d2, th = .61) 
degree20d$data_coords
attr(degree20d, "Xcoords", exact = FALSE)

CNdegree20d <- graph2measure(degree20d)
CNdegree20d$edens

## 30 d graph
load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_30d.rda")
degree30d <- graph_from_Grid(tas_ncep_30d, th = .68) #Check of nu consistent DAG
CNdegree30d<- graph2measure(degree30d)
CNdegree30d$edens

## 45 d graph
load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_45d.rda")
degree45d <- graph_from_Grid(tas_ncep_45d, th = 0.3) #Check of nu consistent DAG
degree45d$graph
CNdegree45d<- graph2measure(degree45d)
CNdegree45d$edens

## 45 / 90 graph
load("/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_9045d.rda")
degree9045d <- graph_from_Grid(tas_ncep_9045d, th = .4) 
CNdegree9045d<- graph2measure(degree9045d)
CNdegree9045d$edens

degree9045data <- TimeCoordsAnom_from_Grid(tas_ncep_9045d)

## [complex Networks properties 20d graph]------
clim20d <- betweenness2clim(degree20, what = "betweenness", ref.grid = tas_ncep_20d2) 
clim20d.close <- betweenness2clim(degree20, what = "closeness", ref.grid = tas_ncep_20d2) 
clim20d.AWcen <- betweenness2clim(degree20, what = "AWcentrality", ref.grid = tas_ncep_20d2)
clim20d.localclus <- betweenness2clim(degree20, what = "localclustering", ref.grid = tas_ncep_20d2) 

plotClimatology(clim20d, backdrop.theme = "coastline")
plotClimatology(clim20d.close, backdrop.theme = "coastline") 
plotClimatology(clim20d.AWcen,backdrop.theme = "coastline")
plotClimatology(clim20d.localclus,backdrop.theme = "coastline")

#[Bayesian Networks properties 20d]
#Inverse Corelation: Conditional independencies between XY given all other variables
library(corpcor)
degree20$corelation

invcor <- cor2pcor(degree20$corelation)
dimnames(invcor) <- dimnames(degree20$corelation)
str(invcor)
invcor[1,2] # Conditional correlation for coordinates 1,2 given all others.

####Learning DAG STRUCTURE from Data_coords ----------
##Conditional independent test
#iamb
#(Variation in test, alpha, amount of nodes result
#in (un)directed, (a)cyclic, extendable....
print(degree20$data_coords)
str(degree20)
tas_position_20d<- data.frame(degree20$data_coords)
save(tas_position_20d, file = "/Users/lisettegraafland/Documents/R_practice/Data/tas_position_20d.rda")
struDegree20Debug <- iamb(tas_position_20d, debug = TRUE)
struDegree20Debug
struDegree20 <- iamb(tas_position_20d, test = "cor") 

struDegree20 <- iamb(data.frame(degree20$data_coords), test = "cor") 
struDegree20op <- iamb(data.frame(degree20$data_coords), test = "cor", optimized = TRUE)
warnings()

skeleton(struDegree20)

directed(struDegree20op)
acyclic(struDegree20op, directed = FALSE)

directed(struDegree20)
acyclic(struDegree20, directed = TRUE)
struDegree20

struDegree20test <- iamb(data.frame(degree20$data_coords), alpha = 0.00000001) # higher a -> higher dependency -> more v-structs
struDegree20un <- iamb(data.frame(degree20$data_coords), test = "cor", undirected = TRUE, debug = TRUE)
length(vstructs(struDegree20un))
length(vstructs(struDegree20test))

struDegree20Marco <- iamb(data.frame(degree20$data_coords), alpha = 1e-10)
struDegree20Marco
directed(struDegree20Marco)
acyclic(struDegree20Marco)
cextend(struDegree20Marco)

#lower alpha, less v-structures, less arcs (dependencies)
degree10d <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
degree10d <- as.data.frame(degree10d)
blub <- pc.stable(degree10d, test = "spearman", alpha = 1*10^-13)

blub_a <- mmpc(degree10d, cluster = NULL, whitelist = NULL, blacklist = NULL, test = NULL,
               alpha = 1*10^-26, B = NULL, debug = FALSE, optimized = FALSE, strict = FALSE,
               undirected = FALSE)
nparams(bn.fit(cextend(blub_a), degree10d))
360/1896
blub2 <- mmhc(degree10d, whitelist = NULL, blacklist = NULL, restrict.args = list(),
                  maximize.args = list(), debug = FALSE)
blub2a <- mmhc(degree10d, whitelist = NULL, blacklist = NULL, restrict.args = list(),
           maximize.args = list(max.iter = 600, start = NULL), debug = FALSE)
blub2b <- mmhc(degree10d, whitelist = NULL, blacklist = NULL, restrict.args = list(alpha = 0.1),
               maximize.args = list(max.iter = 600, start = NULL), debug = FALSE)
blub3 <- hc(degree10d)
blub4 <- rsmax2(degree10d, restrict = "pc.stable", maximize = "hc", restrict.args = list(max.iter = 200), maximize.args = list())
1+1

nparams(bn.fit(blub2b, degree10d))
struDegree10 <- iamb(degree10d, test = "cor") # 50 v-structures
struDegree10a <-iamb(data.frame(degree10d), test = "mi-g", optimized = TRUE, alpha = 1e-100)
struDegree10b <-iamb(data.frame(degree10d), test = "mi-g", optimized = TRUE, alpha = 1e-75) # 225 arcs.
struDegree10c <-iamb(data.frame(degree10d), test = "mi-g", optimized = TRUE, alpha = 1e-60) # 321 arcs.
struDegree10d <-iamb(data.frame(degree10d), test = "mi-g", optimized = TRUE, alpha = 1e-55) # 351 arcs.
cextend(struDegree10d)
struDegree10f  <-iamb(data.frame(degree10d), test = "mi-g", optimized = TRUE, alpha = 1e-53) # 351 arcs.
cextend(struDegree10f) # 363 arcs.
struDegree10e <-iamb(data.frame(degree10d), test = "mi-g", optimized = TRUE, alpha = 1e-50) # 351 arcs.
cextend(struDegree10e) # cycles

struDegreegs_a <- gs(degree10d, alpha =  1*10^-35)
cextend(struDegreegs_a)
struDegree10cor <-iamb(data.frame(degree10d), test = "cor", optimized = TRUE, alpha = 1e-50) # 361 arcs.
cextend(struDegree10cor) # for 1e-50 also cycles!
#vstructure X415 -> X433 <- X469 is not applicable, because one or both arcs introduce cycles in the graph.
#idee verlaag corelationcoef voor de test!! meer dependencies minder v-structures

struDegree30 <- iamb(data.frame(degree30d$data_coords), test = "cor", optimized = FALSE) # 8 v-structures not
struDegree30 <- iamb(data.frame(degree30d$data_coords), test = "cor", optimized = TRUE) # 6 v-structures not
struDegree30 <- iamb(data.frame(degree30d$data_coords), test = "zf") # 8 v-structures
acyclic(struDegree30, directed = TRUE) #FALSE (Both un and directed)
cextend(moral(struDegree30)) #no consistent extension


struDegree30a <- iamb(data.frame(degree30d$data_coords), test = "cor", optimized = TRUE,alpha = 0.009)
acyclic(struDegree30a, directed = FALSE) #No directed cycles, but cycles
directed(struDegree30a)
cextend(struDegree30a) #NOOP
cextend(moral(struDegree30a))


struDegree45 <- iamb(data.frame(degree45d$data_coords), test = "cor") #perfect
struDegree45

struDegree9045 <-iamb(data.frame(degree9045d$data_coords), test = "cor")

str(struDegree20)
amat(struDegree20)

#gs 
stru.gs.Degree20 <- gs(data.frame(degree20$data_coords))
acyclic(stru.gs.Degree20, directed = TRUE) # false; directed cycles
cextend(stru.gs.Degree20)
directed(stru.gs.Degree20) #false

#fast_iamb
stru.fast.iamb.Degree20 <-fast.iamb(data.frame(degree20$data_coords))
acyclic(stru.fast.iamb.Degree20, directed = TRUE) #False
directed(stru.fast.iamb.Degree20) #False
skeleton(stru.fast.iamb.Degree20)
pdag2dag(stru.fast.iamb.Degree20,keepVstruct = FALSE)
cextend(skeleton(stru.fast.iamb.Degree20)) #No consistent extension of skeleton is possible

#inter.iamb.
stru.inter.iamb.Degree20 <- inter.iamb(data.frame(degree20$data_coords))



#Hybrid algorithms: combination of CI -test and Greedysearch
struDegree30mmhc <- mmhc(data.frame(degree30d$data_coords), 
                         restrict.args = list(test = "cor", alpha = 0.1))
acyclic(struDegree30mmhc) # TRUE
directed(struDegree30mmhc) # TRUE
struDegree30hc <- hc(data.frame(degree30d$data_coords))
score(struDegree30mmhc, data.frame(degree30d$data_coords))
score(struDegree30hc, data.frame(degree30d$data_coords))

## characteristics resulting structures
#Edges multiple/simplify no use
#mutual gives twosided igroph
cextend(struDegree20test)
acyclic(struDegree20test, directed = FALSE) # no directed cycles

cextend(struDegree10)
acyclic(struDegree10) #Contains cycles (undirected?)
directed(struDegree10) #Not directed

cextend(struDegree30)
acyclic(struDegree30) #Contains cycles (undirected?)
directed(struDegree30) #Not directed

undirected.arcs(struDegree20)

# Evaluation with Igraph (same conclusions)
iGraphDegree20 <- igraph.from.graphNEL(as.graphNEL(struDegree20))
print(iGraphDegree20)
is.directed(iGraphDegree20) #YES
iGraphDegree20.dir <- as.directed(iGraphDegree20, mode = "arbitrary") #already was
all.equal(iGraphDegree20,iGraphDegree20.dir) #NOT SURE
print(iGraphDegree20.dir)

is.dag(iGraphDegree20) #FALSE: Directed but cyclic, with mutual edges
iGraphDegree20.pdag <- as.undirected(iGraphDegree20, mode = "mutual") #only returns graph with former mutual edges now as undirected
iGraphDegree20.simp <- simplify(iGraphDegree20) #Does nothing because does not find mutuals

all.equal(iGraphDegree20.simp,iGraphDegree20)
is.multiple(iGraphDegree20)# AL FALSE

V(iGraphDegree20)
E(iGraphDegree20)

plot(iGraphDegree20)


#HileClimbing test HILE CLIMBING DEBUG CHECKEN
#20 degree
stru.hc.Degree20 <- hc(data.frame(degree20$data_coords))
stru.hc.Degree20igraph <- igraph.from.graphNEL(as.graphNEL(struDegree20))
str(stru.hc.Degree20)
plot(stru.hc.Degree20)

all.equal(struDegree20,stru.hc.Degree20)
score(struDegree20, data = data.frame(degree20$data_coords), type = "bic-g") #NOT TOTALLY DIRECTED
score(stru.hc.Degree20, data = data.frame(degree20$data_coords), type = "bic-g")


# 9045 degree
stru.hc.Degree9045 <- hc(data.frame(degree9045d$data_coords))
stru.hc.Degree9045igraph <- igraph.from.graphNEL(as.graphNEL(stru.hc.Degree9045))
all.equal(stru.hc.Degree9045, struDegree9045)

#Learning parameters
help(bic)
??big
??k2

bnDegree20 <- bn.fit(stru.hc.Degree20, data.frame(degree20$data_coords), method = "mle")
install.packages("rbmn")
library(rbmn)

## plotting graphstructure ------------


#nodig: Xcoords, Ycoords, grid

TimeCoordsAnom_from_Grid <- function(grid, subind = NULL) {
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
    subsetGrid(grid, season = seas[i]) %>% localScaling(scale = TRUE) %>% redim(drop = TRUE) })
  #localScaling(grid,time.frame = "monthly")
  
  grid <- NULL
  aux <- do.call("bindGrid.time", seas.list) %>% redim(drop = TRUE)
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


TimeCoordsAnom_from_Grid_rms <- function(grid, subind = NULL, rms = FALSE) {
  seas <- getSeason(grid)
  coords <- getCoordinates(grid)
  x <- coords$x
  y <- coords$y
  ref.coords <- expand.grid(y, x)[2:1]
  names(ref.coords) <- c("x", "y")
  ref.dates <- getRefDates(grid)
  seas.list <- lapply(1:length(seas), function(i) {
    subsetGrid(grid, season = seas[i]) %>% localScaling() %>% redim(drop = TRUE) })
    #localScaling(grid,time.frame = "monthly")
 
  grid <- NULL
  aux <- do.call("bindGrid.time", seas.list) %>% redim(drop = TRUE)
  seas.list <- NULL
  time.coords.matrix <- array3Dto2Dmat(aux$Data)
  if (rms == TRUE) {
    time.coords.matrix <- scale(time.coords.matrix, center = FALSE, scale = TRUE)}
  if (!is.null(subind)) {
    time.coords.matrix <- time.coords.matrix[subind,]}
  
  out <- time.coords.matrix
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  attr(out, "VertexCoords") <- ref.coords
  return(out)
}

time.coords.10d.rms <-  TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
time.coords.10d.std <-  TimeCoordsAnom_from_Grid_std(tas_ncep_10d)

var(as.data.frame(time.coords.10d.rms))
var(as.data.frame(time.coords.10d.std))
TimeCoordsAnom_from_Grid(tas_ncep_10d)

cpqu
#Test 10d  
all.equal(TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = ind1981_2010),TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = NULL))
TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = ind1981_2010)
str(TimeCoordsAnom_from_Grid(tas_ncep_10d))
all.equal(data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d)),data.frame(degree10d$data_coords)) #yes

# Plots dags produced by hc or iamb or pc
# Entries: time.coords.matrix from TimeCoordsAnom_from_Grid
#         or obj from graph_from_Grid
#         !!! uses only its attributes : introduce grid could also work !!!
# Outcome: plot of dag with worldmap
plot.Meteodag <- function(time.coords, meteodag){
  if (class(meteodag) == "igraph") 
    igraphDegree <- meteodag
  
  if (class(meteodag) == "graphNEL") 
    igraphDegree <- igraph.from.graphNEL(meteodag)
    
  if (class(meteodag) == "bn") 
    igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))

  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)

  points <- expand.grid(y, x)[2:1]

  plot(wrld)
  plot.igraph(igraphDegree, 
            vertex.size = 100,
            vertex.color = "blue",
            vertex.label = NA,
            edge.color= "green",
            edge.arrow.size = 0.3,
            edge.lty = 2,
            layout = as.matrix(points), add = TRUE, rescale = FALSE)
  
}

#Test
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d), struDegree10)
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_20d2),struDegree20)
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_9045d),struDegree9045)
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_20d2), struDegree20test)


plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_20d2),struDegree20Marco)


#VRAAG JAUCO: ANDERE ONDERGROND METEO DAG
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d), graph_from_Grid(tas_ncep_10d)$graph)

load(file = "/Users/lisettegraafland/Documents/R_practice/Data/tas_ncep_20d2.rda")
tas_ncep_20d2.clim <- climatology(tas_ncep_20d2)
plotClimatology(tas_ncep_20d2.clim, backdrop.theme = "countries")


load(system.file("countries.rda", package = "transformeR"))
class(l1[1])
l1[1]
l1[2]
str(wrld)



#oude plots
x <- attr(degree10, "Xcoords", exact = FALSE)
y <- attr(degree10, "Ycoords", exact = FALSE)

points20D <- expand.grid(y, x)[2:1]

plot(wrld)
edge.attributes(iGraphDegree10, index = )
plot(iGraphDegree20, layout = as.matrix(points20D), add = TRUE, rescale = FALSE)
plot.igraph(iGraphDegree20, 
            vertex.size = 10,
            vertex.color = "blue",
            vertex.label = NA,
            edge.color= "green",
            edge.arrow.size = 0.3,
            edge.lty = 2,
            layout = as.matrix(points20D), add = TRUE, rescale = FALSE)

plot(wrld) #set_vertex_attr(graph = stru.hc.Degree20igraph, name = "size", value = 10)
plot.igraph(stru.hc.Degree20igraph, 
     vertex.size = 10,
     vertex.color = "blue",
     vertex.label = NA,
     edge.color= "green",
     edge.arrow.size = 0.3,
     edge.lty = 2,
     layout = as.matrix(points20D), add = TRUE, rescale = FALSE)


as_data_frame(iGraphDegree20, what = "both")
str(as_data_frame(iGraphDegree20, what = "vertices")$name)

iGraphDegree20.sp <- SpatialPixelsDataFrame(points20D, data.frame(as_data_frame(iGraphDegree20, what = "vertices")$name))
class(iGraphDegree20.sp)
spplot(iGraphDegree20.sp, sp.layout=list(wrld, first =F))

#Plot 45 d
struDegree45
as.graphAM(struDegree45dir)
acyclic(struDegree45) #FALSE
directed(struDegree45) #FALSE
cextend(struDegree45) #CYCLES
ske45 <- skeleton(struDegree45)
cextend(ske45) #no consistent extension possible
cextend(cpdag(ske45))

#undirected arcs
undir45 <- as.matrix(undirected.arcs(struDegree45))
undir45
undir45rm <- undir45[c(-1,-2,-7,-8,-10),]
undir45rm
struDegree45dir <- struDegree45
struDegree45dir <- drop.edge(struDegree45, unname(undir45rm[3,1]), unname(undir45rm[3,2]))
#remove undirected arcs manual
nrow(undir45rm)
for(i in 1:nrow(undir45rm)){
struDegree45dir <- drop.edge(struDegree45dir, unname(undir45rm[i,1]), unname(undir45rm[i,2]))
}

directed(struDegree45dir)
acyclic(struDegree45)

graphviz.plot(ske45)
graphviz.plot(struDegree45, highlight = list("arcs" = as.matrix(undirected.arcs(struDegree45))), layout = "dot", main = NULL, sub = NULL)#X6 - X30 - X2 - X10 - X28 - X29 - X26 - X6
graphviz.plot(struDegree45dir, highlight = list("arcs" = as.matrix(undirected.arcs(struDegree45))), layout = "dot", main = NULL, sub = NULL)
#X6 - X30 - X2 - X10 - X28 - X29 - X26 - X6

#plot 9045 d
struDegree9045
acyclic(struDegree9045) # TRUE!!
directed(struDegree9045) # TRUE !!
graphviz.plot(struDegree9045)

igraphDegree9045 <- igraph.from.graphNEL(as.graphNEL(struDegree9045))
  
x <- attr(degree9045d, "Xcoords", exact = FALSE)
y <- attr(degree9045d, "Ycoords", exact = FALSE)

points9045D <- expand.grid(y, x)[2:1]

plot(wrld)
plot.igraph(igraphDegree9045, 
            vertex.size = 100,
            vertex.color = "blue",
            vertex.label = NA,
            edge.color= "green",
            edge.arrow.size = 0.3,
            edge.lty = 2,
            layout = as.matrix(points9045D), add = TRUE, rescale = FALSE)


#plot 30d a
struDegree30a
acyclic(struDegree30a) # FLASE! But no directed
directed(struDegree30a) # FLASE !!
graphviz.plot(struDegree30a, highlight = list("arcs" = as.matrix(undirected.arcs(struDegree30a))))

igraphDegree30a <- igraph.from.graphNEL(as.graphNEL(struDegree30a))

x <- attr(degree30d, "Xcoords", exact = FALSE)
y <- attr(degree30d, "Ycoords", exact = FALSE)

points30aD <- expand.grid(y, x)[2:1]

plot(wrld)
plot.igraph(igraphDegree30a, 
            vertex.size = 100,
            vertex.color = "blue",
            vertex.label = NA,
            edge.color= "green",
            edge.arrow.size = 0.3,
            edge.lty = 2,
            layout = as.matrix(points30aD), add = TRUE, rescale = FALSE)


## problem 1 ----------------------

graphmaker <- function(data, data.prop = 1, hyp.alpha = 0.05){
  
  if(!is.null(names(data))){nodenames <- names(data)} 
  
  prop <- data.prop
  nobs <- base::nrow(data)
  
  train <- as.data.frame(matrix(data = NA, nrow = prop*nobs, ncol = ncol(data)))
  names(train) <- nodenames
  for(i in 1:ncol(train)){
    train[1:(prop*nobs),i] <- data[1:(prop*nobs),i]
  }
  
  train.pdag <- iamb(train, test = "cor", alpha = hyp.alpha)
  train.dag <- cextend(train.pdag)
  
  return(train.dag)
  
}

#test
class(graphmaker(gaussian.test, 0.5))
graph9045d <- graphmaker(data.frame(degree9045d$data_coords))
nodes(graph9045d)
graph9045d

compare <- function(graph, data, train.prop = 0.5){
  if(!is.null(names(data))){nodenames <- names(data)
  } else {nodenames <- nodes(graph)}
  
  prop <- train.prop
  nobs <- base::nrow(data)
  
  train <- as.data.frame(matrix(data = NA, nrow = prop*nobs, ncol = length(nodenames)))
  names(train) <-nodenames
  
  for(i in names(train)){
    train[1:(prop*nobs),i] <- data[1:(prop*nobs),i]
  }
  
  testdata <- as.data.frame(matrix(data = NA, nrow = ((1-prop)*nobs), ncol = length(nodenames)))
  names(testdata)<- nodenames
  
  for(i in names(testdata)){
    testdata[1:((1-prop)*nobs),i] <- data[(prop*nobs+1):(nobs),i]
  }
  
  train.gbn <- bn.fit(graph, train)
  
  liktrain <- logLik(train.gbn, train, nodenames)
  liktestdata <- logLik(train.gbn, testdata, nodenames)
  
  #scoretrain <- score(graph, train)
  #scoretest <- score(train.gbn, testdata)
  
  print(paste0("Loglikelihood train is: ",liktrain)) 
  print(paste0("Loglikelihood test is: ",liktestdata))
  
  return(c(liktrain,liktestdata))
  
}

#Compare tests
#Gaussian
compare(train.dag, gaussian.test, 0.75) #
compare(graphmaker(gaussian.test),gaussian.test,0.75)

#graph9045d
compare(graph9045d, data.frame(degree9045d$data_coords), 0.5)
compare(graphmaker(data.frame(degree9045d$data_coords), 0.5, 0.01), data.frame(degree9045d$data_coords), 0.5)
class(struDegree9045)
compare(stru.hc.Degree9045, data.frame(degree9045d$data_coords),0.5)


#Background material function Compare
nodenames <- names(gaussian.test)
length(nodenames)
prop<- 0.5
nobs <- base::nrow(gaussian.test)
prop*nobs
(1-prop)*nobs

train <- as.data.frame(matrix(data = NA, nrow = prop*nobs, ncol = length(nodenames)))
names(train) <-nodenames
rownames(train) <- NULL
head(train)
dimnames(train)

for(i in names(train)){
train[1:(prop*nobs),i] <- gaussian.test[1:(prop*nobs),i]
}
train$A

nobs - prop*nobs
testdata <- as.data.frame(matrix(data = NA, nrow = ((1-prop)*nobs), ncol = length(nodenames)))
names(testdata)<- nodenames
head(testdata)

for(i in names(testdata)){
  testdata[1:((1-prop)*nobs),i] <- gaussian.test[(prop*nobs+1):(nobs),i]
}

testdata$A

train.cpdag <- iamb(train)

acyclic(train.cpdag)
directed(train.cpdag)

all.equal(train.cpdag,test.iamb) #
train.dag <- cextend(train.cpdag)
train.gbn <- bn.fit(train.dag, train)

score(train.dag, train, type = "loglik-g")
score(train.dag, testdata, type = "loglik-g")
logLik(train.gbn, train, nodenames)
logLik(train.gbn, testdata, nodenames)

#Function alpha vs. loglik.
#alpha_loglik <- function(data, alpha){
#  if(!class(data) == "data.frame") 
#    data <- as.data.frame(data)
#  eqclass <- iamb(data, test = "cor", alpha = alpha)
#  dag <- cextend(eqclass)
#  gbn <- bn.fit(dag, data)
#  likgbn <- logLik(gbn, data)
#  return(likgbn)
#}

#alpha_loglik_pc <- function(data, alpha){
#  degree30stat <- list(C = data$correlation, n = nrow(data$data_coords))
#  eqclass <- pc(degree30stat, indepTest = gaussCItest, p = ncol(data$data_coords), alpha = alpha)
#  dag <- pcalg::pdag2dag(eqclass@graph)
#  bn.dag <- as.bn(dag$graph)
#  newdata <- as.data.frame(data$data_coords)
#  names(newdata) <- nodes(bn.dag)
#  gbn <- bn.fit(bn.dag, newdata)
#  likgbn <- logLik(gbn, newdata)
#  return(likgbn)
#}

#alpha_loglik_pc2 <- function(data, alpha){
#  degreestat <- list(C = data$correlation, n = nrow(data$data_coords))
#  eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data$data_coords), alpha = alpha)
#  amata <- as(eqclass,"amat")
#  
#  if(isValidGraph(amata, type = "pdag")){
#    dag <- pcalg::pdag2dag(eqclass@graph)
#    bn.dag <- as.bn(dag$graph)
#    newdata <- as.data.frame(data$data_coords)
#    names(newdata) <- nodes(bn.dag)
#    gbn <- bn.fit(bn.dag, newdata)
#    likgbn <- logLik(gbn, newdata)
#  }
#  else likgbn <- NA
#  
#  return(likgbn)
#}


alpha_loglik_pc(degree30d, 0.05)

mapply(alpha_loglik_pc2, list(degree20d), myalphas0)
mapply(alpha_loglik_pc, list(degree30d), myalphas)
mapply(alpha_loglik_pc, list(degree30d), myalphas2) # example alpha mas pequena pero loglik mas pequena
mapply(alpha_loglik_pc, list(degree45d), myalphas)
mapply(alpha_loglik_pc, list(degree9045d), as.vector(seq(0.1,0.9,0.1)))

alpha_loglik_pc2(degree9045data, 0.05)
str(alpha_loglik(degree9045data, 0.14))

## [test gaussian ]------------------------
test.iamb <- iamb(gaussian.test)
gaussian.test2 <- gaussian.test
names(gaussian.test2) <- NULL
str(as.matrix(gaussian.test2))
test.iamb2 <- iamb(as.matrix(gaussian.test2))
test.iamb2
test.cpdag <- cpdag(test.iamb)
undirected.arcs(test.cpdag)
test.cpdagNEL <- as.graphNEL(test.cpdag)
test.cpdag.igraph <- igraph.from.graphNEL(test.cpdagNEL)
plot(test.cpdag.igraph)

test.dag <- cextend(test.cpdag) #no problem
test.dagNEL <- as.graphNEL(test.dag)
test.dag.igraph <- igraph.from.graphNEL(test.dagNEL)
plot(test.dag.igraph)

plot
gaussian.test


data(learning.test)
res = gs(learning.test)

acyclic(res)
directed(res)
res = pdag2dag(res, ordering = LETTERS[1:6])
directed(res)
skeleton(res)

## [pcalg (test)]--------

#klad
degree45d$correlation
degree45stat <- list(C = degree45d$correlation, n = nrow(degree45d$data_coords))
pc.fit.45d <- pc(degree45stat, indepTest = gaussCItest, p = ncol(degree45d$data_coords), alpha = 0.05)
pc.fit.45d
pcdag45d <- pdag2dag(pc.fit.45d@graph)$graph


igraphDegree45d <- igraph.from.graphNEL(pdag2dag(pc.fit.45d@graph)$graph)
igraphDegree45d
acyclic(as.bn(pc.fit.45d@graph))

plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_45d), pc.fit.45d@graph)

degree30d$correlation
degree30stat <- list(C = degree30d$correlation, n = nrow(degree30d$data_coords))
pc.fit.30d <- pc(degree30stat, indepTest = gaussCItest, p = ncol(degree30d$data_coords), alpha = 0.01)
pc.fit.30d
pdag2dag(pc.fit.30d@graph)
pcdag30d <- pdag2dag(pc.fit.30d@graph)$graph
pcdag30d

igraph.from.graphNEL(pc.fit.30d@graph)
acyclic(as.bn(pc.fit.30d@graph))
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_30d), pcdag30d)
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_30d), pc.fit.30d@graph)

#pc alpha 0.5 ANDERE ALPHAS MAKEN.
degree20d$correlation
degree20stat <- list(C = degree20d$correlation, n = nrow(degree20d$data_coords))
degree30stat <- list(C = degree30d$correlation, n = nrow(degree30d$data_coords))
pc.fit.20d <- pc(degree20stat, indepTest = gaussCItest, p = ncol(degree20d$data_coords), alpha = 0.04)

amatpc20da0.4 <-as(pc.fit.20d, "amat")
value <- isValidGraph(amatpc20da0.4, type = "pdag")
value

pc.fit.20d.a0.001 # 20
pc.fit.20d.a0.01 # 10
pc.fit.20d.a0.02 # 13
pc.fit.20d.a0.03 # 4
pc.fit.20d.a0.04 # 7 
pc.fit.20d.a0.05 # 4 
pc.fit.20d.a0.1 # 4

### pcalg ----------

mypc <- function(alpha){
  eqclas<- pc(degree20stat, indepTest = gaussCItest, p = ncol(degree20d$data_coords), alpha = alpha)
  amata <- as(eqclas,"amat")
  value <- isValidGraph(amata, type = "pdag")
  bneqclas <- as.bn(eqclas@graph)
  nund <- length(undirected.arcs(bneqclas))
  nedges <- narcs(bneqclas)
  return(data.frame(nund, nedges, value))
}

#mypc30 <- function(alpha){
#  eqclas<- pc(degree30stat, indepTest = gaussCItest, p = ncol(degree30d$data_coords), alpha = alpha)
#  amata <- as(eqclas,"amat")
#  value <- isValidGraph(amata, type = "pdag")
#  bneqclas <- as.bn(eqclas@graph)
#  nund <- length(undirected.arcs(bneqclas))
#  nedges <- narcs(bneqclas)
#  return(data.frame(nund, nedges, value))
#}
#
#mypc30(0.01)
#
#rep30d <- rep(degree30d,length(myalphas))

mypc2 <- function(data, alpha){
  datastat <- list(C = data$correlation, n = nrow(data$data_coords))
  eqclas<- pc(datastat, indepTest = gaussCItest, p = ncol(data$data_coords), alpha = alpha)
  amata <- as(eqclas,"amat")
  vcpdag <- isValidGraph(amata, type = "cpdag")
  vpdag <- isValidGraph(amata, type = "pdag")
  
  bneqclas <- as.bn(eqclas@graph)
  nund <- length(undirected.arcs(bneqclas))
  nedges <- narcs(bneqclas)
  return(data.frame(nund, nedges, vcpdag, vpdag))
}

myalphas0 <- as.vector(seq(0.0001, 0.001, 0.0001))
myalphas1 <- as.vector(seq(0.001, 0.01, 0.001))
myalphas2 <- as.vector(seq(0.01, 0.1, 0.01))
myalphastot <- as.vector(seq(0.001,0.1,0.001))

mapply(mypc2, list(degree20d), myalphas1) #Ninguno
mapply(mypc2, list(degree20d), myalphas2) #Ninguno
mapply(mypc2, list(degree20d), myalphas0) #Algunos

mapply(mypc2, list(degree30d), myalphas1)
mapply(mypc2, list(degree45d), myalphas1) 
# SEe which Alphas are ok 
#enter in alpha_loglik_pc

apply(0.01,mypc2)

sapply(myalphas, mypc2)
sapply(myalphas2, mypc)
sapply(myalphas, mypc2)
sapply(myalphas, mypc30)




## HIER KIJKEN WAT ER AAN DE HAND IS
undirected.arcs(bn.pc.fit.20d.a0.5)
bn.pc.fit.20d.a0.5.dir <- drop.arc(bn.pc.fit.20d.a0.5,"100","136")
acyclic(bn.pc.fit.20d.a0.5.dir) #FALSE
bn.pc.fit.20d.a0.5.dir <- drop.arc(bn.pc.fit.20d.a0.5.dir,"136","100")
acyclic(bn.pc.fit.20d.a0.5.dir) #FALSE
outgoing.arcs(bn.pc.fit.20d.a0.5.dir, 100) # BEKIJKEN OF DE ARC REMOVED IS

nodes(bn.pc.fit.20d)
pdag2dag(pc.fit.20d@graph)
bn.pdag.pc.fit.20d <- as.bn(pdag2dag(pc.fit.20d@graph)$graph)
pcdag20d <- pdag2dag(pc.fit.20d@graph)$graph

data20d2 <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_20d2))
names(data20d2) <- names(bn.pc.fit.20d$nodes) 
data20d2

logLik(bn.pc.fit.20d,data20d2)
names(bn.pc.fit.20d$nodes)

if (length(setdiff(names(bn.pc.fit.20d$nodes) , names(as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_20d2))))) != 0)
  print("the variables in the data and in the network do not match.")
setdiff(names(bn.pc.fit.20d$nodes) , names(as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_20d2))))

all.equal(pc.fit.20d@graph,pcdag20d)
plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_20d2), pcdag20d)

plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_45d), pc.fit.45d@graph)


