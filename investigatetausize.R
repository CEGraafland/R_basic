# taus that belong to gridgraphs:
taus <- seq(from = 0, to = 1, length.out = 101)
taus

graphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = NULL)
graphssolo  <- lapply(graphs10d, function(m) m$graph)

nedgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(nedgesnetworksCM, length)


load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
graphs <- lapply(gridGraphs, function(x){x$graph})
edges<- lapply(graphs,E)
nedges <- sapply(edges,length)
length(nedges)

tauextra <- c(0.01,0.025,0.026,0.33)
length(tau)
tauchar <- as.character(tau)
labelsCM <- c()

for(i in 1:length(tauchar)){
  labelsCM[i] <- paste0("= ",tauchar[i])
}


graphs10d <- lapply(tauextra, graph_from_Grid, grid = tas_ncep_10d, subind = NULL)
graphssolo  <- lapply(graphs10d, function(m) m$graph)

nedgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(nedgesnetworksCM, length)
nedgesnetworksCM
