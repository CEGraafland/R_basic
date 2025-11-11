################################################################################
# Logarithmic bining CN and Erdos renyi and Lattice and BN
################################################################################
#source("Functions/propagationFunctions.R")
#load("../Data/tas_ncep_10d.rda")
#load("../Data/gridGraphs.rda")

setwd("~/data/Untitled/Trabajo/R_practice/R/")
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")
setwd("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")
##########################################################################################
rm(list = ls())
library(ggplot2)
library(gridExtra)
library(gridGraphics)
library(grid)
library(transformeR)
library(magrittr)
library(igraph)
library(bnlearn)
library(scales)

load("../Data/tas_ncep_10d.rda")
load("../Data/data_hcnetworks10d.rda")
load("../Data/Data_graphsPC.rda")
load("../Data/data_pcnetworks10d.rda")
source("Functions/CN_ConstructionandMeasuresFunctions.R")
source("Functions/BasicNetworkFunctions.R")
#######################################################################################
###############################################################################################
# Combining powerlaw figures bayesian and complex networks
###############################################################################################
# These with spearman in results 3
taucomp <- c(0.8,
             0.7,
             #0.64,# (1512)
             0.62,
             0.52,
             0.49,
             0.41,
             0.345,
             0.32,
             0.31,
             0.3,
             0.25,
             0.2)
# Thesis with pearson 
taucomp <- c(0.49)

grid <- tas_ncep_10d
ggplotsCMcomp <- list()

ggplotsCMcomp <- list()
coefficients <- c()
numberofedges <- c()
trans <- c()

i <- 1
for (i in 1:length(taucomp)){
  graphObj <- graph_from_Grid(grid, th = taucomp[i], method = 'spearman')
  numberofedges[i] <- length(E(graphObj$graph))
  trans[i] <- transitivity(graphObj$graph, type = "global")
  degrees <- igraph::degree(graphObj$graph)
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(graphObj$graph)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  # max amount of counts
  #dist <-hist(degrees, breaks = k_max)
  #max(dist$counts)
  
  # bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
  # Calculate breaks: binlims
  nbins <- ceiling(log2(k_max))
  nbins
  binlims <- 2^(1:nbins)
  binlims <- append(binlims, 0, after = 0)
  binlims
  # Count amount values in bins with hist.
  # replace the counts by their normalized values (divided by binwidth)
  histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
  histdegrees2 <- histdegrees
  binwidth1 <- diff(histdegrees2$breaks, lag = 1)
  histdegrees2$counts <- histdegrees2$counts/binwidth1
  # Convert to adequate dataframe for ggplot.
  # Plot on log2 log2 scale. 
  histdegrees$counts
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks", "counts")
  # form <-lm(datahist$counts ~ datahist$breaks)
  # form$coefficients
  
  histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/2)*2,1), plot =FALSE)
  binwidth2 <- diff(histlinear$breaks, lag = 1)
  histlinear$counts <- histlinear$counts/binwidth2
  datalinear <- data.frame(histlinear$breaks[2:length(histlinear$breaks)], histlinear$counts)
  names(datalinear) <- c("breaks", "counts")
  # plot(log2(datahist))
  # abline(lm(log2(datahist$counts) ~ log2(datahist$breaks)))
  
  xmax <- datahist$breaks[which.max(datahist$counts)]
  taildatax <- datahist$breaks[datahist$breaks >= xmax]
  taildatay <- datahist$counts[datahist$breaks >= xmax]
  taildata <- data.frame(breaks = taildatax, counts = taildatay)
  logtaildata <- log2(taildata)
  logtaildata$breaks[logtaildata$counts == -Inf] <- NA
  logtaildata <- na.omit(logtaildata)
  
  form <- lm(logtaildata$counts ~ logtaildata$breaks)
  form
  # plot(log2(datahist))
  # abline(lm(log2(taildata$counts) ~ log2(taildata$breaks)))
  predicted_df <- data.frame(counts = predict(form), breaks = logtaildata$breaks)
  
  logdata <- log2(datahist)
  loglindata <- log2(datalinear)
  coefficients[i] <- form$coefficients[2]
  plot <- ggplot(logdata, aes(breaks,counts)) + 
    scale_x_continuous("degree (k)",breaks = logdata$breaks, labels = 2^logdata$breaks) +
   scale_y_continuous("frequency of k",breaks = seq(0,ceiling(max(logdata$counts)),2),labels =2^(seq(0,ceiling(max(logdata$counts)),2)))  + 
    ggtitle(paste0("CN: ",numberofedges[i]," tau =",taucomp[i]," coef = ",round(coefficients[i],3),"\n kmax = ",k_max," C =",round(trans[i],4))) +
    
    geom_line(color='red',data = predicted_df) +
    geom_point(data = loglindata, color='orange') +
    geom_point(color='blue') 
  ggplotsCMcomp[[i]]  <- plot
  
}

numberofedges
ggplotsCMcomp[[i]]



##################################################################################
# Load Bayesian Network perm 1
# Create gridGraphsBN
##################################################################################
# IFCA: 
# filesinmap <- list.files("/media/catharina/ubuntu/outofoffice/hc_perm1", full.names = T)
# filesnames <- list.files("/media/catharina/ubuntu/outofoffice/hc_perm1")
# MACBOOK:
# filesinmap <- list.files("/Volumes/ubuntu/outofoffice/hc_perm1", full.names = T)
# filesnames <- list.files("/Volumes/ubuntu/outofoffice/hc_perm1")
# OCEANO
filesinmap <- list.files("../Data/Struct_learn/hciterations/perm1", full.names = T)
filesnames <- list.files("../Data/Struct_learn/hciterations/perm1")



filesnames <- gsub(".rda", "", filesnames)
filesnames

networklist <- list()
names <- c()

for (i in 1:length(filesinmap)){
  variablepos <- get(load(filesinmap[i]))
  networklist[[i]] <- variablepos
}
names(networklist) <- filesnames
rm(list = filesnames)

###########################
load("../Data/interim/tas_interim_10dnew.rda")

  permused <- 1
  for(j in c(permused)){
    pattern <- paste0("int_hc",permused,"_")
    hc_interim_list <- list.files(paste0("../Data/interim_struct/hciterations/perm",permused), full.names = T, pattern = pattern)
    hc_interim_names <- list.files(paste0("../Data/interim_struct/hciterations/perm",permused), pattern = pattern)
    hc_interim_names <- gsub(".rda", "", hc_interim_names)
    
    hc_interim_networks <- list()
    
    for (i in 1:length(hc_interim_list)){
      object <- get(load(hc_interim_list[i]))
      hc_interim_networks[[i]] <- object
    }
  }
  names(hc_interim_networks) <- hc_interim_names
  interimsizes <- sapply(hc_interim_networks,narcs)
  networks <- hc_interim_networks[order(interimsizes)]
  rm(hc_interim_networks)




networklist <- networks

# Convert the graphs to undirected igraphs
graphsNEL <- lapply(networklist, as.graphNEL) 
igraphsdir <- lapply(graphsNEL, igraph.from.graphNEL)
igraphsske<- lapply(igraphsdir, as.undirected)
edgelists <- lapply(igraphsske, E)
nedgeslists <- sapply(edgelists, length)
nedgeslists
indsame <- c()
# Sort the graphs from small to big
for (i in 1:length(nedgeslists)){
  int <- which(nedgeslists == sort(nedgeslists)[i])
  indsame[i] <- as.vector(int[1])
}
igraphsskesort <- igraphsske[indsame]
perm1sort <- networklist[indsame]

# Create the graphObject as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_interim_10dnew)
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObjects <- rep(list(graphObject),length(igraphsskesort))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- igraphsskesort[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsskesort[[i]])
}

gridGraphsBN <- graphObjects
names(gridGraphsBN) <- names(networklist)[indsame]
GraphsBN <- lapply(gridGraphsBN, function(x){x$graph})

# save(perm1sort, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
# save(gridGraphsBN, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsBN.rda")
#################################################################
# Create random erdös-Renyi graph
# Check other possibility in BNLEARN !
#################################################################
edgelistsBN <- lapply(GraphsBN, E)
numberofedgesBN <- lapply(edgelistsBN, length)
part2 <- numberofedgesCN[0:35]
part2 <- rev(part2)

edgesRenyi <- sapply(edgelistsBN, length)
edgesRenyi <- c(numberofedgesBN,part2)

renyigraphs <- lapply(edgesRenyi, sample_gnm, n = 648)
renyiObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
renyiObjects <- rep(list(renyiObject),length(renyigraphs))
for (i in 1:length(renyiObjects)){
  renyiObjects[[i]]$graph <- renyigraphs[[i]]
  renyiObjects[[i]]$adjacency <- as_adjacency_matrix(renyigraphs[[i]])
}

gridGraphsRenyi <- renyiObjects
names(gridGraphsRenyi) <- as.character(edgesRenyi)

# #################################################################
# # Create random regular graph
# # Check other possibility in BNLEARN !
# #################################################################
# edgesRegular <- sapply(edgelistsBN, length)
# degreesRegular <- sapply(edgesRegular, function(x){(2*x)/648})
# degreesRegular <- ceiling(degreesRegular)
# degreesRegular <- 1:25
# regulargraphs <- lapply(degreesRegular, sample_k_regular, no.of.nodes = 648)
# regulargraphs[[1]]
# regularObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
# regularObjects <- rep(list(regularObject),length(regulargraphs))
# for (i in 1:length(regularObjects)){
#   regularObjects[[i]]$graph <- regulargraphs[[i]]
#   regularObjects[[i]]$adjacency <- as_adjacency_matrix(regulargraphs[[i]])
# }
# 
# gridGraphsRegular <- regularObjects
# names(gridGraphsRegular) <- as.character(edgesRegular)
#################################################################
# Create deterministic lattices graph
# Check other possibility in BNLEARN !
#################################################################
g <- graph.lattice( c(18,36), circular = FALSE )
plot(g)
latticesgraphs <- lapply(1:25, connect.neighborhood, graph = g)
# # lapply(latticesgraphs, layout_on_sphere)
# layout(matrix(1:4, nrow=2, byrow=TRUE))
# sapply(graphsLattices[1:4], plot.igraph, vertex.label=NA)
edgelistsLattices <- lapply(latticesgraphs, E)
numberofedgesLattices <- lapply(edgelistsLattices, length)
names(latticesgraphs) <- as.character(numberofedgesLattices)

latticesObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
latticesObjects <- rep(list(latticesObject),length(latticesgraphs))
for (i in 1:length(latticesObjects)){
  latticesObjects[[i]]$graph <- latticesgraphs[[i]]
  latticesObjects[[i]]$adjacency <- as_adjacency_matrix(latticesgraphs[[i]])
}

gridGraphsLattices <- latticesObjects
names(gridGraphsLattices) <- as.character(numberofedgesLattices)

dev.off()

gridGraphsRenyi$`1398`$graph
gridGraphsLattices$`3620`$graph
gridGraphsBN$hc1_1600_1700i$graph

##############################################################
# Bayesian Network
##############################################################
grid <- tas_interim_10dnew
ggplotsCMcomp <- list()

ggplotsCMcomp <- list()
coefficients <- c()
numberofedges <- c()
trans <- c()

graphObj <-gridGraphsRenyi$`1795`
  graphObj <- gridGraphsBN$int_hc1_1700_1800i
  numberofedges[i] <- length(E(graphObj$graph))
  trans[i] <- transitivity(graphObj$graph, type = "global")
  degrees <- igraph::degree(graphObj$graph)
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(graphObj$graph)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  # max amount of counts
  #dist <-hist(degrees, breaks = k_max)
  #max(dist$counts)
  
  # bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
  # Calculate breaks: binlims
  nbins <- ceiling(log2(k_max))
  nbins
  binlims <- 2^(1:nbins)
  binlims <- append(binlims, 0, after = 0)
  binlims
  # Count amount values in bins with hist.
  # replace the counts by their normalized values (divided by binwidth)
  histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
  histdegrees2 <- histdegrees
  binwidth1 <- diff(histdegrees2$breaks, lag = 1)
  histdegrees2$counts <- histdegrees2$counts/binwidth1
  # Convert to adequate dataframe for ggplot.
  # Plot on log2 log2 scale. 
  histdegrees$counts
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks", "counts")
  # form <-lm(datahist$counts ~ datahist$breaks)
  # form$coefficients
  
  histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/2)*2,1), plot =FALSE)
  binwidth2 <- diff(histlinear$breaks, lag = 1)
  histlinear$counts <- histlinear$counts/binwidth2
  datalinear <- data.frame(histlinear$breaks[2:length(histlinear$breaks)], histlinear$counts)
  names(datalinear) <- c("breaks", "counts")
  # plot(log2(datahist))
  # abline(lm(log2(datahist$counts) ~ log2(datahist$breaks)))
  
  xmax <- datahist$breaks[which.max(datahist$counts)]
  taildatax <- datahist$breaks[datahist$breaks >= xmax]
  taildatay <- datahist$counts[datahist$breaks >= xmax]
  taildata <- data.frame(breaks = taildatax, counts = taildatay)
  logtaildata <- log2(taildata)
  logtaildata$breaks[logtaildata$counts == -Inf] <- NA
  logtaildata <- na.omit(logtaildata)
  
  form <- lm(logtaildata$counts ~ logtaildata$breaks)
  form
  # plot(log2(datahist))
  # abline(lm(log2(taildata$counts) ~ log2(taildata$breaks)))
  predicted_df <- data.frame(counts = predict(form), breaks = logtaildata$breaks)
  
  logdata <- log2(datahist)
  loglindata <- log2(datalinear)
  coefficients[i] <- form$coefficients[2]
  plot <- ggplot(logdata, aes(breaks,counts)) + 
    scale_x_continuous("degree (k)",breaks = logdata$breaks, labels = 2^logdata$breaks) +
    scale_y_continuous("frequency of k",breaks = seq(0,ceiling(max(logdata$counts)),2),labels =2^(seq(0,ceiling(max(logdata$counts)),2)))  + 
    ggtitle(paste0("BN: ",numberofedges[i]," coef = ",round(coefficients[i],3),"\n kmax = ",k_max," C =",round(trans[i],4))) +
    
    geom_line(color='red',data = predicted_df) +
    geom_point(data = loglindata, color='orange') +
    geom_point(color='blue') 
  ggplotsCMcomp[[i]]  <- plot
  


numberofedges
ggplotsCMcomp[[i]]



#################################################################
# Erdos renyi
#################################################################
graphObj <-gridGraphsRenyi$`1795`
graphObj <- gridGraphsBN$int_hc1_1700_1800i
ggplotsCMcomp <- list()
coefficients <- c()
numberofedges <- c()
trans <- c()

i <- 1



#graphObj <- graph_from_Grid(grid, th = taucomp[i], method = 'spearman')
numberofedges[i] <- length(E(graphObj$graph))
trans[i] <- transitivity(graphObj$graph, type = "global")
degrees <- igraph::degree(graphObj$graph)

# Calculate degrees from graph 
degrees <- igraph::degree(graphObj$graph)
degrees
# Calculate logarithmic bins:
# maximum degree 
k_max <- max(degrees)
# max amount of counts
#dist <-hist(degrees, breaks = k_max)
#max(dist$counts)

# bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
# Calculate breaks: binlims
nbins <- ceiling(log2(k_max))
nbins
binlims <- 2^(1:nbins)
binlims <- append(binlims, 0, after = 0)
binlims
# Count amount values in bins with hist.
# replace the counts by their normalized values (divided by binwidth)
histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
histdegrees2 <- histdegrees
binwidth1 <- diff(histdegrees2$breaks, lag = 1)
histdegrees2$counts <- histdegrees2$counts/binwidth1
# Convert to adequate dataframe for ggplot.
# Plot on log2 log2 scale. 
histdegrees$counts
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
names(datahist) <- c("breaks", "counts")
# form <-lm(datahist$counts ~ datahist$breaks)
# form$coefficients

histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/2)*2,1), plot =FALSE)
binwidth2 <- diff(histlinear$breaks, lag = 1)
histlinear$counts <- histlinear$counts/binwidth2
datalinear <- data.frame(histlinear$breaks[2:length(histlinear$breaks)], histlinear$counts)
names(datalinear) <- c("breaks", "counts")
# plot(log2(datahist))
# abline(lm(log2(datahist$counts) ~ log2(datahist$breaks)))

xmax <- datahist$breaks[which.max(datahist$counts)]
taildatax <- datahist$breaks[datahist$breaks >= xmax]
taildatay <- datahist$counts[datahist$breaks >= xmax]
taildata <- data.frame(breaks = taildatax, counts = taildatay)
logtaildata <- log2(taildata)
logtaildata$breaks[logtaildata$counts == -Inf] <- NA
logtaildata <- na.omit(logtaildata)

form <- lm(logtaildata$counts ~ logtaildata$breaks)
form
# plot(log2(datahist))
# abline(lm(log2(taildata$counts) ~ log2(taildata$breaks)))
predicted_df <- data.frame(counts = predict(form), breaks = logtaildata$breaks)

logdata <- log2(datahist)
loglindata <- log2(datalinear)
coefficients[i] <- form$coefficients[2]
plota <- ggplot(datalinear, aes(breaks,counts)) + 
  scale_x_continuous("degree (k)") +
  scale_y_continuous("frequency of k")  + 
  ggtitle(paste0("ER: ",numberofedges[i], "\nkmax = ",k_max," C =",round(trans[i],4))) +
  
  #geom_line(color='red',data = predicted_df) +
  #geom_point(data = datahist, color='blue') +
  geom_point(color='orange') 
plota
  
#################################################################
# regular lattice
#################################################################

ggplotsCMcomp <- list()
coefficients <- c()
numberofedges <- c()
trans <- c()

i <- 1

graphObj<-gridGraphsLattices$`3620`


#graphObj <- graph_from_Grid(grid, th = taucomp[i], method = 'spearman')
numberofedges[i] <- length(E(graphObj$graph))
trans[i] <- transitivity(graphObj$graph, type = "global")
degrees <- igraph::degree(graphObj$graph)

# Calculate degrees from graph 
degrees <- igraph::degree(graphObj$graph)
degrees
# Calculate logarithmic bins:
# maximum degree 
k_max <- max(degrees)
# max amount of counts
#dist <-hist(degrees, breaks = k_max)
#max(dist$counts)

# bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
# Calculate breaks: binlims
nbins <- ceiling(log2(k_max))
nbins
binlims <- 2^(1:nbins)
binlims <- append(binlims, 0, after = 0)
binlims
# Count amount values in bins with hist.
# replace the counts by their normalized values (divided by binwidth)
histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
histdegrees2 <- histdegrees
binwidth1 <- diff(histdegrees2$breaks, lag = 1)
histdegrees2$counts <- histdegrees2$counts/binwidth1
# Convert to adequate dataframe for ggplot.
# Plot on log2 log2 scale. 
histdegrees$counts
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
names(datahist) <- c("breaks", "counts")
# form <-lm(datahist$counts ~ datahist$breaks)
# form$coefficients

histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/2)*2,1), plot =FALSE)
binwidth2 <- diff(histlinear$breaks, lag = 1)
histlinear$counts <- histlinear$counts/binwidth2
datalinear <- data.frame(histlinear$breaks[2:length(histlinear$breaks)], histlinear$counts)
names(datalinear) <- c("breaks", "counts")
# plot(log2(datahist))
# abline(lm(log2(datahist$counts) ~ log2(datahist$breaks)))

xmax <- datahist$breaks[which.max(datahist$counts)]
taildatax <- datahist$breaks[datahist$breaks >= xmax]
taildatay <- datahist$counts[datahist$breaks >= xmax]
taildata <- data.frame(breaks = taildatax, counts = taildatay)
logtaildata <- log2(taildata)
logtaildata$breaks[logtaildata$counts == -Inf] <- NA
logtaildata <- na.omit(logtaildata)

form <- lm(logtaildata$counts ~ logtaildata$breaks)
form
# plot(log2(datahist))
# abline(lm(log2(taildata$counts) ~ log2(taildata$breaks)))
predicted_df <- data.frame(counts = predict(form), breaks = logtaildata$breaks)

logdata <- log2(datahist)
loglindata <- log2(datalinear)
coefficients[i] <- form$coefficients[2]
plotb <- ggplot(datalinear, aes(breaks,counts)) + xlab("degree (k)") +xlim(0,13)+

  scale_y_continuous("frequency of k")  + 
  ggtitle(paste0("RL: ",numberofedges[i], "\nkmax = ",k_max," C =",round(trans[i],4))) +
  
  #geom_line(color='red',data = predicted_df) +
  #geom_point(data = datahist, color='blue') +
  geom_point(color='orange') 
plotb

grid.arrange(plot,plota,plotb,ncol = 3)
  
##########################################
# Maak plot RL, ER, en BN
##########################################
graphObj4 <- graph_from_Grid(tas_interim_10dnew, th = 0.615, method = 'pearson')
numberofedges4 <- length(E(graphObj4$graph))

graphObj1 <-gridGraphsRenyi$`1795`
graphObj2 <- gridGraphsBN$int_hc1_1700_1800i
graphObj3 <- gridGraphsLattices$`1242`


numberofedges1<- length(E(graphObj1$graph))
numberofedges2<- length(E(graphObj2$graph))
numberofedges3<- length(E(graphObj3$graph))

trans1 <- transitivity(graphObj1$graph, type = "global")
trans2 <- transitivity(graphObj2$graph, type = "global")
trans3 <- transitivity(graphObj3$graph, type = "global")

degrees1 <- igraph::degree(graphObj1$graph)
degrees2 <- igraph::degree(graphObj2$graph)
degrees3 <- igraph::degree(graphObj3$graph)
degrees4 <- igraph::degree(graphObj4$graph)

k_max1 <- max(degrees1)
k_max2 <- max(degrees2)
k_max3 <- max(degrees3)
k_max4<- max(degrees4)

histdegrees1 <- hist(degrees1, breaks = 1, plot = TRUE)
histdegrees2 <- hist(degrees2, breaks = 1, plot = TRUE)
histdegrees3 <- hist(degrees3, breaks = 1, plot = TRUE)
histdegrees4 <- hist(degrees4,breaks =1,plot = TRUE)

histdegrees$counts
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
names(datahist) <- c("breaks", "counts")
# form <-lm(datahist$counts ~ datahist$breaks)
# form$coefficients

histlinear1 <- hist(degrees1, breaks = seq(0,ceiling(k_max1/2)*2,1), plot =TRUE)
datalinear1 <- data.frame(histlinear1$breaks[2:length(histlinear1$breaks)], histlinear1$counts)
names(datalinear1) <- c("breaks", "counts")
histlinear2 <- hist(degrees2, breaks = seq(0,ceiling(k_max2/2)*2,1), plot =TRUE)
datalinear2 <- data.frame(histlinear2$breaks[2:length(histlinear2$breaks)], histlinear2$counts)
names(datalinear2) <- c("breaks", "counts")
histlinear3 <- hist(degrees3, breaks = seq(0,ceiling(k_max3/2)*2,1), plot =TRUE)
datalinear3 <- data.frame(histlinear2$breaks[2:length(histlinear3$breaks)], histlinear3$counts)
names(datalinear3) <- c("breaks", "counts")
histlinear4 <- hist(degrees4, breaks = seq(0,ceiling(k_max4/2)*2,1), plot =TRUE)
datalinear4 <- data.frame(histlinear4$breaks[2:length(histlinear4$breaks)], histlinear4$counts)
names(datalinear4) <- c("breaks", "counts")


plota <- ggplot(datalinear1, aes(breaks,counts)) + 
  scale_x_continuous("degree (k)") +
  scale_y_continuous("frequency of k")  + 
  ggtitle(paste0("BN: ",numberofedges2, "\nkmax = ",k_max2," C =",round(trans2,4))) +
  
  #geom_line(color='red',data = predicted_df) +
  #geom_point(data = datahist, color='blue') +
  geom_point(color='orange') +
geom_point(data = datalinear2, color='blue') +
geom_point(data = datalinear4, color='red')
plota


