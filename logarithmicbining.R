##########################################################################################
# Powerlaw plot: logarithmic bining
##########################################################################################
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
########################################################################################## 
# working example Complex network
##########################################################################################
# obtain graph
example <- precgraph
example <- graph_from_Grid(tas_ncep_10d, th = 0.5, method = "pearson")
E(example$graph)
# Calculate degrees from graph 
degrees <- igraph::degree(example$graph)
degrees <- igraph::degree(example)
# Calculate logarithmic bins:
# maximum degree 
k_max <- max(degrees)
k_max
# bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
# Calculate breaks: binlims
nbins <- ceiling(log2(k_max))
nbins
binlims <- 2^(1:nbins)
binlims <- append(binlims, 0, after = 0)
binlims
# Count amount values in bins with hist.
# replace the counts by their normalized values (divided by binwidth)
histdegrees <- hist(degrees, breaks = binlims)
histdegrees2 <- histdegrees
histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
plot(histdegrees2, freq = TRUE)
plot(histdegrees, freq = TRUE)
# Convert to adequate dataframe for ggplot.
# Plot on log2 log2 scale. 
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)],histdegrees2$counts)
names(datahist) <- c("breaks","counts")
ggplot(datahist, aes(breaks,counts))+ scale_x_continuous(trans = "log2")+ scale_y_continuous(trans = "log2") + geom_point()
##########################################################################################
# working example complex network list
##########################################################################################
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x, digits = decimals))
}


plplotlist_logbin <- function(grid, tau){

  ggplots = list()
  for (i in 1:length(tau)){
    graphObj <- graph_from_Grid(grid, th = tau[i],method ='pearson')
    numberofedges <- length(E(graphObj$graph))
    
    degrees <- igraph::degree(graphObj$graph)
    k_max <- max(degrees)
    nbins <- ceiling(log2(k_max))
    binlims <- 2^(1:nbins)
    binlims <- append(binlims, 0, after = 0)
    
    histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
    histdegrees2 <- histdegrees
    histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
    
    # Convert to adequate dataframe for ggplot.
    # Plot on log2 log2 scale. 
    datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
    names(datahist) <- c("breaks","counts")
    plot <- ggplot(datahist, aes(breaks,counts)) + 
      scale_x_continuous(trans = "log2") + 
      scale_y_continuous(trans = "log2", labels = fmt_dcimals(2)) + 
      geom_point() +
      ggtitle(paste0("CN: ",numberofedges," tau = ",tau[i])) +
      xlab("degree k") +
      ylab("frequency k /\nbinwidth")
    ggplots[[i]] <- plot
  }
  
  return(ggplots)
}

ggplotsCM <- plplotlist_logbin(tas_ncep_10d, tau = c(#0.87,
                                                     0.8,
                                                     #0.7,
                                                     0.62,
                                                     0.52,
                                                     #0.49,
                                                     #0.41,
                                                     0.345))
ggplotsCM

n <- length(ggplotsCM)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggplotsCM, ncol=nCol))
##########################################################################################
# working example complex network align with bayes
##########################################################################################
fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x, digits = decimals))
}
tau <- c(0.87,
  0.8,
  0.7,
  0.62,
  0.52,
  0.49,
  0.41,
  0.345,
  0.32,
  0.31,
  0.3,
  0.25,
  0.2,
  0.15,
  0.1)

tau <- c(0.87,
         0.8,
         0.75,
         0.7,
         0.66,
         0.62,
         0.52,
         0.49,
         0.41,
         0.345,
         0.32,
         0.3,
         0.25,
         0.2,
         0.1,
         0.15,
         0.05,
         0)

grid <- tas_ncep_10d
ggplotsCM <- list()

ggplotsCM <- list()
coefficientsCM <- c()
numberofedgesCM <- c()
transCM <- c()

for (i in 1:length(tau)){
  graphObj <- graph_from_Grid(grid, th = tau[i],method = 'pearson')
  numberofedgesCM[i] <- length(E(graphObj$graph))
  transCM[i] <- transitivity(graphObj$graph, type = "global")
  degrees <- igraph::degree(graphObj$graph)
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(graphObj$graph)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks", "counts")
  # form <-lm(datahist$counts ~ datahist$breaks)
  # form$coefficients
  
  histlinear <- hist(degrees, breaks = seq(0,floor(k_max)+4,4), plot = FALSE)
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
  coefficientsCM[i] <- form$coefficients[2]
  plot <- ggplot(logdata, aes(breaks,counts)) + 
    # scale_x_continuous(trans = "log2") + 
    # scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    ggtitle(paste0("CN: ",numberofedgesCM[i]," tau =",tau[i]," coef = ",round(coefficientsCM[i],3),"\n kmax = ",k_max," trans =",round(transCM[i],4))) +
    
    geom_line(color='red',data = predicted_df) +
    geom_point(data = loglindata, color='orange') +
    geom_point(color='blue') 
  ggplotsCM[[i]]  <- plot
  
}


n <- length(ggplotsCM)
nCol <- floor(sqrt(n/2))
do.call("grid.arrange", c(ggplotsCM[1:9], ncol=nCol))
do.call("grid.arrange", c(ggplotsCM[10:18], ncol=nCol))

plot(numberofedgesCM,coefficientsCM, xlim = c(0,17000))
text(numberofedgesCM, coefficientsCM, labels = c(numberofedgesCM), cex= 0.7, pos = 3)
plot(numberofedgesCM,transCM, xlab = "|E|",ylab = "C", main = "Transitivity CM")
text(numberofedgesCM, transCM, labels = c(numberofedgesCM), cex= 0.7, pos = 3)

########################################################################################## 
# working example precision network
##########################################################################################
m <- glasso01$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
diag(adj.matrix) <- 0
adj.matrix[abs(adj.matrix) <= th ] <- 0
adj.matrix[abs(adj.matrix) > th ] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
# Calculate degrees from graph 
# degrees <- igraph::degree(covvargraph)
degrees <- igraph::degree(precgraph)
# Calculate logarithmic bins:
# maximum degree 
k_max <- max(degrees)
k_max
# bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
# Calculate breaks: binlims
nbins <- ceiling(log2(k_max))
nbins
binlims <- 2^(1:nbins)
binlims <- append(binlims, 0, after = 0)
binlims
# Count amount values in bins with hist.
# replace the counts by their normalized values (divided by binwidth)
histdegrees <- hist(degrees, breaks = binlims)
histdegrees2 <- histdegrees
histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
plot(histdegrees2, freq = TRUE)
plot(histdegrees, freq = TRUE)
# Convert to adequate dataframe for ggplot.
# Plot on log2 log2 scale. 
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)],histdegrees2$counts)
names(datahist) <- c("breaks","counts")
ggplot(datahist, aes(breaks,counts))+ scale_x_continuous(trans = "log2")+ scale_y_continuous(trans = "log2") + geom_point() + xlab("degree k") +
  ylab("frequency k /\n binwidth")



###################################################################################
# Precision network list glasso01 tresholds
###################################################################################

th <- c(0.11,0.067,0.055,0.0353)
i <- 1
ggplotsPN = list()
for (i in 1:length(th)){
  m <- glasso01$wi
  PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
  adj.matrix <- PartVar
  diag(adj.matrix) <- 0
  adj.matrix[abs(adj.matrix) <= th[i] ] <- 0
  adj.matrix[abs(adj.matrix) > th[i] ] <- 1
  precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  igraphske <- precgraph
  numberofedges <- length(E(igraphske))
  
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(igraphske)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
  
  # Convert to adequate dataframe for ggplot.
  # Plot on log2 log2 scale. 
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks","counts")
  plot <- ggplot(datahist, aes(breaks,counts)) + 
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    geom_point() +
    ggtitle(paste0("PN: ",numberofedges)) +
    xlab("degree k") +
    ylab("frequency k /\n binwidth")
  ggplotsPN[[i]] <- plot
}

n <- length(ggplotsPN)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggplotsPN, ncol=nCol))
##########################################################################################
# Bayesian network: working example logarithmic bining
##########################################################################################
exampleB <- hc_edges_loglik_10d_2800_3000i$networks
igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
igraphske<- as.undirected(igraphdir)
E(igraphske)
# Calculate degrees from graph 
degrees <- igraph::degree(igraphske)
degrees
# Calculate logarithmic bins:
# maximum degree 
k_max <- max(degrees)
k_max
# bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
# Calculate breaks: binlims
nbins <- ceiling(log2(k_max))
nbins
binlims <- 2^(1:nbins)
binlims <- append(binlims, 0, after = 0)
binlims
# Count amount values in bins with hist.
# replace the counts by their normalized values (divided by binwidth)
histdegrees <- hist(degrees, breaks = binlims)
histdegrees2 <- histdegrees
histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
plot(histdegrees2, freq = TRUE)
plot(histdegrees, freq = TRUE)
# Convert to adequate dataframe for ggplot.
# Plot on log2 log2 scale. 
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
names(datahist) <- c("breaks","counts")
ggplot(datahist, aes(breaks,counts))+ 
  scale_x_continuous(trans = "log2")+ 
  scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
  geom_point() +
  xlab("degree k") +
  ylab("frequency k /\n binwidth")










##########################################################################################
# Bayesian network: listplot logarithmic bining
##########################################################################################
#load tabu
for(j in c(1)){
  pattern <- paste0("tabu_",j,"_eBIC_g")
  filestabu1 <- list.files("../Data/Struct_learn/Tabu", full.names = T, pattern = pattern)
  filestabu1names <- list.files("../Data/Struct_learn/Tabu", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1$networks[[i]] <- variablepos$networks[[1]]
    names[i] <- names(variablepos$networks)
    data.frames[i,] <- as.matrix(variablepos$networkdata)
  }
  colnames(data.frames) <- names(variablepos$networkdata)
  data.frames <- as.data.frame(data.frames)
  variablelisttabu1$networkdata<- data.frames
  variablelisttabu1$networkdata
  
  names(variablelisttabu1$networks) <- names
  assign(paste0("tabu_",j,"_eBIC"),variablelisttabu1)
}
class(tabu_1_eBIC_g0[["networks"]][[1]])
hc_list <- list(#tabu_1_eBIC_g3,
                tabu_1_eBIC_g2,
                #tabu_1_eBIC_g1,
                #tabu_1_eBIC_g0.75,
                tabu_1_eBIC_g0.5,
                tabu_1_eBIC_g0.25,
                #tabu_1_eBIC_g0.1,
                tabu_1_eBIC_g0
                )
length(hc_list)


# select only the networksto calculate nedges
networks <- lapply(hc_list, function(m) m[["networks"]][[1]])
networks
nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]


ggplots = list()
for (i in 1:length(hc_list)){
  exampleB <- hc_list[[i]]$networks[[1]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedges <- length(E(igraphske))
 

# Calculate degrees from graph 
degrees <- igraph::degree(igraphske)
degrees
# Calculate logarithmic bins:
# maximum degree 
k_max <- max(degrees)
k_max
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
histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]

# Convert to adequate dataframe for ggplot.
# Plot on log2 log2 scale. 
datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
names(datahist) <- c("breaks","counts")
plot <- ggplot(datahist, aes(breaks,counts)) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
  geom_point() +
  ggtitle(paste0("BN: ",numberofedges)) +
  xlab("degree k") +
  ylab("frequency k /\n binwidth")
ggplots[[i]] <- plot
}
n <- length(ggplots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggplots, ncol=nCol))

dev.off()

##########################################################################################
# Bayesian network: listplot logarithmic bining with eBIC graphs
##########################################################################################

# select only the networksto calculate nedges
hc_list <- list(hc_1_eBIC_g1, 
                hc_1_eBIC_g0.75, 
                hc_1_eBIC_g0.5, 
                hc_1_eBIC_g0.25, 
                hc_1_eBIC_g0.1)
networks <- lapply(hc_list, function(m) m$networks)
networks2 <- lapply(networks, function(m) m[[1]])
nedgesnetworks <- as.character(sapply(networks2, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]

ggplots = list()
for (i in 1:length(hc_list)){
  exampleB <- hc_list[[i]]$networks[[1]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedges <- length(E(igraphske))
  
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(igraphske)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
  
  # Convert to adequate dataframe for ggplot.
  # Plot on log2 log2 scale. 
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks","counts")
  plot <- ggplot(datahist, aes(breaks,counts)) + 
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    geom_point() +
    ggtitle(paste0("BN: ",numberofedges)) +
    xlab("degree k") +
    ylab("frequency k /\n binwidth")
  ggplots[[i]] <- plot
}

n <- length(ggplots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggplots, ncol=nCol))

dev.off()
##########################################################################################
# Bayesian network: listplot logarithmic bining with eBIC graphs in list
##########################################################################################

# select only the networksto calculate nedges

hc_list <- hc_2_eBIC$networks
networks <- lapply(hc_list, function(m) m$networks)
networks2 <- lapply(networks, function(m) m[[1]])
nedgesnetworks <- as.character(sapply(hc_list, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]

ggplots = list()
for (i in 1:length(hc_list)){
  exampleB <- hc_list[[i]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedges <- length(E(igraphske))
  
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(igraphske)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
  
  # Convert to adequate dataframe for ggplot.
  # Plot on log2 log2 scale. 
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks","counts")
  plot <- ggplot(datahist, aes(breaks,counts)) + 
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    geom_point() +
    ggtitle(paste0("BN: ",numberofedges)) +
    xlab("degree k") +
    ylab("frequency k /\n binwidth")
  ggplots[[i]] <- plot
}

n <- length(ggplots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggplots, ncol=nCol))

dev.off()
##########################################################################################
# Bayesian network: listplot logarithmic bining with all hc iteration graphs
##########################################################################################

# select only the networks to calculate nedges
networks3 <- iterationnetworks[1:10]
length(networks3)
nedgesnetworks <- as.character(sapply(networks3, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]
i <- 1
ggplots = list()
for (i in 1:length(networks3)){
  exampleB <- iterationnetworks[[i]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedges <- length(E(igraphske))
  
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(igraphske)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  histdegrees2$counts <- histdegrees2$counts/histdegrees2$breaks[2:length(histdegrees2$breaks)]
  
  # Convert to adequate dataframe for ggplot.
  # Plot on log2 log2 scale. 
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks","counts")
  plot <- ggplot(datahist, aes(breaks,counts)) + 
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    geom_point() +
    ggtitle(paste0("BN: ",numberofedges)) +
    xlab("degree k") +
    ylab("frequency k /\n binwidth")
  ggplots[[i]] <- plot
}

n <- length(ggplots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(ggplots, ncol=nCol))
warnings()
dev.off()
###########################################################################
# Plots and coefficients tails HC ITERATION.
###########################################################################


networks3 <- iterationnetworks

powerlawplots <- list()
coefficientsBN <- c()
numberofedgesBN <- c()
transBN <- c()
i <- 18
for (i in 1:length(networks3)){
  exampleB <- iterationnetworks[[i]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedgesBN[i] <- length(E(igraphske))
  transBN[i] <- transitivity(igraphske, type = "global")
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(igraphske)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks", "counts")
  # form <-lm(datahist$counts ~ datahist$breaks)
  # form$coefficients
  
  histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/2)*2,2), plot = FALSE)
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
  coefficientsBN[i] <- form$coefficients[2]
  plot <- ggplot(logdata, aes(breaks,counts)) + 
    # scale_x_continuous(labels = c(2,4,8,16,32,64)) + 
    # #scale_y_continuous(labels = 2^datahist$counts) + 
    ggtitle(paste0("BN: ",numberofedgesBN[i]," coef = ",round(coefficientsBN[i],3),"\n kmax = ",k_max," trans =",round(transBN[i],4))) +
    geom_line(color='red',data = predicted_df) +
    geom_point(data = loglindata, color='orange') +
    geom_point(color='blue') 
  powerlawplots[[i]]  <- plot
  
}

do.call(grid.arrange, c(powerlawplots[1:9]))
do.call(grid.arrange, c(powerlawplots[10:18]))
do.call(grid.arrange, c(powerlawplots[19:27]))

do.call(grid.arrange, c(powerlawplots[1:20]))
do.call(grid.arrange, c(powerlawplots[21:40]))

plot(numberofedgesBN[4:length(numberofedgesBN)],coefficientsBN[4:length(numberofedgesBN)])
plot(numberofedgesBN,coefficientsBN, col = "blue", xlab = "|E|", ylab = "slope", main = "coefficients BN")
text(numberofedgesBN, coefficientsBN, labels = c(numberofedgesBN), cex= 0.7, pos = 3)
plot(numberofedgesBN,transBN, xlab = "|E|", ylab = "C", main = "Transitivity BN", col = "blue")
text(numberofedgesBN, transBN, labels = c(numberofedgesBN), cex= 0.7, pos = 3)
###########################################################################
# Plots and coefficients tails HC ITERATION. PROBEREN
###########################################################################
# networks3 <- iterationnetworks
# 
# powerlawplots <- list()
# coefficients <- c()
# numberofedges <- c()
# i <- 5
# for (i in 1:length(networks3)){
#   exampleB <- iterationnetworks[[i]]
#   igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
#   igraphske<- as.undirected(igraphdir)
#   
#   
#   numberofedges[i] <- length(E(igraphske))
#   
#   
#   # Calculate degrees from graph 
#   degrees <- igraph::degree(igraphske)
#   degrees
#   # Calculate logarithmic bins:
#   # maximum degree 
#   k_max <- max(degrees)
#   k_max
#   # bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
#   # Calculate breaks: binlims
#   nbins <- ceiling(log2(k_max))
#   nbins
#   binlims <- 2^(1:nbins)
#   binlims <- append(binlims, 0, after = 0)
#   binlims
#   # Count amount values in bins with hist.
#   # replace the counts by their normalized values (divided by binwidth)
#   histdegrees <- hist(degrees, breaks = binlims, plot = FALSE)
#   histdegrees2 <- histdegrees
#   histdegrees2$counts <- histdegrees$counts/c(2,histdegrees$breaks[2:(length(histdegrees$breaks)-1)])
#   
#   # Convert to adequate dataframe for ggplot.
#   # Plot on log2 log2 scale. 
#   datahistsin <- data.frame(histdegrees$breaks[2:length(histdegrees$breaks)], histdegrees$counts)
#   datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
#   names(datahistsin) <- names(datahist) <- c("breaks", "counts")
#   # form <-lm(datahist$counts ~ datahist$breaks)
#   # form$coefficients
#   
#   # plot(log2(datahist))
#   # abline(lm(log2(datahist$counts) ~ log2(datahist$breaks)))
#   
#   xmax <- datahist$breaks[which.max(datahist$counts)]
#   taildatax <- datahist$breaks[datahist$breaks >= xmax]
#   taildatay <- datahist$counts[datahist$breaks >= xmax]
#   taildata <- data.frame(breaks = taildatax, counts = taildatay)
#   
#   form <- lm(log2(taildata$counts) ~ log2(taildata$breaks))
#   # plot(log2(datahist))
#   # abline(lm(log2(taildata$counts) ~ log2(taildata$breaks)))
#   predicted_df <- data.frame(counts = predict(form), breaks = log2(taildata$breaks))
#   
#   logdata <- log2(datahist)
#   coefficients[i] <- form$coefficients[2]
#   plot <- ggplot(logdata, aes(breaks,counts), color = "blue") + 
#     # scale_x_continuous(trans = "log2") + 
#     # scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
#     ggtitle(paste0("BN: ",numberofedges[i],"",coefficients[i],"\n kmax = ",k_max)) +
#     geom_point(color='blue') +
#     geom_point(data = log2(datahistsin))+
#     geom_line(color='red',data = predicted_df)
#   powerlawplots[[i]]  <- plot
#   
# }
# 
# powerlawplots
# do.call(grid.arrange, c(powerlawplots))
# 
# plot(numberofedges,coefficients)
###########################################################################
# Plots and coefficients tails HC ITERATION LOG 3
###########################################################################


networks3 <- iterationnetworks

powerlawplots <- list()
coefficients <- c()
numberofedges <- c()
trans <- c()
i <- 18
for (i in 1:length(networks3)){
  exampleB <- iterationnetworks[[i]]
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  
  
  numberofedges[i] <- length(E(igraphske))
  trans[i] <- transitivity(igraphske, type = "global")
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(igraphske)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
  # bining as in Newman: Every grid 2 times as big as previous. Amount of bins: Where does k_max fall?
  # Calculate breaks: binlims
  nbins <- ceiling(log(k_max, base = 3))
  nbins
  binlims <- 3^(1:nbins)
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
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks", "counts")
  # form <-lm(datahist$counts ~ datahist$breaks)
  # form$coefficients
  
  histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/3)*3,3), plot = FALSE)
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
  logtaildata <- log(taildata, base = 3)
  logtaildata$breaks[logtaildata$counts == -Inf] <- NA
  logtaildata <- na.omit(logtaildata)
  
  form <- lm(logtaildata$counts ~ logtaildata$breaks)
  form
  # plot(log2(datahist))
  # abline(lm(log2(taildata$counts) ~ log2(taildata$breaks)))
  predicted_df <- data.frame(counts = predict(form), breaks = logtaildata$breaks)
  
  logdata <- log(datahist, base = 3)
  loglindata <- log(datalinear, base = 3)
  coefficients[i] <- form$coefficients[2]
  plot <- ggplot(logdata, aes(breaks,counts)) + 
    # scale_x_continuous(trans = "log2") + 
    # scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    ggtitle(paste0("BN: ",numberofedges[i]," coef = ",round(coefficients[i],3),"\n kmax = ",k_max," trans =",round(trans[i],4))) +
    geom_point(color='blue') +
    geom_line(color='red',data = predicted_df) +
    geom_point(data = loglindata, color='orange') 
  powerlawplots[[i]]  <- plot
  
}

powerlawplots
do.call(grid.arrange, c(powerlawplots[1:20]))
do.call(grid.arrange, c(powerlawplots[21:40]))

plot(numberofedges[4:length(numberofedges)],coefficients[4:length(numberofedges)])
plot(numberofedges,coefficients, col = "blue")
text(numberofedges, coefficients, labels = c(numberofedges), cex= 0.7, pos = 3)
plot(numberofedges,trans)
text(numberofedges, trans, labels = c(numberofedges), cex= 0.7, pos = 3)


###############################################################################################
# Combining powerlaw figures bayesian and complex networks
###############################################################################################
# what <- closenesses
# multis = list()
# for (i in 1:length(what)){
#   plot <- plotClimatology(what[[i]], main = list(paste0("BN: ",nedgesnetworks[i]), cex = 0.4),
#                           backdrop.theme = "coastline", 
#                           colorkey = list(width = 0.6, lables = list(cex = 0.4)))
#   multis[[i]] <- plot
# }
# 
# whatCM <- closenessesCM
# multisCM = list()
# for (i in 1:length(whatCM)){
#   plot <- plotClimatology(whatCM[[i]], main = list(paste0("CN: ",nedgesnetworksCM[i]), cex = 0.4, pos = 0.25),
#                           backdrop.theme = "coastline",
#                           colorkey = list(width = 0.6, lables = list(cex = 0.4)))
#   multisCM[[i]] <- plot
# }

# These with spearman in results 3
taucomp <- c(0.8,
         0.7,
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
taucomp <- c(0.50)



grid <- tas_ncep_10d
ggplotsCMcomp <- list()

ggplotsCMcomp <- list()
coefficients <- c()
numberofedges <- c()
trans <- c()

i <- 1
for (i in 1:length(taucomp)){
  graphObj <- graph_from_Grid(grid, th = taucomp[i], method = 'pearson')
  numberofedges[i] <- length(E(graphObj$graph))
  trans[i] <- transitivity(graphObj$graph, type = "global")
  degrees <- igraph::degree(graphObj$graph)
  
  # Calculate degrees from graph 
  degrees <- igraph::degree(graphObj$graph)
  degrees
  # Calculate logarithmic bins:
  # maximum degree 
  k_max <- max(degrees)
  k_max
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
  datahist <- data.frame(histdegrees2$breaks[2:length(histdegrees2$breaks)], histdegrees2$counts)
  names(datahist) <- c("breaks", "counts")
  # form <-lm(datahist$counts ~ datahist$breaks)
  # form$coefficients
  
  histlinear <- hist(degrees, breaks = seq(0,ceiling(k_max/2)*2,2), plot = FALSE)
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
    #scale_x_continuous(trans = "log2") + 
    #scale_y_continuous(trans = "log2",labels = fmt_dcimals(2)) + 
    ggtitle(paste0("CN: ",numberofedges[i]," tau =",taucomp[i]," coef = ",round(coefficients[i],3),"\n kmax = ",k_max," trans =",round(trans[i],4))) +
    
    geom_line(color='red',data = predicted_df) +
    geom_point(data = loglindata, color='orange') +
    geom_point(color='blue') 
  ggplotsCMcomp[[i]]  <- plot
  
}

mat <- rbind(c(1,8),
      c(2,9),
      c(3,10),
      c(4,11),
      c(5,12),
      c(6,13),
      c(7,14),
      c(0,15),
      c(0,16),
      c(0,17),
      c(0,18),
      c(0,19))

do.call(grid.arrange, c(powerlawplots[c(4,5,6,7,8,12,16)], ncol = 1))
do.call(grid.arrange, c(ggplotsCMcomp[1:6], ncol = 1))
do.call(grid.arrange, c(ggplotsCMcomp[7:12], ncol = 1))

ggplotsCMcomp[[1]]

list  <- append(powerlawplots,ggplotsCM)
list <- append(list, ggplotsPN)

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/powerlaw/HCandBNlogbin",firstel,"_",lastel,"_",length(list)/2,".pdf")
plotname
#pdf(plotname, height = 8, width =5) # para 8
pdf(plotname)
do.call(grid.arrange, c(list, nrow = length(ggplots), ncol = 3, as.table = FALSE))
dev.off()

plot(numberofedgesCM,transCM, xlab = "|E|", ylab = "C", main = "Transitivity", ylim = c(0,1), xlim = c(0,11000))
text(numberofedgesCM, transCM, labels = c(numberofedgesCM), cex= 0.7, pos = 3)
points(numberofedgesBN,transBN, col = "blue")

index <- seq(0,floor(length(numberofedgesBN)/5)*5,5)
labels <- rep(NA,length(numberofedgesBN))
labels[index] <- numberofedgesBN[index]
labels
text(numberofedgesBN, transBN, labels = labels, cex= 0.7, pos = 3)
legend("topright",legend = c("CN","BN"), pch = 1,col = c("black","blue"))

plot(numberofedgesCM,coefficientsCM, xlab = "|E|", ylab = "slope", main = "Powerlaw coefficients", ylim = c(-6,-1), xlim = c(0,11000))
text(numberofedgesCM, coefficientsCM, labels = c(numberofedgesCM), cex= 0.7, pos = 3)
points(numberofedgesBN,coefficientsBN, col = "blue")

index <- seq(0,floor(length(numberofedgesBN)/4)*4,4)
labels <- rep(NA,length(numberofedgesBN))
labels[index] <- numberofedgesBN[index]
labels
text(numberofedgesBN, coefficientsBN, labels = labels, cex= 0.7, pos = 3)
legend("topright",legend = c("CN","BN"), pch = 1,col = c("black","blue"))


#######################################################
# Bayesian network: Powerlaw old method
#######################################################
  exampleB <- hc_edges_loglik_10d_3000_3200i$networks
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  ddistB <- degree_distribution(igraphske)
  degree <- igraph::degree(igraphske)
  hist <- hist(degree)
  hist
  
  par(mfrow = c(1,3))
  
  graph1 <- plot(1:length(ddistB),
                 ddistB, 
                 main = paste0("degree_distribution"), 
                 type = "p", 
                 xlab = "degree K" , 
                 ylab = "p(K)")
  graph2 <- plot(log(1:length(ddistB)),
                 log(ddistB), main = paste0("degree_distribution"), 
                 type = "p", 
                 xlab = "log (degree K)" , 
                 ylab = "log (p(K))")
  graph3 <- plot(log2(hist$mids), log2(hist$counts))
  dev.off()
#######################################################################################
# List of old method powerlaw plots Bayesian network
#####################################################################################
  list <- list(hc_edges_loglik_10d_200_400i,
                hc_edges_loglik_10d_600_800i,
                hc_edges_loglik_10d_1400_1600i,
                hc_edges_loglik_10d_3000_3200i,
                hc_edges_loglik_10d_6200_6400i)
  length(list)
  
  list[[1]]$networks
  par(mar=c(4,4,2,2))
  par(mfrow = c(ceiling(length(list)/2), 4))
  
  for (i in 1:length(list)){
    exampleB <- list[[i]]$networks
    igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
    igraphske<- as.undirected(igraphdir)
    ddistB <- degree_distribution(igraphske)
    
    numberofedges <- length(E(igraphske))
    plot(1:length(ddistB),
         ddistB, 
         main = paste0("degree_distribution\n |E| = ",numberofedges), 
         type = "l", 
         xlab = "degree K" , 
         ylab = "p(K)")
    plot(log(1:length(ddistB)),
         log(ddistB), main = paste0("degree_distribution loglog\n |E| = ",numberofedges), 
         type = "l", 
         xlab = "log (degree K)" , 
         ylab = "log (p(K))")
  }

rm(graph10d,graphObj,graphObeject,graphObjects,graphs10d,grid)
rm(hc_list)
#######################################################################################
# List of old method powerlaw plots Bayesian network
#####################################################################################
list <- list(hc_edges_loglik_10d_200_400i,
             hc_edges_loglik_10d_600_800i,
             hc_edges_loglik_10d_1400_1600i,
             hc_edges_loglik_10d_3000_3200i,
             hc_edges_loglik_10d_6200_6400i)
length(list)

list[[1]]$networks
par(mar=c(4,4,2,2))
par(mfrow = c(ceiling(length(list)/2), 4))

for (i in 1:length(list)){
  exampleB <- list[[i]]$networks
  igraphdir <- igraph.from.graphNEL(as.graphNEL(exampleB))
  igraphske<- as.undirected(igraphdir)
  ddistB <- degree_distribution(igraphske)
  
  numberofedges <- length(E(igraphske))
  plot(1:length(ddistB),
       ddistB, 
       main = paste0("degree_distribution\n |E| = ",numberofedges), 
       type = "l", 
       xlab = "degree K" , 
       ylab = "p(K)")
  plot(log(1:length(ddistB)),
       log(ddistB), main = paste0("degree_distribution loglog\n |E| = ",numberofedges), 
       type = "l", 
       xlab = "log (degree K)" , 
       ylab = "log (p(K))")
}

rm(graph10d,graphObj,graphObeject,graphObjects,graphs10d,grid)
rm(hc_list)