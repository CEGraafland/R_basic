##############################################################################
# Complex measures on Bayesian networks and correlation networks
# Para la tesis
# Interim / NCEP?NCAR / JRA55
##############################################################################
setwd("~/data/Untitled/Trabajo/R_practice/R/")
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")
rm(list = ls())
library(gridExtra)
library(gridGraphics)
library(grid)
library("bnlearn")
library("visualizeR")
source("Functions/BasicNetworkFunctions.R")
# source("Functions/HillclimbingFunctions2.R")
source("Functions/CN_ConstructionandMeasuresFunctions.R")
load("../Data/tas_ncep_10d.rda")
load("../Data/interim/tas_interim_10dnew.rda")
load("../Data/Struct_learn/backpermutations.rda")
rean <- "interim"
if(rean == "ncep"){rean.grid <- tas_ncep_10d} else if(rean == "interim"){rean.grid <- tas_interim_10dnew}
############################################################################################
# Correlation networks: 
# Complex measures with lapply and measure2clim
# used in Resumen 2
############################################################################################
# 1 3 6 7 +4
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
         0.2,
         0.1,
         0.05,
         0.0)
tau <- c(0.60,0.55,0.34,0.31)
tau <- c(0.5,0.42,0.35,0.31) # NCEP 
tau <- c(0.6,0.49,0.40,0.32) # Pearson Betweenness interim
graphs10d <- lapply(tau, graph_from_Grid, grid = rean.grid, subind = NULL, method = 'pearson')
measures10d <- lapply(graphs10d, graph2measure)
graphssolo  <- lapply(graphs10d, function(m) m$graph)
edgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(edgesnetworksCM, length)
nedgesnetworksCM

degreesCM <- lapply(measures10d, measure2clim, what = "degree", ref.grid = rean.grid)
closenessesCM <- lapply(measures10d, measure2clim, what = "closeness", ref.grid = rean.grid)
betweennessCM <- lapply(measures10d, measure2clim, what = "betweenness", ref.grid = rean.grid)
awconnectivitiesCM <- lapply(measures10d, measure2clim, what = "awconnectivity", ref.grid = rean.grid)
localclusteringsCM <- lapply(measures10d, measure2clim, what = "localclustering", ref.grid = rean.grid)
#######################################################################################
# Working example measures complex networks for list of bayesian networks * 
# Measures are applicated on non rescaled data (only anomaly)
# (propagation is applicated on rescaled data)
# TMS_as_list is therefore adequate: Graph from non rescaled data in list with non rescaled
# data from function.
#######################################################################################
####################################################################################
# load HC interim iteration data permutation X and make list
####################################################################################
if (rean == "interim"){
permused <- 3
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
} else if (rean=="ncep"){
  permused <- 1
load("../Data/networks_hcnetworks10d.rda")
hc_list <- list(hc_edges_loglik_10d_200_400i,
                hc_edges_loglik_10d_600_800i,
                hc_edges_loglik_10d_800_1000i,
                hc_edges_loglik_10d_1000_1200i,
                hc_edges_loglik_10d_1200_1400i,
                hc_edges_loglik_10d_1400_1600i,
                hc_edges_loglik_10d_1600_1800i,
                hc_edges_loglik_10d_1800_2000i,
                hc_edges_loglik_10d_2400_2600i,
                hc_edges_loglik_10d_3000_3200i,
                hc_edges_loglik_10d_4800_5000i,
                hc_edges_loglik_10d_8600_8695i)
# select only the networks to calculate nedges
networks <- lapply(hc_list, function(m) m[["networks"]])
}


nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]

# Convert the graphs to undirected igraphs
igraphsdir <- lapply(networks, as.graphNEL) 
igraphsdir <- lapply(igraphsdir, igraph.from.graphNEL)
igraphsske<- lapply(igraphsdir, as.undirected)
# Create the graphObject as would have been obtained by graph_from_Grid
graphObject <- TimeCoordsAnom_from_Grid_aslist(rean.grid)
graphObjects <- rep(list(graphObject),length(networks))
for (i in 1:length(graphObjects)){
  graphObjects[[i]]$graph <- igraphsske[[i]]
  graphObjects[[i]]$adjacency <- as_adjacency_matrix(igraphsske[[i]])
}
# apply measure2clim and plotClimatology
# on list with graphObjects
measureslistBay <- lapply(graphObjects,graph2measure)

degrees <- lapply(measureslistBay, function(x) quantity2clim(x$degree, what = "degree", ref.grid = rean.grid,backperm = backpermutations[[permused]]))
closenesses <- lapply(measureslistBay, function(x) quantity2clim(x$closeness, what = "closeness", ref.grid = rean.grid,backperm = backpermutations[[permused]]))
betweenness <- lapply(measureslistBay, function(x) quantity2clim(x$betweenness, what = "betweenness", ref.grid = rean.grid,backperm = backpermutations[[permused]]))
awconnectivities <- lapply(measureslistBay, function(x) quantity2clim(x$awconnectivity, what = "awconnectivity", ref.grid = rean.grid,backperm = backpermutations[[permused]]))
localclusterings <- lapply(measureslistBay, function(x) quantity2clim(x$localclustering, what = "localclustering", ref.grid = rean.grid,backperm = backpermutations[[permused]]))

numberofedgesBN <- sapply(igraphsske, function(x)length(E(x)))


#bn <- 18
#whats <- c("degree","closeness","betweenness","awconnectivity","localclustering")
#clims <- lapply(whats, function(x) lapply(measureslistBay[[bn]],measure2clim,what = x,ref.grid =rean.grid))
#measureslistBay[[bn]]$
#clims[[1]]
#clims[['degree']]
#movs <- lapply(clims[whats], function(x) movingmedias(x[[bn]]))
#measurenames <- 
#clims

#################################################################
# Create random erdös-Renyi graph
# Check other possibility in BNLEARN !
#################################################################
edgesRenyi <- numberofedgesBN

renyigraphs <- lapply(numberofedgesBN, sample_gnm, n = 648)
renyiObject <- TimeCoordsAnom_from_Grid_aslist(rean.grid)
renyiObjects <- rep(list(renyiObject),length(renyigraphs))
for (i in 1:length(renyiObjects)){
  renyiObjects[[i]]$graph <- renyigraphs[[i]]
  renyiObjects[[i]]$adjacency <- as_adjacency_matrix(renyigraphs[[i]])
}

gridGraphsRenyi <- renyiObjects
names(gridGraphsRenyi) <- as.character(edgesRenyi)

measureslistER <- lapply(gridGraphsRenyi,graph2measure)
degreesER <- lapply(measureslistER, measure2clim,what = "degree", ref.grid = rean.grid)

###############################################################################################
# combining measures bayesian and complex networks and precision networks
###############################################################################################
##################################################################
# All measures CN
##################################################################
CN <- 2

plot1e<- spatialPlot(closenessesCM[[CN]], main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8, pos = 0.25),
                     backdrop.theme = "coastline", rev.colors = TRUE,
                     colorkey = list(width = 0.6, lables = list(cex = 0.5)))

movclos1.cm<-movingmedias(measures10d[[CN]]$closeness)
movclos1matrixmean.cm <- apply(movclos1.cm, MARGIN = 3, FUN = mean,na.rm = TRUE)
difclos.cm <- (max(movclos1matrixmean.cm) - min(movclos1matrixmean.cm))/14
minclos.cm <- min(movclos1matrixmean.cm) - difclos.cm
maxclos.cm <- max(movclos1matrixmean.cm) + difclos.cm

memClim.cm <- quantity2clim(movclos1matrixmean.cm, "clos mov", rean.grid)
plot1f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), 
                      rev.colors = TRUE,  at = seq(minclos.cm,maxclos.cm,difclos.cm),
                      colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1f

plot2e <- spatialPlot(betweennessCM[[CN]], main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8, pos = 0.25),
                      backdrop.theme = "coastline", lonCenter = 0, color.theme = "Reds", set.max = 10,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,10,1))

movbet1.cm<-movingmedias(measures10d[[CN]]$betweenness)
movbet1matrixmean.cm <- apply(movbet1.cm, MARGIN = 3, FUN = mean,na.rm = TRUE)
difbet.cm <- (max(movbet1matrixmean.cm) - min(movbet1matrixmean.cm))/14
minbet.cm <- min(movbet1matrixmean.cm) - difbet.cm
maxbet.cm <- max(movbet1matrixmean.cm) + difbet.cm
memClim.cm <- quantity2clim(movbet1matrixmean.cm, "bet mov", rean.grid)
plot2f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                      set.max = 10, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), 
                      color.theme = "Reds",  at = seq(minbet.cm,maxbet.cm,difbet.cm),
                      colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot2f


plot3e <- spatialPlot(grid = awconnectivitiesCM[[CN]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = 0.020, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8),
                      color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))

movawc1.cm<-movingmedias(measures10d[[CN]]$awconnectivity)
movawc1matrixmean.cm <- apply(movawc1.cm, MARGIN = 3, FUN = mean,na.rm=TRUE)

memClim.cm <- quantity2clim(movawc1matrixmean.cm, "awc mov", rean.grid)
plot3f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot3f


plot4e <- spatialPlot(grid = localclusteringsCM[[CN]], backdrop.theme = "coastline", set.min = NULL,# at = seq(0,1,0.1),
                      set.max = NULL, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8),
                      color.theme = "Blues",colorkey = list(width = 0.6, lables = list(cex = 0.6)))


movlcl1.cm<-movingmedias(measures10d[[CN]]$localclustering)
movlcl1matrixmean.cm <- apply(movlcl1.cm, MARGIN = 3, FUN = mean, na.rm = TRUE)
diflcl.cm <- (max(movlcl1matrixmean.cm) - 0)/15
minlcl.cm <- 0
maxlcl.cm <- max(movlcl1matrixmean.cm) + diflcl.cm

memClim.cm <- quantity2clim(movlcl1matrixmean.cm, "lcl mov", rean.grid)
plot4f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), 
                      color.theme = "Blues", at = seq(minlcl.cm,maxlcl.cm,diflcl.cm),
                      colorkey = list(width = 0.6, lables = list(cex = 0.6)))

plot4f

plot5e <- spatialPlot(grid = degreesCM[[CN]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter = 0, main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8),
                      color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))

movdeg1.cm<-movingmedias(measures10d[[CN]]$degree)
movdeg1matrixmean.cm <- apply(movdeg1.cm, MARGIN = 3, FUN = mean,na.rm = TRUE)

memClim.cm <- quantity2clim(movdeg1matrixmean.cm, "deg mov", rean.grid)
plot5f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot5f



grid.arrange(plot1e,plot2e,plot3e,plot4e,plot5e)
grid.arrange(plot1f,plot2f,plot3f,plot4f,plot5f)


plotname <- paste0("../plots/movingmedias_thesis/movs_CN",nedgesnetworksCM[CN],"_",rean,"_all.pdf")
pdf(plotname, width = 25, height = 7)

grid.arrange(plot1e,plot2e,plot3e,plot4e,plot5e,plot1f,plot2f,plot3f,plot4f,plot5f,nrow =2)
dev.off()
##################################################################
# All measures BN
##################################################################
bn <-30
plot1a <- spatialPlot(closenesses[[bn]], main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                    backdrop.theme = "coastline", rev.colors = TRUE,
                    colorkey = list(width = 0.6, lables = list(cex = 0.5)))

movclos1<-movingmedias(measureslistBay[[bn]]$closeness[backpermutations[[permused]]])
movclos1matrixmean <- apply(movclos1, MARGIN = 3, FUN = mean,na.rm = TRUE)
memClim <- quantity2clim(movclos1matrixmean, "clos mov", rean.grid)
plot1b <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]), cex = 0.8), 
                     rev.colors = TRUE, at = seq(minclos.cm,maxclos.cm,difclos.cm),
                     colorkey = list(width = 0.6, lables = list(cex = 0.6)))
plot1b

plot2a <- spatialPlot(betweenness[[bn]], main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                      backdrop.theme = "coastline", lonCenter = 0, color.theme = "Reds", set.max = 10,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,12,2))
plot2a

                      
movbet1<-movingmedias(measureslistBay[[bn]]$betweenness[backpermutations[[permused]]])
movbet1matrixmean <- apply(movbet1, MARGIN = 3, FUN = mean, na.rm = TRUE)
memClim <- quantity2clim(movbet1matrixmean, "bet mov", rean.grid)
plot2b <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                      set.max = 8.5, lonCenter =0, color.theme ="Reds", at = seq(minbet.cm,maxbet.cm,difbet.cm),
                      main = list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]), cex = 0.8),
                      colorkey = list(width = 0.6, lables = list(cex = 0.6)))
plot2b

plot3a <- spatialPlot(grid = awconnectivities[[bn]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = 0.016, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                      color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
movawc1<-movingmedias(measureslistBay[[bn]]$awconnectivity[backpermutations[[permused]]])
movawc1matrixmean <- apply(movawc1, MARGIN = 3, FUN = mean, na.rm = TRUE)

memClim <- quantity2clim(movawc1matrixmean, "awc mov", rean.grid)
plot3b <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]), cex = 0.8), color.theme = "Greens",colorkey = list(width = 0.6, lables = list(cex = 0.6)))

plot3b

plot4a <-spatialPlot(grid = localclusterings[[bn]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                      color.theme = "Blues",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
  
movlcl1<-movingmedias(measureslistBay[[bn]]$localclustering[backpermutations[[permused]]])
movlcl1matrixmean <- apply(movlcl1, MARGIN = 3, FUN = mean, na.rm = TRUE)

memClim <- quantity2clim(movlcl1matrixmean, "lcl mov", rean.grid)
plot4b <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main =  list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]),cex =0.8), 
                     color.theme = "Blues", at = seq(minlcl.cm,maxlcl.cm,diflcl.cm),
                     colorkey = list(width = 0.6, lables = list(cex = 0.6)))

plot4b
min(movlcl1matrixmean)
max(movlcl1matrixmean)

plot5a <- spatialPlot(grid = degrees[[bn]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = 9, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                      color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))

movdeg1<-movingmedias(measureslistBay[[bn]]$degree[backpermutations[[permused]]])
movdeg1matrixmean <- apply(movdeg1, MARGIN = 3, FUN = mean, na.rm = TRUE)

memClim <- quantity2clim(movdeg1matrixmean, "deg mov", rean.grid)
plot5b <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main =  list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]),cex =0.8), color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
plot5b




grid.arrange(plot1a,plot2a,plot3a,plot4a,plot5a)
grid.arrange(plot1b,plot2b,plot3b,plot4b,plot5b)


plotname <- paste0("../plots/movingmedias_thesis/movs_BN",nedgesnetworks[bn],"_perm_",permused,"_",rean,"_all.pdf")
pdf(plotname, width = 25, height = 7)
grid.arrange(plot1a,plot2a,plot3a,plot4a,plot5a,plot1b,plot2b,plot3b,plot4b,plot5b,nrow =2)
dev.off()
##################################################################
# All measures ER
##################################################################
ER <-18
plot5c <- spatialPlot(grid = degreesER[[ER]], backdrop.theme = "coastline", set.min = NULL,
                      set.max = 9, lonCenter = 0, main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                      color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
plot5c
plot5a
movdeg1.ER<-movingmedias(measureslistER[[ER]]$degree)
movdeg1matrixmean.ER <- apply(movdeg1.ER, MARGIN = 3, FUN = mean)

memClim.ER <- quantity2clim(movdeg1matrixmean.ER, "deg mov", rean.grid)
plot5d <- spatialPlot(grid = memClim.ER, backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter =0, main =  list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]),cex =0.8), color.theme = "Oranges",colorkey = list(width = 0.6, lables = list(cex = 0.6)))
plot5b
plot5d

c(plotsbn, nrow = 2, as.table = TRUE)

grid.arrange(plot5a,plot5b,plot5c,plot5d)

######################################################
# Measures CN
######################################################
CN <- 3
plot1e <- spatialPlot(closenessesCM[[CN]], main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8, pos = 0.25),
                    backdrop.theme = "coastline", rev.colors = TRUE,
                    colorkey = list(width = 0.6, lables = list(cex = 0.5)))

plot1e


movclos1.cm<-movingmedias(measures10d[[CN]]$closeness)
movclos1matrixmean.cm <- apply(movclos1.cm, MARGIN = 3, FUN = mean)

memClim.cm <- quantity2clim(movclos1matrixmean.cm, "clos mov", rean.grid)
plot1f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), rev.colors = TRUE,colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1f

grid.arrange(plot1a,plot1b,plot1e,plot1f)
grid.arrange(plot1a,plot1e, ncol = 2)
##################################################################
# measures BN
##################################################################
bn <-8
at1 <- seq(2.49*10^-5,2.5*10^-5,0.0025*10^-5)
at2 <- seq(0.00020,0.00050,0.000033)
plot1a <- spatialPlot(closenesses[[bn]], main = list(paste0("BN: ",nedgesnetworks[bn]), cex = 0.8),
                      backdrop.theme = "coastline", rev.colors = TRUE, at= at2,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5)), set.max = 0.0005)

movclos1<-movingmedias(measureslistBay[[bn]]$closeness[backpermutations[[permused]]])
movclos1matrixmean <- apply(movclos1, MARGIN = 3, FUN = mean,na.rm = TRUE)
memClim <- quantity2clim(movclos1matrixmean, "clos mov", rean.grid)
plot1b <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                      set.max = 0.0005, lonCenter =0, main = list(paste0("Neighbor's average BN: ",nedgesnetworks[bn]), cex = 0.8), rev.colors = TRUE,colorkey = list(width = 0.6, lables = list(cex = 0.6)))
plot1b
plot1a
######################################################
# Measures CN
######################################################
CN <- 1
at1 <- seq(2.49*10^-5,2.5*10^-5,0.0025*10^-5)
at2 <- seq(0.00020,0.00050,0.000033)
plot1e <- spatialPlot(closenessesCM[[CN]], main = list(paste0("CN: ",nedgesnetworksCM[CN]), cex = 0.8, pos = 0.25),
                      backdrop.theme = "coastline", rev.colors = TRUE,
                      colorkey = list(width = 0.6, lables = list(cex = 0.5))
                      #,
                      #at = at2,set.max = 0.000025
                      )
max(closenessesCM[[CN]]$Data)
plot1e


movclos1.cm<-movingmedias(measures10d[[CN]]$closeness)
movclos1matrixmean.cm <- apply(movclos1.cm, MARGIN = 3, FUN = mean,na.rm = TRUE)

memClim.cm <- quantity2clim(movclos1matrixmean.cm, "clos mov", rean.grid)
plot1f <- spatialPlot(grid = memClim.cm, backdrop.theme = "coastline", set.min = NULL,
                      set.max = NULL, lonCenter =0, main = list(paste0("Neighbor's average CN: ",nedgesnetworksCM[CN]), cex = 0.8), rev.colors = TRUE,colorkey = list(width = 0.6, lables = list(cex = 0.6)))
# set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1f


grid.arrange(plot1a,plot1b,plot1e,plot1f)






