library(gridExtra)
library(gridGraphics)
library(grid)
library(glasso)
setwd("~/data/Untitled/Trabajo/R_practice/R/")
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/")
###########################################################################################################
#Complex Networks 10d
###########################################################################################################
#load the data 10 degrees
load("../Data/tas_ncep_10d.rda")
# Calculate graph / adjacency matrix for specific tau. Insert possible subindex in Time. 
tau = 0.5
graph10d <- graph_from_Grid(tas_ncep_10d, th = tau, subind = NULL, method = "pearson")

graph10d$data_coords
graph10d$correlation
graph10d$VertexCoords
graph10d$adjacency
graph10d$graph

# Calculate measures from adjacency matrix and graph properties.
measures10d <- graph2measure(graph10d)

measures10d$edens
measures10d$ddist
measures10d$betweenness
measures10d$closeness
measures10d$awconnectivity
measures10d$localclustering

# visualize graph: with plot.Meteodag (time.coords only uses X and Y coords attributes)
dev.off()
numberofnodes <- length(V(graph10d$graph))
numberofedges <- length(E(graph10d$graph))
edgedensity <- round(measures10d$edens, digits = 4)
plot.Meteodag(graph10d$graph,time.coords = graph10d, lis = TRUE)
title(main = paste0("Complex graph \n number of edges = ", numberofedges,"\n egde density = ", edgedensity))


# create climatologyobjects of measures
clim.awcon.10d <- measure2clim(measures10d, what = "awconnectivity", ref.grid = tas_ncep_10d)
clim.betwe.10d <- measure2clim(measures10d, what = "betweenness", ref.grid = tas_ncep_10d)
clim.close.10d <- measure2clim(measures10d, what = "closeness", ref.grid = tas_ncep_10d)
clim.lclus.10d <- measure2clim(measures10d, what = "localclustering", ref.grid = tas_ncep_10d)

# visualize with plotClimatology
c <- spatialPlot(clim.awcon.10d, backdrop.theme = "coastline", main = "area weighted connectivity")
d <- spatialPlot(clim.betwe.10d, backdrop.theme = "coastline", main = "betweenness")#, at = seq(0,6,0.5))
e <- spatialPlot(clim.close.10d, backdrop.theme = "coastline", main = "closeness")#, at = seq(2.38e-6,2.53e-6,0.01e-6))
f <- spatialPlot(clim.lclus.10d, backdrop.theme = "coastline", main = "local clustering")
#,at = seq(0.01,0.25,0.015)
#, at = seq(2.38e-0.6,2.50e-0.6,0.02e-0.6) at = seq(2.38e-6,2.50e-6,0.02e-6)

# visualize measures per tau
#plotname <- paste0("/home/catharina/Documents/PlotsResumen/MeasuresAndDegree10dtau",tau,".pdf")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/MeasuresAndDegree10dtau",tau,".pdf")
plplot(tas_ncep_10d, tau = tau)
plplot
grid <- tas_ncep_10d

grid.echo()
a <- grid.grab()
text <- textGrob(paste0("Complex graph \n nodes = ", numberofnodes,"\n egde density = ", edgedensity))
pdf(plotname)
grid.arrange(text,a,c,d,e,f, heights=c(3,1,1), widths=c(1,1))
dev.off()

# visualize graphs per tau
#plotname <- paste0("/home/catharina/Documents/PlotsResumen/ComplexGraph10dtau",tau,".pdf")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/ComplexGraph10dtau",tau,".pdf")
pdf(plotname)
plot.Meteodag(graph10d$graph,time.coords = graph10d, lis = TRUE)
title(main = paste0("Complex graph \n number of nodes = ", numberofnodes,"\n egde density = ", edgedensity))
dev.off()

#visualize degree distribution (ddist) with respect to the tresholds
#plotname <- "/home/catharina/Documents/PlotsResumen/Powerlaw10dtau01_08.pdf"
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Powerlaw10dtau01_08.pdf")
pdf(plotname)
plplotlist(tas_ncep_10d, tau = c(0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
dev.off()

#################################################################################################
# Bayesian networks
#################################################################################################
load("../Data/data_hcnetworks10d.rda")
# visualize 1 graph: with plot.Meteodag (time.coords only uses X and Y coords attributes)
dag <- hc_edges_loglik_10d_1400_1600i$networks
numberofnodes <- length(nodes(dag))
numberofedges <- narcs(dag)

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Bayesiangraph",narcs(dag),".pdf")
plotname
pdf(plotname)

edgedensity <- round(numberofedges/(choose(numberofnodes,2)), digits = 4)
TimeCoordsAnom_from_Grid(tas_ncep_10d)
plot.Meteodag(dag, time.coords = TimeCoordsAnom_from_Grid(tas_ncep_10d), lis = TRUE)
title(main = paste0("Bayesian graph \n number of edges = ", numberofedges,"\n egde density = ", edgedensity))
dev.off()

#################################################################################################
# Get information of graph: mean distance HC
#################################################################################################
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
hc_list <- list(hc_edges_loglik_10d_2800_3000i,
                hc_edges_loglik_10d_3000_3200i,
                hc_edges_loglik_10d_3200_3400i,
                hc_edges_loglik_10d_3400_3600i,
                hc_edges_loglik_10d_3600_3800i,
                hc_edges_loglik_10d_3800_4000i,
                hc_edges_loglik_10d_4000_4200i,
                hc_edges_loglik_10d_4200_4400i,
                hc_edges_loglik_10d_4400_4600i)
# select networks
networks <- lapply(hc_list, function(m) m[["networks"]])
networks
nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]
what <- "both"
# Calculate mean distance
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
grids <- lapply(networks, mean_distance_edges, grid = tas_ncep_10d, data.dag = data10d, what = what)
# Plot
plotname <- paste0("/home/catharina/Documents/PlotsResumen/HCmeandists",firstel,"_",lastel,what,".pdf")
pdf(plotname)
multis = list()
for (i in 1:length(grids)){
  plot <- spatialPlot(grids[[i]], main = paste0("number of edges = ",nedgesnetworks[i]) ,backdrop.theme = "coastline")
  multis[[i]] <- plot
}

n <- length(multis)

nCol <- floor(sqrt(n))
do.call("grid.arrange", c(multis, ncol=nCol))
dev.off()
#################################################################################################
# Get information of graph: Long distance HC
#################################################################################################
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
#hc_list <- list(hc_edges_loglik_10d_2800_3000i,
#                hc_edges_loglik_10d_3000_3200i,
#                hc_edges_loglik_10d_3200_3400i,
#                hc_edges_loglik_10d_3400_3600i,
#                hc_edges_loglik_10d_3600_3800i)

# select networks
networks <- lapply(hc_list, function(m) m[["networks"]])
nedgesnetworks <- as.character(sapply(networks, narcs))
firstel <- nedgesnetworks[1]
lastel <- nedgesnetworks[length(nedgesnetworks)]
minimdist = 10000
#calculate area weighted dag degree for all hillclimbing dags
plotname <- paste0("/home/catharina/Documents/PlotsResumen/HClongdists",firstel,"_",lastel,"_",minimdist,"km.pdf")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HClongdists",firstel,"_",lastel,"_",minimdist,"km.pdf")

pdf(plotname)
par(mfrow = c(4,2),cex = 0.4)
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
lapply(networks, plot_long_distances, data.dag = data10d, minimdist = minimdist)

dev.off()



###############################################################################################
# Compare HC and PC
# with loglikelihood function
###############################################################################################
# PC algorithm. Spearman correlation. DAG scores if DAG is produced to number of edges
# List of alphas
myalphas7 <- as.vector(seq(1e-15,1e-13,1e-15))
myalphas8 <- as.vector(seq(1e-17,1e-15,1e-17))
myalphas9 <- as.vector(seq(1e-19,1e-17,1e-19))
myalphas10 <- as.vector(seq(1e-21,1e-19,1e-21))
myalphas10
# data 
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
corrspear10d <- cor(data10d, method = "spearman")
# Dataframes include: loglik(dag|data), number of undirected edges in eqclass (times 2), number of directed edges. 


loglik10d_21_19 <- sapply(myalphas10, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_19_17 <- sapply(myalphas9, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_17_15 <- sapply(myalphas8, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_15_13 <- sapply(myalphas7, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_19_13 <-cbind(loglik10d_19_17,loglik10d_17_15, loglik10d_15_13)
loglik10d_21_13 <- cbind(loglik10d_21_19,loglik10d_19_13)

# bind some extra individuals
vars <- load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Data_graphsPC.rda")
as.array(vars)
loglik10d_21_13_mas <- cbind(loglik10d_21_13, t(test10d_12_e_42),t(test10d_4_e_70))
loglik10d_21_13_mas
loglik10d_21_13
t(test10d_12_e_42)

plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Plotloglikpc10000hc.pdf"
pdf(plotname)
plot.nedge.vs.loglik_CNvsBN(pc3obj = loglik10d_21_13_mas,it_hc_obj = 
                              hc_edges_loglik_10d_8695i, xlim = c(0,10076))
abline(v = c(648, 800,2200,3059,10076), col = c("red","blue","blue" ,"orange","orange"))
legend("topright",legend = c("information range BN","information range CN"),
       col = c("blue" ,"orange"),
       lty = c(1,1),
       cex = 1,
       bty = "n")
dev.off()

# The extended version with a range of alphas can be found in DocumentScutariScript
plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Plotloglikpc1000hc.pdf"
pdf(plotname)
plot.nedge.vs.loglik_combi(loglik10d_21_13_mas,hc_edges_loglik_10d_1000i)
abline(v = 648, col = "red")
dev.off()

plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/Plotloglikpc8695hc.pdf"
pdf(plotname)
plot.nedge.vs.loglik_combi(loglik10d_21_13_mas,hc_edges_loglik_10d_8695i)
abline(v = 648, col = "red")
dev.off()


save(loglik10d_15_13,loglik10d_17_15,loglik10d_19_17,loglik10d_21_19, 
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_pcnetworks10d.rda")


###############################################################################################
# Compare HC and PC
# with strongarcs function. 25 procent strongest arcs.
###############################################################################################
#plotname <- "/home/catharina/Documents/PlotsResumen/HCvsPCquantile25.pdf"
#plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HCvsPCquantile25.pdf"
plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HCvsPCpercent25.pdf"
pdf(plotname, width = , hei)
par(mfrow = c(2,1), cex = 0.6)
# HC:
data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
data.frame.dag <- as.data.frame(data.dag)
quantile2pc <-signif(quantile(arc.strength(hc_edges_loglik_10d_400_600i$networks, data.frame.dag)[,3]), digits =3)
quantile2pc
nedges2 <- narcs(hc_edges_loglik_10d_400_600i$networks)
nedges2
plot_strongweak_edges(hc_edges_loglik_10d_400_600i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d),
                      treshold = quantile2pc[2], which = "strong", cex = 1)

title(main = paste0("dag from HC 25 percent strongest edges \n number of edges = ",nedges2))

# PC:
#data.dag <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
#names(data.dag) <- nodes(pc_10d_1e_13$networks$dag)
nedges <- narcs(pc_10d_1e_13$networks$dag)
nedges


#quantile(duration, c(.32, .57, .98))
names(data.frame.dag) <- nodes(pc_10d_1e_13$networks$dag)
quantile1pc <-signif(quantile(arc.strength(pc_10d_1e_13$networks$dag, data.frame.dag)[,3]), digits =3)
quantile1pc[2]
plot_strongweak_edges(pc_10d_1e_13$networks$dag, data.dag,treshold = quantile1pc[2], which = "strong", cex = 1)
title(main = paste0("dag from PC 25 percent strongest \n number of edges = ",nedges))

dev.off()
##############################################################################################
# Compare HC and PC 
# With long distances.
###############################################################################################
plotname <- "/home/catharina/Documents/PlotsResumen/HCvsPCminimdistances.pdf"
pdf(plotname)
10 * pi*6000/180
par(mfrow = c(3,2))
data10d <- data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plot_long_distances(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 1150)
plot_long_distances(pc_10d_1e_13$networks$dag, data10d, minimdist = 1150)
plot_long_distances(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 1500)
plot_long_distances(pc_10d_1e_13$networks$dag, data10d, minimdist = 1500)
plot_long_distances(hc_edges_loglik_10d_400_600i$networks, data10d, minimdist = 2000)
plot_long_distances(pc_10d_1e_13$networks$dag, data10d, minimdist = 2000)
dev.off()
##############################################################################################
# Compare HC and PC 
# With plot meteodag.
###############################################################################################
load(file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Data_graphsPC.rda")
# In the following also include the graph of 400 edges of PC
data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plotname <- "/home/catharina/Documents/PlotsResumen/HCvsPCnormalgraphs.pdf"
plotname <- "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/HCvsPCnormalgraphs.pdf"
pdf(plotname)
par(mfcol = c(3,2), cex = 0.5)
plot.Meteodag(data10d, hc_edges_loglik_10d_200i$networks)
title(main = paste0("dag from HC \n number of edges = ",narcs(hc_edges_loglik_10d_200i$networks)))
plot.Meteodag(data10d, hc_edges_loglik_10d_200_400i$networks)
title(main = paste0("dag from HC \n number of edges = ",narcs(hc_edges_loglik_10d_200_400i$networks)))
plot.Meteodag(data10d, hc_edges_loglik_10d_400_600i$networks)
title(main = paste0("dag from HC \n number of edges = ",narcs(hc_edges_loglik_10d_400_600i$networks)))
plot.Meteodag(data10d, pc_10d_4e_70$networks$dag)
title(main = paste0("dag from PC \n number of edges = ",narcs(pc_10d_4e_70$networks$dag)))
plot.Meteodag(data10d, pc_10d_12e_42$networks$dag)
title(main = paste0("dag from PC \n number of edges = ",narcs(pc_10d_12e_42$networks$dag)))
plot.Meteodag(data10d, pc_10d_1e_13$networks$dag)
title(main = paste0("dag from PC \n number of edges = ",narcs(pc_10d_1e_13$networks$dag)))
dev.off()
# Ideas for functions with spplot instead of plot
########################################################################################
install.packages("ggplot2")
ggplot2::
  library(transformeR)
ras <- grid2mopa(climatology(NCEP_Iberia_psl))

wpols <- wrld@polygons
w <- (SpatialPolygons(wpols))

spplot(ras, sp.layout=list(w, first=F))
#https://blogs.oii.ox.ac.uk/bright/2015/12/07/getting-ggplot2-to-work-with-igraph/


class(a)
b <- plot(1:10)
a <- grob(a)
b
a
grid.echo()
a <- grid.grab()
grid.newpage()
a <- editGrob(a, vp=viewport(width=unit(2,"in")), gp=gpar(fontsize=10))
a
plplot(tas_ncep_10d,tau = 0.4)
grid.echo()
b <- grid.grab()
####################################################################################
# Try out to combine par / grid.arrange / 
####################################################################################

install.packages("cowplot")
library(cowplot)
matrix
par(mfrow = c(2,2))
plot(1:10) 
b
grid.arrange(c,d,e,f)
a <- recordPlot()
replayPlot(a)


b <- recordPlot()
b
dev.off()

10*pi/180 * 6378
pi*6378^2

plot_grid(a,b)

p <- recordPlot()
plot.new() ## clean up device
p # redraw

## grab the scene as a grid object


grid.echo()
a <- grid.grab()

## draw it, changes optional
grid.newpage()
a <- editGrob(a, vp=viewport(width=unit(2,"in")), gp=gpar(fontsize=10))
grid.draw(a)
############################################################################