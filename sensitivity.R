#############################################################################
# Precision matrix with Glasso method
#############################################################################
library(glasso)
data
# precisionr05 <- glasso(cormatrix, rho = 0.5)
# precisionr04 <- glasso(cormatrix, rho = 0.4)
# precisionr025 <- glasso(cormatrix, rho = 0.25) 677 opnieuw
 precisionr0225 <- glasso(cormatrix, rho = 0.225) # 847
# precisionr022 <- glasso(cormatrix, rho = 0.22) 894
# precisionr020<-  glasso(cormatrix, rho = 0.20) 1122
# precisionr0175 <- glasso(cormatrix, rho = 0.175) 1451
# precisionr017 <- glasso(cormatrix, rho = 0.17) 1556
precisionr0165 <- glasso(cormatrix, rho = 0.165) # 1643
precisionr013 <- glasso(cormatrix, rho = 0.13) # 2584
#precisionr011 <- glasso(cormatrix, rho = 0.11) # 3432
#precisionr01 <- glasso(cormatrix, rho = 0.1) # 3981 
#precisionr009 <- glasso(cormatrix, rho = 0.09) #4686
#precisionr0075 <- glasso(cormatrix, rho = 0.075) # 6136
#precisionr007<- glasso(cormatrix, rho = 0.07) #6765
#precisionr0065 <- glasso(cormatrix, rho = 0.065) #7476
precisionr0063 <- glasso(cormatrix, rho = 0.063) #7789
#precisionr006 <- glasso(cormatrix, rho = 0.06) #8239
#precisionr005 <- glasso(cormatrix, rho = 0.05, trace = TRUE) # 10204
# precision <- glasso(cormatrix, rho = 0.01) 36356 

save(precision, precisionr005,precisionr006,precisionr0065,precisionr009, precisionr0075,precisionr007, precisionr01,precisionr011,precisionr013,precisionr0165,precisionr017,precisionr0175, precisionr020,precisionr022, precisionr0225,precisionr025, precisionr04, precisionr05, file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/HC10dprecision.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/HC10dprecision.rda")
precisionr0063$wi
which(precision$wi != 0)
adj.matrix <- precisionr0063$wi
# m <- matrix(data = 1:9, nrow =3, ncol = 3)
# m
#-m[2,3]/sqrt(m[2,2]*m[3,3])
f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
# f(2,3,m)
m <- precisionr0063$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
adj.matrix
diag(adj.matrix) <- 0
adj.matrix[adj.matrix != 0 ] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
dev.off()
plot.Meteodag(data,precgraph)
graph_from_adjacency_matrix
################################################################################################
# False but interesting method:
################################################################################################
m <- precisionr0065$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
diag(adj.matrix) <- 0
adj.matrix[adj.matrix > .15 ] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
dev.off()
plot.Meteodag(data,precgraph)


adj.matrix <- PartVar

################################################################################################
# or obtain precision graph with fixed rho and treshold tau as in CN: 
################################################################################################
th = 0.05
# Adjacency matrix
f <- function(r,c,m){-m[r,c]/sqrt(m[r,r]*m[c,c])}
Vecf <- Vectorize(f,vectorize.args = c('r','c'))
# f(2,3,m)
m <- precision$wi
PartVar <- outer(1:nrow(m),1:ncol(m),Vecf,m)
adj.matrix <- PartVar
# adj.matrix <- abs(precision$wi) # beta = 0.01
diag(adj.matrix) <- 0
adj.matrix[adj.matrix <= th ] <- 0
adj.matrix[abs(adj.matrix) > th ] <- 1
precgraph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
E(precgraph)
dev.off()
plot.Meteodag(data,precgraph)


################################################################################################
# Create measuress precision graph
################################################################################################
igraphske<- precgraph
# create graph_from_Grid object with TimeCoordsAnom_from_Grid 
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
#check
measuresPrec <- graph2measure(graphObject)
clim.awcon.Prec10d <- measure2clim(measuresPrec, what = "awconnectivity", ref.grid = tas_ncep_10d)
clim.betw.Prec10d<- measure2clim(measuresPrec, what = "betweenness", ref.grid = tas_ncep_10d)
clim.close.Prec10d <- measure2clim(measuresPrec, what = "closeness", ref.grid = tas_ncep_10d)
clima.lclus.Prec10d <- measure2clim(measuresPrec, what = "localclustering", ref.grid = tas_ncep_10d)
plotClimatology(clim.betw.Prec10d, backdrop.theme = "coastline")
plotClimatology(clim.close.Prec10d, backdrop.theme = "coastline")
plotClimatology(clim.awcon.Prec10d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.Prec10d,backdrop.theme = "coastline")

# Compare with bayesian 1400 
plotClimatology(clim.betw.Bay10d, backdrop.theme = "coastline")
plotClimatology(clim.close.Bay10d, backdrop.theme = "coastline")
plotClimatology(clim.awcon.Bay10d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.Bay10d,backdrop.theme = "coastline")

##################################################################################
# Investigate locations of vertices on worldmap
##################################################################################
example <- hc_edges_loglik_10d_1400_1600i$networks
str(example)
igraphdir <- igraph.from.graphNEL(as.graphNEL(example))
igraphske<- as.undirected(igraphdir)
# create graph_from_Grid object with TimeCoordsAnom_from_Grid
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
check
plot.Meteodag(graphObject,igraphdir)
dev.off()
# Make bn 
fitted <- bn.fit(example, as.data.frame(graphObject$data_coords))
# mutilate <- mutilated(fitted, evidence = list("V81" > 0.5))
# dagfitted <- bn.net(fitted)
# dagfitted
# dagfittedigraph <- igraph.from.graphNEL(as.graphNEL(dagfitted))
# dagfittedigraph
# dagmutilateigraph <- igraph.from.graphNEL(as.graphNEL(mutilate))
# dagmutilateigraph
# V(dagfittedigraph)[c("V338","V229")]$color <- "purple"
# V(dagfittedigraph)[V(dagfittedigraph)$name != "V338"]$color <- "green"
# V(dagfittedigraph)$color <- "yellow"
# plot
# V(dagfittedigraph)!="V338"
# V(dagfittedigraph)["V338"]$size <- 300
# V(dagfittedigraph)["V338"]$color
# V(dagfittedigraph)["V229"]$color
# dagmutilate$nodes$V338



meteodag <- dagfittedigraph
class(meteodag)

vertex_attr(meteodag, "color", index = V(meteodag))
vertex_attr(meteodag, "name", index = V(meteodag))

  if (class(meteodag) == "igraph") 
    igraphDegree <- meteodag
  
  if (class(meteodag) == "graphNEL") 
    igraphDegree <- igraph.from.graphNEL(meteodag)
  
  if (class(meteodag) == "bn") 
    igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))
  
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  points
  
  plot(wrld)
  plot.igraph(meteodag, 
              vertex.label = vertex_attr(meteodag, "name", index = V(meteodag)),
              vertex.label.cex = 0.5,
              vertex.label.degree = -pi/2,
              vertex.color = vertex_attr(meteodag, "color", index = V(meteodag)),
              vertex.size = 100,
              edge.arrow.size = 0,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  



pdf("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/Nodesgraph.pdf",
    width = 10, height = 6 )
dev.off()

##################################################################################
# Sensitivity Analysis Bayesian networks
# Exploring example
##################################################################################
example <- hc_edges_loglik_10d_400_600i$networks
str(example)
# Make bn 
fitted <- bn.fit(example, as.data.frame(graphObject$data_coords))

with <- cpquery(fitted, event = (V99 >= 0.5), evidence = (V81 >= 0.5))
withcomplement <- cpquery(fitted, event = (V299 >= 0.5), evidence = (V81 < 0.5))
without <- cpquery(fitted, event = (V99 >= 0.5), evidence = TRUE)

probs600V99 <- c(with,withcomplement,without)
probs600V286 <- c(with,withcomplement,without)
probs1000V99 <- c(with,withcomplement,without)
probs1600V99 <- c(with,withcomplement,without)
probs1000V286 <- c(with,withcomplement,without)
probs1600V286 <- c(with,withcomplement,without)

probs600V99
probs600V286
probs1000V99
probs1600V99
probs1000V286
probs1600V286


igraphdir <- igraph.from.graphNEL(as.graphNEL(example))
igraphske<- as.undirected(igraphdir)
# create graph_from_Grid object with TimeCoordsAnom_from_Grid
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
check
plot.Meteodag(graphObject,igraphdir)
dev.off()
##################################################################################
# Sensitivity Analysis Bayesian networks
# Loop 
##################################################################################
example <- hc_edges_loglik_10d_400_600i$networks
str(example)
# Make bn 

fitted <- bn.fit(example, as.data.frame(graphObject$data_coords))
df <- as.data.frame(graphObject$data_coords)
str(df)

with <- list()
names(fitted)[1]
str <- paste0("(", names(fitted)[i], ">=", 0.5, ")")
parse(text = str)
for(i in 1:2) {
  str <- paste0("(", names(fitted)[i], ">=", 0.5, ")")
  with[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 >= 0.5))
}
with

######## ya está hecho Lis ;)     ----
withcomplement <- list()
for(i in 1:2) {
  str <- paste0("(", names(fitted)[i], ">=", 0.5, ")")
  withcomplement[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 < 0.5))
}
withcomplement
########

without <- list()
for(i in 1:2) {
  str <- paste0("(", names(fitted)[i], ">=", 0.5, ")")
  without[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = TRUE)
}
without
data.frame(names,with,withcomplement,without)
##################################################################################
# Sensitivity Analysis Bayesian networks
# Loops all in one
##################################################################################

Propagation <- function(BN, startnode, endnode, evidencevalue){
 
  a <- evidencevalue
  str2 <- paste0("(V81>=", a, ")")
  str3 <- paste0("(V81<", a, ")")
  
  with <- c()
  withcomplement <- c()
  without <- c()
  nod <- c()
  
  for(i in startnode:endnode) {
  str1 <- paste0("(", names(BN)[i], ">=", a, ")")
  with[i] <- cpquery(BN, event = eval(parse(text = str1)), evidence = eval(parse(text = str2)))
  withcomplement[i] <- cpquery(BN, event = eval(parse(text = str1)), evidence = eval(parse(text = str3)))
  without[i] <- cpquery(BN, event = eval(parse(text = str1)), evidence = TRUE)
  nod[i] <- names(BN)[i]
  }

  out <- data.frame(nodes = nod,
                    with, 
                    withcomplement,
                    without)
  
  return(out)
}

#############################################################################3
# Estimate propagtion time
#############################################################################
example <- hc_edges_loglik_10d_1400_1600i$networks
df <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
fitted <- bn.fit(example, df)

length(nodes(fitted))
test <- Propagation(example,284,285,1)
test[c(284,285),]
system.time(Propagation(fitted,284,285,1), gcFirst = TRUE)

#################################################################################
# Comando propagation
#
#################################################################################
df <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
example1 <- hc_edges_loglik_10d_200_400i$networks
example2 <- hc_edges_loglik_10d_600_800i$networks
example3 <- hc_edges_loglik_10d_1400_1600i$networks
bn1 <- bn.fit(example1, df)
bn2 <- bn.fit(example2, df)
bn3 <- bn.fit(example3, df)

length(nodes(BN))

propagationbn1 <- Propagation(bn1,1,length(nodes(bn1)),0.5)
propagationbn2 <- Propagation(bn2,1,length(nodes(bn2)),0.5)
propagationbn3 <- Propagation(bn3,1,length(nodes(bn3)),0.5)


##################################################################################
# Correlation matrixes 
# Indices come from Nino script
##################################################################################
library(corrplot)
# Select data for correlations
data <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
datanino <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = ind1981_2010)
dataColdbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = coldmonths)
dataWarmbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = warmmonths)
# make different correlation matrices
cormatrix <- cor(data, method = "spearman")
cormatrixnino <- cor(datanino, method = "spearman")
cormatrixwarm <- cor(dataWarmbasin, method ="spearman")
# Give variable names (Vertex Names)
# names as in Nino script
class(cormatrix)
rownames(cormatrix) <- colnames(cormatrix) <- names
rownames(cormatrixnino) <- colnames(cormatrixnino) <- names
rownames(cormatrixwarm) <- colnames(cormatrixwarm) <- names
# Plot some scatterplot
par("mar")
par(mar = c(1,1,1,1))
par(mfrow = c(1,3))
plot(cormatrix)
plot(cormatrixnino)
plot(cormatrixwarm)
# Corrplot from library corrplot
# Corrplots: "heatmaps" 
par("mar")
par(mar = c(1,1,1,1))
par(mfrow = c(1,3))
corrplot(cormatrix, method = "color", type = "full", is.corr = FALSE)
corrplot(cormatrixnino, method = "color", type = "full", is.corr = FALSE)
corrplot(cormatrixwarm, method = "color", type = "full", is.corr = FALSE)


pairs(cormatrix[1:5,1:5]) 



dev.off()


########################################################################################
# Correlation matrix with ggplot
########################################################################################
# Convert correlation matrices to ggplot format
library(reshape2)
library(ggplot2)
str(cormatrix)
melted_cormat <- melt(cormatrix, as.is = TRUE) # GOED: met variable namen!!! NU GROOT UITPRINTEN
str(melted_cormat)
melted_cormatabs <- melt(abs(cormatrix), as.is = TRUE)
melted_cormatnino <- melt(cormatrixnino, as.is = TRUE)
melted_cormatabsnino <- melt(abs(cormatrixnino), as.is = TRUE)
melted_cormatwarm <- melt(cormatrixwarm, as.is = TRUE)
melted_cormatabswarm <- melt(abs(cormatrixwarm), as.is = TRUE)

melted_cormatdif <- melted_cormat
melted_cormatdif[,3] <- melted_cormat[,3]-melted_cormatnino[,3]
melted_cormatdifabs <- melted_cormatabs
melted_cormatdifabs[,3] <- abs(melted_cormat[,3]-melted_cormatnino[,3])
melted_cormatdifabs

# Plot correlation matrices with geom_raster() 
# and highlight / select variables with scale _ discrete breaks /limits.
# For example Nino basin or various regions
somevar <- c(81,286,297,205,604)
somevar <- indNino
allvar <- as.vector(seq(1:ncol(data)))
offvar <- setdiff(allvar,somevar)

# normal correlation matrices
plotcor <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill= value))+ 
  geom_raster()+
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor all data\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])        
plotcor
plotcornino <- ggplot(data = melted_cormatnino, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster()+
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor nino exact\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrixnino)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrixnino)[somevar]) 
plotcornino
plotcorwarm <- ggplot(data = melted_cormatwarm, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor warm basin\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrixwarm)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrixwarm)[somevar]) 

grid.arrange(plotcor,plotcornino, plotcorwarm)

# Absolute correlation matrices 
plotcorabs <- ggplot(data = melted_cormatabs, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = c("white","orange","red"), 
                       #midpoint = 0, 
                       limit = c(0,.8), space = "Lab", 
                       name="|cor| all data\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])   

plotcorabsnino <- ggplot(data = melted_cormatabsnino, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = c("white","orange","red"), 
                       #midpoint = 0, 
                       limit = c(0,.8), space = "Lab", 
                       name="|cor| nino exact\nSpearman") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])   

plotcorabswarm <- ggplot(data = melted_cormatabswarm, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = c("white","orange","red"), 
                       #midpoint = 0, 
                       limit = c(0,.8), space = "Lab", 
                       name="|cor| warm basin\nSpearman") + 
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar])   
grid.arrange(plotcorabs,plotcorabsnino, plotcorabswarm)

# differencias between all data and nino data or warm months.
plotcordif <- ggplot(data = melted_cormatdif, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-0.5,0.5), space = "Lab", 
                       name="Spearman\n cor all - cor nino") +
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar]) 

plotcordifabs <- ggplot(data = melted_cormatdifabs, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("white","yellow","green","blue"), 
                       #midpoint = 0, 
                       limit = c(0,0.5), space = "Lab", 
                       name="Spearman\n |cor all - cor nino|") + 
  scale_x_discrete(breaks = rownames(cormatrix)[somevar]) +
  scale_y_discrete(breaks = rownames(cormatrix)[somevar]) 

grid.arrange(plotcordifabs, plotcordif)


####################################################################################
# Try to order by magnitude correlation. (NOT YET: VARIABLE NAMES UNCLEAR)
#####################################################################################
order.AOE <- corrMatOrder(cormatrix, order = "AOE")
corrAOE <- cormatrix[order.AOE,order.AOE]
rownames(corrAOE)
colnames(corrAOE)
melted_cormatAOE <- melt(corrAOE, as.is = TRUE)
melted_cormatAOE

orderedcordifabs <- melted_cormatAOE
rownames(corrAOE)[somevar,]
rownames(corrAOE)[527]
indvector <- c()
for (i in 1:length(somevar)){
indvector <- append(indvector, which(rownames(corrAOE) == names[somevar][i]))
       }
rownames(corrAOE)[indvector]
plotorderedcordifabs<- ggplot(data = melted_cormatAOE, aes(x=Var1, y=Var2, fill= value)) + 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  scale_x_discrete(breaks = rownames(corrAOE)[indvector]) +
  scale_y_discrete(breaks = rownames(corrAOE)[indvector])
plotorderedcordifabs
ord <- hclust(melted_cormatdifabs, method = "complete" )
dist(melted_cormatdifabs[,3])
melted_cormatdifabs
########################################################################################
# Correlation matrix with ggcorrplot: Cluster option
########################################################################################
install.packages("ggcorrplot")
library("ggcorrplot")
hclust(cormatrix)
names <- c()
names
for(i in 1:648) names <- append(names,paste0("V",i))
names
matrix
rownames(cormatrix) <- names
colnames(cormatrix) <- names
rownames(cormatrixnino) <- names
colnames(cormatrixnino) <- names
rownames(cormatrixwarm) <- names
colnames(cormatrixwarm) <- names

cormatrix
?hclust
x <- hclust(as.dist(cormatrix), method ="complete")
ggcorrplot(cormatrix, hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6) + scale_x_discrete(limits = c(81,286,315,205,604))
ggcorrplot(cormatrixnino,hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6) # HOW ORDER?
plotcorabsO <- ggcorrplot(abs(cormatrix), outline.color = NA, tl.cex = 2, tl.srt = pi/6)
plotcorabsninoO <- ggcorrplot(abs(cormatrixnino), outline.color = NA, tl.cex = 2, tl.srt = pi/6)
ggcorrplot(cormatrixnino - cormatrix, hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6)
plotcorabsdifO <- ggcorrplot(abs(cormatrixnino - cormatrix), hc.order = TRUE, outline.color = NA, tl.cex = 2, tl.srt = pi/6)

grid.arrange(plotcorabsO,plotcorabsninoO,plotcorabsdifO)
###################################################################################
# Clustering from data correlation matrix.
###################################################################################
x <- hclust(as.dist(cormatrix), method ="single")
x$order
x$labels
cormatrix
cormatrixPer <- cormatrix[x$order,x$order]
melt(cormatrixPer, as.is = TRUE)
plotcorPer <- ggplot(data = melt(cormatrixPer), aes(x=Var1, y=Var2, fill= value))+ 
  geom_raster() +
  scale_fill_gradientn(colours = c("red","orange","white","green","blue"), 
                       #midpoint = 0, 
                       limit = c(-1,1), space = "Lab", 
                       name="cor all data\nSpearman") +
  scale_x_discrete(breaks = c("V81","V286","V315","V205","V604")) +
  scale_y_discrete(breaks = c("V81","V286","V315","V205","V604"))

plotcorPer
####################################################################################
# Clustering by mean temperature ?
####################################################################################
data <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
plotClimatology(climatology(tas_ncep_10d))
tas_ncep_10d$Data
time.coords.matrix <- array3Dto2Dmat(tas_ncep_10d$Data)
colMeans(time.coords.matrix)

datanino <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = ind1981_2010)
dataColdbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = coldmonths)
dataWarmbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = warmmonths)
# Different correlation matrices
str(data)

cormatrix <- cor(data, method = "spearman")
cormatrixnino <- cor(datanino, method = "spearman")
cormatrixwarm <- cor(dataWarmbasin, method ="spearman")


####################################################################################
# Investigate detection long range arcs.
####################################################################################
example <- hc_edges_loglik_10d_800_1000i$networks
igraphdir <- igraph.from.graphNEL(as.graphNEL(example))
igraphske<- as.undirected(igraphdir)
is_connected(igraphdir)
compdir <- component_distribution(igraphdir)
plot(compdir)
components(igraphdir)
subcomponent(igraphdir, "V244", mode = "all")



##############################################################################
# Try cpquery and mutilate on gaussian.test
#############################################################################
data(gaussian.test)
dag <- hc(gaussian.test)
fitted <- bn.fit(dag, gaussian.test)
# the result should be around 0.04.
cpquery(fitted,
        event = ((A >= 0) & (A <= 1)) & ((B >= 0) & (B <= 3)),
        evidence = (C + D < 10))
muti <- mutilated(fitted, evidence = list(C = 10))
bn.net(muti)
bn.net(fitted)
hamming(muti,fitted)
as.bn(fitted)

