### BF for HC and PC -------------
# compare
# arc. strength
# costum.strength
# BF

# make networks PC - 600 nodes
degree10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
test10d_1_e_200 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-200)
test10d_1_e_160 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-160) 
test10d_1_e_150 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-150) 
test10d_1_e_113 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-113) 
test10d_1_e_100 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-100) 
test10d_1_e_81 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-81) 
test10d_1_e_71 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-71) 
test10d_1_e_70 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-70) 
test10d_4_e_70 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 4e-70) 
test10d_5_e_70 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 5e-70) 
test10d_1_e_69 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-69) 
test10d_1_e_67 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-67) 
test10d_1_e_65 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-65) 
test10d_1_e_63 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-63) 
test10d_1_e_61 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-61) 
test10d_1_e_45 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-45) 
test10d_1_e_43 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-43)
test10d_1_e_42 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-42)
test10d_1_e_41 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-41)
test10d_12_e_42 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 12e-42)
test10d_15_e_42 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 15e-42)
test10d_2_e_41 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 2e-41)
test10d_1_e_40 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-40)
test10d_1_e_35 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-35)
test10d_1_e_30 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-30)
test10d_1_e_25 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-25)
test10d_1_e_22 <- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-22)
test10d_1_e_20<- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 1e-20) #
test10d_2_e_13<- alpha_loglik_pc3_simple(degree10d, method = "spearman",alpha = 2e-13) # DOES NOT WORK ANYMORE

test10d_1_e_200
test10d_1_e_160
test10d_4_e_70

pc_10d_12e_42 <- alpha_loglik_pc3_withdag(degree10d, method = "spearman", alpha = 1.2e-41)
pc_10d_3.6e_18 <- alpha_loglik_pc3_withdag(degree10d, method = "spearman",alpha = 3.6e-18)
pc_10d_1e_15 <- alpha_loglik_pc3_withdag(degree10d, method = "spearman", alpha = 1e-15)
pc_10d_1e_13 <- alpha_loglik_pc3_withdag(degree10d, method = "spearman", alpha = 1e-13)
pc_10d_4e_70 <- alpha_loglik_pc3_withdag(degree10d, method = "spearman", alpha = 4e-70)
pc_10d_4e_70

dev.off()
plot.Meteodag(degree10d,pc_10d_12e_42$networks$dag)
# Save PCs with particular edges densities above.
save(test10d_12_e_42, test10d_15_e_42, test10d_1_e_22, test10d_1_e_25, test10d_1_e_30, test10d_1_e_35,test10d_1_e_40, test10d_1_e_45, test10d_2_e_41,test10d_4_e_70,test10d_1_e_200, test10d_1_e_160, test10d_1_e_150, test10d_1_e_113,test10d_1_e_100,test10d_1_e_81 ,test10d_1_e_71 , test10d_1_e_70, test10d_4_e_70, test10d_5_e_70, test10d_1_e_67, test10d_1_e_65, test10d_1_e_63, test10d_1_e_61 ,test10d_1_e_45, test10d_1_e_43,
  pc_10d_1e_13,pc_10d_1e_15,pc_10d_3.6e_18,pc_10d_12e_42, pc_10d_4e_70,file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Data_graphsPC.rda")
save(pc_10d_1e_13,pc_10d_1e_15,pc_10d_3.6e_18,pc_10d_12e_42, pc_10d_4e_70,file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_graphsPC.rda")
rm(test10d_12_e_42, test10d_15_e_42, test10d_1_e_22, test10d_1_e_25, test10d_1_e_30, test10d_1_e_35,test10d_1_e_40, test10d_1_e_45, test10d_2_e_41,test10d_4_e_70,
   pc_10d_1e_13,pc_10d_1e_15,pc_10d_3.6e_18,pc_10d_12e_42, pc_10d_4e_70)

# Use BayesFactor to compare PC-PC , HC-HC and PC-HC
BF(pc_10d_1e_13$networks$dag, pc_10d_1e_15$networks$dag, data10d)
BF(hc_edges_loglik_10d_400_600i$networks,
   pc_10d_1e_13$networks$dag, 
   as.data.frame(degree10d$data_coords), log = TRUE)
nodes(pc_10d_1e_13$networks$dag)  <- nodes(hc_edges_loglik_10d_1000_1200i$networks)
nodes(pc_10d_1e_15$networks$dag) <- nodes(hc_edges_loglik_10d_1000_1200i$networks)


#what is difference between BF and loglikelihoods difference?
hc_edges_loglik_10d_400_600i$networkdata$logliks[112]- pc_10d_1e_13$networkdata$likgbn

BF(hc_edges_loglik_10d_2000_2200i$networks,hc_edges_loglik_10d_1800_2000i$networks, data10d, log = TRUE)

#other comparisons
arc.strength(hc_edges_loglik_10d_2000_2200i$networks,data10d)
custom.strength(list(hc_edges_loglik_10d_200i$networks, 
                     hc_edges_loglik_10d_200_400i$networks,
                     hc_edges_loglik_10d_400_600i$networks), nodes = nodes(hc_edges_loglik_10d_200i$networks) )
nodes(hc_edges_loglik_10d_200i$networks)

narcs(hc_edges_loglik_10d_400_600i$networks)
vstructs(hc_edges_loglik_10d_400_600i$networks)
hc_edges_loglik_10d_400_600i$networkdata$logliks[200]

narcs(pc_10d_1e_13$networks$dag)
vstructs(pc_10d_1e_13$networks$dag)
vstructs(hc_edges_loglik_10d_1400_1600i$networks, moral = FALSE)
pc_10d_1e_13$networkdata$likgbn

as.matrix(hc_edges_loglik_10d_400_600i$networks$arcs)[1,]
tp[1,]
tp <- bnlearn::compare(hc_edges_loglik_10d_400_600i$networks,pc_10d_1e_13$networks$dag,arcs = TRUE)$tp
which(hc_edges_loglik_10d_400_600i$networks$arcs == tp[1,])

### arc.strength visualizing -------------------------------------------
#arc.strength binnen package bnlearn
data10d <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
strength10d <- arc.strength(hc_edges_loglik_10d_400_600i$networks, data10d)
strength10d
strength10dmat <- as.matrix(strength10d)
strength10dmat

#sterkste link:
strength10d[strength10d["strength"] == min(strength10d["strength"]),,]
#zwakste link
strength10d[strength10d["strength"] == max(strength10d["strength"]),,]
#kleiner dan treshold:
strength10d[strength10d["strength"] <= -230.00000,,]
strongarcs <- as.matrix(strength10d[strength10d["strength"] <= -230.00000,,])[,1:2]
strongarcs


#create igraph 600 nodes
igraph600 <- igraph.from.graphNEL(as.graphNEL(hc_edges_loglik_10d_400_600i$networks))
E(igraph600)

#Set some colors to sequence#identify indices sequence with 'strong arcs' in igraph600
#identify indices sequence in igraph600
indsame <- c()
nrow(strength10d)

for (i in 1:nrow(strength10d)){
  int <- intersect(which(as_edgelist(igraph600)[,1] == strength10dmat[i,1]),
                   which(as_edgelist(igraph600)[,2] == strength10dmat[i,2]))
  indsame[i] <- int
}
indsame

#test voor same indices igraph600
strength10d[2,]
as_edgelist(igraph600)[165,]

#make weights
#empty weight vector
weights <- numeric(nrow(as_edgelist(igraph600)))
# fill with strengths
for (i in 1:nrow(as_edgelist(igraph600))){
  ind <- indsame[i]
  as_edgelist(igraph600)[ind,]
  weights[ind] <- strength10d[i,3]
}

weights[165]

#identify indices sequence with 'strong arcs' in igraph600
indstrong <- c()
nrow(strongarcs)
strongarcs
for (i in 1:nrow(strongarcs)){
  int <- intersect(which(as_edgelist(igraph600)[,1] == unname(strongarcs[i,1])),
                   which(as_edgelist(igraph600)[,2] == unname(strongarcs[i,2])))
  indstrong[i] <- int
}
indstrong
sort(indstrong)
as_edgelist(igraph600)[sort(indstrong),]

#Zwakke links:
no.indstrong <- as.vector(1:nrow(as_edgelist(igraph600)))
no.indstrong <- no.indstrong[-sort(indstrong)]
no.indstrong

#STel de kleuren in:
E(igraph600)[ sort(indstrong) ]$color <- "red"
E(igraph600)[ sort(no.indstrong) ]$color <- "blue"
E(igraph600)$color

#stel de weights in:
E(igraph600)$weights <- weights
max(weights)
normalize <- max(weights)-min(weights)
E(igraph600)$normweights <- (weights - min(weights))/(max(weights) - min(weights))
c_scale <- heat.colors(5)
E(igraph600)$color <- as.factor(E(igraph600)$weights)
E(igraph600)$color <- E(igraph600)$normweights*255
E(igraph600)[normweights <0.5 ]$color <- c_scale[1]
E(igraph600)[(normweights >= 0.50) & (normweights < 0.60) ]$color <- c_scale[2]
E(igraph600)[(normweights >= 0.60) & (normweights < 0.70) ]$color <- c_scale[3]
E(igraph600)[(normweights >= 0.70) & (normweights < 0.80) ]$color <- c_scale[4]
E(igraph600)[normweights >=0.90 ]$color <- c_scale[5]
as.factor(E(igraph600)$weights)
#maak het plaatje:
meteodag <- igraph600
time.coords <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
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
  plot.igraph(igraph600, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              edge.lty = 1+(E(igraph600)$normweights)*2,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)

### Visualizing fomrat -----------------------------------------------
#arc.strength binnen package bnlearn
#Lower arc.strength -> Large difference in removal! Strong arc.
  #data10d <- as.data.frame(time.coords)
  dag <- hc_edges_loglik_10d_200_400i$networks
  dag <- pc_10d_1e_13$networks$dag
  time.coords <- TimeCoordsAnom_from_Grid(tas_ncep_10d) # behoudt grid atrributes
  dag.data <- as.data.frame(time.coords)
  
  linkstrength <- arc.strength(dag, dag.data)
  linkstrength
  
  #groter dan treshold:
  strongarcs <- as.matrix(linkstrength[linkstrength["strength"] <= -230.00000,,])[,1:2]
  strongarcs
  #alle strengths in matrix format
  linkstrengthmat <- as.matrix(linkstrength)
  linkstrengthmat
  
  #create igraph from bn graph for plotting purpose
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))
  E(igraph)
  
  #identify which indices igraph correspond to indices in bn
  indsame <- c()
  nrow(linkstrength)
  
  for (i in 1:nrow(linkstrength)){
    int <- intersect(which(as_edgelist(igraph)[,1] == linkstrengthmat[i,1]),
                     which(as_edgelist(igraph)[,2] == linkstrengthmat[i,2]))
    indsame[i] <- int
  }
  indsame
  
  #create strengths vector: 
  #low strength corresponds to high weight
  strengths <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    strengths[ind] <- linkstrength[i,3]
  }
  
  strengths[165]
  
  #identify indices sequence with 'strong arcs' in igraph600
  indstrong <- c()
  nrow(strongarcs)
  strongarcs
  for (i in 1:nrow(strongarcs)){
    int <- intersect(which(as_edgelist(igraph)[,1] == unname(strongarcs[i,1])),
                     which(as_edgelist(igraph)[,2] == unname(strongarcs[i,2])))
    indstrong[i] <- int
  }
  indstrong
  sort(indstrong)
  as_edgelist(igraph600)[sort(indstrong),]
  
  #Zwakke links:
  no.indstrong <- as.vector(1:nrow(as_edgelist(igraph)))
  no.indstrong <- no.indstrong[-sort(indstrong)]
  no.indstrong
  
  #STel de kleuren in voor sterk en zwak:
  E(igraph)[ sort(indstrong) ]$color <- "red"
  E(igraph)[ sort(no.indstrong) ]$color <- "blue"
  E(igraph)$color
  
  #stel de strengths en weights in:
  E(igraph)$strengths <- strengths
  max(strengths)
  min(strengths)
  normalize <- max(strengths)-min(strengths)
  E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  E(igraph)$weights
  log((max(strengths) - strengths)/(max(strengths) - min(strengths)))
  c_scale <- rev(topo.colors(5))
  
  E(igraph)[weights <0.5]
  quantile( E(igraph)$weights)
  
  
  #for HC
  E(igraph)[weights <0.1 ]$color <- c_scale[1]
  E(igraph)[(weights >= 0.1) & (weights < 0.2) ]$color <- c_scale[2]
  E(igraph)[(weights >= 0.2) & (weights < 0.3) ]$color <- c_scale[3]
  E(igraph)[(weights >= 0.3) & (weights < 0.4) ]$color <- c_scale[4]
  E(igraph)[weights >=0.4 ]$color <- c_scale[5]
  E(igraph)[weights >= 0.5]$color <- "red"
  
  #quantiles:
  quantile(E(igraph)$weights)[2]

  E(igraph)[weights <= quantile(E(igraph)$weights)[2]]$color <- c_scale[1]
  E(igraph)[weights > quantile(E(igraph)$weights)[2] & weights <=quantile(E(igraph)$weights)[3]]$color <- c_scale[2]
  E(igraph)[weights > quantile(E(igraph)$weights)[3] & weights <=quantile(E(igraph)$weights)[4]]$color <- c_scale[3]
  E(igraph)[weights > quantile(E(igraph)$weights)[4] & weights <=quantile(E(igraph)$weights)[5]]$color <- c_scale[4]
  
  
  #maak het plaatje:
  meteodag <- igraph
  
  if (class(meteodag) == "igraph") 
    igraphDegree <- meteodag
  
  if (class(meteodag) == "graphNEL") 
    igraphDegree <- igraph.from.graphNEL(meteodag)
  
  if (class(meteodag) == "bn") 
    igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))
  
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)
  
  points <- expand.grid(y, x)[2:1]
  data(wrld)
  plot(wrld)
  plot.igraph(igraphDegree, 
              vertex.size = 100,
              vertex.color = "blue",
              vertex.label = NA,
              edge.arrow.size = 0.2,
              edge.lty = 1+(E(igraph)$normweights)*2,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  
  
  legend(-250,90,legend=c(paste0("<",round(quantile(E(igraph)$weights)[2],digits = 2)),
                           paste0("<",round(quantile(E(igraph)$weights)[3],digits = 2)),
                           paste0("<",round(quantile(E(igraph)$weights)[4],digits = 2)),
                           paste0("<",round(quantile(E(igraph)$weights)[5],digits = 2))),
         fill=c_scale[1:4])

  
  
  
### visualizing format PC---------------
  #for PC HELE ANDERE SCHAAL!!!! DENK AAN CRITERIOM STRENGTH
  E(igraph)$strengths <- strengths
  max(strengths)
  min(strengths)
  normalize <- max(strengths)-min(strengths)
  E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  c_scale <- rev(topo.colors(5))
  E(igraph)[weights <0.2 ]$color <- c_scale[1]
  E(igraph)[(weights >= 0.2) & (weights < 0.4) ]$color <- c_scale[2]
  E(igraph)[(weights >= 0.4) & (weights < 0.6) ]$color <- c_scale[3]
  E(igraph)[(weights >= 0.6) & (weights < 0.8) ]$color <- c_scale[4]
  E(igraph)[weights >=0.8 ]$color <- c_scale[5]
  E(igraph)[weights >= 0.9]$color <- "red"
  
  #arc.strength binnen package bnlearn
  #Lower arc.strength -> Large difference in removal! Strong arc.
  dag <- pc_10d_1e_13$networks$dag
  dag.data <- data10d
  time.coords <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
  
  linkstrength <- arc.strength(dag, dag.data, criterion = "zf")
  linkstrength
  
  #groter dan treshold:
  strongarcs <- as.matrix(linkstrength[linkstrength["strength"] <= 1e-110,,])[,1:2]
  strongarcs
  #alle strengths in matrix format
  linkstrengthmat <- as.matrix(linkstrength)
  linkstrengthmat
  
  #create igraph from bn graph
  igraph <- igraph.from.graphNEL(as.graphNEL(dag))
  E(igraph)
  
  #identify which indices igraph correspond to indices in bn
  indsame <- c()
  nrow(linkstrength)
  
  for (i in 1:nrow(linkstrength)){
    int <- intersect(which(as_edgelist(igraph)[,1] == linkstrengthmat[i,1]),
                     which(as_edgelist(igraph)[,2] == linkstrengthmat[i,2]))
    indsame[i] <- int
  }
  indsame
  
  #create strengths vector: 
  #low strength corresponds to high weight
  strengths <- numeric(nrow(as_edgelist(igraph)))
  for (i in 1:nrow(as_edgelist(igraph))){
    ind <- indsame[i]
    as_edgelist(igraph)[ind,]
    strengths[ind] <- linkstrength[i,3]
  }
  
  strengths[165]
  
  #identify indices sequence with 'strong arcs' in igraph600
  indstrong <- c()
  nrow(strongarcs)
  strongarcs
  for (i in 1:nrow(strongarcs)){
    int <- intersect(which(as_edgelist(igraph)[,1] == unname(strongarcs[i,1])),
                     which(as_edgelist(igraph)[,2] == unname(strongarcs[i,2])))
    indstrong[i] <- int
  }
  indstrong
  sort(indstrong)
  as_edgelist(igraph600)[sort(indstrong),]
  
  #Zwakke links:
  no.indstrong <- as.vector(1:nrow(as_edgelist(igraph)))
  no.indstrong <- no.indstrong[-sort(indstrong)]
  no.indstrong
  
  #STel de kleuren in voor sterk en zwak:
  E(igraph)[ sort(indstrong) ]$color <- "red"
  E(igraph)[ sort(no.indstrong) ]$color <- "blue"
  E(igraph)$color
  
  #stel de strengths en weights in:
  E(igraph)$strengths <- strengths
  max(strengths)
  min(strengths)
  normalize <- max(strengths)-min(strengths)
  E(igraph)$weights <- (max(strengths) - strengths)/(max(strengths) - min(strengths))
  c_scale <- rev(topo.colors(5))
  
  E(igraph)[weights <0.5]
  
  
  
  #for HC
  E(igraph)[weights <0.5]$color <- c_scale[1]
  E(igraph)[(weights >= 0.5) & (weights < 0.7) ]$color <- c_scale[2]
  E(igraph)[(weights >= 0.7) & (weights < 0.8) ]$color <- c_scale[3]
  E(igraph)[(weights >= 0.8) & (weights < 0.9) ]$color <- c_scale[4]
  E(igraph)[weights >=0.9 ]$color <- c_scale[5]
  E(igraph)[weights >= 0.5]$color <- "red"
  
  #maak het plaatje:
  meteodag <- igraph
  
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
              edge.arrow.size = 0.2,
              edge.lty = 1+(E(igraph)$normweights)*2,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  
  
  legend(-300,150,legend=c(" <0.1","0.1-0.2","0.2-0.3","0.3-0.4",">0.4"),fill=c_scale)
  
  
  
#opslaan Show strengths
plotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/example_linkstrength.pdf"
pdf(plotname)
par(mfrow = c(2,1))
#nu de twee plaatjes implementeren
dev.off()



#blabla tecniek
#create graph with only strongarcs
E(graph_from_edgelist(unname(strongarcs)))

#TEST: INdeed they are the same
as_edgelist(igraph600)[396,]
unname(strongarcs[1,])
all.equal(as_edgelist(igraph600)[396,],
          unname(strongarcs[1,])) #TRUE

########## bn.fit for HC compared with Gauassian denstiy from CN
grid =tas_ncep_10d
th = 0.4

graph_from_Gridb <- function(grid, th = 0.8) {
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
  mean.vector <- colMeans(time.coords.matrix)
  
  dmvnorm
  
  # Adjacency matrix
  diag(adj.matrix) <- 0
  adj.matrix[adj.matrix <= th ] <- 0
  adj.matrix[adj.matrix > th ] <- 1
  # Graph
  graph <- graph_from_adjacency_matrix(adj.matrix, mode = "undirected")
  
  out <- list("graph" = graph,
              "data_coords" = time.coords.matrix,
              "correlation" = cor.matrix,
              "VertexCoords" = ref.coords ,
              "adjacency" = adj.matrix)
  
  attr(out, "Xcoords") <- x
  attr(out, "Ycoords") <- y
  attr(out, "ref.dates") <- ref.dates
  return(out)
}
grid =tas_ncep_10d
th = 0.4

fit <- bn.fit(hc_edges_loglik_10d_1000_1200i$networks,data10d)
bn.fit.qqplot(fit$V1)
bn.fit.histogram(fit$V1)
bn.fit.xyplot(fit$V1)
bn.fit.barchart(fit$V1)
bn.fit()
fitted(fit$V1)
residuals(fit$V1)
coefficients(fit$V1)



test10d_1_e_200, test10d_1_e_160, test10d_1_e_150, test10d_1_e_113,test10d_1_e_100,test10d_1_e_81 ,test10d_1_e_71 , test10d_1_e_70, test10d_4_e_70, test10d_5_e_70, test10d_1_e_67, test10d_1_e_65, test10d_1_e_63, test10d_1_e_61 ,test10d_1_e_45, test10d_1_e_43

