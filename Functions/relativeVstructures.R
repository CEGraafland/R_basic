###########################
# Codigo para un solo red 
###########################
library(bnlearn)
# Select graph
g <- perm1sort$hc1_1400_1500i
# Give names (To avoid problems)
names <- c()
for(i in 1:648) names <- append(names,paste0("V",i))
# obtain vstructs immoral (with BNLEARN)
voefimmor <- vstructs(g, moral = FALSE)
# Which of those vstructs leads to dseperation?
dsepvimmor <- logical()
for(i in 1:nrow(voefimmor)){
  answer <- dsep(g,voefimmor[i,1],voefimmor[i,3])
  dsepvimmor[i] <- answer
}
blockedim <- which(dsepvimmor)
# Matrix with Dseperating Vstructs
permv <- voefimmor[blockedim,]

###########################################################################################
# Now calculate all paths of length 2 in the graph 
degarray <- array(dim = c(length(nodes(g))))
twopathsarray <- array(dim = c(length(nodes(g))))

# per node 
deg <- c()
twopaths <- c()
for(i in 1:length(nodes(g))){
  node <- nodes(g)[i]
  deg[i] <- bnlearn::degree(g,node)
  if(deg[i]>1){
    twopaths[i] <- ncol(combn(deg[i],2))
  } else{twopaths[i] <- 0}
}
degarray <- deg
twopathsarray <- twopaths

# All possible two paths in the graph are the sum of possible twopaths per node. 
possible2pathsmatrix <-  sum(twopathsarray)

# relative d seperating vstructs
nvstructs <- nrow(permv)
relnumbervstructs <- nvstructs/possible2pathsmatrix



#################################################################################
# Function:
###################################################################################
library(bnlearn)
dag <- h2pc_1_eBIC$networks$h2pc_10d_g0.1
dag <- hc_1_eBIC$networks$hc_10d_g1
dag <- pc_1$networks
dag <- h2pc_1_eBIC$networks[[7]]
dag <- eqclass
vstructs(dag)
vstructsprop(dag)
narcs(dag)
wgt

dags <- h2pc_1_eBIC$networks

vlist <- list()
for (i in 1:length(h2pc_1_eBIC$networks)){
  vlist[[i]] <- vstructsprop(dags[[i]])
}

vstructsprop <- function(dag){
  # Select graph
  g <- dag

  # Give names (To avoid problems)
  names <- c()
  for(i in 1:648) names <- append(names,paste0("V",i))
  # obtain vstructs immoral (with BNLEARN)
  voefimmor <- vstructs(g, moral = FALSE)
  nvoefimmor <- nrow(voefimmor)
  # Which of those vstructs leads to dseperation?
  if (!nvoefimmor == 0){
    dsepvimmor <- logical()
    for(i in 1:nrow(voefimmor)){
      answer <- bnlearn::dsep(g,voefimmor[i,1],voefimmor[i,3])
      dsepvimmor[i] <- answer
    }
  blockedim <- which(dsepvimmor)
  # Matrix with Dseperating Vstructs
  permv <- voefimmor[blockedim,]
  } else {permv <- blockedim <- dsepvimmor <- 0}
  ###########################################################################################
  # Now calculate all paths of length 2 in the graph 
  degarray <- array(dim = c(length(nodes(g))))
  twopathsarray <- array(dim = c(length(nodes(g))))

  # per node 
  deg <- c()
  twopaths <- c()
  for(i in 1:length(nodes(g))){
    node <- nodes(g)[i]
    deg[i] <- bnlearn::degree(g,node)
    if(deg[i]>1){
      twopaths[i] <- ncol(combn(deg[i],2))
    } else{twopaths[i] <- 0}
  }
  degarray <- deg
  twopathsarray <- twopaths

  # All possible two paths in the graph are the sum of possible twopaths per node. 
  possible2pathsmatrix <-  sum(twopathsarray)

  # relative d seperating vstructs
  if(length(blockedim) ==0){
    nvstructs <- 0
  } else if(length(blockedim) ==1){
      nvstructs <- 1
      } else {nvstructs <- nrow(permv)}

  relnumbervstructs <- nvstructs/possible2pathsmatrix

  return(c(nvoefimmor, nvstructs,relnumbervstructs))
}

lapply(h2pc_1_eBIC$networks,vstructsprop)


vstructssimpprop <- function(dag, moral = FALSE){
  # Select graph
  g <- dag
  
  # Give names (To avoid problems)
  names <- c()
  for(i in 1:648) names <- append(names,paste0("V",i))
  # obtain vstructs immoral (with BNLEARN)
  voefimmor <- vstructs(g, moral = moral)
  nvoefimmor <- nrow(voefimmor)
  # # Which of those vstructs leads to dseperation?
  # if (!nvoefimmor == 0){
  #   dsepvimmor <- logical()
  #   for(i in 1:nrow(voefimmor)){
  #     answer <- bnlearn::dsep(g,voefimmor[i,1],voefimmor[i,3])
  #     dsepvimmor[i] <- answer
  #   }
  #   blockedim <- which(dsepvimmor)
  #   # Matrix with Dseperating Vstructs
  #   permv <- voefimmor[blockedim,]
  # } else {permv <- blockedim <- dsepvimmor <- 0}
  ###########################################################################################
  # Now calculate all paths of length 2 in the graph 
  degarray <- array(dim = c(length(nodes(g))))
  twopathsarray <- array(dim = c(length(nodes(g))))
  
  # per node 
  deg <- c()
  twopaths <- c()
  for(i in 1:length(nodes(g))){
    node <- nodes(g)[i]
    deg[i] <- bnlearn::degree(g,node)
    if(deg[i]>1){
      twopaths[i] <- ncol(combn(deg[i],2))
    } else{twopaths[i] <- 0}
  }
  degarray <- deg
  twopathsarray <- twopaths
  
  # All possible two paths in the graph are the sum of possible twopaths per node. 
  possible2pathsmatrix <-  sum(twopathsarray)
  
  # # relative d seperating vstructs
  # if(length(blockedim) ==0){
  #   nvstructs <- 0
  # } else if(length(blockedim) ==1){
  #   nvstructs <- 1
  # } else {nvstructs <- nrow(permv)}
  
  relnumbersimpvstructs <- nvoefimmor/possible2pathsmatrix
  
  return(c(nvoefimmor,relnumbersimpvstructs))
}
