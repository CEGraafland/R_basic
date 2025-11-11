#####################################################################################
# Betweenness moving medias. 
#####################################################################################
test <- runif(648, 0,1)
movingmedias <- function(x){
  vect <- x
  names <- character(length = 648)
  for (i in 1:648){names[i] <- paste0("V",i)}
  
  dimnames <-list("row" = c("r1","r2","r3"),"column"=c("c1",'c2',"c3"),"coord" = names)
  matrixlist <- array(data =NA, dim = c(3,3,648), dimnames = dimnames)
  
  for (i in 1:648){
  if (i<=18|i%%18 == 0) {
    if (i%%18 == 0){ # j1 in matrix
      matrixlist[1,1,i] <- 0
    } else {matrixlist[1,1,i] <- vect[i-17 + 648]} 
  } else {matrixlist[1,1,i] <- vect[i-17]} 
  
  if (i%%18 == 0){ # j2 in matrix
    matrixlist[1,2,i] <- 0
  } else {matrixlist[1,2,i] <-vect[i+1]}
  
  if (i>= 631|i%%18 == 0){ # j3 in matrix
    if (i%%18 == 0){
      matrixlist[1,3,i] <- 0 
    } else  {matrixlist[1,3,i] <- vect[i+19 - 648]}
  } else {matrixlist[1,3,i] <- vect[i+19]}
  
  if (i<=18){ # j4 in matrix
    matrixlist[2,1,i] <- vect[i-18 + 648]
  } else {matrixlist[2,1,i] <- vect[i-18]}
  
  matrixlist[2,2,i] <- vect[i] # j5 in matrix
  
  if (i>= 631){ # j6 in matrix
    matrixlist[2,3,i] <- vect[i+18 -648]
  } else matrixlist[2,3,i] <- vect[i+18]
    
  if (i<=18|i%%18 == 1){ # j7 in matrix
    if (i%%18 == 1) {matrixlist[3,1,i] <- 0
    } else {matrixlist[3,1,i] <- vect[i-19 + 648]}
  } else {matrixlist[3,1,i] <- vect[i-19]}
  
  if (i%%18 == 1){ # j8 in matrix
    matrixlist[3,2,i] <- 0
  } else {matrixlist[3,2,i] <-vect[i-1]}
  
  if (i%%18 == 1 |i>= 631){ # j9 in matrix
    if (i%%18 == 1){
      matrixlist[3,3,i] <- 0
    } else {matrixlist[3,3,i] <- vect[i+17-648]}
  } else {matrixlist[3,3,i] <- vect[i+17]}
  }
  return(matrixlist)
}

all.equal(movingmedias(test),movingmedias(betwmeans25))
movbetw25 <- movingmedias(betwmeans25)
movbetw25matrixmean <- apply(movbetw25, MARGIN = 3, FUN = mean)

memClim <- quantity2clim(movbetw25matrixmean, "betw mov", tas_ncep_10d)
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800\n unweighted", color.theme = "Reds")
                      # set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 25 permutations BN 1800", color.theme = "Reds")
plot1

movbetw1 <- movingmedias(betwmeans1)
movbetw1matrixmean <- apply(movbetw1, MARGIN = 3, FUN = mean)

memClim <- quantity2clim(movbetw1matrixmean, "betw mov", tas_ncep_10d)
plot1 <- spatialPlot(grid = memClim, backdrop.theme = "coastline", set.min = NULL,
                     # set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 1st permutation BN 1800", color.theme = "Reds")
                     set.max = NULL, lonCenter = 180, main = "Moving medias log(1 +Betweenness) 1st permutation BN 1800\n large weight = short distance", color.theme = "Reds")
plot1

# memClim2 <- quantity2clim(log(1+movbetw25matrixmean), "betw mov", tas_ncep_10d)
# plot1 <- spatialPlot(grid = memClim2, backdrop.theme = "coastline", set.min = NULL,
#                      set.max = NULL, lonCenter = 180, main = "Moving medias Betweenness",color.theme = "Reds")
# 
# plot1
#########################################################################################
# To do: apply Backpermutations to betweennesslist Average Networking. 
# Then sum and average over 25 permutations.
# Then do moving media.
#########################################################################################
################################################
# Backpermutations BN 1- 20
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations20.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations20.rda")


permutations25<- c(permutations,permutations20)
datapermutations25 <- c(datapermutations,datapermutations20)
backpermutations25 <- list()
for (j in 1:length(datapermutations25)){
  indback <- c()
  for (i in 1:ncol(datapermutations25[[1]])){
    intback <- which(colnames(datapermutations25[[j]]) == colnames(datapermutations25[[1]][i]))
    indback[i] <- intback
  }
  backpermutations25[[j]] <- indback
}
backpermutations25

backbetwlist <- list()
i<-1
for(i in 1:length(betwlist)){
  backbetwlist[[i]] <- betwlist[[i]][backpermutations25[[i]]]
}

dfbetw <- data.frame(backbetwlist[[1]])
for (i in 2:length(backbetwlist)){
  dfbetw[,i] <- backbetwlist[[i]]
}

betwmeans25 <- rowMeans(dfbetw)
betwmeans1 <- dfbetw[,1]
betwmeans2 <- dfbetw[,2]
###################################################################
# New moving window now in Basicnetworkfunctions
####################################################################
movingmedias <- function(x){
  vect <- x
  names <- character(length = 648)
  for (i in 1:648){names[i] <- paste0("V",i)}
  
  dimnames <-list("row" = c("r1","r2","r3"),"column"=c("c1",'c2',"c3"),"coord" = names)
  matrixlist <- array(data =NA, dim = c(3,3,648), dimnames = dimnames)
  
  for (i in 1:648){
    if (i<=18|i%%18 == 0) {
      if (i%%18 == 0){ # j1 in matrix
        matrixlist[1,1,i] <- NA
      } else {matrixlist[1,1,i] <- vect[i-17 + 648]} 
    } else {matrixlist[1,1,i] <- vect[i-17]} 
    
    if (i%%18 == 0){ # j2 in matrix
      matrixlist[1,2,i] <- NA
    } else {matrixlist[1,2,i] <-vect[i+1]}
    
    if (i>= 631|i%%18 == 0){ # j3 in matrix
      if (i%%18 == 0){
        matrixlist[1,3,i] <- NA 
      } else  {matrixlist[1,3,i] <- vect[i+19 - 648]}
    } else {matrixlist[1,3,i] <- vect[i+19]}
    
    if (i<=18){ # j4 in matrix
      matrixlist[2,1,i] <- vect[i-18 + 648]
    } else {matrixlist[2,1,i] <- vect[i-18]}
    
    matrixlist[2,2,i] <- vect[i] # j5 in matrix
    
    if (i>= 631){ # j6 in matrix
      matrixlist[2,3,i] <- vect[i+18 -648]
    } else matrixlist[2,3,i] <- vect[i+18]
    
    if (i<=18|i%%18 == 1){ # j7 in matrix
      if (i%%18 == 1) {matrixlist[3,1,i] <- NA
      } else {matrixlist[3,1,i] <- vect[i-19 + 648]}
    } else {matrixlist[3,1,i] <- vect[i-19]}
    
    if (i%%18 == 1){ # j8 in matrix
      matrixlist[3,2,i] <- NA
    } else {matrixlist[3,2,i] <-vect[i-1]}
    
    if (i%%18 == 1 |i>= 631){ # j9 in matrix
      if (i%%18 == 1){
        matrixlist[3,3,i] <- NA
      } else {matrixlist[3,3,i] <- vect[i+17-648]}
    } else {matrixlist[3,3,i] <- vect[i+17]}
  }
  return(matrixlist)
}
