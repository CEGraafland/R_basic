rm(list = ls())
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworld.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")

library(visualizeR)
library(transformeR)
library(igraph)
library(bnlearn)
library(gridExtra)

library(gridGraphics)
library(grid)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
time.coords <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
meteodag <- perm1sort[[16]]
igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))

        
plot.Meteodag.shift <- function(time.coords, meteodag, lis = FALSE){
  if (class(meteodag) == "igraph") 
    igraphDegree <- meteodag
  
  if (class(meteodag) == "graphNEL") 
    igraphDegree <- igraph.from.graphNEL(meteodag)
  
  if (class(meteodag) == "bn") 
    igraphDegree <- igraph.from.graphNEL(as.graphNEL(meteodag))
  
  x <- attr(time.coords, "Xcoords", exact = FALSE)
  y <- attr(time.coords, "Ycoords", exact = FALSE)
  
  x <- (x+360)%%360
  length(x)
  
  points <- expand.grid(y, x)[2:1]
  nodes(meteodag)
  if (lis == TRUE){raster::plot(lisworld)}
  else {raster::plot(wrld)}
  dev.off()
  raster::plot(m)

  plot.igraph(igraphDegree, 
              vertex.size = 1,
              vertex.color = "blue",
              vertex.label = NA,
              vertex.label.cex = 0.7,
              edge.color= "green",
              edge.arrow.size = 0.3,
              edge.lty = 2,
              layout = as.matrix(points), add = TRUE, rescale = FALSE)
  
  plot.igraph(igraphDegree, 
              vertex.size = 1,
              vertex.color = "blue",
              vertex.label = NA,
              vertex.label.cex = 0.7,
              edge.color= "green",
              edge.arrow.size = 0.3,
              edge.lty = 2,
              layout = as.matrix(points), add = FALSE, rescale = TRUE)
  
}


raster::plot(m)


# for(i in 1:length(m)) {
#   ag <- rep(FALSE, length(m@lines[[i]]@Lines[[1]]@coords[,1]))
#   for(j in 1:(length(m@lines[[i]]@Lines[[1]]@coords[,1]) - 1)){
#   if(abs(m@lines[[i]]@Lines[[1]]@coords[j,1] - m@lines[[i]]@Lines[[1]]@coords[j+1,1]) > 300)
#     ag[j] <- TRUE
#   }
#   if(any(ag)) m@lines[[i]]@Lines[[1]]@coords <- m@lines[[i]]@Lines[[1]]@coords[1:which(ag),]
# }
# raster::plot(m)


##################################################################################
# New lisworld
##################################################################################
m <- lisworld
for(i in 1:length(m)) {
  m@lines[[i]]@Lines[[1]]@coords[,1] <- (m@lines[[i]]@Lines[[1]]@coords[,1]+360) %% 360
}
m@bbox[1,] <- c(0, 360) 


for(i in 1:length(m)) {
  ag <- rep(FALSE, length(m@lines[[i]]@Lines[[1]]@coords[,1]))
  for(j in 1:(length(m@lines[[i]]@Lines[[1]]@coords[,1]))){
    if (j == length(m@lines[[i]]@Lines[[1]]@coords[,1])) {
      z <- abs(m@lines[[i]]@Lines[[1]]@coords[j,1] - m@lines[[i]]@Lines[[1]]@coords[1,1]) 
    } else {
      z <- abs(m@lines[[i]]@Lines[[1]]@coords[j,1] - m@lines[[i]]@Lines[[1]]@coords[j+1,1]) 
    }
    if(z > 300) ag[j] <- TRUE
    }
    if(any(ag)) {
      m@lines[[i]]@Lines[[1]]@coords <- m@lines[[i]]@Lines[[1]]@coords[1:which(ag),]
    }
}

raster::plot(m)
save(m, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/lisworldm.rda")


# for(i in 1:length(m)) {
#   ind <- order(m@lines[[i]]@Lines[[1]]@coords[,1])
#     ind2 <- order(m@lines[[i]]@Lines[[1]]@coords[,1], m@lines[[i]]@Lines[[1]]@coords[,2])
#   m@lines[[i]]@Lines[[1]]@coords <- m@lines[[i]]@Lines[[1]]@coords[ind,] 
# }
######################################################################################################




dev.off()
points$Var2
begin <- 1
cutend <- length(points$Var2)/2
cutbegin <- cutend +1
end <- length(points$Var2)
Var2b <- points$Var2[c(cutbegin:end,begin:cutend)]
Var2b

points

points$Var2 <- Var2b

37%/%36
plot.Meteodag(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d), perm1sort[[3]], lis = TRUE)

warnings()

clim <- climatology(tas_ncep_10d)
spatialPlot(clim, backdrop.theme =  "coastline",lonCenter = 180)
spatialPlot(climatology(tas.ncep), backdrop.theme =  "coastline", lonCenter = 180)

getDim <- function(obj) {
  attr(obj[["Data"]], "dimensions")
}
grid <- tas_ncep_10d
lonCenter = 0
if (!is.null(lonCenter)) {
  indcenter <- which(abs(grid$xyCoords$x - lonCenter) == min(abs(grid$xyCoords$x - lonCenter)))
  indcenter <- abs(length(grid$xyCoords$x)/2 - indcenter)
  newlon <- if (indcenter == 0) grid$xyCoords$x else c(tail(grid$xyCoords$x, -indcenter), 
                                                       head(grid$xyCoords$x, indcenter))
  dimNames <-  attr(grid$Data, "dimensions")
  climfun <-  attr(grid$Data, "climatology:fun")
  ind <- match(newlon, grid$xyCoords$x)
  grid$Data <- asub(grid$Data, idx = ind, 
                    dims = which(getDim(grid) == "lon"), drop = FALSE)
  attr(grid$Data, "dimensions") <- dimNames
  attr(grid$Data, "climatology:fun") <- climfun
  grid$xyCoords$x <- newlon
}

