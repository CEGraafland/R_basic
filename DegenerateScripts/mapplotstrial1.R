install.packages("rgdal")
install.packages("mapproj")
library(rgdal)

coastline_pol <- readOGR("path", layer = "filename without extension")

map('world', fill = TRUE, col = 1:10, wrap=c(-180,180) )
map('world', fill = TRUE, col = 1:10, wrap=c(-180,180) )


map("world", projection="rectangular", parameter=0, 
    orientation=c(90,0,180), wrap=c(90,), fill=T, resolution=0,col=0)

plot.map<- function(database,center,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  Obj[[1]] <- polygons[,1]
  Obj[[2]] <- polygons[,2]
  
  map(Obj,...)
}

plot.map("world", center=180, col="white",bg="gray",
         fill=TRUE,ylim=c(-90,90),mar=c(0,0,0,0))
