load(file = "/Users/lisettegraafland/Documents/R_practice/Data/ncep/tas_ncep.rda")
tas_ncep

dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
attr(dataRMS,"VertexCoords")[205,]
plot(dataRMS[205,],type = "l")
plot(dataRMS[205,],type = "l")
jan <- seq(1,by = 12, length.out = 30)
feb <- seq(2, by = 12, length.out = 30)
jul <- seq(7, by = 12, length.out = 30)
lines(dataRMS[205,][feb],type = "l", col= "red")
lines(dataRMS[205,][jul],type = "l", col= "orange")
lines(dataRMS[286,],type = "l", col= "red")

grid <- tas_ncep_10d
grid <- tas_ncep

seas <- getSeason(grid)
coords <- getCoordinates(grid)
x <- coords$x
y <- coords$y
ref.coords <- expand.grid(y, x)[2:1]
names(ref.coords) <- c("x", "y")
ref.dates <- getRefDates(grid)
seas.list <- lapply(1:length(seas), function(i) {
  subsetGrid(grid, season = seas[i]) %>% scaleGrid(type = "center") %>% redim(drop = TRUE) })
aux <-  bindGrid(seas.list, dimension = "time")
aux <-  redim(aux, drop = TRUE)

time.coords.matrix <- array3Dto2Dmat(aux$Data)
rms <- TRUE
if (rms == TRUE) {
  time.coords.matrix <- scale(time.coords.matrix, center = FALSE, scale = TRUE)}
if (!is.null(subind)) {
  time.coords.matrix <- time.coords.matrix[subind,]}

temporalPlot(seas.list, latLim = c(-30,-28), lonLim = c(-69,-67))
temporalPlot(aux, latLim = c(-30,-28), lonLim = c(-69,-67))
seas.list

aux$xyCoords
aux$Data

plot(time.coords.matrix[,5578],type = "l")

databigres <- TimeCoordsAnom_from_Grid_rms(tas_ncep, rms = TRUE)
xy_databigres$x <- attr(databigres,"VertexCoords")
which(-69 <xy_databigres$x& xy_databigres$x< -67 & -30 <xy_databigres$y& xy_databigres$y< -28) 
xy_databigres[5578,,]
databigres[,5578]
str(tas_ncep_10d)
str(tas_ncep)
# resolution tas_ncep
tas_ncep$xyCoords$x -> z
a <- numeric()
for(i in 2:length(z)){
  a[i] <- z[i] - z[i-1]
  
}
mean(a,na.rm = TRUE)

temporalPlot(tas_ncep, latLim = c(-30,-28), lonLim = c(-69,-67))
temporalPlot(tas_ncep_10d,latLim = c(-30,-28), lonLim = c(-69,-67))

V205data <- subsetGrid(tas_ncep_10d,latLim = c(-30,-28), lonLim = c(-69,-67))
dim(V205data$Data)
V205data2 <- subsetGrid(tas_ncep,latLim = c(-30,-28), lonLim = c(-69,-67))
V205data2$xyCoords
dim(V205data2$Data)
#########################################################################################
# Load data second time
#########################################################################################
library(loadeR)
tas_ncep_orig <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/ncepReanalysis1/ncepReanalysis1_4xDaily.ncml",
             years = 1981:2010,
             var = "tas",
             time = "DD",
             aggr.d = "mean",
             aggr.m = "mean")
spatialPlot(climatology(tas_ncep))
str(tas_ncep)



######################################################################################
# Load NCEP data.
######################################################################################rm(list = ls())
library(loadeR)
tas_ncep_orig <- loadGridData(dataset = "http://meteo.unican.es/tds5/dodsC/ncepReanalysis1/ncepReanalysis1_4xDaily.ncml",
                              years = 1981:2010,
                              var = "tas",
                              time = "DD",
                              aggr.d = "mean",
                              aggr.m = "mean")
save(tas_ncep_orig, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/ncep/tas_ncep_orig.rda")
# Invastigate V205
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/ncep/tas_ncep_orig.rda")
temporalPlot(tas_ncep_orig, latLim = c(-30,-28), lonLim = c(-69,-67))

# Subset
a <- subsetGrid(tas_ncep_orig, years = 1981:1999, latLim = c(-40, 0), lonLim = c(-80, -60))
b <- subsetGrid(tas_ncep_orig, years = 2000:2010, latLim = c(-40, 0), lonLim = c(-80, -60))

pt <- SpatialPoints(t(as.matrix(c(-68, -29))))

aa <- spatialPlot(climatology(a, clim.fun = list(FUN = "mean")), 
                  at = seq(265, 300, 2), 
                  sp.layout = list(list(pt, first = F, col = "black")),
                  backdrop.theme = "coastline",
                  main = "1981:1999 NCEP")



bb <- spatialPlot(climatology(b, clim.fun = list(FUN = "mean")), 
                  at = seq(265, 300, 2), 
                  sp.layout = list(list(pt, first = F, col = "black")),
                  backdrop.theme = "coastline",
                  main = "2000:2010 NCEP")



grid.arrange(aa, bb, ncol = 1)
#######################################################################################
# Load interim y JRA55
#######################################################################################
library(loadeR)
# Which datasets are available with name interim
UDG.datasets(pattern = "Interim")
# First datasets are reanlysis datasets
UDG.datasets()$name[1:10]
# Asign urls
datasets <- UDG.datasets()$url[1:2]
# asign names
datasetnames <- UDG.datasets()$name[1:2] 
# What is in it?
di <- dataInventory(datasets[1])
# load reanlysis datasets
rean <- lapply(datasets, function(x) loadGridData(x, var = "tas", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean"))
loginUDG("Lisette", "UDG01!")
rean1 <- loadGridData(datasets[1], var = "tas", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean")
rean2 <- loadGridData(datasets[2], var = "tas", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean")
# Method 2
rean <- lapply(1:length(datasets), function(x) loadGridData(datasets[x], var = c("tas", "t2m")[x], years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean", lonLim = c(-80, -60), latLim = c(-40, 0)))
names(rean) <- datasetnames

# Asign 
aa2 <- spatialPlot(climatology(subsetGrid(rean$`ECMWF_ERA-Interim-ESD`, years = 1981:1999)), backdrop.theme = "coastline", main = "1981:1999 Interim", at = seq(265, 300, 2),sp.layout = list(list(pt, first = F, col = "black")))
bb2 <- spatialPlot(climatology(subsetGrid(rean$`ECMWF_ERA-Interim-ESD`, years = 2000:2010)), backdrop.theme = "coastline", main = "2000:2010 Interim", at = seq(265, 300, 2),sp.layout = list(list(pt, first = F, col = "black")))

aa3 <- spatialPlot(climatology(subsetGrid(rean$JRA55, years = 1981:1999), clim.fun = list(FUN = "mean")),
                   main = "1981:1999 JRA55", 
                   at = seq(265, 300, 2), 
                   sp.layout = list(list(pt, first = F, col = "black")),
                   backdrop.theme = "coastline")



bb3 <- spatialPlot(climatology(subsetGrid(rean$JRA55, years = 2000:2010), clim.fun = list(FUN = "mean")), main = "2000:2010 JRA55", 
                   at = seq(265, 300, 2), 
                   sp.layout = list(list(pt, first = F, col = "black")),
                   backdrop.theme = "coastline")

# SpatialPlots All in one
grid.arrange(aa,bb,aa2,bb2,aa3,bb3, nrow = 3)
# TemporalPlots all in one
t1 <- temporalPlot("NCEP"= tas_ncep_orig, latLim = c(-30,-28), lonLim = c(-69,-67), xyplot.custom = list(ylim = c(270,290)))
t2 <- temporalPlot("Interim"=rean$`ECMWF_ERA-Interim-ESD`, latLim = c(-30,-28), lonLim = c(-69,-67), xyplot.custom = list(ylim = c(270,290)))
t3 <- temporalPlot("JRA55"=rean$JRA55, latLim = c(-30,-28), lonLim = c(-69,-67), xyplot.custom = list(ylim = c(270,290)))
grid.arrange(t1,t2,t3)

#######################################################################################
# Correlation V205 with rest. 
#######################################################################################
data10 <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
all.equal(cor(data10),cov(data10))
cor10 <- cor(data10)
corV205 <- quantity2clim(cor10[205,],what = "correlation V205",tas_ncep_10d)
spatialPlot(corV205, backdrop.theme = 'coastline', lonCenter = 180, main = "Correlation with V205", color.theme = "RdBu", at = seq(-1,1,0.1), rev.colors = TRUE)

##########################################################################################
# Load Interim data
##########################################################################################
library(loadeR)
loginUDG("Lisette", "UDG01!")
# Which datasets are available with name interim
UDG.datasets(pattern = "Interim")
# First datasets are reanlysis datasets
UDG.datasets()$name[1:10]
# Asign urls
JRA55data <- UDG.datasets()$url[1]
Interimdata <- UDG.datasets()$url[2]
# asign names
JRA55name <- UDG.datasets()$name[1] 
Interimname <- UDG.datasets()$name[2] 
# What is in it?
di <- dataInventory(Interimname)
# load reanlysis datasets
Interim <- loadGridData(Interimdata, var = "2T", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean")
JRA55 <- loadGridData(JRA55data, var = "tas", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean")

rean1 <- loadGridData(datasets[1], var = "tas", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean")
rean2 <- loadGridData(datasets[2], var = "tas", years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean")
# Method 2
rean <- lapply(1:length(datasets), function(x) loadGridData(datasets[x], var = c("tas", "t2m")[x], years = 1981:2010, time = "DD", aggr.d = "mean", aggr.m = "mean", lonLim = c(-80, -60), latLim = c(-40, 0)))
names(rean) <- datasetnames