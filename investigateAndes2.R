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
