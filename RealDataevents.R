#########################################################################
#Real data events tas_ncep_10d
#########################################################################
setwd("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/exp_GCMs/")
setwd("~/data/Untitled/Trabajo/R_practice/exp_GCMs/")
source("../R/Functions/BasicNetworkFunctions.R")
source("../R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("../R/Functions/propagationFunctions.R")
load("../Data/Struct_learn/datapermutations.rda")
load("../Data/Struct_learn/permutations.rda")
load("../Data/tas_ncep_10d.rda")
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
library(visualizeR)
library('RColorBrewer')
library(gridExtra)
col.blue <- rev(brewer.pal(4,"Blues"))
col.blue
col.red <- brewer.pal(4,"Reds")

col.div <- c(col.blue, col.red)
col.div[4:5] <- "white"
col.b <- c(col.blue,col.red)
length(col.b)

##########################################################################
# select simple / combine with paperScript single data evaluation. 
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

# worst 1:.....  or best .....:360
whichdates <- 345:360
# From CN / BN or differences
sizeCN <- 1
nedgesnetworksCM[sizeCN]
whichorder <- orderCNs[,sizeCN]
# no
sizeBN <- 9
whichorder <- orderBNs[,sizeBN]
whichorder <- orderdif
# The event indices:
whichorder[whichdates]

plots <- list()
for(i in 1:length(whichdates)){
anomaliesdate <- dataRMS2[whichorder[whichdates[i]],]

climBig2 <- quantity2clim(as.matrix(anomaliesdate), what = "Big 2",tas_ncep_10d)
plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                          col.regions = col.div,set.min = -2.5, set.max = 2.5,
                          at = seq(-3,3,1))
}

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/worst16_CN",nedgesnetworksCM[sizeCN],".pdf")
pdf(plotname,9,5)
do.call(grid.arrange,c(plots, top = paste0("worst CN",nedgesnetworksCM[sizeCN])))
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/worst16_BN",nedgesIT[sizeBN],".pdf")
pdf(plotname,9,5)
do.call(grid.arrange,c(plots, top = paste0("worst BN",nedgesIT[sizeBN])))
dev.off()

nedgesIT


where <- c()
for(i in 1:length(orderBN)){
  where[i] <-  which(orderBN == orderCM[i])
}
orderBN[where]
##############################################################################
# Difference in "apreciation" full CN and 1800 BN
##############################################################################
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/16bestCNfullaboveBN1800.pdf")
pdf(plotname,9,5)
do.call(grid.arrange,c(plots, top = "16 best CN full above BN 1800")) #1795
dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/worstDataevents/16bestBN1800aboveCNfull.pdf")
pdf(plotname,9,5)
do.call(grid.arrange,c(plots, top = "16 best BN 1800 above CN full")) #1795
dev.off()


##########################################################################
# V81 V280
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)
  
Big2ind <-  which(dataRMS2$V81>=1.5 & dataRMS2$V280 >=1.5)
Big2ind
event1 <- 1:5
event2 <- 6:7
Big2ind <- Big2ind[event1]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])

  
climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))
##########################################################################
# V81 V227
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

Big2ind <-  which(dataRMS2$V81>=1.8 & dataRMS2$V227 >=1.8)
Big2ind
event1 <- 1
event2 <- 2:3
Big2ind <- Big2ind[event1]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))

##########################################################################
# V81 V424
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

Big2ind <-  which(dataRMS2$V81>=1.8 & dataRMS2$V424 >=1.8)
Big2ind
event1 <- 2:6
event2 <- 1
event3 <- 7
event4 <- 12
Big2ind <- Big2ind[event1]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))

##########################################################################
# V280 V424
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

Big2ind <-  which(dataRMS2$V280>=1.8 & dataRMS2$V424>= 1.8)
Big2ind
event1 <- 1:3
event2 <- 4

Big2ind <- Big2ind[event1]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))

##########################################################################
# V280 V497
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

Big2ind <-  which(dataRMS2$V280>=2 & dataRMS2$V497>= 2)
Big2ind
event1 <- 1
event2 <- 2:4

Big2ind <- Big2ind[event2]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))
##########################################################################
# V81
##########################################################################
dataRMS2 <- as.data.frame(dataRMS)

Big2ind <-  which(dataRMS2$V81>=2)
Big2ind
event1 <- 1:4
event2 <- 6:9
event3 <- 11:15

Big2ind <- Big2ind[event3]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))
############################################################
# V85
############################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V85 >=2)
Big2ind
Big2ind <- Big2ind[5:10]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))
############################################################
# V227 
############################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V227 >=2)
Big2ind
event1 <- 1
event2 <- 2:8
event3 <- 9
event4 <- 10:15
Big2ind <- Big2ind[event1]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))

############################################################
# V280
############################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V280 >=2)
Big2ind
event1 <- 1:3
event2 <- 6:15
event3 <- 4
event4 <- 5
Big2ind <- Big2ind[event2]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))


###########################################################################
# V459 
###########################################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind <- which(dataRMS2$V81 >=0.5& dataRMS2$V171 >= 2)
Big2ind
event0 <- 12
event1 <- 13:17
event2 <- 3:8
event3 <- 10:12
Big2indevent <- Big2ind
dataRMS2[Big2indevent,]
extremes <- colMeans(dataRMS2[Big2indevent,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))
###########################################################################
# V424
###########################################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind <- which(dataRMS2$V81 <=-1.5& dataRMS2$V227 >= 1.5)
dataRMS2$V81[Big2ind]
dataRMS2$V227[Big2ind]
Big2ind <- which(dataRMS2$V28 >=1.8& dataRMS2$V280 >= 1.8)
Big2ind <- which(dataRMS2$V153 >=1.5& dataRMS2$V639 <=-1)
Big2ind <- which(dataRMS2$V424 >= 1.8)
Big2ind <- which(dataRMS2$V205 >= 1.5)
Big2ind <- which(dataRMS2$V532 >= 1.8)
Big2ind <- which(dataRMS2$V28>= 1.5)
Big2ind <- which(dataRMS2$V81 >=1.5& dataRMS2$V171 >= -1 $ dataRMS2$V171 <= 1)
Big2ind <- which(dataRMS2$V81 >=0.5& dataRMS2$V171 >= 2) # WP
Big2ind <- which(dataRMS2$V81 >=2& dataRMS2$V171 <= 1) # Cold Tong
Big2ind <- which(dataRMS2$V171 >= 2)
Big2ind <- which(dataRMS2$V227 >= 2)
Big2ind <- which(dataRMS2$V568 >= 1.8)
Big2ind <- which(dataRMS2$V621 >= 1)
Big2ind <- which(dataRMS2$V621 <= -2)
Big2ind <- which(dataRMS2$V205 <= -1.5)
Big2ind <- which(dataRMS2$V568 <= -1.8)
Big2ind <- which(dataRMS2$V459 >= 1.8)
Big2ind <- which(dataRMS2$V459 <= -1.8)
Big2ind <- which(dataRMS2$V81 <= -1.8)

events <- list()
events[[1]] <- c(Big2ind[1])
m <- 1
if (!length(Big2ind)==1){
  for (k in 2:length(Big2ind)){
    if ((Big2ind[k]-1) == Big2ind[k-1]){events[[m]] <- append(events[[m]], Big2ind[k])
    } else {events[[m+1]] <- c(Big2ind[k])
            m <- m+1}
  }
} 
plots <- list()
for(i in 1:length(events)){
extremes <- colMeans(dataRMS2[events[[i]],])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = list(paste0(length(events[[i]]))))
}

do.call(grid.arrange,c(plots))
350/12
###########################################################################
# V424
###########################################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind <- which(dataRMS2$V424 >=1.5)
Big2ind
event0 <- 2:3
event1 <- 5:10
event2 <- 14:18
event4 <- 19:20

Big2indevent <- Big2ind
dataRMS2[Big2indevent,]
extremes <- colMeans(dataRMS2[Big2indevent,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))



##########################################################################
# V568
###########################################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind <- which(dataRMS2$V568 >=1.8)
Big2ind
event0 <- 4:13
event1 <- 46:50
event2 <- 5
Big2indevent <- Big2ind
dataRMS2[Big2indevent,]
extremes <- colMeans(dataRMS2[Big2indevent,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))
############################################################
# V205
############################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V205 >=  1.5)
Big2ind
Big2ind <- Big2ind[3]
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))

############################################################
# V316
############################################################
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V316 >=1.5)
Big2ind <- Big2ind[14:17]

dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1))

############################################################
# - and + Simple All data
############################################################
par(mfrow = c(1,2))
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V81 >=1.8)
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
warm <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = paste0("V81 >= 1.8: ",length(Big2ind),"-month mean in 1981 - 2010"))


Big2ind<- which(dataRMS2$V81 <=-1.8)
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
cold <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = paste0("V81 <= -1.8: ",length(Big2ind),"-month mean in 1981 - 2010"))

do.call(grid.arrange, c(list(warm,cold)))
############################################################
# - and + Simple All data V568 V550 V532
############################################################

par(mfrow = c(1,2))
dataRMS2 <- as.data.frame(dataRMS)
Big2ind<- which(dataRMS2$V550 >=1.8)
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
warm <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                    col.regions = col.div,set.min = -2.5, set.max = 2.5,
                    at = seq(-3,3,1),
                    main = paste0("V550 >= 1.8: ",length(Big2ind),"-month mean in 1981 - 2010"))


Big2ind<- which(dataRMS2$V550 <=-1.8)
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])

climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
cold <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                    col.regions = col.div,set.min = -2.5, set.max = 2.5,
                    at = seq(-3,3,1),
                    main = paste0("V550 <= -1.8: ",length(Big2ind),"-month mean in 1981 - 2010"))

do.call(grid.arrange, c(list(warm,cold)))
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropVFavourBN/DataV550.pdf"
pdf(plotname)
dev.off()
############################################################
# With half of data
############################################################
dataRMS3 <- dataRMS2[1:180,]
dataRMS4 <- dataRMS2[181:360,]

Big2ind <- which(dataRMS4$V81 >=1.4& dataRMS4$V227 >= 1.4)
Big2ind <- which(dataRMS3$V81 <=-1.8)
Big2ind <- which(dataRMS3$V81 >=1.5)
Big2ind <- which(dataRMS3$V227 >=1)
dataRMS4[,c(81,227)][c(170,171),]

events <- list()
events[[1]] <- c(Big2ind[1])
m <- 1
if (!length(Big2ind)==1){
  for (k in 2:length(Big2ind)){
    if ((Big2ind[k]-1) == Big2ind[k-1]){events[[m]] <- append(events[[m]], Big2ind[k])
    } else {events[[m+1]] <- c(Big2ind[k])
    m <- m+1}
  }
} 
plots <- list()
for(i in 1:length(events)){
  extremes <- colMeans(dataRMS3[events[[i]],])
  
  climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
  plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                            col.regions = col.div,set.min = -2.5, set.max = 2.5,
                            at = seq(-3,3,1),
                            main = list(paste0(length(events[[i]]))))
}

do.call(grid.arrange,c(plots))


############################################################
# With half of data V81: Simple
############################################################
# First 15 year
dataRMS2 <- as.data.frame(dataRMS)
dataRMS3 <- dataRMS2[1:180,]
dataRMS4 <- dataRMS2[181:360,]

Big2ind <- which(dataRMS3$V81 >=1.5)

dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS3[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = paste0("V81 >= 1.5: ",length(Big2ind),"-month mean in 1981 - 1995"))


Big2ind <- which(dataRMS4$V81 >=1.5)
dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS4[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = paste0("V81 >= 1.5: ",length(Big2ind),"-month mean in 1996 - 2010"))

Big2ind <- which(dataRMS4$V81 >=1.4& dataRMS4$V227 >= 1.4)
Big2ind <- which(dataRMS3$V81 >=1.5)
Big2ind <- which(dataRMS3$V227 >=1)
dataRMS4[,c(81,227)][c(170,171),]

events <- list()
events[[1]] <- c(Big2ind[1])
m <- 1
if (!length(Big2ind)==1){
  for (k in 2:length(Big2ind)){
    if ((Big2ind[k]-1) == Big2ind[k-1]){events[[m]] <- append(events[[m]], Big2ind[k])
    } else {events[[m+1]] <- c(Big2ind[k])
    m <- m+1}
  }
} 
plots <- list()
for(i in 1:length(events)){
  extremes <- colMeans(dataRMS3[events[[i]],])
  
  climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
  plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                            col.regions = col.div,set.min = -2.5, set.max = 2.5,
                            at = seq(-3,3,1),
                            main = list(paste0(length(events[[i]]))))
}

do.call(grid.arrange,c(plots))


############################################################
# With all data: Simple
############################################################
node <- "V205"
value <- 1.5

Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," >=",value))))

dataRMS2[Big2ind,]
extremes <- colMeans(dataRMS2[Big2ind,])


climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = paste0(node," >= ",value, ": ",length(Big2ind),"-month mean in 1981 - 2010"))


Big2indneg <- which(eval(parse(text = paste0("dataRMS2$",node," <=",-value))))
extremesneg <- colMeans(dataRMS2[Big2indneg,])


climBig2neg <- quantity2clim(extremesneg, what = "Big 2",tas_ncep_10d)
spatialPlot(climBig2neg, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
            col.regions = col.div,set.min = -2.5, set.max = 2.5,
            at = seq(-3,3,1),
            main = paste0(node," <= ", -value,": ",length(Big2ind),"-month mean in 1981 - 2010"))


Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," >=",value))))
Big2ind <- which(eval(parse(text = paste0("dataRMS2$",node," <=",-value))))
events <- list()
events[[1]] <- c(Big2ind[1])
m <- 1
if (!length(Big2ind)==1){
  for (k in 2:length(Big2ind)){
    if ((Big2ind[k]-1) == Big2ind[k-1]){events[[m]] <- append(events[[m]], Big2ind[k])
    } else {events[[m+1]] <- c(Big2ind[k])
    m <- m+1}
  }
} 
plots <- list()
for(i in 1:length(events)){
  extremes <- colMeans(dataRMS2[events[[i]],])
  
  climBig2 <- quantity2clim(extremes, what = "Big 2",tas_ncep_10d)
  plots[[i]] <- spatialPlot(climBig2, lonCenter = 180, backdrop.theme = "coastline", rev.colors = TRUE, 
                            col.regions = col.div,set.min = -2.5, set.max = 2.5,
                            at = seq(-3,3,1),
                            main = list(paste0(length(events[[i]]))))
}

do.call(grid.arrange,c(plots))

cluster_edge_betweenness()