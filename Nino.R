##############################################################################
# Niño complex networks
##############################################################################
# script niño data
#####  
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/ncep/tas_ncep.rda")
str(tas_ncep)
tas_ncep$Dates$end

#import el nino data
url <- url(description = "https://www.esrl.noaa.gov/psd/data/correlation/nina34.data")
data <- readLines(con = url)
typeof(data)
data[1:5]


#select years, build table with names.
period1 <- c(1966,1995)  
yearsperiod1 <- c(1981,1985)
period2 <- c(1971,2000)  
yearsperiod2 <- c(1986,1990)
period3 <- c(1976,2005)  
yearsperiod3 <- c(1991,1995)
period4 <- c(1981,2010)
yearsperiod4 <- c(1996,2000)
period5 <- c(1986,2015)
yearsperiod5 <- c(2001,2005)
period6 <- c(1986,2016)
yearsperiod6 <-c(2006,2016)

#Timeframe
periods <- rbind(period1,period2,period3,period4,period5,period6)
rownames(periods) <- NULL
nrow(periods)
yearsperiods <- rbind(yearsperiod1, yearsperiod2, yearsperiod3, yearsperiod4, yearsperiod5,yearsperiod6)
rownames(yearsperiods) <- NULL

nrowlastperiod <- yearsperiods[nrow(yearsperiods),2] - yearsperiods[nrow(yearsperiods),1] +1

#Make empty run3month matrix 
run3month <- matrix(data = NA, nrow = (nrow(yearsperiods)-1)*5+nrowlastperiod, ncol = 13, dimnames = list(NULL,c("year","DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ")))
run3month[,"year"] <- (yearsperiods[1,1]):(yearsperiods[nrow(yearsperiods),ncol(yearsperiods)])
class(run3month)
run3month <- as.data.frame(run3month)
run3month

#Obtain data run3month in periods of 5 year (with exception of the last)
for(i in 1:nrow(periods)){
  #Grep de 30 yearperiod that is used to produce anamolies
  #for the middle 5 years
  start.indic <- as.character(periods[i,1])
  end.indic <- as.character(periods[i,2])
  start <- grep(start.indic, data)
  end <- grep(end.indic, data)
  data30y = data[start:end]
  data30y = read.table(text = data30y, sep = "", colClasses = "numeric", col.names = c("year","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"))
  

  # Compute anomalies data
  # vector with means:
  means <- colMeans(data30y[2:ncol(data30y)])
  
  # anomalie matrix
  anom <- data30y
  for (j in 1:12){
  mean <- means[j]
  anom[,j+1] <- anom[,j+1] - mean
  }


  # 3 month running matrix: Calculate the mean over
  # every 3 overlapping months for all years 
  run3month5y <- anom
  #Data van december 1979 niet beschikbaar voor DJF
  run3month5y[1,1+1] <- NA
  #DJF: Altijd de waarde van het jaar ervoor van December
  for(k in 2:nrow(run3month5y)){
    run3month5y[k,1+1] = (anom[k-1,1+12] + anom[k,1+1] + anom[k,1+2])/3
  }
  # JFM - OND: Simple equation
  for (l in 2:11){
    run3month5y[1+l] = (anom[,1+(l-1)] + anom[,1+l] + anom[,1+(l+1)])/3}
  # NDJ: Altijd de waarde van het jaar erna van Januari
  for (m in 1:(nrow(run3month5y)-1)){
    run3month5y[m,1 + 12] = (anom[m,1+11] + anom[m,1+12] + anom[m+1,1+1])/3
  }
  # NDJ: Data van Januari 2012 niet beschikbaar
  run3month5y[nrow(run3month5y),1+12] <- NA
  # correcte namen
  colnames(run3month5y) <- c("year","DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ")
  #yearsperiods[2,1]
  #i
  #as.character(yearsperiods[i,1])
  #as.character(yearsperiods[i,2])
  #run3month5y[,"year"]
  begin5y <- which(run3month5y[,"year"] == yearsperiods[i,1])
  end5y <- which(run3month5y[,"year"] == yearsperiods[i,2])
  run3month5y[begin5y:end5y,]
  periodstart <- which(run3month[,"year"] == yearsperiods[i,1])
  periodend <- which(run3month[,"year"] == yearsperiods[i,2])
  run3month[periodstart:periodend,] <- run3month5y[begin5y:end5y,]
}

# run 3 month matrix for 2 decimals
run3monthrounded <- round(run3month,digits = 1)
str(run3monthrounded)
#Saves (rounded) run 3 month matrices 
library(gridExtra)
pdf("Nino.pdf", height=8, width=23)
grid.table(run3month)
dev.off()


pdf("NinoRounded.pdf", height=11, width=13)
grid.table(round(run3month,digits = 1))
dev.off()

#identifying El nino/La nina
#Make one large vector of the anamolies.
test <- run3monthrounded[2:length(run3monthrounded)]
datamonths <- as.vector(t(test))
totalmonths <- length(as.vector(t(test)))

# El Nino: Check if degree > 0.5 for 5 consecutive months 
i<-1
ind <- c()

while (i<totalmonths){
  if (datamonths[i] < 0.5){
    i <- i+1
  } else {
    j <- 1
    while (datamonths[i+j] >= 0.5) {
      j<- j+1}
    
    if (j >= 5){
      ind <- append(ind,i:(i+j-1))
      i = i+j+1}
    else {i <-i+1}
  }
}

# indices vector /matrix: 1 voor el nino month 0 for normal month:
ind
#vectorind <- vector(mode = "numeric", totalmonths)
#vectorind[ind] <- 1
#vectorind
#matrixind<- matrix(data = vectorind, byrow = TRUE, nrow = nrow(run3monthrounded),ncol = 12) 
#matrixind

# and a reverse:
no.ind <- as.vector(1:totalmonths)
no.ind
no.ind <- no.ind[-ind]
no.ind

# Cut indices for specific period of 1981_2010 (Data of ncep)
ind1981_2010 <- ind[ind <= 360]
ind1981_2010
no.ind1981_2010 <- no.ind[no.ind<= 360]
no.ind1981_2010
length(ind1981_2010)
vectorind1981_2010 <- vectorind[1:360]

############################################################
# Measure plots with Niño data
############################################################
# complex measures Nino data 
# 10 d
tau <- 0.5
tauNino <- 0.5
tas10dNino <- graph_from_Grid(tas_ncep_10d, th = tauNino, subind = ind1981_2010)
measures10dNino <- graph2measure(tas10dNino)
rhoNino <- round(measures10dNino$edens,digits = 5)
measure2clim(measures10dNino, what = "betweenness", ref.grid = tas_ncep_10d)
plotClimatology(measure2clim(measures10dNino, what = "localclustering", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dNino, what = "betweenness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dNino, what = "awconnectivity", ref.grid = tas_ncep_10d), backdrop.theme = "coastline", at = seq(0.01,0.25,0.015))
plotClimatology(measure2clim(measures10dNino, what = "closeness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline",at = seq(0.00010,0.0009,0.00001))

tas10d <- graph_from_Grid(tas_ncep_10d, th = tau)
measures10d <- graph2measure(tas10d)
rhoAll <- round(measures10d$edens, digits = 5)
measure2clim(measures10d, what = "betweenness", ref.grid = tas_ncep_10d)
plotClimatology(measure2clim(measures10d, what = "localclustering", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10d, what = "betweenness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10d, what = "awconnectivity", ref.grid = tas_ncep_10d), backdrop.theme = "coastline",at = seq(0.01,0.25,0.015))
plotClimatology(measure2clim(measures10d, what = "closeness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline", at = seq(0.00010,0.00018,0.00001))

tas10dnoNino <- graph_from_Grid(tas_ncep_10d, th = tau, subind = no.ind1981_2010)
measures10dnoNino <- graph2measure(tas10dnoNino)
measures10dnoNino$edens
measure2clim(measures10dnoNino, what = "betweenness", ref.grid = tas_ncep_10d)
plotClimatology(measure2clim(measures10dnoNino, what = "localclustering", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dnoNino, what = "betweenness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline", col.regions = colores2(69))
plotClimatology(measure2clim(measures10dnoNino, what = "awconnectivity", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dnoNino, what = "closeness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")



measure2clim(measures10dNino, what = "localclustering", ref.grid = tas_ncep_10d)
measure2clim(measures10d, what = "localclustering", ref.grid = tas_ncep_10d)
measure2clim(measures10dNino, what = "betweenness", ref.grid = tas_ncep_10d)
measure2clim(measures10d, what = "betweenness", ref.grid = tas_ncep_10d)
measure2clim(measures10dNino, what = "awconnectivity", ref.grid = tas_ncep_10d)
measure2clim(measures10d, what = "awconnectivity", ref.grid = tas_ncep_10d)
measure2clim(measures10dNino, what = "closeness", ref.grid = tas_ncep_10d)
measure2clim(measures10d, what = "closeness", ref.grid = tas_ncep_10d)

#############################################################################
# Compare Niño años with all años 
#############################################################################
# nino vs total : localclustering: "are my friends also friends"
locclus.clim <- makeMultiGrid(measure2clim(measures10dNino, what = "localclustering", ref.grid = tas_ncep_10d),
                              measure2clim(measures10d, what = "localclustering", ref.grid = tas_ncep_10d))
a <- plotClimatology(locclus.clim, backdrop.theme = "coastline", 
                     names.attr = c(paste0("El Niño months rho = ",rhoNino," tau = ",tauNino),paste0("All months rho = ",rhoAll," tau = ",tau)), 
                     main = "Local clustering coefficient",
                     par.strip.text = list(cex = 0.75))
a
# nino vs total: awconnectivity: "proporcion of friends I have" 
awcon.clim <- makeMultiGrid(measure2clim(measures10dNino, what = "awconnectivity", ref.grid = tas_ncep_10d),
                            measure2clim(measures10d, what = "awconnectivity", ref.grid = tas_ncep_10d))
b <- plotClimatology(awcon.clim, backdrop.theme = "coastline",
                     names.attr = c("El Niño months","All months"), 
                     main = "Area weighted connectivity")
b
# nino vs total: betweenness: "proporcion of shortes paths that cross me"
betw.clim <- makeMultiGrid(measure2clim(measures10dNino, what = "betweenness", ref.grid = tas_ncep_10d),
                           measure2clim(measures10d, what = "betweenness", ref.grid = tas_ncep_10d))
c<- plotClimatology(betw.clim, backdrop.theme = "coastline",
                    #col.regions = colores2(69),
                    names.attr = c("El Niño months","All months"),
                    main = "Betweenness centrality")
c
# nino vs. total: closeness; "
close.clim <- makeMultiGrid(measure2clim(measures10dNino, what = "closeness", ref.grid = tas_ncep_10d),
                            measure2clim(measures10d, what = "closeness", ref.grid = tas_ncep_10d))
d <- plotClimatology(close.clim, backdrop.theme = "coastline",
                     names.attr = c("El Niño months","All months"), 
                     main = "Closeness centrality")
d

# opslaan tau = 0,4, tauNino = 0.4425
# tau = 0,5, tauNino = 0.52795
library(gridExtra)
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/complexnetworksnino/meas_nino_total_10d_tau",tau,"_tauNino",tauNino,".pdf")
pdf(plotname)
grid.arrange(a,c,b,d)
dev.off()

logLik(hc_edges_loglik_10d_1000_1200i$networks,data10d)
logLik(bn.fit(hc_edges_loglik_10d_1000_1200i$networks,data10d),data10d)
# To compare: ncep 5d
tas5d <- graph_from_Grid(tas_ncep_5d2, th = .68)
measures5d <- graph2measure(tas5d)
measures5d$edens
plotClimatology(measure2clim(measures5d, what = "localclustering", ref.grid = tas_ncep_5d2), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures5d, what = "betweenness", ref.grid = tas_ncep_5d2), backdrop.theme = "coastline",at = seq(0,15,0.5), col.regions = colores2(69))
plotClimatology(measure2clim(measures5d, what = "awconnectivity", ref.grid = tas_ncep_5d2), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures5d, what = "closeness", ref.grid = tas_ncep_5d2), backdrop.theme = "coastline",at = seq(0.0000025,0.000003,0.00000005))



plotClimatology(clim.betw.test5d, backdrop.theme = "coastline", at = seq(0,13,.5), col.regions = topo.colors(69))
plotClimatology(clim.betw.test5d, backdrop.theme = "coastline", at = seq(0,15,0.5), col.regions = colores2(69))
plotClimatology(clim.close.test5d, backdrop.theme = "coastline",at = seq(0.0000025,0.000003,0.00000005))
plotClimatology(clim.awcon.test5d,backdrop.theme = "coastline")
plotClimatology(clima.lclus.test5d,backdrop.theme = "coastline")

tas10dnoNino <- graph_from_Grid(tas_ncep_10d, th = 0.4, subind = no.ind1981_2010)
measures10dnoNino <- graph2measure(tas10dnoNino)
measures10dnoNino$edens
measure2clim(measures10dnoNino, what = "betweenness", ref.grid = tas_ncep_10d)
plotClimatology(measure2clim(measures10dnoNino, what = "localclustering", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dnoNino, what = "betweenness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dnoNino, what = "awconnectivity", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")
plotClimatology(measure2clim(measures10dnoNino, what = "closeness", ref.grid = tas_ncep_10d), backdrop.theme = "coastline")





###### [Functie maken om totale data te maken <- climAdjust] -------------------------------------
#import el nino data kant en klaar
url2 <- url(description = "http://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/detrend.nino34.ascii.txt")
data2 <- readLines(con = url2)
typeof(data2)
data2[1:5]


#delete head, select years, build table with names.
start.indic <- "1980"
end.indic <- "2011"
start <- grep(start.indic, data2)
start
end <- grep(end.indic, data2)
data2 = data2[start[1]:end[12]]
data2
data2 = read.table(text = data2, sep = "", colClasses = "numeric", col.names = c("year","month","total","climAdjust","anom"))
data2

##anomalies data
#vector with means:
months <- list()

for(i in 1:12){
  monthname <- which(data2['month'] == i)
  months[[i]] = data2[monthname,]
}

meanclimAdjust <- vector()

for (i in 1:12){
meanclimAdjust[i] <- colMeans(months[[i]]['climAdjust'])
}
meanclimAdjust

anom<- matrix(data = NA , nrow = 32, ncol = 12)

for (i in 1:12){
  colum <- months[[i]]['climAdjust'] - meanclimAdjust[i]
  names(colum) <- NULL
  rownames(colum) <- NULL
  as.matrix(colum)
  colum <- as.matrix(colum)
  anom[,i] <- colum
}


# 3 month running matrix
run3month <- anom
#Data van december 1979 niet beschikbaar voor DJF
run3month[1,1] <- NA
#DJF: Altijd de waarde van het jaar ervoor van December
for(j in 2:nrow(run3month)){
  run3month[j,1] = (anom[j-1,12] + anom[j,1] + anom[j,2])/3
}
# JFM - OND: Simple equation
for (i in 2:11){
  run3month[,i] = (anom[,i-1] + anom[,i] + anom[,i+1])/3
  }
# NDJ: Altijd de waarde van het jaar erna van Januari
for (j in 1:(nrow(run3month)-1)){
  run3month[j,12] = (anom[j,11] + anom[j,12] + anom[j+1,1])/3
}
# NDJ: Data van Januari 2012 niet beschikbaar
run3month[nrow(run3month),12] <- NA
# correcte namen
colnames(run3month) <- c("DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ")
run3month


### [Old script Nino] ---------------
#tas_ncep$Dates
str(tas_ncep)
which(is.na(tas_ncep_10d$Data))

#import el nino data
url <- url(description = "https://www.esrl.noaa.gov/psd/data/correlation/nina34.data")
data <- readLines(con = url)
typeof(data)
data[1:5]


#delete head, select years, build table with names.
start.indic <- "1980"
end.indic <- "2011"
start <- grep(start.indic, data)
end <- grep(end.indic, data)
data = data[start:end]
data = read.table(text = data, sep = "", colClasses = "numeric", col.names = c("year","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"))


##anomalies data
#vector with means:
means <- colMeans(data[2:ncol(data)])
means

#anomalie matrix
anom <- data
for (i in 1:12){
  mean <- means[i]
  anom[,i+1] <- anom[,i+1] - mean
}


# 3 month running matrix
run3month <- anom

#Data van december 1979 niet beschikbaar voor DJF
run3month[1,1+1] <- NA
#DJF: Altijd de waarde van het jaar ervoor van December
for(j in 2:nrow(run3month)){
  run3month[j,1+1] = (anom[j-1,1+12] + anom[j,1+1] + anom[j,1+2])/3
}
# JFM - OND: Simple equation
for (i in 2:11){
  run3month[1+i] = (anom[,1+(i-1)] + anom[,1+i] + anom[,1+(i+1)])/3}
# NDJ: Altijd de waarde van het jaar erna van Januari
for (j in 1:(nrow(run3month)-1)){
  run3month[j,1 + 12] = (anom[j,1+11] + anom[j,1+12] + anom[j+1,1+1])/3
}
# NDJ: Data van Januari 2012 niet beschikbaar
run3month[nrow(run3month),1+12] <- NA
# correcte namen
colnames(run3month) <- c("year","DJF","JFM","FMA","MAM","AMJ","MJJ","JJA","JAS","ASO","SON","OND","NDJ")
run3month


###################################################################################
# Identify Niño months using terciles of TAS (coldmonths/warmmonths data)
###################################################################################
data <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
# Add variable names to data
names <- c()
names
1:ncol(data)
for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
colnames(data) <- names
str(data)
# Get VertexCoords from data
VertexCoords <- attr(data,"VertexCoords")
VertexCoords$x
# Define which variables are in Nino basin. (Lat: 120 -- 172, Lon: -10 -- 12)
indx <- which(findInterval(VertexCoords$x, vec = c(-172,-120), left.open = TRUE)==1)
indy <- which(findInterval(VertexCoords$y, vec = c(-10,12), left.open = TRUE)==1)
indNino <- which((findInterval(VertexCoords$x, vec = c(-172,-120), left.open = TRUE)==1) &
(findInterval(VertexCoords$y, vec = c(-10,12), left.open = TRUE)==1))
indNino
# Make new data with only Ninomonths
TimeCoordsAnomNino <- data[,indNino]
# Calculate mean, terciles.
# Define indices coldmomths, medmonths, warmmonths: to insert in TimeCoordsAnomFromGrid function.
meanmonthNino <- rowMeans(TimeCoordsAnomNino)
terciles <-quantile(meanmonthNino, probs = c(1/3,2/3))
terciles
coldmonths <- which(meanmonthNino <= terciles[1])
medmonths <- which(meanmonthNino > terciles[1] & meanmonthNino <= terciles[2])
warmmonths <- which(meanmonthNino > terciles[2])
# Check
dataColdbasin2 <- data[coldmonths,]
dataColdbasin <- TimeCoordsAnom_from_Grid(tas_ncep_10d, subind = coldmonths)
all.equal(as.data.frame(dataColdbasin), as.data.frame(dataColdbasin2)) # TRUE 

##################################################################################

