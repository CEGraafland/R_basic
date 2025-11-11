#########################################################################################
# Visualize evidence propagation. 
#########################################################################################
rm(list = ls())
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
library(corpcor)
library(bnlearn)
library(transformeR)
library(visualizeR)
library(condMVNorm)
library(RColorBrewer)
library(grid)
library(gridExtra)
################################################
# Backpermutations BN
################################################
backpermutations <- list()
for (j in 1:length(datapermutations)){
  indback <- c()
  for (i in 1:ncol(datapermutations[[1]])){
    intback <- which(colnames(datapermutations[[j]]) == colnames(datapermutations[[1]][i]))
    indback[i] <- intback
  }
  backpermutations[[j]] <- indback
}
##################################################################################
# Complex Network treshold
##################################################################################
tau<- seq(0, 0.99, 0.01)
tau <- c(tau,0.026)
tau <- c(0.87,
         0.8,
         0.7,
         0.66,
         0.62,
         0.52,
         0.49,
         0.41,
         0.345,
         0.32,
         0.3,
         0.2,
         0.1,
         0.06,
         0.05,
         0.0375,
         0.029,
         0.026,
         0)


length(tau)
tauchar <- as.character(tau)
labelsCM <- c()

for(i in 1:length(tauchar)){
  labelsCM[i] <- paste0("= ",tauchar[i])
}



graphs10d <- lapply(tau, graph_from_Grid, grid = tas_ncep_10d, subind = NULL)
graphssolo  <- lapply(graphs10d, function(m) m$graph)

nedgesnetworksCM <- lapply(graphssolo, E)
nedgesnetworksCM <- sapply(nedgesnetworksCM, length)


corssolo <- lapply(graphs10d, function(m) m$correlation)
corcutsolo <- list()
for (i in 1:length(corssolo)){
  corcutsolo[[i]] <- corssolo[[i]]
  corcutsolo[[i]][abs(corcutsolo[[i]]) <= tau[i] ] <- 0
}

logliksCM <- c()
data <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
cormatrix <- cor(data, method = "spearman")
for(i in 1:length(tau)){
  nncor <- make.positive.definite(corcutsolo[[i]], tol = 0.06)
  density <- mvtnorm::dmvnorm(dataRMS, sigma = nncor,log = TRUE)
  logliksCM[i] <- sum(density)
}
####################################################################################
# load HC iteration data and make list
####################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")
load(file ="/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/hc_edges_loglik_10d_200i.rda")
x <- seq(200,8600,200)
y <- seq(400,8800,200)
i <- 200
iterationnetworks <- list()

iterationnetworks[[1]] <- hc_edges_loglik_10d_200i$networks
for(i in 1:length(x)){
  iterationnetworks[[i+1]] <- eval(parse(text = paste0("hc_edges_loglik_10d_",x[i],"_",y[i],"i$networks")))
}
####################################################################################
# hc iteration analyse
####################################################################################3
dataRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
iterationfits <- lapply(iterationnetworks, bn.fit, data = dataRMS)

logliksIT <- c()
for (i in 1:length(iterationfits)){
  logliksIT[i] <- logLik(iterationfits[[i]], data = dataRMS)
}
logliksIT
nedgesIT <- sapply(iterationnetworks, narcs)

data.frame(nedges = nedgesIT, logliks = logliksIT)
plot(nedgesIT,logliksIT, col = "blue", xlim = c(0,3000))



###############################################################################
# Propagation part Integrate functions V205pos
# Bayesian network
###############################################################################
# Load HC propagation iterations V205
fileshciterationV205<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V205pos", full.names = T)
hcV205pos<- list()
nameshcV205pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV205)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV205[i], temp.space),envir = temp.space)
  nameshcV205pos[i] <- ls(envir = temp.space)
  hcV205pos[[i]] <- variable
  rm(temp.space)
}
names(hcV205pos) <- nameshcV205pos
hcV205pos$prop_hc_1000_1200i_V205_equal2

x <- seq(200,2400,200)
y <- seq(400,2600,200)

hcV205posord <- list()
nameshcV205posord <- c()
for(i in 1:length(x)){
  hcV205posord[[i]] <- eval(parse(text =paste0("hcV205pos$prop_hc_",x[i],"_",y[i],"i_V205_equal2")))
  nameshcV205posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V205_equal2")
}
names(hcV205posord) <- nameshcV205posord


# Load HC propagation iterations V227
fileshciterationV227<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227pos", full.names = T)
hcV227pos<- list()
nameshcV227pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV227)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV227[i], temp.space),envir = temp.space)
  nameshcV227pos[i] <- ls(envir = temp.space)
  hcV227pos[[i]] <- variable
  rm(temp.space)
}
names(hcV227pos) <- nameshcV227pos
hcV227pos$prop_hc_1000_1200i_V227_equal2

x <- seq(200,2400,200)
y <- seq(400,2600,200)

hcV227posord <- list()
nameshcV227posord <- c()
for(i in 1:length(x)){
  hcV227posord[[i]] <- eval(parse(text =paste0("hcV227pos$prop_hc_",x[i],"_",y[i],"i_V227_equal2")))
  nameshcV227posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V227_equal2")
}
names(hcV227posord) <- nameshcV227posord

# Load HC propagation iterations V227 neg
pattern <- paste0("propneg_")
fileshciterationV227neg <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227neg", full.names = T, pattern = pattern)
fileshciterationV227neg
hcpropnegV227list <- list()
namesnegV227 <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV227neg )){
  temp.space <- new.env()
  variableneg1 <- get(load(fileshciterationV227neg[i], temp.space),envir = temp.space)
  namesnegV227[i] <- ls(envir = temp.space)
  hcpropnegV227list[[i]] <- variableneg1
  rm(temp.space)
}
names(hcpropnegV227list) <- namesnegV227

x <- seq(200,2200,200)
y <- seq(400,2400,200)

hcpropnegV227listord <- list()
namesnegV227ord <- c()
for(i in 1:length(x)){
  hcpropnegV227listord[[i]] <- eval(parse(text =paste0("hcpropnegV227list$propneg_hc_",x[i],"_",y[i],"i_V227_equal2")))
  namesnegV227ord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V227_equal2")
}
names(hcpropnegV227listord) <- namesnegV227ord

# Load HC propagation iterations V280 
fileshciterationV280<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280pos", full.names = T)
hcV280pos<- list()
nameshcV280pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV280)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV280[i], temp.space),envir = temp.space)
  nameshcV280pos[i] <- ls(envir = temp.space)
  hcV280pos[[i]] <- variable
  rm(temp.space)
}
names(hcV280pos) <- nameshcV280pos
hcV280pos$prop_hc_1000_1200i_V280_equal2

x <- seq(200,1600,200)
y <- seq(400,1800,200)

hcV280posord <- list()
nameshcV280posord <- c()
for(i in 1:length(x)){
  hcV280posord[[i]] <- eval(parse(text =paste0("hcV280pos$prop_hc_",x[i],"_",y[i],"i_V280_equal2")))
  nameshcV280posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V280_equal2")
}
names(hcV280posord) <- nameshcV280posord

# Load HC propagation iterations V280 neg
pattern <- paste0("propneg_")
fileshciterationV280neg <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V280neg", full.names = T, pattern = pattern)
fileshciterationV280neg
hcpropnegV280list <- list()
namesnegV280 <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV280neg )){
  temp.space <- new.env()
  variableneg1 <- get(load(fileshciterationV280neg[i], temp.space),envir = temp.space)
  namesnegV280[i] <- ls(envir = temp.space)
  hcpropnegV280list[[i]] <- variableneg1
  rm(temp.space)
}
names(hcpropnegV280list) <- namesnegV280

x <- seq(200,1800,200)
y <- seq(400,2000,200)

hcpropnegV280listord <- list()
namesnegV280ord <- c()
for(i in 1:length(x)){
  hcpropnegV280listord[[i]] <- eval(parse(text =paste0("hcpropnegV280list$propneg_hc_",x[i],"_",y[i],"i_V280_equal2")))
  namesnegV280ord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V280_equal2")
}
names(hcpropnegV280listord) <- namesnegV280ord

# Load HC propagation iterations V81
fileshciterationV81<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos", full.names = T)
hcV81pos<- list()
nameshcV81pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV81[i], temp.space),envir = temp.space)
  nameshcV81pos[i] <- ls(envir = temp.space)
  hcV81pos[[i]] <- variable
  rm(temp.space)
}
names(hcV81pos) <- nameshcV81pos
hcV81pos$prop_hc_1000_1200i_V81_equal2

x <- seq(200,2800,200)
y <- seq(400,3000,200)

hcV81posord <- list()
nameshcV81posord <- c()
for(i in 1:length(x)){
  hcV81posord[[i]] <- eval(parse(text =paste0("hcV81pos$prop_hc_",x[i],"_",y[i],"i_V81_equal2")))
  nameshcV81posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81_equal2")
}
names(hcV81posord) <- nameshcV81posord

# Load HC propagation iterations V81 neg
pattern <- paste0("propneg_")
fileshciterationV81neg <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81neg", full.names = T, pattern = pattern)
fileshciterationV81neg
hcpropnegV81list <- list()
namesnegV81 <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81neg )){
  temp.space <- new.env()
  variableneg1 <- get(load(fileshciterationV81neg[i], temp.space),envir = temp.space)
  namesnegV81[i] <- ls(envir = temp.space)
  hcpropnegV81list[[i]] <- variableneg1
  rm(temp.space)
}
names(hcpropnegV81list) <- namesnegV81

x <- seq(200,3000,200)
y <- seq(400,3200,200)

hcpropnegV81listord <- list()
namesnegV81ord <- c()
for(i in 1:length(x)){
  hcpropnegV81listord[[i]] <- eval(parse(text =paste0("hcpropnegV81list$propneg_hc_",x[i],"_",y[i],"i_V81_equal2")))
  namesnegV81ord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V81_equal2")
}
names(hcpropnegV81listord) <- namesnegV81ord


# Load HC propagation iterations V568
fileshciterationV568<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V568pos", full.names = T)
hcV568pos<- list()
nameshcV568pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV568)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV568[i], temp.space),envir = temp.space)
  nameshcV568pos[i] <- ls(envir = temp.space)
  hcV568pos[[i]] <- variable
  rm(temp.space)
}
names(hcV568pos) <- nameshcV568pos
hcV568pos$prop_hc_1000_1200i_V568_equal2

x <- seq(200,2000,200)
y <- seq(400,2200,200)

hcV568posord <- list()
nameshcV568posord <- c()
for(i in 1:length(x)){
  hcV568posord[[i]] <- eval(parse(text =paste0("hcV568pos$prop_hc_",x[i],"_",y[i],"i_V568_equal2")))
  nameshcV568posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V568_equal2")
}
names(hcV568posord) <- nameshcV568posord
###############################################################################
# Load HC propagation iterations V81V280
fileshciterationV81V280<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V280pos", full.names = T)
hcV81V280pos<- list()
nameshcV81V280pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81V280)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV81V280[i], temp.space),envir = temp.space)
  nameshcV81V280pos[i] <- ls(envir = temp.space)
  hcV81V280pos[[i]] <- variable
  rm(temp.space)
}
names(hcV81V280pos) <- nameshcV81V280pos
hcV81V280pos$prop_hc_1000_1200i_V81V280_equal2

x <- seq(200,3000,200)
y <- seq(400,3200,200)

hcV81V280posord <- list()
nameshcV81V280posord <- c()
for(i in 1:length(x)){
  hcV81V280posord[[i]] <- eval(parse(text =paste0("hcV81V280pos$prop_hc_",x[i],"_",y[i],"i_V81V280_equal22")))
  nameshcV81V280posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81V280_equal22")
}
names(hcV81V280posord) <- nameshcV81V280posord


# Load HC propagation iterations V81V280 negative
pattern <- paste0("propneg_")
fileshciterationV81V280neg <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V280neg", full.names = T, pattern = pattern)
fileshciterationV81V280neg
hcpropV81V280neglist <- list()
namesV81V280neg <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81V280neg)){
  temp.space <- new.env()
  variableneg <- get(load(fileshciterationV81V280neg[i], temp.space),envir = temp.space)
  namesV81V280neg[i] <- ls(envir = temp.space)
  hcpropV81V280neglist[[i]] <- variableneg
  rm(temp.space)
}
names(hcpropV81V280neglist) <- namesV81V280neg

x <- seq(200,2800,200)
y <- seq(400,3000,200)

hcpropV81V280neglistord <- list()
namesV81V280negord <- c()
for(i in 1:length(x)){
  hcpropV81V280neglistord[[i]] <- eval(parse(text =paste0("hcpropV81V280neglist$propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22")))
  namesV81V280negord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V81V280_equal22")
}
names(hcpropV81V280neglistord) <- namesV81V280negord

nedgesord <- sapply(iterationnetworks, narcs)

###############################################################################
# Load HC propagation iterations V81V227
fileshciterationV81V227<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V227pos", full.names = T)
hcV81V227pos<- list()
nameshcV81V227pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81V227)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV81V227[i], temp.space),envir = temp.space)
  nameshcV81V227pos[i] <- ls(envir = temp.space)
  hcV81V227pos[[i]] <- variable
  rm(temp.space)
}
names(hcV81V227pos) <- nameshcV81V227pos
hcV81V227pos$prop_hc_1000_1200i_V81V227_equal2

x <- seq(200,2200,200)
y <- seq(400,2400,200)

hcV81V227posord <- list()
nameshcV81V227posord <- c()
for(i in 1:length(x)){
  hcV81V227posord[[i]] <- eval(parse(text =paste0("hcV81V227pos$prop_hc_",x[i],"_",y[i],"i_V81V227_equal22")))
  nameshcV81V227posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81V227_equal22")
}
names(hcV81V227posord) <- nameshcV81V227posord

# Load HC propagation iterations V81V227 negative
pattern <- paste0("propneg_")
fileshciterationV81V227neg <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V227neg", full.names = T, pattern = pattern)
fileshciterationV81V227neg
hcpropV81V227neglist <- list()
namesneg <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81V227neg)){
  temp.space <- new.env()
  variableneg <- get(load(fileshciterationV81V227neg[i], temp.space),envir = temp.space)
  namesneg[i] <- ls(envir = temp.space)
  hcpropV81V227neglist[[i]] <- variableneg
  rm(temp.space)
}
names(hcpropV81V227neglist) <- namesneg

x <- seq(200,2200,200)
y <- seq(400,2400,200)

hcpropV81V227neglistord <- list()
namesV81V227negord <- c()
for(i in 1:length(x)){
  hcpropV81V227neglistord[[i]] <- eval(parse(text =paste0("hcpropV81V227neglist$propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22")))
  namesV81V227negord[i] <- paste0("propneg_hc_",x[i],"_",y[i],"i_V81V227_equal22")
}
names(hcpropV81V227neglistord) <- namesV81V227negord

nedgesord <- sapply(iterationnetworks, narcs)


###############################################################################
# Load HC propagation iterations V205V227
fileshciterationV205V227<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V205V227pos", full.names = T)
hcV205V227pos<- list()
nameshcV205V227pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV205V227)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV205V227[i], temp.space),envir = temp.space)
  nameshcV205V227pos[i] <- ls(envir = temp.space)
  hcV205V227pos[[i]] <- variable
  rm(temp.space)
}
names(hcV205V227pos) <- nameshcV205V227pos
hcV205V227pos$prop_hc_1000_1200i_V205V227_equal2

x <- seq(200,2200,200)
y <- seq(400,2400,200)

hcV205V227posord <- list()
nameshcV205V227posord <- c()
for(i in 1:length(x)){
  hcV205V227posord[[i]] <- eval(parse(text =paste0("hcV205V227pos$prop_hc_",x[i],"_",y[i],"i_V205V227_equal22")))
  nameshcV205V227posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V205V227_equal22")
}
names(hcV205V227posord) <- nameshcV205V227posord


###############################################################################
# Load HC propagation iterations V227V568pos
fileshciterationV227V568<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227V568pos", full.names = T)
hcV227V568pos<- list()
nameshcV227V568pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV227V568)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV227V568[i], temp.space),envir = temp.space)
  nameshcV227V568pos[i] <- ls(envir = temp.space)
  hcV227V568pos[[i]] <- variable
  rm(temp.space)
}
names(hcV227V568pos) <- nameshcV227V568pos
hcV227V568pos$prop_hc_1000_1200i_V227V568_equal2

x <- seq(200,2000,200)
y <- seq(400,2200,200)

hcV227V568posord <- list()
nameshcV227V568posord <- c()
for(i in 1:length(x)){
  hcV227V568posord[[i]] <- eval(parse(text =paste0("hcV227V568pos$prop_hc_",x[i],"_",y[i],"i_V227V568_equal22")))
  nameshcV227V568posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V227V568_equal22")
}
names(hcV227V568posord) <- nameshcV227V568posord




##########################################################################################
# Visualize V227V568 pos
##########################################################################################
nedgesord <- sapply(iterationnetworks, narcs)
nedgesord

props <- hcV227V568posord
Center3 <- 180
plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  library(RColorBrewer)
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  colsindex <- rev(brewer.pal(n = 9, "RdBu"))
  cb2 <- colorRampPalette(colsindex) 
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center3, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center3, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center3, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n |E| = ",nedgesord[i]," log = ",round(logliksIT[i])),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  b
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, ncol = 4, nrow = 4))

###############################################################################
# Load HC propagation iterations V81V171 pos
fileshciterationV81V171<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171pos", full.names = T)
hcV81V171pos<- list()
nameshcV81V171pos <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81V171)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV81V171[i], temp.space),envir = temp.space)
  nameshcV81V171pos[i] <- ls(envir = temp.space)
  hcV81V171pos[[i]] <- variable
  rm(temp.space)
}
names(hcV81V171pos) <- nameshcV81V171pos
hcV81V171pos$prop_hc_1200_1400i_V81V171_equal22

x <- seq(1200,1600,200)
y <- seq(1400,1800,200)

hcV81V171posord <- list()
nameshcV81V171posord <- c()
for(i in 1:length(x)){
  hcV81V171posord[[i]] <- eval(parse(text =paste0("hcV81V171pos$prop_hc_",x[i],"_",y[i],"i_V81V171_equal22")))
  nameshcV81V171posord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal22")
}
names(hcV81V171posord) <- nameshcV81V171posord
###############################################################################
# Load HC propagation iterations V81V171 posneg 2 - 2
fileshciterationV81V171posneg<- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171posneg", full.names = T)
hcV81V171posneg<- list()
nameshcV81V171posneg <- c()
temp.space <- new.env()

for (i in 1:length(fileshciterationV81V171posneg)){
  temp.space <- new.env()
  variable <- get(load(fileshciterationV81V171posneg[i], temp.space),envir = temp.space)
  nameshcV81V171posneg[i] <- ls(envir = temp.space)
  hcV81V171posneg[[i]] <- variable
  rm(temp.space)
}
names(hcV81V171posneg) <- nameshcV81V171posneg
hcV81V171posneg$
x <- seq(1200,1600,200)
y <- seq(1400,1800,200)

hcV81V171posnegord <- list()
nameshcV81V171posnegord <- c()
for(i in 1:length(x)){
  hcV81V171posnegord[[i]] <- eval(parse(text =paste0("hcV81V171posneg$prop_hc_",x[i],"_",y[i],"i_V81V171_equal22")))
  nameshcV81V171posnegord[i] <- paste0("prop_hc_",x[i],"_",y[i],"i_V81V171_equal2min2")
}
names(hcV81V171posnegord) <- nameshcV81V171posnegord

##########################################################################################
# Visualize V81V171 posneg 
##########################################################################################
nedgesord <- sapply(iterationnetworks, narcs)
nedgesord
numberBN <- 8
numberCM <-4
nedgesnetworksCM[numberCM]
logliksCM[numberCM]

test1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,171), valueEvidence = c(2,2))
test2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,171), valueEvidence = c(2,-2))
all.equal(test1,test2)

props <- list(hcV81V171posord$prop_hc_1600_1800i_V81V171_equal22,
              hcV81V171posnegord$prop_hc_1600_1800i_V81V171_equal2min2, 
              test1, 
              test2)
nedgesord <-c(nedgesIT[numberBN+1],nedgesIT[numberBN+1],nedgesnetworksCM[numberCM],nedgesnetworksCM[numberCM])
nedgesord
logliks <- c(logliksIT[numberBN+1],logliksIT[numberBN+1],logliksCM[numberCM],logliksCM[numberCM])
props[[3]]
Center3 <- 180
plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  library(RColorBrewer)
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  colsindex <- rev(brewer.pal(n = 9, "RdBu"))
  cb2 <- colorRampPalette(colsindex) 
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  cb <- colorRampPalette(brewer.pal(9, "OrRd"))(80)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center3, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center3, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center3, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n |E| = ",nedgesord[i]," log = ",round(logliks[i])),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  b
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, ncol = 2, nrow = 2))


##################################################################################
# V81V280 pos
##################################################################################


numberCM <- length(corcutsolo)-1
numberBN <- 6


test1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(280,81), valueEvidence = c(2,2))
test2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(280), valueEvidence = c(2))
test3 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))

test4 <- test3
test4$with <- test3$with + test2$with - test2$without
# test4$with <- (test3$with + test2$with)/2
all.equal(test4,test3)

test5 <- hcV81V280posord[[numberBN]]
test6 <- hcV280posord[[numberBN]]
test7<- hcV81posord[[numberBN]]


test8 <- test7
test8$with <- (test7$with +test6$with - test6$without)
# test8$with <- (test7$with +test6$with)/2

props <- list(test1,test2,test3,test4,test5,test6,test7,test8)
logliksplot <- c(rep(logliksCM[numberCM],4),rep(logliksIT[[numberBN+1]],4))
nedgesplot <- c(rep(nedgesnetworksCM[numberCM],4) ,rep(nedgesIT[numberBN+1],4))
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center9 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}
all <- cbind(plotdifferences)
all
do.call("grid.arrange", c(all, nrow = 2))


  prop_dif_plus <- quantity2clim(test5$with - test8$with, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)), rev.colors = TRUE)
  all.equal(test5$with-test5$without,test8$with-test8$without)
  all.equal(test1$with-test1$without,test4$with-test4$without)  

##############################################################################
# 227 568
##############################################################################  
  numberCM <- length(corcutsolo)-1
  numberCM <- 35
  numberBN <- 8
  
  
  test1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227,568), valueEvidence = c(2,2))
  test2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  test3 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(568), valueEvidence = c(2))
  
  test4 <- test3
  test4$with <- test3$with + test2$with - test2$without
  # test4$with <- (test3$with + test2$with)/2
  all.equal(test4,test3)
  
  test5 <- hcV227V568posord[[numberBN]]
  test6 <- hcV227posord[[numberBN]]
  test7 <- hcV568posord[[numberBN]]
  
  test8 <- test7
  test8$with <- (test7$with +test6$with - test6$without)
  # test8$with <- (test7$with +test6$with)/2
  
  props <- list(test1,test2,test3,test4,test5,test6,test7,test8)
  logliksplot <- c(rep(logliksCM[numberCM],4),rep(logliksIT[[numberBN+1]],4))
  nedgesplot <- c(rep(nedgesnetworksCM[numberCM],4) ,rep(nedgesIT[numberBN+1],4))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center9 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  all <- cbind(plotdifferences)
  all
  do.call("grid.arrange", c(all, nrow = 2))
  
  
  ##############################################################################
  # 205 227
  ##############################################################################  
  numberCM <- length(corcutsolo)-2
  numberBN <- 7
  
  
  test1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(205,227), valueEvidence = c(2,2))
  test2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(205), valueEvidence = c(2))
  test3 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  
  test4 <- test3
  test4$with <- test3$with + test2$with - test2$without
  # test4$with <- (test3$with + test2$with)/2
  all.equal(test4,test3)
  
  test5 <- hcV205V227posord[[numberBN]]
  test6 <- hcV205posord[[numberBN]]
  test7 <- hcV227posord[[numberBN]]
  
  test8 <- test7
  test8$with <- (test7$with +test6$with - test6$without)
  # test8$with <- (test7$with +test6$with)/2
  
  props <- list(test1,test2,test3,test4,test5,test6,test7,test8)
  logliksplot <- c(rep(logliksCM[numberCM],4),rep(logliksIT[[numberBN+1]],4))
  nedgesplot <- c(rep(nedgesnetworksCM[numberCM],4) ,rep(nedgesIT[numberBN+1],4))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center9 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
    
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  all <- cbind(plotdifferences)
  all
  do.call("grid.arrange", c(all, nrow = 2))
  do.call("grid.arrange", c(all[c(2,3,1,6,7,5)], nrow = 2))
  
  prop_dif_plus <- quantity2clim(test5$with - test8$with, paste0(attr(test1$with, "probability"),"-", attr(test1$without, "probability")), tas_ncep_10d)
  spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)), rev.colors = TRUE)
  all.equal(test5$with-test5$without,test8$with-test8$without)
  all.equal(test1$with-test1$without,test4$with-test4$without)  
  
###############################################################################
# V81V227
###############################################################################
  numberCM <- length(corcutsolo)-2
  numberCM <- 3
  numberBN <- 8
  
  
  test1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  test2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  test3 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  
  test4 <- test3
  test4$with <- test3$with + test2$with - test2$without
  # test4$with <- (test3$with + test2$with)/2
  all.equal(test4,test3)
  
  test5 <- hcV81V227posord[[numberBN]]
  test6 <- hcV81posord[[numberBN]]
  test7 <- hcV227posord[[numberBN]]
  
  test8 <- test7
  test8$with <- (test7$with +test6$with - test6$without)
  # test8$with <- (test7$with +test6$with)/2
  
  props <- list(test1,test2,test3,test4,test5,test6,test7,test8)
  logliksplot <- c(rep(logliksCM[numberCM],4),rep(logliksIT[[numberBN+1]],4))
  nedgesplot <- c(rep(nedgesnetworksCM[numberCM],4) ,rep(nedgesIT[numberBN+1],4))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center9 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
    
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  all <- cbind(plotdifferences)
  all
  do.call("grid.arrange", c(all, nrow = 2))
  do.call("grid.arrange", c(all[c(2,3,1,6,7,5)], nrow = 2))
  
  prop_dif_plus <- quantity2clim(test5$with - test8$with, paste0(attr(test1$with, "probability"),"-", attr(test1$without, "probability")), tas_ncep_10d)
  spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)), rev.colors = TRUE)
  all.equal(test5$with-test5$without,test8$with-test8$without)
  all.equal(test1$with-test1$without,test4$with-test4$without)  
  
  ###############################################################################
  # V81V227
  ###############################################################################
  numberCM <- length(corcutsolo)-2
  numberBN <- 7
  
  
  test1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  test2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  test3 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  
  test4 <- test3
  test4$with <- test3$with + test2$with - test2$without
  # test4$with <- (test3$with + test2$with)/2
  all.equal(test4,test3)
  
  test5 <- hcV81V227posord[[numberBN]]
  test6 <- hcV81posord[[numberBN]]
  test7 <- hcV227posord[[numberBN]]
  
  test8 <- test7
  test8$with <- (test7$with +test6$with - test6$without)
  # test8$with <- (test7$with +test6$with)/2
  
  props <- list(test1,test2,test3,test4,test5,test6,test7,test8)
  logliksplot <- c(rep(logliksCM[numberCM],4),rep(logliksIT[[numberBN+1]],4))
  nedgesplot <- c(rep(nedgesnetworksCM[numberCM],4) ,rep(nedgesIT[numberBN+1],4))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center9 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(80)
    
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  all <- cbind(plotdifferences)
  all
  do.call("grid.arrange", c(all, nrow = 2))
  do.call("grid.arrange", c(all[c(2,3,1,6,7,5)], nrow = 2))
  
  prop_dif_plus <- quantity2clim(test5$with - test8$with, paste0(attr(test1$with, "probability"),"-", attr(test1$without, "probability")), tas_ncep_10d)
  spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)), rev.colors = TRUE)
  all.equal(test5$with-test5$without,test8$with-test8$without)
  all.equal(test1$with-test1$without,test4$with-test4$without)  
  
  
  
##############################################################################
# Compare from CM and BN 1800 V81V280
##############################################################################
  numberCM <- length(corcutsolo)-1
  numberCM <- 
  numberBN <- 9
  
  posdoubleCN14 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  posdoubleBN <- hcproplistord[[numberBN]]
  negdoubleCN14 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  negdoubleBN <- hcpropneglistord[[numberBN]]
  possingleCN14 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleBN <- hcproppos1listord[[numberBN]]
  negsingleCN14 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleBN <- hcpropneg1listord[[numberBN]]
  
  props <- list(possingleBN, possingleCN14,negsingleBN,negsingleCN14,posdoubleBN,posdoubleCN14,negdoubleBN,negdoubleCN14)
  logCM <- logliksCM[numberCM]
  logBN <- logliksIT[[numberBN]]
  logliksplot <- c(rep(c(logBN,logCM),4))
  nCM <- nedgesnetworksCM[numberCM]
  nBM <- nedgesIT[numberBN]
  nedgesplot <- c(rep(c(nBM,nCM),4))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center7 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.03,0.85,0.01), set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    a
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  
  all <- cbind(plotdifferences)
  all
  
  plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper2/figures/probabilitiesComp1800.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all, ncol=2))
  dev.off()

  
  
  ##############################################################################
  # Compare from CM and BN 1800 V81V227
  ##############################################################################
  numberCM <- length(corcutsolo)-1
  numberBNedge <- 9
  numberBN <- 8
  
  
  possingleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleBN1 <- hcV81posord[[numberBN]]
  negsingleCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleBN1 <- hcpropnegV81listord[[numberBN]]
  possingleCN2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  possingleBN2 <- hcV227posord[[numberBN]]
  negsingleCN2 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(227), valueEvidence = c(2))
  negsingleBN2 <- hcpropnegV227listord[[numberBN]]
  posdoubleCN12 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  posdoubleBN12 <- hcV81V227posord[[numberBN]]
  negdoubleCN12 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  negdoubleBN12 <- hcpropV81V227neglistord[[numberBN]]
  
  attr(possingleCN1, "sign") <- attr(possingleCN2, "sign") <- attr(possingleBN1, "sign") <-  attr(possingleBN2, "sign") <-  attr(posdoubleBN12, "sign") <- attr(posdoubleCN12, "sign") <- "pos"
  attr(negsingleCN1, "sign") <- attr(negsingleCN2, "sign") <- attr(negsingleBN1, "sign") <-  attr(negsingleBN2, "sign") <-  attr(negdoubleBN12, "sign") <- attr(negdoubleCN12, "sign") <- "neg"
  
  props <- list(possingleBN1, possingleCN1,negsingleBN1,negsingleCN1,
                possingleBN2, possingleCN2,negsingleBN2,negsingleCN2,
                posdoubleBN12,posdoubleCN12,negdoubleBN12,negdoubleCN12)
  
  
  logCM <- logliksCM[numberCM]
  logBN <- logliksIT[[numberBNedge]]
  logliksplot <- c(rep(c(logBN,logCM),6))
  nCM <- nedgesnetworksCM[numberCM]
  nBM <- nedgesIT[numberBNedge]
  nedgesplot <- c(rep(c(nBM,nCM),6))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center7 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
    col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
    if(attr(prop, "sign") == "pos"){cb <- col.r}
    if(attr(prop, "sign") == "neg"){cb <- col.b}
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    a
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  
  all <- cbind(plotdifferences)
  all
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V227BNCN1800.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all, ncol=2))
  dev.off()

  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V227BNCN1800warm.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(1,2,5,6,9,10)], ncol=2))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V227BNCN1800cold.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(3,4,7,8,11,12)], ncol=2))
  dev.off()
  
  
  ##############################################################################
  # Compare from CM and BN 1800 V81V280
  ##############################################################################
  numberCM <- length(corcutsolo)-1
  numberBNedge <- 9
  numberBN <- 8
  
  
  possingleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleBN1 <- hcV81posord[[numberBN]]
  negsingleCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleBN1 <- hcpropnegV81listord[[numberBN]]
  possingleCN2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(280), valueEvidence = c(2))
  possingleBN2 <- hcV280posord[[numberBN]]
  negsingleCN2 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(280), valueEvidence = c(2))
  negsingleBN2 <- hcpropnegV280listord[[numberBN]]
  posdoubleCN12 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  posdoubleBN12 <- hcV81V280posord[[numberBN]]
  negdoubleCN12 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  negdoubleBN12 <- hcpropV81V280neglistord[[numberBN]]
  
  attr(possingleCN1, "sign") <- attr(possingleCN2, "sign") <- attr(possingleBN1, "sign") <-  attr(possingleBN2, "sign") <-  attr(posdoubleBN12, "sign") <- attr(posdoubleCN12, "sign") <- "pos"
  attr(negsingleCN1, "sign") <- attr(negsingleCN2, "sign") <- attr(negsingleBN1, "sign") <-  attr(negsingleBN2, "sign") <-  attr(negdoubleBN12, "sign") <- attr(negdoubleCN12, "sign") <- "neg"
  
  props <- list(possingleBN1, possingleCN1,negsingleBN1,negsingleCN1,
                possingleBN2, possingleCN2,negsingleBN2,negsingleCN2,
                posdoubleBN12,posdoubleCN12,negdoubleBN12,negdoubleCN12)
  
  
  logCM <- logliksCM[numberCM]
  logBN <- logliksIT[[numberBNedge]]
  logliksplot <- c(rep(c(logBN,logCM),6))
  nCM <- nedgesnetworksCM[numberCM]
  nBM <- nedgesIT[numberBNedge]
  nedgesplot <- c(rep(c(nBM,nCM),6))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center7 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
    col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
    if(attr(prop, "sign") == "pos"){cb <- col.r}
    if(attr(prop, "sign") == "neg"){cb <- col.b}
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    a
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  
  all <- cbind(plotdifferences)
  all
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V280BNCN1800.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all, ncol=2))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V280BNCN1800warm.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(1,2,5,6,9,10)], ncol=2))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V280BNCN1800cold.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(3,4,7,8,11,12)], ncol=2))
  dev.off()
  
  ##############################################################################
  # Compare from CM and BN 1800 V81V280 CN 8000 170000
  ##############################################################################
  numberCM <- length(corcutsolo)-1
  numberCM2 <- 8
  numberBNedge <- 9
  numberBN <- 8
  
  
  possingleCN1a <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleCN1b <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleBN1 <- hcV81posord[[numberBN]]
  
  negsingleCN1a <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleCN1b <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleBN1 <- hcpropnegV81listord[[numberBN]]
  
  possingleCN2a <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(280), valueEvidence = c(2))
  possingleCN2b <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(280), valueEvidence = c(2))
  possingleBN2 <- hcV280posord[[numberBN]]
  
  negsingleCN2a <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(280), valueEvidence = c(2))
  negsingleCN2b <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(280), valueEvidence = c(2))
  negsingleBN2 <- hcpropnegV280listord[[numberBN]]
  
  posdoubleCN12a <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  posdoubleCN12b <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  posdoubleBN12 <- hcV81V280posord[[numberBN]]
  
  negdoubleCN12a <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  negdoubleCN12b <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
  negdoubleBN12 <- hcpropV81V280neglistord[[numberBN]]
  
  attr(possingleCN1a, "sign") <- attr(possingleCN1b, "sign") <- attr(possingleCN2a, "sign") <- attr(possingleCN2b, "sign") <- attr(possingleBN1, "sign") <-  attr(possingleBN2, "sign") <-  attr(posdoubleBN12, "sign") <- attr(posdoubleCN12a, "sign")  <- attr(posdoubleCN12b, "sign") <- "pos"
  attr(negsingleCN1a, "sign") <-  attr(negsingleCN1b, "sign") <-attr(negsingleCN2a, "sign") <-  attr(negsingleCN2b, "sign") <- attr(negsingleBN1, "sign") <-  attr(negsingleBN2, "sign") <-  attr(negdoubleBN12, "sign") <- attr(negdoubleCN12a, "sign") <- attr(negdoubleCN12b, "sign") <- "neg"
  
  props <- list(possingleBN1, possingleCN1a, possingleCN1b, negsingleBN1, negsingleCN1a, negsingleCN1b,
                possingleBN2, possingleCN2a, possingleCN2b, negsingleBN2, negsingleCN2a, negsingleCN2b,
                posdoubleBN12,posdoubleCN12a, posdoubleCN12b, negdoubleBN12, negdoubleCN12a, negdoubleCN12b)
  
  
  logCM <- logliksCM[numberCM]
  logCM2 <- logliksCM[numberCM2]
  logBN <- logliksIT[[numberBNedge]]
  logliksplot <- c(rep(c(logBN,logCM,logCM2),6))
  nCM <- nedgesnetworksCM[numberCM]
  nCM2 <- nedgesnetworksCM[numberCM2]
  nBM <- nedgesIT[numberBNedge]
  nedgesplot <- c(rep(c(nBM,nCM,nCM2),6))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center7 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
    col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
    if(attr(prop, "sign") == "pos"){cb <- col.r}
    if(attr(prop, "sign") == "neg"){cb <- col.b}
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    a
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  
  all <- cbind(plotdifferences)
  all
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V280BNCN1CN21800.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all, ncol=3))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V280BNCN1CN21800warm.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(1,2,5,6,9,10)], ncol=2))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V280BNCN1CN21800cold.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(3,4,7,8,11,12)], ncol=2))
  dev.off()

##############################################################################
# Compare from CM and BN 1800 V81V227 CN 2200 170000
##############################################################################
  numberBNedge <- 9
  numberBN <- 8
  numberCM <- 56
  numberCM2 <- 101
  
  
  possingleCN1a <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleCN1b <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  possingleBN1 <- hcV81posord[[numberBN]]
  
  negsingleCN1a <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleCN1b <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
  negsingleBN1 <- hcpropnegV81listord[[numberBN]]
  
  possingleCN2a <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  possingleCN2b <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
  possingleBN2 <- hcV227posord[[numberBN]]
  
  negsingleCN2a <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(227), valueEvidence = c(2))
  negsingleCN2b <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(227), valueEvidence = c(2))
  negsingleBN2 <- hcpropnegV227listord[[numberBN]]
  
  posdoubleCN12a <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  posdoubleCN12b <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  posdoubleBN12 <- hcV81V227posord[[numberBN]]
  
  negdoubleCN12a <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  negdoubleCN12b <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  negdoubleBN12 <- hcpropV81V227neglistord[[numberBN]]
  
  attr(possingleCN1a, "sign") <- attr(possingleCN1b, "sign") <- attr(possingleCN2a, "sign") <- attr(possingleCN2b, "sign") <- attr(possingleBN1, "sign") <-  attr(possingleBN2, "sign") <-  attr(posdoubleBN12, "sign") <- attr(posdoubleCN12a, "sign")  <- attr(posdoubleCN12b, "sign") <- "pos"
  attr(negsingleCN1a, "sign") <-  attr(negsingleCN1b, "sign") <-attr(negsingleCN2a, "sign") <-  attr(negsingleCN2b, "sign") <- attr(negsingleBN1, "sign") <-  attr(negsingleBN2, "sign") <-  attr(negdoubleBN12, "sign") <- attr(negdoubleCN12a, "sign") <- attr(negdoubleCN12b, "sign") <- "neg"
  
  props <- list(possingleBN1, possingleCN1a, possingleCN1b, negsingleBN1, negsingleCN1a, negsingleCN1b,
                possingleBN2, possingleCN2a, possingleCN2b, negsingleBN2, negsingleCN2a, negsingleCN2b,
                posdoubleBN12,posdoubleCN12a, posdoubleCN12b, negdoubleBN12, negdoubleCN12a, negdoubleCN12b)
  
  
  logCM <- logliksCM[numberCM]
  logCM2 <- logliksCM[numberCM2]
  logBN <- logliksIT[[numberBNedge]]
  logliksplot <- c(rep(c(logBN,logCM,logCM2),6))
  nCM <- nedgesnetworksCM[numberCM]
  nCM2 <- nedgesnetworksCM[numberCM2]
  nBM <- nedgesIT[numberBNedge]
  nedgesplot <- c(rep(c(nBM,nCM,nCM2),6))
  nedgesplot
  logliksplot
  
  plotwiths <- list()
  plotwithouts <- list()
  plotdifferences <- list()
  Center7 <- 180
  
  for(i in 1:length(props)){
    props[[1]]
    prop <- props[[i]]
    
    col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
    col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
    if(attr(prop, "sign") == "pos"){cb <- col.r}
    if(attr(prop, "sign") == "neg"){cb <- col.b}
    
    prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
    prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
    prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
    
    str(prop_with)
    
    col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    a
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  
  all <- cbind(plotdifferences)
  all
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V227BNCN1CN21800.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all, ncol=3))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V227BNCN1CN21800warm.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(1,2,3,7,8,9,13,14,15)], ncol=3))
  dev.off()
  
  plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/V81V227BNCN1CN21800cold.pdf")
  pdf(plotname, height = 7, width = 5)
  do.call("grid.arrange", c(all[c(4,5,6,10,11,12,16,17,18)], ncol=3))
  dev.off()

##############################################################################
#Evidences points chosen from comunities
##############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V85pos/prop_hc_1600_1800i_V85_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V205pos/prop_hc_1600_1800i_V205_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227pos/prop_hc_1600_1800i_V227_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V424pos/prop_hc_1600_1800i_V424_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V459pos/prop_hc_1600_1800i_V459_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V459pos/prop_hc_1600_1800i_V459_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V532pos/prop_hc_1600_1800i_V532_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V568pos/prop_hc_1600_1800i_V568_equal2.rda")

  
numberCM <- length(corcutsolo)-1
numberBNedge <- 9
numberBN <- 8

possingleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(85), valueEvidence = c(2))
possingleBN1 <- prop_hc_1600_1800i_V85_equal2
possingleCN2 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(205), valueEvidence = c(2))
possingleBN2 <- prop_hc_1600_1800i_V205_equal2
possingleCN3 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
possingleBN3 <- prop_hc_1600_1800i_V227_equal2
possingleCN4 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(424), valueEvidence = c(2))
possingleBN4 <- prop_hc_1600_1800i_V424_equal2
possingleCN5 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(459), valueEvidence = c(2))
possingleBN5 <- prop_hc_1600_1800i_V459_equal2
possingleCN6 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(532), valueEvidence = c(2))
possingleBN6 <- prop_hc_1600_1800i_V532_equal2
possingleCN7 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(568), valueEvidence = c(2))
possingleBN7 <- prop_hc_1600_1800i_V568_equal2

attr(possingleCN1, "sign") <- attr(possingleCN2, "sign") <- attr(possingleBN1, "sign") <-  attr(possingleBN2, "sign") <-  attr(possingleBN3, "sign") <- attr(possingleCN3, "sign") <- "pos"
attr(possingleCN4, "sign") <- attr(possingleCN5, "sign") <- attr(possingleBN4, "sign") <-  attr(possingleBN5, "sign") <-  attr(possingleBN6, "sign") <- attr(possingleCN6, "sign") <- "pos"
attr(possingleBN7, "sign") <- attr(possingleCN7, "sign") <- "pos"

props <- list(possingleBN1, possingleCN1,possingleBN2,possingleCN2,
              possingleBN3, possingleCN3,possingleBN4,possingleCN4,
              possingleBN5, possingleCN5,possingleBN6,possingleCN6,
              possingleBN7, possingleCN7)

logCM <- logliksCM[numberCM]
logBN <- logliksIT[[numberBNedge]]
logliksplot <- c(rep(c(logBN,logCM),7))
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBNedge]
nedgesplot <- c(rep(c(nBM,nCM),7))
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180

for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
do.call("grid.arrange", c(all, ncol=2))
do.call("grid.arrange", c(all[c(13,14)], ncol=2))
do.call("grid.arrange", c(all[c(3,4,7,8,9,10)], ncol=2))
do.call("grid.arrange", c(all[c(1,2,5,6,11,12)], ncol=2))



##############################################################################
# Evidence Warm Pool Cold Tong Niño.
##############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos/prop_hc_1600_1800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos/prop_hc_2600_2800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171CT/prop_hc_1600_1800i_V81V171_equal12.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171CT/prop_hc_2600_2800i_V81V171_equal12.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171WP/prop_hc_8600_8800i_V81V171_equal20.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171WP/prop_hc_1600_1800i_V81V171_equal20.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V171WP/prop_hc_2600_2800i_V81V171_equal20.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V171pos/prop_hc_1600_1800i_V171_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V171neg/propneg_hc_1600_1800i_V171_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28pos/prop_hc_1600_1800i_V28_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V28V171WP2pos/prop_hc_1600_1800i_V28V171_equal20.rda")


x <- corcutsolo[[numberCM]]
sigma <- x
chkcorr(x)
chkcorr <- function(x) {
  if (!is.matrix(x) || (d <- dim(x))[1] != d[2])
    return(FALSE)
  rownames(x) <- colnames(x) <- NULL
  storage.mode(x) <- "numeric"
  ONE <- 1 + sqrt(.Machine$double.eps)
  
  ## return
  -ONE <= min(x) && max(x) <= ONE && isTRUE(all.equal(diag(x), rep(1, d[1])))
}

if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) < 0))
  stop(sQuote("sigma"), " is not a covariance matrix")

nedgesnetworksCM
numberCM <- 61
numberCM2 <- 56
numberCM3 <- 54
numberCM4 <- 35
numberCM5 <- 1
nedgesnetworksCM[c(numberCM,numberCM2,numberCM3,numberCM4,numberCM5)]

numberBN <- 8
numberBNedge <- 9
numberBN2 <- 12
numberBNedge2 <- 13
numberBN3 <- 43
numberBNedge3 <- 44
nedgesIT[c(numberBNedge,numberBNedge2,numberBNedge3)]


posmixCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN4 <- propagationCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN5 <- propagationCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))

posmixBN1 <- prop_hc_1600_1800i_V81_equal2
posmixBN2 <- prop_hc_2600_2800i_V81_equal2

posCTCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(171), valueEvidence = c(2))
posCTCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(171), valueEvidence = c(2))
posCTCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(171), valueEvidence = c(2))
posCTCN4 <- propagationCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(171), valueEvidence = c(2))
posCTCN5 <- propagationCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(171), valueEvidence = c(2))

posCTBN1 <- prop_hc_1600_1800i_V171_equal2

posWPCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28,171), valueEvidence = c(2,0))
posWPCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28,171), valueEvidence = c(2,0))
posWPCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28,171), valueEvidence = c(2,0))
posWPCN4 <- propagationCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28,171), valueEvidence = c(2,0))
posWPCN5 <- propagationCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28,171), valueEvidence = c(2,0))

posWPBN1 <- prop_hc_1600_1800i_V81V171_equal20
posWPBN2 <- prop_hc_2600_2800i_V81V171_equal20
posWPBN3 <- prop_hc_8600_8800i_V81V171_equal20

posmixCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28), valueEvidence = c(2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28), valueEvidence = c(2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28), valueEvidence = c(2))
posmixCN4 <- propagationCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28), valueEvidence = c(2))
posmixCN5 <- propagationCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(28), valueEvidence = c(2))

posmix2BN1 <- prop_hc_1600_1800i_V28_equal2
posWP2BN1 <- prop_hc_1600_1800i_V28V171_equal20

attr(posmix2BN1, "sign") <- "pos"
attr(posWP2BN1, "sign") <- "pos"
attr(posmixCN1, "sign") <- attr(posmixCN2, "sign") <- attr(posmixCN3, "sign") <- "pos"
attr(posmixCN4, "sign") <- attr(posmixCN5, "sign") <- "pos"
attr(posmixBN1, "sign") <- attr(posmixBN2, "sign") <- "pos"
attr(posCTCN1, "sign") <- attr(posCTCN2, "sign") <- attr(posCTCN3, "sign") <- "pos"
attr(posCTCN4, "sign") <- attr(posCTCN5, "sign") <- "pos"
attr(posCTBN1, "sign")  <- "pos"
attr(posWPBN1, "sign") <- attr(posWPBN2, "sign") <- attr(posWPBN3, "sign")<- "pos"
attr(posWPCN1, "sign") <- attr(posWPCN2, "sign") <- attr(posWPCN3, "sign") <- "pos"
attr(posWPCN4, "sign") <- attr(posWPCN5, "sign") <- "pos"
  # <- attr(possingleCN3, "sign") <- "pos"
# attr(possingleCN4, "sign") <- attr(possingleCN5, "sign") <- attr(possingleBN4, "sign") <-  attr(possingleBN5, "sign") <-  attr(possingleBN6, "sign") <- attr(possingleCN6, "sign") <- "pos"
# attr(possingleBN7, "sign") <- attr(possingleCN7, "sign") <- "pos"

props <- list(posmix2BN1)
props <- list(posWP2BN1)
props <- list(posmixCN1, posmixCN2, posmixCN3) 
props <- list(posmixCN4, posmixCN5)
props <- list(posmixCN1, posmixCN2, posmixCN3,posmixCN4, posmixCN5)
props <- list(posmixBN1, posmixBN2)
props <- list(posCTCN1, posCTCN2,posCTCN3,posCTCN4,posCTCN5)
props <- list(posCTBN1)
props <- list(posWPBN1, posWPBN2,posWPBN3)
props <- list(posWPCN4,posWPCN5)
props <- list(posWPCN1,posWPCN2,posWPCN3,posWPCN4,posWPCN5)
props <- list(posmixCN3,posCTCN3,posWPCN3)
props <- list(posmixCN2,posCTCN2)
props <- list(posmixBN1,posCTBN1,posWPBN1)
props <- list(posmixBN2,posWPBN2)

lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3,numberCM4,numberCM5)]
logBN <- c(logliksIT[c(numberBNedge,numberBNedge2,numberBNedge3)])
logliksplot <- logBN
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3,numberCM4,numberCM5)]
nBM <- nedgesIT[c(numberBNedge,numberBNedge2,numberBNedge3)]
nedgesplot <- nBM
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
do.call("grid.arrange", c(all, ncol=1))

##############################################################################
# Evidence Warm Pool Cold Tong Niño. NEgative
##############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V171neg/propneg_hc_1600_1800i_V171_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V171neg/propneg_hc_2600_2800i_V171_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81neg/propneg_hc_1600_1800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81neg/propneg_hc_2600_2800i_V81_equal2.rda")



nedgesnetworksCM
numberCM <- 61
numberCM2 <- 56
numberCM3 <- 54
numberCM4 <- 35
numberCM5 <- 1
nedgesnetworksCM[c(numberCM,numberCM2,numberCM3,numberCM4,numberCM5)]

numberBN <- 8
numberBNedge <- 9
numberBN2 <- 12
numberBNedge2 <- 13
numberBN3 <- 43
numberBNedge3 <- 44
nedgesIT[c(numberBNedge,numberBNedge2,numberBNedge3)]


negmixCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN4 <- propagationnegCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN5 <- propagationnegCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))

negmixBN1 <- propneg_hc_1600_1800i_V81_equal2
negmixBN2 <- propneg_hc_2600_2800i_V81_equal2

negCTCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(171), valueEvidence = c(2))
negCTCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(171), valueEvidence = c(2))
negCTCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(171), valueEvidence = c(2))
negCTCN4 <- propagationnegCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(171), valueEvidence = c(2))
negCTCN5 <- propagationnegCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(171), valueEvidence = c(2))

negCTBN1 <- propneg_hc_1600_1800i_V171_equal2
negCTBN2 <- propneg_hc_2600_2800i_V171_equal2

negWPCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,171), valueEvidence = c(2,0))
negWPCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,171), valueEvidence = c(2,0))
negWPCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,171), valueEvidence = c(2,0))
negWPCN4 <- propagationnegCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,171), valueEvidence = c(2,0))
negWPCN5 <- propagationnegCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,171), valueEvidence = c(2,0))

negWPBN1 <- propneg_hc_1600_1800i_V81V171_equal20
negWPBN2 <- propneg_hc_2600_2800i_V81V171_equal20

negmix2CN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(28), valueEvidence = c(2))
negmix2CN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(28), valueEvidence = c(2))
negmix2CN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(28), valueEvidence = c(2))
negmix2CN4 <- propagationCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(28), valueEvidence = c(2))
negmix2CN5 <- propagationCorr(corcutsolo[[numberCM5]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(28), valueEvidence = c(2))

negmix2BN1 <- propneg_hc_1600_1800i_V28_equal2

attr(negmixBN1, "sign") <- attr(negmixBN2, "sign") <- "neg"
attr(negmixCN1, "sign") <- attr(negmixCN2, "sign") <- attr(negmixCN3, "sign") <- "neg"
attr(negmixCN4, "sign") <- attr(negmixCN5, "sign") <- "neg"

attr(negCTCN1, "sign") <- attr(negCTCN2, "sign") <- attr(negCTCN3, "sign") <- "neg"
attr(negCTCN4, "sign") <- attr(negCTCN5, "sign") <- "neg"
attr(negCTBN1, "sign")  <- attr(negCTBN2, "sign") <- "neg"

attr(negWPCN1, "sign") <- attr(negWPCN2, "sign") <- attr(negWPCN3, "sign") <- "neg"
attr(negWPCN4, "sign") <- attr(negWPCN5, "sign") <- "neg"
attr(negWPBN1, "sign") <- attr(negWPBN2, "sign") <- "neg"

attr(negmix2BN1, "sign") <- "neg"

# <- attr(possingleCN3, "sign") <- "pos"
# attr(possingleCN4, "sign") <- attr(possingleCN5, "sign") <- attr(possingleBN4, "sign") <-  attr(possingleBN5, "sign") <-  attr(possingleBN6, "sign") <- attr(possingleCN6, "sign") <- "pos"
# attr(possingleBN7, "sign") <- attr(possingleCN7, "sign") <- "pos"

props <- list(negmixCN1, negmixCN2, negmixCN3,negmixCN4, negmixCN5)
props <- list(negmixBN1, negmixBN2)
props <- list(negCTCN1, negCTCN2,negCTCN3,negCTCN4,negCTCN5)
props <- list(negCTBN1,negCTBN2)
props <- list(negWPCN1,negWPCN2,negWPCN3,negWPCN4,negWPCN5)
props <- list(negWPCN4,negWPCN5)
props <- list(negWPBN1, negWPBN2)
props <- list(negmix2BN1)



lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3,numberCM4,numberCM5)]
logBN <- c(logliksIT[c(numberBNedge,numberBNedge2,numberBNedge3)])
logliksplot <- logBN
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3,numberCM4,numberCM5)]
nBM <- nedgesIT[c(numberBNedge,numberBNedge2,numberBNedge3)]
nedgesplot <- nBM
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
do.call("grid.arrange", c(all, ncol=1))
#######################################################################################
# Full BN Network 7862 links. V280 V81 
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V81V280_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V280_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_8600_8800i_V81_equal2.rda")
numberBNedge <- 44
numberBN <- 43
numberCM <- 1


possingleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
possingleBN1 <- prop_hc_8600_8800i_V81_equal2

possingleCN2<- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(280), valueEvidence = c(2))
possingleBN2 <- prop_hc_8600_8800i_V280_equal2

posdoubleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,280), valueEvidence = c(2,2))
posdoubleBN1 <- prop_hc_8600_8800i_V81V280_equal22


attr(possingleCN1, "sign") <- "pos"
attr(possingleBN1, "sign") <- "pos"
attr(possingleCN2, "sign") <- "pos"
attr(possingleBN2, "sign") <- "pos"
attr(posdoubleCN1, "sign") <- "pos"
attr(posdoubleBN1, "sign") <- "pos"

props <- list(possingleBN1,possingleCN1, possingleBN2,possingleCN2 ,posdoubleBN1,posdoubleCN1) 

logCM <- logliksCM[numberCM]
logBN <- logliksIT[[numberBNedge]]
logliksplot <- rep(c(logBN,logCM),3)
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBNedge]
nedgesplot <- rep(c(nBM,nCM),3)
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  prop <- props[[i]]
  attr(possingleCN1, "sign")
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/overestimate.pdf"
pdf(plotname)
do.call("grid.arrange", c(all, ncol=2))
dev.off()

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/FULLBNFULLCNV81.pdf"
pdf(plotname)
do.call("grid.arrange", c(all[1:2], ncol=2))
dev.off()

#######################################################################################
# Full BN Network 7862 links. V227 V81 
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V227_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_8600_8800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_8600_8800i_V227_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_8600_8800i_V81V227_equal22.rda")

numberBNedge <- 44
numberBN <- 43
numberCM <- 1

possingleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
possingleBN1 <- prop_hc_8600_8800i_V81_equal2
negsingleCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negsingleBN1 <- propneg_hc_8600_8800i_V81_equal2

possingleCN2<- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(227), valueEvidence = c(2))
possingleBN2 <- prop_hc_8600_8800i_V227_equal2
negsingleCN2<- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(227), valueEvidence = c(2))
negsingleBN2 <- propneg_hc_8600_8800i_V227_equal2

posdoubleCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
posdoubleBN1 <- prop_hc_8600_8800i_V81V227_equal22
negdoubleCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
negdoubleBN1 <- propneg_hc_8600_8800i_V81V227_equal22

attr(possingleCN1, "sign") <- "pos"
attr(possingleBN1, "sign") <- "pos"
attr(possingleCN2, "sign") <- "pos"
attr(possingleBN2, "sign") <- "pos"
attr(posdoubleCN1, "sign") <- "pos"
attr(posdoubleBN1, "sign") <- "pos"
attr(negsingleCN1, "sign") <- "neg"
attr(negsingleBN1, "sign") <- "neg"
attr(negsingleCN2, "sign") <- "neg"
attr(negsingleBN2, "sign") <- "neg"
attr(negdoubleCN1, "sign") <- "neg"
attr(negdoubleBN1, "sign") <- "neg"

propspos <- list(possingleBN1,possingleCN1, possingleBN2,possingleCN2 ,posdoubleBN1,posdoubleCN1) 
propsneg <- list(negsingleBN1,negsingleCN1, negsingleBN2,negsingleCN2 ,negdoubleBN1,negdoubleCN1) 
props <- propspos

logCM <- logliksCM[numberCM]
logBN <- logliksIT[[numberBNedge]]
logliksplot <- rep(c(logBN,logCM),3)
nCM <- nedgesnetworksCM[numberCM]
nBM <- nedgesIT[numberBNedge]
nedgesplot <- rep(c(nBM,nCM),3)
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  prop <- props[[i]]
  attr(possingleCN1, "sign")
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/FullV227V81all.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/FullV227V81allneg.pdf"
pdf(plotname)
do.call("grid.arrange", c(all, ncol=2))
dev.off()

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/FULLBNFULLCNV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/FULLBNFULLCNV81neg.pdf"
pdf(plotname)
do.call("grid.arrange", c(all[1:2], ncol=2))
dev.off()


##############################################################################
# V81 BN Full 1800 CN Full 2200 17000  
##############################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81pos/prop_hc_1600_1800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81neg/propneg_hc_1600_1800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/prop_hc_8600_8800i_V81_equal2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/lastBN/propneg_hc_8600_8800i_V81_equal2.rda")
nedgesnetworksCM
numberCM1 <- 56
numberCM2 <- 43
numberCM3 <- 35
numberCM4 <- 101
nedgesnetworksCM[c(numberCM1,numberCM2, numberCM3,numberCM4)]

numberBN <- 8
numberBNedge <- 9

nedgesIT[c(numberBNedge)]

posmixCN1 <- propagationCorr(corcutsolo[[numberCM1]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
posmixCN4 <- propagationCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN1 <- propagationnegCorr(corcutsolo[[numberCM1]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))
negmixCN4 <- propagationnegCorr(corcutsolo[[numberCM4]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(2))

posmixBN1 <- prop_hc_1600_1800i_V81_equal2
negmixBN1 <- propneg_hc_1600_1800i_V81_equal2

attr(posmixCN1, "sign") <- attr(posmixCN2, "sign") <- attr(posmixCN3, "sign") <- attr(posmixCN4, "sign") <- "pos"
attr(posmixBN1, "sign") <- "pos"
attr(negmixCN1, "sign") <- attr(negmixCN2, "sign") <- attr(negmixCN3, "sign") <- attr(negmixCN4, "sign") <- "neg"
attr(negmixBN1, "sign") <- "neg"

propspos <- list(posmixBN1, posmixCN1, posmixCN2, posmixCN3, posmixCN4) 
propsneg <- list(negmixBN1, negmixCN1, negmixCN2, negmixCN3, negmixCN4) 
props <- propsneg

lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3)]
logBN <- c(logliksIT[c(numberBNedge)])
logliksplot <- c(logBN,logCM)
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]
nBM <- nedgesIT[c(numberBNedge)]
nedgesplot <- c(nBM,nCM)
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
blank <- grid.rect(gp=gpar(col="white"))
plotwindow <- list(all[[1]],all[[2]],blank,all[[3]],blank,all[[4]],blank,all[[5]])
plotwindow <- list(all[[1]],all[[3]],blank,all[[4]])
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_V81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_V81neg.pdf"

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_V81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_V81neg.pdf"
pdf(plotname)

do.call("grid.arrange", c(plotwindow, ncol=2))
dev.off()


#######################################################################################
# BN 1800 V227 and V81V227  pos and neg
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V227pos/prop_hc_1600_1800i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V81V227neg/propneg_hc_1600_1800i_V81V227_equal22.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227pos/prop_hc_1600_1800i_V227_equal2.rda")
oad("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/V227neg/propneg_hc_1600_1800i_V227_equal2.rda")

numberBNedge <- 9
numberBN <- 8

possingleBN1 <- prop_hc_1600_1800i_V227_equal2
negsingleBN1 <- propneg_hc_1600_1800i_V227_equal2
posdoubleBN1 <- prop_hc_1600_1800i_V81V227_equal22
negdoubleBN1 <- propneg_hc_1600_1800i_V81V227_equal22

attr(possingleBN1, "sign") <- "pos"
attr(negsingleBN1, "sign") <- "neg"
attr(posdoubleBN1, "sign") <- "pos"
attr(negdoubleBN1, "sign") <- "neg"

props <- list(possingleBN1,negsingleBN1,posdoubleBN1,negdoubleBN1) 

logBN <- logliksIT[[numberBNedge]]
logliksplot <- rep(c(logBN),4)
nBM <- nedgesIT[numberBNedge]
nedgesplot <- rep(c(nBM),4)
nedgesplot
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  prop <- props[[i]]
  attr(possingleCN1, "sign")
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

all <- cbind(plotdifferences)
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropV81V227/BNV227.pdf"
pdf(plotname)
do.call("grid.arrange", c(all, ncol=2))
dev.off()

#######################################################################################
# CNBN negative evidence V81
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81pos/prop_hc_1600_1800i_V81_equalmin2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81neg/propneg_hc_1600_1800i_V81_equalmin2.rda")


nedgesnetworksCM
numberCM <- 56
numberCM2 <- 35
numberCM3 <- 101
nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]

numberBN <- 8
numberBNedge <- 9

nedgesIT[c(numberBNedge)]

posmixCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(-2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(-2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(-2))
negmixCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(-2))
negmixCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(-2))
negmixCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81), valueEvidence = c(-2))

posmixBN1 <- prop_hc_1600_1800i_V81_equalmin2
negmixBN1 <- propneg_hc_1600_1800i_V81_equalmin2

attr(posmixCN1, "sign") <- attr(posmixCN2, "sign") <- attr(posmixCN3, "sign") <- "pos"
attr(posmixBN1, "sign") <- "pos"
attr(negmixCN1, "sign") <-  attr(negmixCN2, "sign") <- attr(negmixCN3, "sign") <- "neg"
attr(negmixBN1, "sign") <- "neg"

propspos <- list(posmixBN1, posmixCN1, posmixCN2, posmixCN3) 
propspos <- list(posmixCN1, posmixCN2, posmixCN3) 
propsneg <- list(negmixBN1, negmixCN1, negmixCN2, negmixCN3) 
propsneg <- list(negmixCN1,negmixCN2, negmixCN3) 
props <- propspos

lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3)]
logBN <- c(logliksIT[c(numberBNedge)])
logliksplot <- c(logBN,logCM)
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]
nBM <- nedgesIT[c(numberBNedge)]
nedgesplot <- c(nBM,nCM)
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  b
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}


posmixCN1$with[81]

all <- cbind(plotdifferences)
blank <- grid.rect(gp=gpar(col="white"))
plotwindow <- list(all[[1]],all[[2]],blank,all[[3]],blank,all[[4]])
plotwindow <- list(all[[1]],all[[3]],blank,all[[4]])
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81neg.pdf"

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81neg.pdf"
pdf(plotname)
do.call("grid.arrange", c(plotwindow, ncol=2))
dev.off()


#######################################################################################
# CN  negative evidence V568 
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV568pos/prop_hc_1600_1800i_V568_equalmin2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV568neg/propneg_hc_1600_1800i_V568_equalmin2.rda")


nedgesnetworksCM
numberCM <- 56
numberCM2 <- 35
numberCM3 <- 101
nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]

numberBN <- 8
numberBNedge <- 9

nedgesIT[c(numberBNedge)]

posmixCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(568), valueEvidence = c(-2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(568), valueEvidence = c(-2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(568), valueEvidence = c(-2))
negmixCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(568), valueEvidence = c(-2))
negmixCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(568), valueEvidence = c(-2))
negmixCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(568), valueEvidence = c(-2))

posmixBN1 <- prop_hc_1600_1800i_V568_equalmin2
negmixBN1 <- propneg_hc_1600_1800i_V568_equalmin2

attr(posmixCN1, "sign") <- attr(posmixCN2, "sign") <- attr(posmixCN3, "sign") <- "pos"
attr(posmixBN1, "sign") <- "pos"
attr(negmixCN1, "sign") <-  attr(negmixCN2, "sign") <- attr(negmixCN3, "sign") <- "neg"
attr(negmixBN1, "sign") <- "neg"

propspos <- list(posmixBN1, posmixCN1, posmixCN2, posmixCN3) 
# propspos <- list(posmixCN1, posmixCN2, posmixCN3) 
propsneg <- list(negmixBN1, negmixCN1, negmixCN2, negmixCN3) 
# propsneg <- list(negmixCN1,negmixCN2, negmixCN3) 
props <- propsneg

lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3)]
logBN <- c(logliksIT[c(numberBNedge)])
logliksplot <- c(logBN,logCM)
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]
nBM <- nedgesIT[c(numberBNedge)]
nedgesplot <- c(nBM,nCM)
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  b
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}


all <- cbind(plotdifferences)
blank <- grid.rect(gp=gpar(col="white"))
plotwindow <- list(all[[1]],all[[2]],blank,all[[3]],blank,all[[4]])
plotwindow <- list(all[[1]],all[[3]],blank,all[[4]])
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81neg.pdf"

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81neg.pdf"
pdf(plotname)
do.call("grid.arrange", c(plotwindow, ncol=2))
dev.off()

#######################################################################################
# CN  negative evidence V532
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV532pos/prop_hc_1600_1800i_V532_equalmin2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV532neg/propneg_hc_1600_1800i_V532_equalmin2.rda")


nedgesnetworksCM
numberCM <- 56
numberCM2 <- 35
numberCM3 <- 101
nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]

numberBN <- 8
numberBNedge <- 9

nedgesIT[c(numberBNedge)]

posmixCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(532), valueEvidence = c(-2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(532), valueEvidence = c(-2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(532), valueEvidence = c(-2))
negmixCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(532), valueEvidence = c(-2))
negmixCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(532), valueEvidence = c(-2))
negmixCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(532), valueEvidence = c(-2))

posmixBN1 <- prop_hc_1600_1800i_V532_equalmin2
negmixBN1 <- propneg_hc_1600_1800i_V532_equalmin2
# negmixBN1 <- propneg_hc_1600_1800i_V227_equal2

attr(posmixCN1, "sign") <- attr(posmixCN2, "sign") <- attr(posmixCN3, "sign") <- "pos"
attr(posmixBN1, "sign") <- "pos"
attr(negmixCN1, "sign") <-  attr(negmixCN2, "sign") <- attr(negmixCN3, "sign") <- "neg"
attr(negmixBN1, "sign") <- "neg"

propspos <- list(posmixBN1, posmixCN1, posmixCN2, posmixCN3) 
# propspos <- list(posmixCN1, posmixCN2, posmixCN3) 
propsneg <- list(negmixBN1, negmixCN1, negmixCN2, negmixCN3) 
# propsneg <- list(negmixCN1,negmixCN2, negmixCN3) 
props <- propsneg

lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3)]
logBN <- c(logliksIT[c(numberBNedge)])
logliksplot <- c(logBN,logCM)
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]
nBM <- nedgesIT[c(numberBNedge)]
nedgesplot <- c(nBM,nCM)
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  b
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}


posmixCN1$with[81]

all <- cbind(plotdifferences)
blank <- grid.rect(gp=gpar(col="white"))
plotwindow <- list(all[[1]],all[[2]],blank,all[[3]],blank,all[[4]])
plotwindow <- list(all[[1]],all[[3]],blank,all[[4]])
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81neg.pdf"

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81neg.pdf"
pdf(plotname)
do.call("grid.arrange", c(plotwindow, ncol=2))

dev.off()


#######################################################################################
# CN  negative evidence V81V227
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV532pos/prop_hc_1600_1800i_V532_equalmin2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV532neg/propneg_hc_1600_1800i_V532_equalmin2.rda")


nedgesnetworksCM
numberCM <- 56
numberCM2 <- 35
numberCM3 <- 101
nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]

numberBN <- 8
numberBNedge <- 9

nedgesIT[c(numberBNedge)]

posmixCN1 <- propagationCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(-2,2))
posmixCN2 <- propagationCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence =c(81,227), valueEvidence = c(-2,2))
posmixCN3 <- propagationCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(-2,2))
negmixCN1 <- propagationnegCorr(corcutsolo[[numberCM]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(-2,2))
negmixCN2 <- propagationnegCorr(corcutsolo[[numberCM2]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(-2,2))
negmixCN3 <- propagationnegCorr(corcutsolo[[numberCM3]], nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(81,227), valueEvidence = c(-2,2))

posmixBN1 <- prop_hc_1600_1800i_V81V227_equalmin2plus2
negmixBN1 <- propneg_hc_1600_1800i_V81V227_equalmin2plus2
# negmixBN1 <- propneg_hc_1600_1800i_V227_equal2

attr(posmixCN1, "sign") <- attr(posmixCN2, "sign") <- attr(posmixCN3, "sign") <- "pos"
attr(posmixBN1, "sign") <- "pos"
attr(negmixCN1, "sign") <-  attr(negmixCN2, "sign") <- attr(negmixCN3, "sign") <- "neg"
attr(negmixBN1, "sign") <- "neg"

propspos <- list(posmixBN1, posmixCN1, posmixCN2, posmixCN3) 
# propspos <- list(posmixCN1, posmixCN2, posmixCN3) 
propsneg <- list(negmixBN1, negmixCN1, negmixCN2, negmixCN3) 
# propsneg <- list(negmixCN1,negmixCN2, negmixCN3) 
props <- propspos

lapply(props, function(x) attr(x,"sign") )

logCM <- logliksCM[c(numberCM,numberCM2,numberCM3)]
logBN <- c(logliksIT[c(numberBNedge)])
logliksplot <- c(logBN,logCM)
nCM <- nedgesnetworksCM[c(numberCM,numberCM2,numberCM3)]
nBM <- nedgesIT[c(numberBNedge)]
nedgesplot <- c(nBM,nCM)
logliksplot

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,0.85,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  a
  b
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}


posmixCN1$with[81]

all <- cbind(plotdifferences)
blank <- grid.rect(gp=gpar(col="white"))
plotwindow <- list(all[[1]],all[[2]],blank,all[[3]],blank,all[[4]])
plotwindow <- list(all[[1]],all[[3]],blank,all[[4]])
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_170719_CN_minV81neg.pdf"

plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81.pdf"
plotname <- "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propagation/1800BN_2260_8066_170719_CN_minV81neg.pdf"
pdf(plotname)
do.call("grid.arrange", c(plotwindow, ncol=2))

dev.off()


diffCN1y2 <- quantity2clim((posmixCN2 - posmixCN3)$with, what = "diff", ref.grid = tas_ncep_10d)
n1 <- 478
posmixCN1[n1,]
posmixCN2[n1,]
posmixCN3[n1,]
spatialPlot(diffCN1y2, backdrop.theme = "coastline", lonCenter = 180)
negmixCN1[587,]
negmixCN2[587,]
negmixCN3[587,]


#######################################################################################
# Automatic
#######################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81plusV227neg/propneg_hc_1600_1800i_V81V227_equalmin2plus2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/minV81plusV227pos/prop_hc_1600_1800i_V81V227_equalmin2plus2.rda")

nodenumber <- 532
nodenumber <- c(81,227)

whichnode <- paste0("V",nodenumber)
if (length(nodenumber) >=1){
 whichnode <- paste0(whichnode, collapse = "")}

numberBN <- 8
numberBNedge <- 9
nedgesselectedBN <- nedgesIT[c(numberBNedge)]



files <- c(paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/min",whichnode,"pos/prop_hc_1600_1800i_",whichnode,"_equalmin2.rda"),
          paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/min",whichnode,"neg/propneg_hc_1600_1800i_",whichnode,"_equalmin2.rda"))
for (i in 1:length(files)){load(file = files[i])}


posmixBN1 <- eval(parse(text = paste0("prop_hc_1600_1800i_",whichnode,"_equalmin2plus2")))
negmixBN1 <- eval(parse(text = paste0("propneg_hc_1600_1800i_",whichnode,"_equalmin2plus2")))

posmixBN1 <- eval(parse(text = paste0("prop_hc_1600_1800i_",whichnode,"_equalmin2")))
negmixBN1 <- eval(parse(text = paste0("propneg_hc_1600_1800i_",whichnode,"_equalmin2")))
attr(posmixBN1, "sign") <- "pos"
attr(negmixBN1, "sign") <- "neg"
assign(paste0("negposBN_",whichnode),posmixBN1)
assign(paste0("negnegBN_",whichnode),negmixBN1)

nedgesnetworksCM
numberCM <- 56
numberCM2 <- 35
numberCM3 <- 101
numbersCM <- c(numberCM,numberCM2,numberCM3)
nedgesselectedCN <- nedgesnetworksCM[numbersCM]
dim(nedgesselectedCN) <- c(1,length(nedgesselectedCN))
nodenumber
all.equal(nodenumber,c(81,227))

posCNs <- lapply(numbersCM, function(x) {y <- propagationCorr(corcutsolo[[x]],nodesEvents = 1:648, valueEvent = 1, nodesEvidence = nodenumber, valueEvidence = c(-2,2))
                  attr(y, "sign") <- "pos" 
                  return(y)})
negCNs <- lapply(numbersCM, function(x) {y <- propagationnegCorr(corcutsolo[[x]],nodesEvents = 1:648, valueEvent = -1, nodesEvidence = nodenumber, valueEvidence = c(-2,2))
                  attr(y, "sign") <- "neg" 
                  return(y)})
posCNs <- lapply(numbersCM, function(x) {y <- propagationCorr(corcutsolo[[x]],nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(nodenumber), valueEvidence = c(-2))
                                              attr(y, "sign") <- "pos" 
                                              return(y)})
negCNs <- lapply(numbersCM, function(x) {y <- propagationnegCorr(corcutsolo[[x]],nodesEvents = 1:648, valueEvent = -1, nodesEvidence = c(nodenumber), valueEvidence = c(-2))
                                              attr(y, "sign") <- "neg" 
                                              return(y)})
attr(posCNs[[1]],"sign")
names(posCNs) <- apply(nedgesselectedCN, MARGIN = 1, FUN = function(x) paste0("CN",x))
names(negCNs) <- apply(nedgesselectedCN, MARGIN = 1, FUN = function(x) paste0("CN",x))

assign(paste0("negposCN_",whichnode),posCNs)
assign(paste0("negnegCN_",whichnode),negCNs)

negnegCN_V81$CN2260[81,]

propspos <- c(list(negposBN_V81V227), negposCN_V81V227) 
propsneg <- c(list(negnegBN_V81V227), negnegCN_V81V227) 

propspos <- c(list(negposBN_V81), negposCN_V81) 
propsneg <- c(list(negnegBN_V81), negnegCN_V81) 
propspos <- c(list(negposBN_V532), negposCN_V532) 
propsneg <- c(list(negnegBN_V532), negnegCN_V532) 
propspos <- c(list(negposBN_V568), negposCN_V568) 
propsneg <- c(list(negnegBN_V568), negnegCN_V568) 
props <- propspos
lapply(propspos, function(x)attr(x,"sign"))

logCM <- logliksCM[numbersCM]
logBN <- c(logliksIT[c(numberBNedge)])
logliksplot <- c(logBN,logCM)
nCM <- nedgesnetworksCM[numbersCM]
nBM <- nedgesIT[c(numberBNedge)]
nedgesplot <- c(nBM,nCM)

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()
Center7 <- 180
i <- 1
for(i in 1:length(props)){
  props[[1]]
  prop <- props[[i]]
  
  col.r <- colorRampPalette(brewer.pal(9, "Reds"))(100)
  col.b <-  colorRampPalette(brewer.pal(9, "Blues"))(100)
  if(attr(prop, "sign") == "pos"){cb <- col.r}
  if(attr(prop, "sign") == "neg"){cb <- col.b}
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.14,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center7, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center7, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",round(logliksplot[i])," |E| = ",nedgesplot[i]),cex = 0.6), at = seq(0.05,1,0.01),region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))

  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}



all <- cbind(plotdifferences)
blank <- grid.rect(gp=gpar(col="white"))
plotwindow <- list(all[[1]],all[[2]],blank,all[[3]],blank,all[[4]])
do.call("grid.arrange", c(plotwindow, ncol=2))

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropVFavourBN/",nedgesselectedBN,"BN_",nedgesselectedCN[1],"_",nedgesselectedCN[2],"_",nedgesselectedCN[3],"CN_min",whichnode,".pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropVFavourBN/",nedgesselectedBN,"BN_",nedgesselectedCN[1],"_",nedgesselectedCN[2],"_",nedgesselectedCN[3],"CN_min",whichnode,"neg.pdf")
plotname
pdf(plotname)
do.call("grid.arrange", c(plotwindow, ncol=2))
dev.off()

# specific colorscale
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropVFavourBN/",nedgesselectedBN,"BN_",nedgesselectedCN[1],"_",nedgesselectedCN[2],"_",nedgesselectedCN[3],"CN_min",whichnode,"_cut3.pdf")
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/PropVFavourBN/",nedgesselectedBN,"BN_",nedgesselectedCN[1],"_",nedgesselectedCN[2],"_",nedgesselectedCN[3],"CN_min",whichnode,"neg_cut3.pdf")
plotname
pdf(plotname)
do.call("grid.arrange", c(plotwindow, ncol=2))
dev.off()
