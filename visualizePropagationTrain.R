#######################################################################################
# Propagation train
#######################################################################################
rm(list = ls())
# traintest.
library(mvtnorm)
library(corpcor)
library(gridExtra)
library(condMVNorm)
library(transformeR)
library(visualizeR)
library(magrittr)
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
library(RColorBrewer)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
#########################################################################################
# Load data train
#########################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1train/hcTrain1_1500_1600i.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/indTRAIN1.rda")
indTEST1 <- (1:360)[-indTRAIN1]
##########################################################################################
# load data train
filestrain1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1train", full.names = T)
filestrain1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1train")
filestrain1names <- gsub(".rda", "", filestrain1names)
filestrain1names

variablelisttrain1 <- list(networks = list())
names <- c()

for (i in 1:length(filestrain1)){
  variablepos <- get(load(filestrain1[i]))
  variablelisttrain1$networks[[i]] <- variablepos
}

names(variablelisttrain1$networks) <- filestrain1names
assign(paste0("hc_1_train1"),variablelisttrain1)
nedgestrain1 <- sapply(hc_1_train1$networks, narcs)
sorted <- sort(nedgestrain1, index.return = TRUE)
hc_1_train1$networks <- hc_1_train1$networks[sorted$ix]
sapply(hc_1_train1$networks, narcs)
##############################################################################################
# load data test
#############################################################################################
filestest1 <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1test", full.names = T)
filestest1names <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1test")
filestest1names <- gsub(".rda", "", filestest1names)
filestest1names

variablelisttest1 <- list(networks = list())
names <- c()

for (i in 1:length(filestest1)){
  variablepos <- get(load(filestest1[i]))
  variablelisttest1$networks[[i]] <- variablepos
}

names(variablelisttest1$networks) <- filestest1names
assign(paste0("hc_1_test1"),variablelisttest1)
nedgestest1 <- sapply(hc_1_test1$networks, narcs)
sorted <- sort(nedgestest1, index.return = TRUE)
hc_1_test1$networks <- hc_1_test1$networks[sorted$ix]
sapply(hc_1_test1$networks, narcs)



###########################################################################
# Make fits on structure. 
###########################################################################
dataRMStrain <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTRAIN1))
dataRMStest <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE, subind = indTEST1))
fitstrain <- lapply(hc_1_train1$networks, bn.fit, data = dataRMStrain)
fitstest <- lapply(hc_1_test1$networks, bn.fit, data = dataRMStest)

loglikstrain <- c()
loglikstest <- c()
loglikstesttest <- c()
for (i in 1:length(fitstrain)){
  loglikstrain[i] <- logLik(fitstrain[[i]], data = dataRMStrain)
  loglikstesttest[i] <- logLik(fitstest[[i]], data = dataRMStest)
  loglikstest[i] <- logLik(fitstrain[[i]], data = dataRMStest)
}
nedgestrain1 <- sapply(hc_1_train1$networks, narcs)
nedgestest1 <- sapply(hc_1_test1$networks,narcs)
x <- seq(0,8800,100)
y <- seq(100,8900,100)

############################################################################################
# HC TRAIN propagation V81 V227 extra evidence (positive)
############################################################################################
proprange <- 10:18  
for (i in proprange){
  
  assign(paste0("prop_hc_perm1train_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1train_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81V227pos/prop_hc_perm1train_",x[i],"_",y[i],"i_V81V227_equal22.rda"))

  assign(paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81V227pos/prop_hc_perm1test_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
}

############################################################################################
# HC TRAIN propagation V81 V227 extra evidence (negative)
############################################################################################

for (i in proprange){
  
  assign(paste0("propneg_hc_perm1train_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_perm1train_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81V227neg/propneg_hc_perm1train_",x[i],"_",y[i],"i_V81V227_equal22.rda"))

  assign(paste0("propneg_hc_perm1test_",x[i],"_",y[i],"i_V81V227_equal22"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81,227),
                                     valueEvidence = c(2,2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_perm1test_",x[i],"_",y[i],"i_V81V227_equal22"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81V227neg/propneg_hc_perm1test_",x[i],"_",y[i],"i_V81V227_equal22.rda"))
}

#################################################################################
# Single evidence. V227 postive
#################################################################################
for (i in proprange){
  
  assign(paste0("prop_hc_perm1train_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1train_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V227pos/prop_hc_perm1train_",x[i],"_",y[i],"i_V227_equal2.rda"))

  assign(paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V227pos/prop_hc_perm1test_",x[i],"_",y[i],"i_V227_equal2.rda"))
}

#################################################################################
# Single evidence. V227 neg
#################################################################################
for (i in proprange){
  
  assign(paste0("propneg_hc_perm1train_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_perm1train_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V227neg/propneg_hc_perm1train_",x[i],"_",y[i],"i_V227_equal2.rda"))

  assign(paste0("propneg_hc_perm1test_",x[i],"_",y[i],"i_V227_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(227),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_perm1test_",x[i],"_",y[i],"i_V227_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V227neg/propneg_hc_perm1test_",x[i],"_",y[i],"i_V227_equal2.rda"))
}


#################################################################################
# Single evidence. V81 postive
#################################################################################
for (i in proprange){
  
  assign(paste0("prop_hc_perm1train_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1train_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81pos/prop_hc_perm1train_",x[i],"_",y[i],"i_V81_equal2.rda"))

  assign(paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81pos/prop_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2.rda"))
}

#################################################################################
# Single evidence. V81 neg
#################################################################################
for (i in proprange){
  
  assign(paste0("propneg_hc_perm1train_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_perm1train_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81neg/propneg_hc_perm1train_",x[i],"_",y[i],"i_V81_equal2.rda"))

  assign(paste0("propneg_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = "<= -1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("propneg_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81neg/propneg_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2.rda"))
}

###############################################################################
#
###############################################################################
assign(paste0("prop_hc_perm1train1500_1600i_V81_equal2"),
         PropagationExactGeneralPerm(baysnet = hcfits[[i]],
                                     nodesEvents = 1:648,
                                     valueEvent = ">= 1",
                                     nodesEvidence = c(81),
                                     valueEvidence = c(2),
                                     perm = permutations[[1]]))
  save(list = paste0("prop_hc_perm1train1500_1600i_V81_equal2"),
       file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/prop_hc_perm1train1500_1600i_V81_equal2.rda"))
  

  
  i <- 1
  proptrain16 <- PropagationExactGeneralPerm(baysnet = fitstrain[[i]],
                                             nodesEvents = 1:648,
                                             valueEvent = ">= 1",
                                             nodesEvidence = c(81),
                                             valueEvidence = c(2),
                                             perm = permutations[[1]])
  proptest16 <- PropagationExactGeneralPerm(baysnet = fitstest[[i]],
                                            nodesEvents = 1:648,
                                            valueEvent = ">= 1",
                                            nodesEvidence = c(81),
                                            valueEvidence = c(2),
                                            perm = permutations[[1]])

###############################################################################
# Load HC propagation iterations V81
###############################################################################
  files <- list.files("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81pos", full.names = T)
  fileslist <- list()
  names <- c()
  temp.space <- new.env()
  
  for (i in 1:length(files)){
    temp.space <- new.env()
    variable <- get(load(files[i], temp.space),envir = temp.space)
    names[i] <- ls(envir = temp.space)
    fileslist[[i]] <- variable
    rm(temp.space)
  }
  names(fileslist) <- names

  
  x <- seq(900,1700,100)
  y <- seq(1000,1800,100)
  
  fileslistord <- list()
  namesord <- c()
  for(i in 1:length(x)){
    fileslistord[[i]] <- eval(parse(text =paste0("fileslist$prop_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2")))
    namesord[i] <- paste0("prop_hc_perm1test_",x[i],"_",y[i],"i_V81_equal2")
  }
  names(fileslistord) <- namesord
  
  
###############################################################################
# Visualize V81 MAKE CORCUTSOLO TRAIN and TEST In traintest.r
###############################################################################
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81pos/prop_hc_perm1train_1600_1700i_V81_equal2.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81pos/prop_hc_perm1test_1600_1700i_V81_equal2.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81pos/prop_hc_perm1train_1700_1800i_V81_equal2.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81pos/prop_hc_perm1test_1700_1800i_V81_equal2.rda")
  numberCM <- 1
  numberBN <- 18
  
  propcortrain1 <- propagationCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  propcortest1 <- propagationCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81), valueEvidence = c(2))
  diffcor1 <- propcortrain1
  diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
  diffcor1$without <- 0
  
  prophctrain1 <- prop_hc_perm1train_1600_1700i_V81_equal2
  prophctest1 <- prop_hc_perm1test_1600_1700i_V81_equal2
  diffhc <- prophctrain1 
  diffhc$with <- abs(prophctrain1$with - prophctest1$with)
  diffhc$without <- 0
  
  props <- list(propcortrain1,propcortest1,diffcor1, prophctrain1,prophctest1,diffhc)

  logliksplot <- c(loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA)
  nedgesplot <- c(nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA)
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
    cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
    
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.025,0.875,0.01), region = TRUE, col.regions= cb,set.max = 0.875, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  all <- cbind(plotdifferences)
  all
  do.call("grid.arrange", c(all, nrow = 2))
  
  prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
  dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
  dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))
  
  do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))
  
  
  ###############################################################################
  # Visualize V81V227 MAKE CORCUTSOLO TRAIN and TEST. 
  ###############################################################################
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81V227pos/prop_hc_perm1train_1400_1500i_V81V227_equal22.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81V227pos/prop_hc_perm1test_1400_1500i_V81V227_equal22.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81V227pos/prop_hc_perm1train_1500_1600i_V81V227_equal22.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81V227pos/prop_hc_perm1test_1500_1600i_V81V227_equal22.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/train1/tr1_V81V227pos/prop_hc_perm1train_1700_1800i_V81V227_equal22.rda")
  load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/hc_iteration/test1/te1_V81V227pos/prop_hc_perm1test_1700_1800i_V81V227_equal22.rda")
  numberCM <- 1
  numberBN <- 18
  
  propcortrain1 <- propagationCorr(corcutsolotrain[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  propcortest1 <- propagationCorr(corcutsolotest[[numberCM]], nodesEvents = 1:648, valueEvent = 1, nodesEvidence = c(81,227), valueEvidence = c(2,2))
  diffcor1 <- propcortrain1
  diffcor1$with <- abs(propcortrain1$with - propcortest1$with)
  diffcor1$without <- 0
  
  prophctrain1 <- prop_hc_perm1train_1500_1600i_V81V227_equal22
  prophctest1 <- prop_hc_perm1test_1500_1600i_V81V227_equal22
  prophctrain1 <- prop_hc_perm1train_1400_1500i_V81V227_equal22
  prophctest1 <- prop_hc_perm1test_1400_1500i_V81V227_equal22
  prophctrain1 <- prop_hc_perm1train_1700_1800i_V81V227_equal22
  prophctest1 <- prop_hc_perm1test_1700_1800i_V81V227_equal22
  diffhc <- prophctrain1 
  diffhc$with <- abs(prophctrain1$with - prophctest1$with)
  diffhc$without <- 0
  
  props <- list(propcortrain1,propcortest1,diffcor1, prophctrain1,prophctest1,diffhc)
  
  logliksplot <- c(loglikstrain1CN[[numberCM]],loglikstesttest1CN[[numberCM]],NA,loglikstrain[[numberBN]],loglikstesttest[[numberBN]],NA)
  nedgesplot <- c(nedgesnetworkstrainCM[[numberCM]],nedgesnetworkstestCM[[numberCM]],NA, nedgestrain1[[numberBN]],nedgestest1[[numberBN]],NA)
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
    cb <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
    
    
    a <- spatialPlot(prop_without,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )), rev.colors = TRUE)
    b <- spatialPlot(prop_with,backdrop.theme = "coastline", lonCenter = Center9, main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05), rev.colors = TRUE)
    c <- spatialPlot(prop_dif,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
    
    plotwithouts[[i]] <- a
    plotwiths[[i]] <- b
    plotdifferences[[i]] <- c
  }
  all <- cbind(plotdifferences)
  all
  do.call("grid.arrange", c(all, nrow = 2, top = "train grow, test grow, abs(train-test)"))
 
  prop_dif_plus <- quantity2clim(abs(prophctrain1$with - prophctrain1$without - prophctest1$with + prophctest1$without), paste0(attr(prophctrain1$with, "probability"),"-", attr(prophctrain1$without, "probability")), tas_ncep_10d)
  dif1 <- spatialPlot(prop_dif_plus,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  prop_dif_plus2 <- quantity2clim(abs(propcortrain1$with - propcortrain1$without - propcortest1$with + propcortest1$without), paste0(attr(propcortrain1$with, "probability"),"-", attr(propcortrain1$without, "probability")), tas_ncep_10d)
  dif2 <- spatialPlot(prop_dif_plus2,backdrop.theme = "coastline", lonCenter = Center9, main = list(paste0(attr(prop_dif_plus$Data,"climatology:fun"),"\n loglik = ",logliksplot[i]," |E| = ",nedgesplot[i]),cex = 0.5), at = seq(0.05,0.85,0.01), region = TRUE, col.regions= cb,set.max = 0.85, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  do.call("grid.arrange", c(list(dif1,dif2), nrow = 1, top = "abs(train grow - test grow)"))
  # all.equal(test5$with-test5$without,test8$with-test8$without)
  # all.equal(test1$with-test1$without,test4$with-test4$without)  
  
  do.call("grid.arrange", c(append(all,list(dif1,dif2)),nrow = 3))
  