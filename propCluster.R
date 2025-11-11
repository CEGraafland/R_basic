#####################################################################################
# Sensitivity Analysis Bayesian Networks
#####################################################################################
library(transformeR, lib.loc = "vols/....")
library(magrittr)
library(bnlearn, lib.loc = "vols/....")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/utilnetworks_hcnetworks10d.rda")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")

tryPropagationExactGeneral <- function(baysnet, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # baysnet <- fitted3
  # nodesEvents <- 284:285
  # valueEvent <- ">= 1"
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V ", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
  }
  
  
  if (length(nodesEvidence) > 1){
    proves <- c()
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0("V",nodesEvidence[j]," = ", valueEvidence[j])
    }
    
    text <- "list("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],",")
    }
    text <- paste0(text,proves[length(nodesEvidence)],")")
    str2 <- text
    
    
    probname <- paste0("P(V ", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  # str2
  # i <- 1
  
  for(i in 1:length(nodesEvents)) {
    l <- nodesEvents[i]
    l
    # str <- paste0("(", names(baysnet)[l], ">=", valueEvent, ")")
    str <- paste0("(", names(baysnet)[l], valueEvent, ")")
    str
    nparams(baysnet)
    cmd1 = paste0("cpquery(baysnet, ", str, ", ", str2, ", method = ","'lw'",")")
    cmd3 = paste0("cpquery(baysnet, ", str, ", ", "TRUE", ", method = ","'lw'",")")
    cmd1
    cmd3
    
    with[i] <- eval(parse(text = cmd1))
    with[i]
    
    without[i] <- eval(parse(text = cmd3))
    without[i]
    
    
    # with[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
    # withcomplement[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
    # without[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = TRUE)
    
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V ", valueEvent,")")
  df <- data.frame(names = names(baysnet)[nodesEvents], with = with, without = without)
  return(df)
}

dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
dfRMS <- as.data.frame(dataRMS)
example3 <- hc_edges_loglik_10d_1400_1600i$networks
fitted3 <- bn.fit(example3, dfRMS)

propHC1600_Vpos_V280_equal2 <- tryPropagationExactGeneral(fitted3, 1:648,">= 1", c(280), c(2))
propHC1600_Vneg_V280_equal2 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(280), c(2))
propHC1600_Vneg_V81V280_equal22 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(81, 280), c(2,2))

save(propHC1600_Vpos_V280_equal2,
     propHC1600_Vneg_V280_equal2,
     propHC1600_Vneg_V81V280_equal22,
     file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_propagationV81V280cluster.rda")