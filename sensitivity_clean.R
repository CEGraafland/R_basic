#####################################################################################
# Sensitivity Analysis Bayesian Networks
#####################################################################################
library(transformeR)
library(visualizeR)
library(magrittr)
library(bnlearn)
library(ggplot2)
library(grid)
library(gridExtra)
library(gridGraphics)
library(lattice)
rm(list = ls())
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/utilnetworks_hcnetworks10d.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_graphsPC.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
#####################################################################################
# Evidence propagation function
# necesario timeCoordsAnom_from_Grid_rms              (same as FINAL FUNCTION in "old)
# -> functions in propagationFunctions.
#####################################################################################

tryPropagation <- function(baysnet, nodesEvents, valueEvent, valueEvidence){
  with <- numeric(length = length(nodesEvents))
  withcomplement <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  str2 <- paste0("(V81 >=", valueEvidence, ")")
  str3 <- paste0("(V81 <", valueEvidence, ")")
  
  
  for(i in 1:length(nodesEvents)) {
    l <- nodesEvents[i]
    str <- paste0("(", names(baysnet)[l], ">=", valueEvent, ")")
    cmd1 = paste("cpquery(baysnet, ", str, ", ", str2, ")", sep = "")
    cmd2 = paste("cpquery(baysnet, ", str, ", ", str3, ")", sep = "")
    cmd3 = paste("cpquery(baysnet, ", str, ", ", "TRUE", ")", sep = "")
    with[i] <- eval(parse(text = cmd1))
    withcomplement[i] <- eval(parse(text = cmd2))
    without[i] <- eval(parse(text = cmd3))
    # with[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
    # withcomplement[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
    # without[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = TRUE)
    
  }
  
  return(data.frame(names = names(baysnet)[nodesEvents], with = with, withcomplement = withcomplement, without = without))
}
tryPropagationV459 <- function(baysnet, nodesEvents, valueEvent, valueEvidence){
  with <- numeric(length = length(nodesEvents))
  withcomplement <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  str2 <- paste0("(V459 >=", valueEvidence, ")")
  str3 <- paste0("(V459 <", valueEvidence, ")")
  
  
  for(i in 1:length(nodesEvents)) {
    l <- nodesEvents[i]
    str <- paste0("(", names(baysnet)[l], ">=", valueEvent, ")")
    cmd1 = paste("cpquery(baysnet, ", str, ", ", str2, ")", sep = "")
    cmd1
    cmd2 = paste("cpquery(baysnet, ", str, ", ", str3, ")", sep = "")
    cmd3 = paste("cpquery(baysnet, ", str, ", ", "TRUE", ")", sep = "")
    with[i] <- eval(parse(text = cmd1))
    withcomplement[i] <- eval(parse(text = cmd2))
    without[i] <- eval(parse(text = cmd3))
    # with[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
    # withcomplement[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
    # without[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = TRUE)
    
  }
  
  return(data.frame(names = names(baysnet)[nodesEvents], with = with, withcomplement = withcomplement, without = without))
}

tryPropagationExact <- function(baysnet, nodesEvents, valueEvent, valueEvidence){
  
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  str2 <- paste0("list(V81 = ", valueEvidence, ")")
  
  i <- 99
  for(i in 1:length(nodesEvents)) {
    l <- nodesEvents[i]
    l
    str <- paste0("(", names(baysnet)[l], ">=", valueEvent, ")")
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
  
  return(data.frame(names = names(baysnet)[nodesEvents], with = with, without = without))
}

tryPropagationExactGeneral <- function(baysnet, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # fitted <- bn.fit(hc_edges_loglik_10d_1400_1600i$networks, as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)))
  # baysnet <- fitted
  # nodesEvents <- c(81,280)
  # valueEvent <- ">= 1"
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
 
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V ", valueEvent,"|V",nodesEvidence[1]," = ", valueEvidence[1], ")")
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

# tryPropagationExactGeneral(bn.fit(hc_edges_loglik_10d_1400_1600i$networks, as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))),
# c(298,299),">=1",c(81,280),c(2,2))

PropagationExactGeneralPerm <- function(baysnet, nodesEvents, valueEvent, nodesEvidence, valueEvidence, perm){
  
  # baysnet <- fitted
  # nodesEvents <- c(81,280)
  # valueEvent <- ">= 1"
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(2,2)
  # dataperm <- datapermutations[[1]]
  # perm <- permutations[[1]]
  
  # baysnet = fitted
  # nodesEvents = c(298,299)
  # valueEvent = ">=1"
  # nodesEvidence = c(81,280)
  # valueEvidence = c(2,2)
  # perm = permutations[[1]]
  
  
  
  nodesEventsRef <- c()
  for (i in 1:length(nodesEvents)){
    nodesEventsRef[i] <- which(perm == nodesEvents[i])
  }
  
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
  # i <- 2
  
  for(i in 1:length(nodesEvents)) {
    # l <- nodesEvents[i]
    # l
    l <- nodesEventsRef[i]
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
  df <- data.frame(names = names(baysnet)[nodesEventsRef], with = with, without = without)
  return(df)
  
}

# initialdata <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE))
# fitted <- bn.fit(hc_edges_loglik_10d_1400_1600i$networks, as.data.frame(initialdata[permutations[[1]]]))
# test1 <- PropagationExactGeneralPerm(baysnet = fitted,
#                             nodesEvents = c(298,299),
#                             valueEvent = ">=1",
#                             nodesEvidence = c(81,280),
#                             valueEvidence = c(2,2),
#                             perm = permutations[[1]])


# tryPropagationAprox <- function(baysnet, nodesEvents, valueEvent, valueEvidence){
#   with <- numeric(length = length(nodesEvents))
#   withcomplement <- numeric(length = length(nodesEvents))
#   without <- numeric(length = length(nodesEvents))
#   
#   str2 <- paste0("(V81 >=", valueEvidence, ")")
#   str3 <- paste0("(V81 <", valueEvidence, ")")
#   
#   
#   for(i in 1:length(nodesEvents)) {
#     l <- nodesEvents[i]
#     str <- paste0("(", names(baysnet)[l], ">=", valueEvent, ")")
#     cmd1 = paste("cpquery(baysnet, ", str, ", ", str2, ")", sep = "")
#     cmd2 = paste("cpquery(baysnet, ", str, ", ", str3, ")", sep = "")
#     cmd3 = paste("cpquery(baysnet, ", str, ", ", "TRUE", ")", sep = "")
#     with[i] <- eval(parse(text = cmd1))
#     withcomplement[i] <- eval(parse(text = cmd2))
#     without[i] <- eval(parse(text = cmd3))
#     # with[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
#     # withcomplement[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
#     # without[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = TRUE)
#     
#   }
#   
#   return(data.frame(names = names(baysnet)[nodesEvents], with = with, withcomplement = withcomplement, without = without))
# }
 
# cpquery(fitted3, event = (V285 >= 1), evidence = list(V81 =2,V280 = 2), method = "lw")
# 
# propHC1600_V81_equal1 <- tryPropagationExact(fitted3, 1:5,1,1)
# test <- tryPropagationExactGeneral(fitted3, c(266),1, c(81, 280), c(2,2))
# test2 <- tryPropagationExactGeneral(fitted3, c(266), "<=-1", c(81), c(5))
# test3 <- tryPropagationExactGeneral(fitted3, c(266),1, c(280), c(2))
# test
# test2
# test3
#################################################################################
# Comando propagation
# Exact evidence, condtional V81 V280
#################################################################################
# test_Vneg_V81V280_equal22 <- tryPropagationExactGeneral(fitted3, 283:285,"<= -1", c(81, 280), c(2,2))
# test_Vpos_V81V280_equal22 <- tryPropagationExactGeneral(fitted3, 1:648,">= 1", c(81, 280), c(2,2))
# testHC1600_Vneg_V81V280V424_equal222 <- tryPropagationExactGeneral(fitted3, 265:266,"<= -1", c(81, 280,424), c(2,2,2))
# str(testHC1600_Vneg_V81V280V424_equal222)
# test_Vneg_V81V280_equal22
# propHC1600_Vpos_V81V280_equal22
# propHC1600_Vpos_V81_equal2
# 
# a <- test_Vneg_V81V280_equal22[,2]
# b <- propHC1600_Vneg_V81_equal2[c(283,284,285),2]
# c <- propHC1600_Vneg_V81_equal2[c(283,284,285),3]
# 
# 
# propHC1600_Vneg_V81_equal2[c(283,284,285),]
# propHC1600_Vpos_V81_equal1 <- tryPropagationExactGeneral(fitted3, 1:648,">= 1", c(81), c(1))
# propHC1600_Vneg_V81_equal1 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(81), c(1))
# propHC1600_Vpos_V81_equal2 <- tryPropagationExactGeneral(fitted3, 1:648,">= 1", c(81), c(2))
# propHC1600_Vneg_V81_equal2 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(81), c(2))

# propHC1600_Vpos_V280_equal2 <- tryPropagationExactGeneral(fitted3, 1:648,">= 1", c(280), c(2))
# # propHC1600_Vneg_V280_equal2 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(280), c(2))

# propHC1600_Vneg_V81V280_equal22 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(81, 280), c(2,2))
# propHC1600_Vpos_V81V280_equal22 <- tryPropagationExactGeneral(fitted3, 1:648,">= 1", c(81, 280), c(2,2))
# # propHC1600_Vpos_V81V280_equal52 <- tryPropagationExactGeneral(fitted3, 1:648,"<= -1", c(81, 280), c(5,2))
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280_t1.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280_t2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280_t3.rda")
str(propHC1600_Vneg_V81V280_equal22)
str(propHC1600_Vneg_V280_equal2)
# attr(propHC1600_Vpos_V81V280_equal22$with, "probability") <- "P(V > 1 |V81 = 2, V280 = 2)"
# attr(propHC1600_Vpos_V81V280_equal22$without, "probability") <- "P(V > 1)"
# attr(propHC1600_Vpos_condV280V81_22_condV81_2$with, "probability") <- "P(V > 1 |V81 = 2, V280 = 2)" 
# attr(propHC1600_Vpos_condV280V81_22_condV81_2$without, "probability") <- "P(V > 1 |V81 = 2)"
# attr(propHC1600_Vpos_condV280V81_22_condV280_2$with, "probability") <- "P(V > 1 |V81 = 2, V280 = 2)" 
# attr(propHC1600_Vpos_condV280V81_22_condV280_2$without, "probability") <- "P(V > 1 |V280 = 2)"
# attr(propHC1600_Vpos_V81_equal2$with, "probability") <- "P(V > 1 |V81 = 2)" 
# attr(propHC1600_Vpos_V81_equal2$without, "probability") <- "P(V > 1)" 
# attr(propHC1600_Vneg_V81_equal2$with, "probability") <- "P(V < -1 |V81 = 2)" 
# attr(propHC1600_Vneg_V81_equal2$without, "probability") <- "P(V < -1)"
# attr(propHC1600_Vneg_V81_equal1$with, "probability") <- "P(V < -1 |V81 = 1)" 
# attr(propHC1600_Vneg_V81_equal1$without, "probability") <- "P(V < -1)"
# attr(propHC1600_Vpos_V81_equal1$with, "probability") <- "P(V > 1 |V81 = 1)" 
# attr(propHC1600_Vpos_V81_equal1$without, "probability") <- "P(V > 1)" 


# save(propHC1600_Vpos_V81_equal1,
#      propHC1600_Vneg_V81_equal1,
#      propHC1600_Vpos_V81_equal2,
#      propHC1600_Vneg_V81_equal2,
#      propHC1600_Vpos_V280_equal2,
#      propHC1600_Vpos_V81V280_equal22,
#      propHC1600_Vneg_V81V280_equal22,
#      file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280.rda")
#################################################################################
# Add together:
# P(V > 1 |V81V280)       P(V > 1 |V81V280) 
# P(V > 1 |V81)           P(V > 1 |V280)  
# P(V > 1)
#################################################################################

# # replace without of propHC1600_Vpos_V81V280_equal22 with P(V > 1 |V81).
# propHC1600_Vpos_condV280V81_22_condV81_2 <- propHC1600_Vpos_V81V280_equal22
# propHC1600_Vpos_condV280V81_22_condV81_2[3] <- propHC1600_Vpos_V81_equal2[2]
# # replace without of propHC1600_Vpos_V81V280_equal22 with P(V > 1 |V280).
# propHC1600_Vpos_condV280V81_22_condV280_2 <- propHC1600_Vpos_V81V280_equal22
# propHC1600_Vpos_condV280V81_22_condV280_2[3] <- propHC1600_Vpos_V280_equal2[2]
# # check
# all.equal(propHC1600_Vpos_condV280V81_22_condV81_2,propHC1600_Vpos_V81V280_equal22)

# # replace without of propHC1600_Vneg_V81V280_equal22 with P(V < -1 |V81).
# propHC1600_Vneg_condV280V81_22_condV81_2 <- propHC1600_Vneg_V81V280_equal22
# propHC1600_Vneg_condV280V81_22_condV81_2[3] <- propHC1600_Vneg_V81_equal2[2]
# str(propHC1600_Vneg_condV280V81_22_condV81_2)

# # replace without of propHC1600_Vneg_V81V280_equal22 with P(V <- 1 |V280).
# propHC1600_Vneg_condV280V81_22_condV280_2 <- propHC1600_Vneg_V81V280_equal22
# propHC1600_Vneg_condV280V81_22_condV280_2[3] <- propHC1600_Vneg_V280_equal2[2]

save(propHC1600_Vpos_condV280V81_22_condV81_2,
     propHC1600_Vpos_condV280V81_22_condV280_2,
     propHC1600_Vneg_condV280V81_22_condV81_2,
     propHC1600_Vneg_condV280V81_22_condV280_2,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280_cond.rda")


# replace without of propHC1600_Vpos_V280V424_equal22 with P(V > 1 |V280).
propHC1600_Vpos_condV280V424_22_condV280_2 <- propHC1600_Vpos_V280V424_equal22
propHC1600_Vpos_condV280V424_22_condV280_2[3] <- propHC1600_Vpos_V280_equal2[2]

# replace without of propHC1600_Vpos_V280V424_equal22 with P(V > 1 |V424).
propHC1600_Vpos_condV280V424_22_condV424_2 <- propHC1600_Vpos_V280V424_equal22
propHC1600_Vpos_condV280V424_22_condV424_2[3] <- propHC1600_Vpos_V424_equal2[2]

# replace without of propHC1600_Vneg_V280V424_equal22 with P(V < -1 |V280).
propHC1600_Vneg_condV280V424_22_condV280_2 <- propHC1600_Vneg_V280V424_equal22
propHC1600_Vneg_condV280V424_22_condV280_2[3] <- propHC1600_Vneg_V280_equal2[2]

# replace without of propHC1600_Vneg_V280V424_equal22 with P(V < -1 |V424).
propHC1600_Vneg_condV280V424_22_condV424_2 <- propHC1600_Vneg_V280V424_equal22
propHC1600_Vneg_condV280V424_22_condV424_2[3] <- propHC1600_Vneg_V424_equal2[2]

save(propHC1600_Vpos_condV280V424_22_condV280_2,
     propHC1600_Vpos_condV280V424_22_condV424_2,
     propHC1600_Vneg_condV280V424_22_condV280_2,
     propHC1600_Vneg_condV280V424_22_condV424_2,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV280V424_cond.rda")

propHC1600_Vpos_V81V424_equal22
# replace without of propHC1600_Vpos_V81V424_equal22 with P(V > 1 |V81).
propHC1600_Vpos_condV81V424_22_condV81_2 <- propHC1600_Vpos_V81V424_equal22
propHC1600_Vpos_condV81V424_22_condV81_2[3] <- propHC1600_Vpos_V81_equal2[2]

# replace without of propHC1600_Vpos_V81V424_equal22 with P(V > 1 |V424).
propHC1600_Vpos_condV81V424_22_condV424_2 <- propHC1600_Vpos_V81V424_equal22
propHC1600_Vpos_condV81V424_22_condV424_2[3] <- propHC1600_Vpos_V424_equal2[2]

# replace without of propHC1600_Vneg_V81V424_equal22 with P(V < -1 |V81).
propHC1600_Vneg_condV81V424_22_condV81_2 <- propHC1600_Vneg_V81V424_equal22
propHC1600_Vneg_condV81V424_22_condV81_2[3] <- propHC1600_Vneg_V81_equal2[2]

# replace without of propHC1600_Vneg_V81V424_equal22 with P(V < -1 |V424).
propHC1600_Vneg_condV81V424_22_condV424_2 <- propHC1600_Vneg_V81V424_equal22
propHC1600_Vneg_condV81V424_22_condV424_2[3] <- propHC1600_Vneg_V424_equal2[2]

save(propHC1600_Vpos_condV81V424_22_condV81_2,
     propHC1600_Vpos_condV81V424_22_condV424_2,
     propHC1600_Vneg_condV81V424_22_condV81_2,
     propHC1600_Vneg_condV81V424_22_condV424_2,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V424_cond.rda")


#################################################################################
# Add together: withouts by propagations that are result of propagation Simple
#################################################################################
# laad de pos files en neg files in een list
filespos <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation", full.names = T, pattern = "HC1600.Vpos.")
filenamespos <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation", pattern = "HC1600.Vpos.")
filenamespos <- gsub(".rda", "", filenamespos)
filesneg <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation", full.names = T, pattern = "HC1600.Vneg.")
filenamesneg <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation", pattern = "HC1600.Vneg.")
filenamesneg <- gsub(".rda", "", filenamesneg)

# pattern1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/", pattern = "[HC1600]{2}")
# pattern2 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/", pattern = "[equal2]{2}")

# Add without to each positivie file
withoutpos <- propHC1600_Vpos_V81_equal2$without

variablelistpos <- list()
for (i in 1:length(filespos)){
  variablepos <- get(load(filespos[i]))
  variablelistpos[[i]] <- variablepos
}
names(variablelistpos) <- filenamespos
str(variablelistpos)

for (i in 1:length(variablelistpos)){
  if (length(variablelistpos[[i]]) <= 3){
    variablelistpos[[i]]$without <- withoutpos
    assign(names(variablelistpos[i]),variablelistpos[[i]])
  }
}
str(variablelistpos)

# Add without to each neg file.
withoutneg <- propHC1600_Vneg_V81_equal2$without
variablelistneg <- list()
for (i in 1:length(filesneg)){
  variableneg <- get(load(filesneg[i]))
  variablelistneg[[i]] <- variableneg
}
names(variablelistneg) <- filenamesneg
str(variablelistneg)

for (i in 1:length(variablelistneg)){
  if (length(variablelistneg[[i]]) <= 3){
    variablelistneg[[i]]$without <- withoutneg
    assign(names(variablelistneg[i]),variablelistneg[[i]])
  }
}

# Save pos en neg files (overwrite)
for (i in 1:length(filespos)){
  save(list = filenamespos[i], file = filespos[i])
}

for (i in 1:length(filesneg)){
  save(list = filenamesneg[i], file = filesneg[i])
}
#################################################################################
# Load Files propagation
#################################################################################
files<- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation", full.names = T)
files
for (i in 1:length(files)){
  load(files[i])
}

#################################################################################
# Comando propagation
# data with Root Mean Square.  (script BayesianGraph)
#################################################################################
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
str(dataRMS)
dfRMS <- as.data.frame(dataRMS)
dfRMSPC <- dfRMS
colnames(dfRMSPC) <- as.character(1:648)
colnames(dfRMSPC)


example <- hc_edges_loglik_10d_200_400i$networks
example1 <- hc_edges_loglik_10d_400_600i$networks
example2 <- hc_edges_loglik_10d_600_800i$networks
example3 <- hc_edges_loglik_10d_1400_1600i$networks
example4 <- hc_edges_loglik_10d_2400_2600i$networks
examplePC <- pc_10d_1e_13$networks$dag
examplePC

# save(example, example1, example2,example3,example4,examplePC, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_for_sens10d.rda")

fitted <- bn.fit(example, dfRMS)
fitted1 <- bn.fit(example1, dfRMS)
fitted2 <- bn.fit(example2, dfRMS)
fitted3 <- bn.fit(example3, dfRMS)
fitted4 <- bn.fit(example4, dfRMS)

fittedPC <- bn.fit(examplePC, dfRMSPC) 
nodes(fittedPC) <- nodes(example1) # recuperacion of node names
nodes(fittedPC)
nparams(fitted3)*500
# propHC400_V81_1  <- tryPropagation(fitted, 1:648,1,1)
# propHC800_V81_1  <- tryPropagation(fitted2, 1:648,1,1)
# propHC1600_V81_1 <- tryPropagation(fitted3, 1:648,1,1)
# propHC600_V81_1  <- tryPropagation(fitted1, 1:648,1,1)
# propPC600_V81_1  <- tryPropagation(fittedPC, 1:648,1,1)
# propHC2600_V81_1 <- tryPropagation(fitted4, 1:648,1,1)
# 
# try2 <- tryPropagation(fitted,80:90,1,1)
# try
# try2
# save(propHC400_V81_1, propHC800_V81_1, propHC1600_V81_1, file = 
#        "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81.rda")

# propHC1600_V459_1 <- tryPropagationV459(fitted3, 1:648,1,1)
# save(propHC1600_V459_1, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV459.rda")

#################################################################################
# Comando propagation automatic
# data with Root Mean Square. 
#################################################################################
str(tas_ncep_10d)
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
dfRMS <- as.data.frame(dataRMS) 

# network <- hc_edges_loglik_10d_2400_2600i$networks
# model <- bn.fit(network, dfRMS)

# propagation  <- tryPropagation(model, c(81,99,286,604,371),1,1)
# propagation2 <- tryPropagation(model, c(323,244,248,540,387,77),1,1)
# rbind(propagation,propagation2)
# propHC400_V81_1[c(81,99,286,604,371,323,244,248,540,387,77),]
# propHC800_V81_1[c(81,99,286,604,371,323,244,248,540,387,77),]
# propHC1600_V81_1[c(81,99,286,604,371,323,244,248,540,387,77),]

# save(propagation, propagation2, propHC400_V81_1, propHC600_V81_1, propHC800_V81_1, propHC1600_V81_1, propHC2600_V81_1, file = 
#        "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81.rda")
# save(propHC600_V81_1,propHC2600_V81_1, propPC600_V81_1,file = 
#        "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81_2.rda")
#################################################################################
# Visualize Bayesian propagation with plotClimatology
#################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81neg.rda")
prop <- propHC1600_V81_equal1
prop <- propHC1600_V81_neg1

matwithout <- matrix(prop$without, nrow = 1)
matwith <- matrix(prop$with, nrow = 1)
matcomplement <- matrix(prop$withcomplement, nrow = 1)
matdif <- matwith - matwithout
matdifcompl <- matcomplement - matwithout

prop_without <- tas_ncep_10d
prop_with <- tas_ncep_10d
prop_dif <- tas_ncep_10d
prop_complement <- tas_ncep_10d
prop_difcompl <- tas_ncep_10d

prop_without$Data <- mat2Dto3Darray(matwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_with$Data <- mat2Dto3Darray(matwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_dif$Data <- mat2Dto3Darray(matdif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_complement$Data <- mat2Dto3Darray(matcomplement, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
prop_difcompl$Data <- mat2Dto3Darray(matdifcompl, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))

attr(prop_without$Data, "climatology:fun") <- "without evidence"
attr(prop_with$Data, "climatology:fun") <- "with evidence"
attr(prop_dif$Data, "climatology:fun") <- "with - without evidence"
attr(prop_complement$Data, "climatology:fun") <- "complement of evidence"
attr(prop_difcompl$Data, "climatology:fun") <- "withcomplement - without evidence"
paste0("(V81 >=",1, ")")


plotClimatology(prop_without,backdrop.theme = "coastline", main = list(paste0("P(V >=1)")))
plotClimatology(prop_with,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ")")),at = seq(0.0,1,0.01))
plotClimatology(prop_dif,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ") -P(V >=1)")), at = seq(0.1,0.3,0.01))
plotClimatology(prop_complement,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81<",1, ")")))
plotClimatology(prop_difcompl,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81<",1, ") -P(V >=1)")))

#####################################################################################
# Methode 2:
# Visualize Byasian propagation using quantity2clim and PlotClimatology
#####################################################################################
test <- quantity2clim(propHC1600_Vpos_condV280V81_22_condV81_2$without, "P(V |V81 = 2)", tas_ncep_10d)
plotClimatology(test, backdrop.theme = "coastline", main = list(attr(test$Data,"climatology:fun")))
test
str(test)
attr(climatology(tas_ncep_10d)$Dates, "season")
#################################################################################
# Visualize Bayesian propagation with plotClimatology
# Set difference plotting parameter
# #################################################################################
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(gridGraphics)
# library(lattice)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81_2.rda")

props <- list(propHC1600_V459_1)
props <- list(propHC400_V81_1,propHC800_V81_1, propHC1600_V81_1, propHC2600_V81_1)
props <- list(propHC600_V81_1, propPC600_V81_1)
props <- list(propPC600_V81_1,propHC600_V81_1,propHC800_V81_1,propHC1600_V81_1)
props <- list(propHC1600_Vpos_condV280V81_22_condV81_2, propHC1600_Vpos_V81_equal2, propHC1600_Vpos_V280_equal2)

plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  prop <- props[[i]]

  matwithout <- matrix(prop$without, nrow = 1)
  matwith <- matrix(prop$with, nrow = 1)
  # matcomplement <- matrix(prop$withcomplement, nrow = 1)
  matdif <- matwith - matwithout
  # matdifcompl <- matwithout - matcomplement

  prop_without <- tas_ncep_10d
  prop_with <- tas_ncep_10d
  prop_dif <- tas_ncep_10d
  # prop_complement <- tas_ncep_10d
  # prop_difcompl <- tas_ncep_10d

  prop_without$Data <- mat2Dto3Darray(matwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  prop_with$Data <- mat2Dto3Darray(matwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  prop_dif$Data <- mat2Dto3Darray(matdif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # prop_complement$Data <- mat2Dto3Darray(matcomplement, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # prop_difcompl$Data <- mat2Dto3Darray(matdifcompl, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))

  attr(prop_without$Data, "climatology:fun") <- "without evidence"
  attr(prop_with$Data, "climatology:fun") <- "with evidence"
  attr(prop_dif$Data, "climatology:fun") <- "with - without evidence"
  # attr(prop_complement$Data, "climatology:fun") <- "complement of evidence"
  # attr(prop_difcompl$Data, "climatology:fun") <- "withcomplement - without evidence"

  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))

  a <- plotClimatology(prop_without,backdrop.theme = "coastline",main = list(paste0("P(V >=1)"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- plotClimatology(prop_with,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ")"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- plotClimatology(prop_dif,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81>=",1, ") -P(V >=1)"),cex = 0.5), at = seq(0.05,0.5,0.005), set.max = 0.5, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  # d <- plotClimatology(prop_complement,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81<",1, ")"), cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  # e <- plotClimatology(prop_difcompl,backdrop.theme = "coastline",colorkey = list(width = 0.6, lables = list(cex = 0.5)))

  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}


sizes <- list(textGrob("HC 1600"),textGrob("HC 1600"))
sizes <- list(textGrob("HC 400"), textGrob("HC 800"), textGrob("HC 1600"), textGrob("HC 2600"))
sizes <- list(textGrob("HC 600"), textGrob("PC 600"))
sizes <- list(textGrob("PC 600"), textGrob("HC 600"), textGrob("HC 800"), textGrob("HC 1600"))
all <- c(sizes,plotwithouts, plotwiths, plotdifferences)
all <- c(plotwithouts, plotwiths, plotdifferences)
length(all)

n <- length(plotwithouts)
m <- length(all)/n

nCol <- n
nRow <- m

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/V81HC400_2600_dif03.pdf")
pdf(plotname, height = 6, width = 10)
do.call("grid.arrange", c(all, ncol=n, nrow=m))
dev.off()

plotname2 <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/V81HCPC_600_dif03.pdf")
pdf(plotname2, height = 7, width = 5)
do.call("grid.arrange", c(all, ncol=n, nrow=m))
dev.off()

#################################################################################
# Visualize Bayesian propagation with plotClimatology
# Set difference plotting parameter
# method 2:
# #################################################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81_2.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81V280_cond.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV280V424_cond.rda")

source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")

props <- list(propHC1600_V459_1)
props <- list(propHC400_V81_1,propHC800_V81_1, propHC1600_V81_1, propHC2600_V81_1)
props <- list(propHC600_V81_1, propPC600_V81_1)
props <- list(propPC600_V81_1,propHC600_V81_1,propHC800_V81_1,propHC1600_V81_1)
props1 <- list(propHC1600_Vpos_V81V280_equal22, propHC1600_Vpos_condV280V81_22_condV81_2, propHC1600_Vpos_condV280V81_22_condV280_2, propHC1600_Vpos_V81_equal2, propHC1600_Vpos_V280_equal2)
props2 <- list(propHC1600_Vneg_V81V280_equal22, propHC1600_Vneg_condV280V81_22_condV81_2, 
              propHC1600_Vneg_condV280V81_22_condV280_2, 
              propHC1600_Vneg_V81_equal2,
              propHC1600_Vneg_V280_equal2
              )
props1 <- list(propHC1600_Vpos_V280V424_equal22, propHC1600_Vpos_condV280V424_22_condV280_2, propHC1600_Vpos_condV280V424_22_condV424_2, propHC1600_Vpos_V280_equal2, propHC1600_Vpos_V424_equal2)
props2 <- list(propHC1600_Vneg_V280V424_equal22, propHC1600_Vneg_condV280V424_22_condV280_2, 
               propHC1600_Vneg_condV280V424_22_condV424_2, 
               propHC1600_Vneg_V280_equal2,
               propHC1600_Vneg_V424_equal2
               )

props1 <- list(propHC1600_Vpos_V81V424_equal22, propHC1600_Vpos_condV81V424_22_condV81_2, propHC1600_Vpos_condV81V424_22_condV424_2, propHC1600_Vpos_V81_equal2, propHC1600_Vpos_V424_equal2)
props2 <- list(propHC1600_Vneg_V81V424_equal22, propHC1600_Vneg_condV81V424_22_condV81_2, 
               propHC1600_Vneg_condV81V424_22_condV424_2, 
               propHC1600_Vneg_V81_equal2,
               propHC1600_Vneg_V424_equal2
)
props <- list(propHC1600_Vpos_V280V424_equal22, propHC1600_Vpos_V424_equal2, propHC1600_Vpos_V280_equal2, propHC1600_Vneg_V280_equal2)
props <- rbind(props1,props2)


plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  prop <- props[[i]]
  
  # matwithout <- matrix(prop$without, nrow = 1)
  # matwith <- matrix(prop$with, nrow = 1)
  # # matcomplement <- matrix(prop$withcomplement, nrow = 1)
  # matdif <- matwith - matwithout
  # # matdifcompl <- matwithout - matcomplement
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  # prop_complement <- tas_ncep_10d
  # prop_difcompl <- tas_ncep_10d
  
  # prop_without$Data <- mat2Dto3Darray(matwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # prop_with$Data <- mat2Dto3Darray(matwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # prop_dif$Data <- mat2Dto3Darray(matdif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # # prop_complement$Data <- mat2Dto3Darray(matcomplement, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # # prop_difcompl$Data <- mat2Dto3Darray(matdifcompl, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # 
  # attr(prop_without$Data, "climatology:fun") <- "without evidence"
  # attr(prop_with$Data, "climatology:fun") <- "with evidence"
  # attr(prop_dif$Data, "climatology:fun") <- "with - without evidence"
  # # attr(prop_complement$Data, "climatology:fun") <- "complement of evidence"
  # # attr(prop_difcompl$Data, "climatology:fun") <- "withcomplement - without evidence"
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- plotClimatology(prop_without,backdrop.theme = "coastline",main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- plotClimatology(prop_with,backdrop.theme = "coastline", main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- plotClimatology(prop_dif,backdrop.theme = "coastline", main = list(attr(prop_dif$Data,"climatology:fun"),cex = 0.5), at = seq(0.05,0.8,0.0125), set.max = 0.8, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  # d <- plotClimatology(prop_complement,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81<",1, ")"), cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  # e <- plotClimatology(prop_difcompl,backdrop.theme = "coastline",colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}

sizes <- list(textGrob("HC 1600"),textGrob("HC 1600"), textGrob("HC 1600"))
all <- cbind(plotdifferences)
all
length(all)

n <- length(plotwithouts)
m <- length(all)/n

nCol <- n
nRow <- m
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/propHC1600_V81V424_dif05.pdf")
pdf(plotname, height = 10, width = 6)
do.call("grid.arrange", c(all, ncol=2, nrow=5))
dev.off()

#################################################################################
# Visualize Bayesian propagation with plotClimatology
# HIER AAN HET WERK!!!!
# Compare eBIC obtained propagation
# #################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC/prop_1_tabu_10d_g0.1_V81_equal2.rda")

source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
for(j in c(1,2)){
  pattern <- paste0("prop_",j,"_tabu_10d_g")
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/tabu_eBIC", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  
  propagationlisttabu1 <- list()
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    propagationlisttabu1[[i]] <- variablepos
  }
  
  
  names(propagationlisttabu1) <- filestabu1names
  assign(paste0("proplist_",j,"_tabu_10d"),propagationlisttabu1)
}

props <-list(prop_1_tabu_10d_g0_V81V280_equal22,prop_2_tabu_10d_g0.75_V81V280_equal22)
props <- proplist_2_tabu_10d


plotwiths <- list()
plotwithouts <- list()
plotdifferences <- list()

for(i in 1:length(props)){
  prop <- props[[i]]
  
  # matwithout <- matrix(prop$without, nrow = 1)
  # matwith <- matrix(prop$with, nrow = 1)
  # # matcomplement <- matrix(prop$withcomplement, nrow = 1)
  # matdif <- matwith - matwithout
  # # matdifcompl <- matwithout - matcomplement
  
  prop_with <- quantity2clim(prop$with, attr(prop$with, "probability"), tas_ncep_10d)
  prop_without <- quantity2clim(prop$without, attr(prop$without, "probability"), tas_ncep_10d)
  prop_dif <- quantity2clim(prop$with - prop$without, paste0(attr(prop$with, "probability"),"-", attr(prop$without, "probability")), tas_ncep_10d)
  
  str(prop_with)
  
  # prop_complement <- tas_ncep_10d
  # prop_difcompl <- tas_ncep_10d
  
  # prop_without$Data <- mat2Dto3Darray(matwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # prop_with$Data <- mat2Dto3Darray(matwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # prop_dif$Data <- mat2Dto3Darray(matdif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # # prop_complement$Data <- mat2Dto3Darray(matcomplement, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # # prop_difcompl$Data <- mat2Dto3Darray(matdifcompl, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  # 
  # attr(prop_without$Data, "climatology:fun") <- "without evidence"
  # attr(prop_with$Data, "climatology:fun") <- "with evidence"
  # attr(prop_dif$Data, "climatology:fun") <- "with - without evidence"
  # # attr(prop_complement$Data, "climatology:fun") <- "complement of evidence"
  # # attr(prop_difcompl$Data, "climatology:fun") <- "withcomplement - without evidence"
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- plotClimatology(prop_without,backdrop.theme = "coastline",main = list(attr(prop_without$Data,"climatology:fun"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- plotClimatology(prop_with,backdrop.theme = "coastline", main = list(attr(prop_with$Data,"climatology:fun"),cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)),at = seq(0,1,0.05))
  c <- plotClimatology(prop_dif,backdrop.theme = "coastline", main = list(attr(prop_dif$Data,"climatology:fun"),cex = 0.5), at = seq(0.05,0.8,0.0125), set.max = 0.8, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  # d <- plotClimatology(prop_complement,backdrop.theme = "coastline", main = list(paste0("P(V >=1|V81<",1, ")"), cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  # e <- plotClimatology(prop_difcompl,backdrop.theme = "coastline",colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  
  plotwithouts[[i]] <- a
  plotwiths[[i]] <- b
  plotdifferences[[i]] <- c
}


all <- cbind(plotdifferences)
all
length(all)

n <- length(plotwithouts)
m <- length(all)/n

nCol <- n
nRow <- m

do.call("grid.arrange", c(all))
dev.off()





############################################################################################
# Prepare session P(V < -1|V81>=1)
# 
############################################################################################
negPropagation <- function(baysnet, nodesEvents, valueEvent, valueEvidence){
  negwith <- numeric(length = length(nodesEvents))
  negwithcomplement <- numeric(length = length(nodesEvents))
  negwithout <- numeric(length = length(nodesEvents))
  
  str2 <- paste0("(V81 >=", valueEvidence, ")")
  str3 <- paste0("(V81 <", valueEvidence, ")")
  
  # baysnet <- fitted
  # i <- 1
  # nodesEvents <- 80:90
  # valueEvent <- -1
  # valueEvidence <- 1
   for(i in 1:length(nodesEvents)) {
    l <- nodesEvents[i]
    str <- paste0("(", names(baysnet)[l], "<=",valueEvent, ")")
    cmd1 = paste("cpquery(baysnet, ", str, ", ", str2, ")", sep = "")
    cmd2 = paste("cpquery(baysnet, ", str, ", ", str3, ")", sep = "")
    cmd3 = paste("cpquery(baysnet, ", str, ", ", "TRUE", ")", sep = "")
    cmd1
    negwith[i] <- eval(parse(text = cmd1))
    negwithcomplement[i] <- eval(parse(text = cmd2))
    negwithout[i] <- eval(parse(text = cmd3))
    # with[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
    # withcomplement[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
    # without[i] <- cpquery(baysnet, event = eval(parse(text = str)), evidence = TRUE)
    
  }
  
  return(data.frame(names = names(baysnet)[nodesEvents], negwith = negwith, negwithcomplement = negwithcomplement, negwithout = negwithout))
}

###############################################################################################
# System time tests and commandos 
###############################################################################################
negPropagation(fitted3,80:90,-1,1)
system.time(negPropagation(fitted3,274:284,-1,1))
846/60
65*15/60
system.time(negPropagation(fitted4,274:284,-1,1))
2137/60/11*648/60

# propHC400_V81_neg1  <- negPropagation(fitted, 1:648,-1,1)
# propHC600_V81_neg1 <- negPropagation(fitted1, 1:648,-1,1)
# propHC800_V81_neg1 <- negPropagation(fitted2, 1:648,-1,1)
# propHC1600_V81_neg1 <- negPropagation(fitted3, 1:648,-1,1)
# propHC2600_V81_neg1 <- negPropagation(fitted4, 1:648,-1,1)
# propPC600_V81_neg1 <- negPropagation(fittedPC, 1:648,-1,1)

load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81neg.rda")
save(propHC400_V81_neg1, propHC600_V81_neg1, propPC600_V81_neg1, propHC800_V81_neg1, propHC1600_V81_neg1, propHC2600_V81_neg1, file = 
       "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81neg.rda")
################################################################################################
# visualize negative possibilties
################################################################################################
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propagation/data_propagationV81neg.rda")

props <- list(propHC400_V81_neg1, propHC800_V81_neg1, propHC1600_V81_neg1, propHC2600_V81_neg1)
props <- list(propHC600_V81_neg1, propPC600_V81_neg1)
props <- list(propPC600_V81_neg1, propHC600_V81_neg1, propHC800_V81_neg1, propHC1600_V81_neg1)

plotnegwiths <- list()
plotnegwithouts <- list()
plotnegdifferences <- list()

for(i in 1:length(props)){
  prop <- props[[i]]
  
  matnegwithout <- matrix(prop$negwithout, nrow = 1)
  matnegwith <- matrix(prop$negwith, nrow = 1)
  matnegcomplement <- matrix(prop$negwithcomplement, nrow = 1)
  matnegdif <- matnegwith - matnegwithout
  matnegdifcompl <- matnegwithout - matnegcomplement
  
  prop_negwithout <- tas_ncep_10d
  prop_negwith <- tas_ncep_10d
  prop_negdif <- tas_ncep_10d
  prop_negcomplement <- tas_ncep_10d
  prop_negdifcompl <- tas_ncep_10d
  
  prop_negwithout$Data <- mat2Dto3Darray(matnegwithout, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  prop_negwith$Data <- mat2Dto3Darray(matnegwith, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  prop_negdif$Data <- mat2Dto3Darray(matnegdif, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  prop_negcomplement$Data <- mat2Dto3Darray(matnegcomplement, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  prop_negdifcompl$Data <- mat2Dto3Darray(matnegdifcompl, x = attr(dataRMS, "Xcoords"), y = attr(dataRMS, "Ycoords"))
  
  attr(prop_negwithout$Data, "climatology:fun") <- "without evidence"
  attr(prop_negwith$Data, "climatology:fun") <- "with evidence"
  attr(prop_negdif$Data, "climatology:fun") <- "with - without evidence"
  attr(prop_negcomplement$Data, "climatology:fun") <- "complement of evidence"
  attr(prop_negdifcompl$Data, "climatology:fun") <- "withcomplement - without evidence"
  
  
  col.l<- colorRampPalette(c('purple', 'blue', 'cyan'))
  
  a <- plotClimatology(prop_negwithout,backdrop.theme = "coastline",main = list(paste0("P(V <=-1)"),cex = 0.5), at = c(0,0.13,0.18,1) ,region = TRUE, col.regions= col.l,colorkey = list(col = col.l, width = 0.6, at = c(0,0.13,0.18,1),lables = list(cex = 0.5, labels =c("0","0.13","0.18","1"),at = c(0,0.13,0.18,1) )))
  b <- plotClimatology(prop_negwith,backdrop.theme = "coastline", main = list(paste0("P(V <=-1|V81>=",1, ")"),cex = 0.5), at = seq(0,1,0.05), colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  c <- plotClimatology(prop_negdif,backdrop.theme = "coastline", main = list(paste0("P(V <=-1|V81>=",1, ") -P(V <=-1)"),cex = 0.5), at = seq(0,0.3,0.05), set.min = 0.03, set.max = 0.3, colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  d <- plotClimatology(prop_negcomplement,backdrop.theme = "coastline", main = list(paste0("P(V <=-1|V81<",1, ")"), cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  e <- plotClimatology(prop_negdifcompl,backdrop.theme = "coastline", main = list(paste0("P(V <=-1) - P(V <=-1|V81<",1, ")"), cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))
  e
  plotnegwithouts[[i]] <- a
  plotnegwiths[[i]] <- b
  plotnegdifferences[[i]] <- c
}

sizes <- list(textGrob("HC 400"), textGrob("HC 800"), textGrob("HC 1600"), textGrob("HC 2600"))
sizes <- list(textGrob("HC 600"), textGrob("PC 600"))
all <- c(sizes,plotnegwithouts, plotnegwiths, plotnegdifferences)
sizes

n <- length(plotnegwithouts)
m <- length(all)/n

nCol <- n
nRow <- m

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/V81negHC400_2600_dif03.pdf")
pdf(plotname, height = 6, width = 10)
do.call("grid.arrange", c(all, ncol=m, nrow=n))
dev.off()

plotname2 <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/plots/V81negHCPC_600_dif03.pdf")
pdf(plotname2, height = 7, width = 5)
do.call("grid.arrange", c(all, ncol=n, nrow=m))
dev.off()

################################################################################################
# visualize possibilties for Marco
################################################################################################
props <- list(propPC600_V81_1,propHC600_V81_1,propHC1600_V81_1)
propsneg <- list(propPC600_V81_neg1, propHC600_V81_neg1, propHC800_V81_neg1, propHC1600_V81_neg1)


sizes <- list(textGrob("PC 603"), textGrob("HC 600"), textGrob("HC 1596"))
all <- c(sizes, plotwiths, plotdifferences)
do.call("grid.arrange", c(all, ncol=3, nrow=3))

plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Resumen/V81PCHC6001600dif03.pdf")
pdf(plotname, height = 10, width = 10)
do.call("grid.arrange", c(all, ncol=3, nrow=3))
dev.off()

max(propHC600_V81_1$without)
min(propHC600_V81_1$without)

hist(propPC600_V81_1$without)
hist(propHC600_V81_1$without)
hist(propHC1600_V81_1$without)


propHC1600_V81_1$with[117]

quantiles()
############################################################################################
#
############################################################################################
set.seed(1)
cpquery(fitted3, event = (V99>= 1), evidence = (V81 >=1))
cpquery(fitted3, event = (V99 >= 1), evidence = TRUE)
cpquery(fitted3, event = (V99>= 1), evidence = (V81 >=1))
cpquery(fitted3, event = (V99 >= 1), evidence = TRUE)
