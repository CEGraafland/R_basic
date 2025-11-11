##################################################################################
# Sensitivity Analysis Bayesian networks
# Exploring example                                   TEST BEFORE DELETE
##################################################################################
example <- hc_edges_loglik_10d_400_600i$networks
str(example)
# Make bn 
fitted <- bn.fit(example, as.data.frame(graphObject$data_coords))

with <- cpquery(fitted, event = (V99 >= 0.5), evidence = (V81 >= 0.5))
withcomplement <- cpquery(fitted, event = (V299 >= 0.5), evidence = (V81 < 0.5))
without <- cpquery(fitted, event = (V99 >= 0.5), evidence = TRUE)

probs600V99 <- c(with,withcomplement,without)
probs600V286 <- c(with,withcomplement,without)
probs1000V99 <- c(with,withcomplement,without)
probs1600V99 <- c(with,withcomplement,without)
probs1000V286 <- c(with,withcomplement,without)
probs1600V286 <- c(with,withcomplement,without)

probs600V99
probs600V286
probs1000V99
probs1600V99
probs1000V286
probs1600V286


igraphdir <- igraph.from.graphNEL(as.graphNEL(example))
igraphske<- as.undirected(igraphdir)
# create graph_from_Grid object with TimeCoordsAnom_from_Grid
# and add graph and adjacency separately.
graphObject <- TimeCoordsAnom_from_Grid_aslist(tas_ncep_10d)
graphObject$graph <- igraphske
graphObject$adjacency <- as_adjacency_matrix(igraphske)
check
plot.Meteodag(graphObject,igraphdir)
dev.off()
##################################################################################
# Sensitivity Analysis Bayesian networks
# First seperated loops                                         
##################################################################################
example <- hc_edges_loglik_10d_200_400i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
str(example)
# Make bn 

fitted <- bn.fit(example, dfRMS)

with <- list()

names(fitted)[1]
str <- paste0("(", names(fitted)[i], ">=", 1, ")")
parse(text = str)
for(i in c(1,81,99,229)) {
  str <- paste0("(", names(fitted)[i], ">=", 1, ")")
  with[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 >= 1))
}
with[c(1,81,99,229)]

######## ya está hecho Lis ;)     ----
withcomplement <- list()
for(i in c(1,81,99,229)) {
  str <- paste0("(", names(fitted)[i], ">=", 1, ")")
  withcomplement[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 < 1))
}
withcomplement[c(1,81,99,229)]
########

without <- list()
for(i in c(1,81,99,229)) {
  str <- paste0("(", names(fitted)[i], ">=", 1, ")")
  without[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = TRUE)
}
system.time()
without[c(1,81,99,229)]
data.frame(names,with,withcomplement,without)
#####################################################################################
# Combine GREAT
# loops with fixed evidence
#####################################################################################
example <- hc_edges_loglik_10d_1400_1600i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
str(example)
# Make bn 

fitted <- bn.fit(example, dfRMS)

with <- list()
withcomplement <- list()
without <- list()
names(fitted)[1]
str <- paste0("(", names(fitted)[i], ">=", 1, ")")
parse(text = str)
for(i in c(1,81,99,229)) {
  str <- paste0("(", names(fitted)[i], ">=", 1, ")")
  with[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 >= 1))
  withcomplement[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 < 1))
  without[[i]] <- cpquery(fitted, event = eval(parse(text = str)), evidence = TRUE)
}


names(fitted)[c(1,81,99,229)]
with[c(1,81,99,229)]
withcomplement[c(1,81,99,229)]
without[c(1,81,99,229)]
#####################################################################################
# Combine GREAT loops
# Working Script!!                SCRIPT WITH FALSE OUTPUT (FULL ENVIRONMENT)
#####################################################################################
example <- hc_edges_loglik_10d_600_800i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
str(example)
# Make bn 

fitted <- bn.fit(example, dfRMS)

with <- c()
withcomplement <- c()
without <- c()
names(fitted)[1]
str <- paste0("(", names(fitted)[i], ">=", 1, ")")
parse(text = str)
for(i in 1:length(names(fitted))) {
  str <- paste0("(", names(fitted)[i], ">=", 1, ")")
  with <- append(with, cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 >= 1)))
  withcomplement <- append(withcomplement, cpquery(fitted, event = eval(parse(text = str)), evidence = (V81 < 1)))
  without <- append(without, cpquery(fitted, event = eval(parse(text = str)), evidence = TRUE))
}

names(fitted)
length(with)
withcomplement
without

data.frame(names = names(fitted)[1:2], with = with, withcomplement = withcomplement, without = without)
#####################################################################################
# Combine 
# loops in function
# necesario timeCoordsAnom_from_Grid_rms              FINAL FUNCTION
#####################################################################################
library(transformeR)
library(magrittr)
library(bnlearn)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda")
rm(str2)
example <- hc_edges_loglik_10d_200_400i$networks
example2 <- hc_edges_loglik_10d_600_800i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
fitted <- bn.fit(example, dfRMS)
baysnet <- fitted
fitted2 <- bn.fit(example2, dfRMS)


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


propHC400_V81_1  <- tryPropagation(fitted, 1:648,1,1)
propHC800_V81_1  <- tryPropagation(fitted2, 1:648,1,1)
propHC800_V81_1
try2 <- tryPropagation(fitted,80:90,1,1)
try
try2
save(propHC400_V81_1, file = 
       "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_propagationV81.rda")
#####################################################################################
# Test of seperated loops on gaussian.test GREAT
#####################################################################################
data(gaussian.test)
dagnormal <- hc(gaussian.test)
stdgaussian.test <- scale(gaussian.test)
stdgaussian.test <- as.data.frame(stdgaussian.test)
dagstd <- hc(stdgaussian.test)
fittednor <- bn.fit(dagnormal,gaussian.test)
fittedstd<- bn.fit(dagstd,stdgaussian.test)
fittedcross<- bn.fit(dagstd, gaussian.test)

cpquery(fitted,
        event = ((A >= 0) & (A <= 1)) & ((B >= 0) & (B <= 3)),
        evidence = (C + D < 10))


with <- list()
names(fittednor)[1]
str <- paste0("(", names(fittednor)[i], ">=", 0, ") & (", names(fittednor)[i], "<=", 1, ")")
parse(text = str)
for(i in 1:length(names(fittednor))) {
  str <- paste0("(", names(fittednor)[i], ">=", 0, ") & (", names(fittednor)[i], "<=", 1, ")")
  with[[i]] <- cpquery(fittednor, event = eval(parse(text = str)), evidence = (C + D < 10))
}
with

with <- list()
names(fittedstd)[1]
str <- paste0("(", names(fittedstd)[i], ">=", 1, ")")
parse(text = str)
for(i in 1:length(names(fittedstd))) {
  str <- paste0("(", names(fittedstd)[i], ">=", 1, ")")
  with[[i]] <- cpquery(fittedstd, event = eval(parse(text = str)), evidence = (C > 1))
}
with

without <- list()
names(fittedstd)[1]
str <- paste0("(", names(fittedstd)[i], ">=", 1, ")")
parse(text = str)
for(i in 1:length(names(fittedstd))) {
  str <- paste0("(", names(fittedstd)[i], ">=", 1, ")")
  without[[i]] <- cpquery(fittedstd, event = eval(parse(text = str)), evidence = TRUE)
}
without
##################################################################################
# Sensitivity Analysis Bayesian networks
# First propagation function                    GAVE SIMILAR RESULTS WITH FULL ENVIRONMENT
# Loops all in one                              GIVES ERROR BECAUSE STR2 is not founded with empty environment
##################################################################################
BN <- fitted
evidencevalue <- 0.5
eventvalue <- 1
startnode <- 640
endnode <- 640
rm(BN,evidencevalue,eventvalue,startnode,endnode)
startnode:endnode
rm(i)
eventnodes <- c(81,99)

Propagation <- function(BN, eventnodes, eventvalue, evidencevalue){
  
  a <- evidencevalue
  b <- eventvalue
  str2 <- paste0("(V81 >= ", a, ")")
  str3 <- paste0("(V81 < ", a, ")")
  str2
  str3
  with <- c()
  withcomplement <- c()
  without <- c()
  nod <- c()
  # str1 <- c()
  #for(i in eventnodes) {
  #  str1[i] <- paste0("(", names(BN)[i], ">=", b, ")")
  #}
  # str1
  
  for(i in eventnodes) {
    str <- paste0("(", names(BN)[i], ">=", b, ")")
    with[i] <- cpquery(BN, event = eval(parse(text = str)), evidence = eval(parse(text = str2)))
    withcomplement[i] <- cpquery(BN, event = eval(parse(text = str)), evidence = eval(parse(text = str3)))
    without[i] <- cpquery(BN, event = eval(parse(text = str)), evidence = TRUE)
    nod[i] <- names(BN)[i]
    str <- NULL
  }
  
  
  return(list(nodes, with, withcomplement, without))
}


example <- hc_edges_loglik_10d_200_400i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
fitted <- bn.fit(example, dfRMS)
str(example)

test<- Propagation(BN = fitted, eventnodes = c(81,99), eventvalue = 1, evidencevalue = 1)
test
scale(df[,81], center = FALSE, scale = TRUE)

########################################################################################
# Second Propagation function 
# different way: with list and elementwiseadding of one vector per node
########################################################################################
Propagation2 <- function(BN, eventnodes, eventvalue, evidencevalue){
  
  a <- evidencevalue
  b <- eventvalue
  str2 <- paste0("(V81 >= ", a, ")")
  str3 <- paste0("(V81 < ", a, ")")
  
  
  list <- list()
  
  for(i in eventnodes) {
    str1 <- paste0("(", names(BN)[i], ">=", b, ")")
    with <- cpquery(BN, event = eval(parse(text = str1)), evidence = eval(parse(text = str2)))
    withcomplement <- cpquery(BN, event = eval(parse(text = str1)), evidence = eval(parse(text = str3)))
    without <- cpquery(BN, event = eval(parse(text = str1)), evidence = TRUE)
    nod <- names(BN)[i]
    vec <- c(nod,with,withcomplement,without)
    list[[i]] <- vec
    rm(str1,vec)
  }
  
  
  
  return(list)
}

example <- hc_edges_loglik_10d_200_400i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
fitted <- bn.fit(example, dfRMS)


testseq2 <- Propagation2(BN = fitted, eventnodes = c(280,267,322,323,371,556,604), eventvalue = 1, evidencevalue = 0.5)
testseq1 <- Propagation(BN = fitted, startnode = 308, endnode = 323, eventvalue = 1, evidencevalue = 0.5)


testseq1[250:500,]
testseq2[c(280,267,322,323,371,556,604)]

testseq3 <- Propagation2(BN = fitted, eventnodes = c(1,81,99,229), eventvalue = 1, evidencevalue = 1)
testseq3[c(1,81,99,229)]

######################################################################################
# without SAVE
#######################################################################################


Propagation <- function(BN, startnode, endnode, eventvalue, evidencevalue){
  
  a <- evidencevalue
  b <- eventvalue
  str2 <- paste0("(V81 >= ", a, ")")
  str3 <- paste0("(V81 < ", a, ")")
  str2
  str3
  with <- c()
  withcomplement <- c()
  without <- c()
  nod <- c()
  
  for(i in startnode:endnode) {
    str1 <- paste0("(", names(BN)[i], ">=", b, ")")
    with[i] <- cpquery(BN, event = eval(parse(text = str1)), evidence = eval(parse(text = str2)))
    withcomplement[i] <- cpquery(BN, event = eval(parse(text = str1)), evidence = eval(parse(text = str3)))
    without[i] <- cpquery(BN, event = eval(parse(text = str1)), evidence = TRUE)
    nod[i] <- names(BN)[i]
  }
  
  out <- data.frame(nodes = nod,
                    with = with, 
                    withcomplement = withcomplement,
                    without = without)
  
  return(out)
}
















#############################################################################3
# Estimate propagtion time
#############################################################################
example <- hc_edges_loglik_10d_1400_1600i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)) 
df <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
fitted <- bn.fit(example, dfRMS)

length(nodes(fitted))
test <- Propagation(BN= fitted,startnode= 284,endnode = 285,eventvalue = 0.5,evidencevalue = 1)

testV206 <- Propagation(BN= fitted,startnode= 206,endnode = 206,eventvalue = 0.5,evidencevalue = 1)
testV280 <- Propagation(BN= fitted,startnode= 280,endnode = 280,eventvalue = 0.5,evidencevalue = 1)
testV273 <- Propagation(BN= fitted,startnode= 273,endnode = 273,eventvalue = 0.5,evidencevalue = 1)
testV371 <- Propagation(BN= fitted,startnode= 371,endnode = 371,eventvalue = 0.5,evidencevalue = 1)
testV481 <- Propagation(BN= fitted,startnode= 481,endnode = 481,eventvalue = 0.5,evidencevalue = 1)
testV556 <- Propagation(BN= fitted,startnode= 556,endnode = 556,eventvalue = 0.5,evidencevalue = 1)
testV528 <- Propagation(BN= fitted,startnode= 528,endnode = 528,eventvalue = 0.5,evidencevalue = 1)
testV601 <- Propagation(BN= fitted,startnode= 601,endnode = 601,eventvalue = 0.5,evidencevalue = 1)

testframe <- data.frame(nodes = vector(mode = "factor", length = 648), 
                        with = vector(mode = "numeric", length = 648),
                        withcomplement = vector(mode = "numeric", length = 648),
                        without = vector(mode = "numeric", length = 648))
testframe[206,] <- testV206[206,]
testframe[280,] <- testV280[280,]
testframe[273,] <- testV273[273,]
testframe[371,] <- testV371[371,]
testframe[481,] <- testV481[481,]
testframe[528,] <- testV528[528,]
testframe[556,] <- testV556[556,]
testframe[601,] <- testV601[601,]

testframe

rbind(testV206[206,],
      testV280[280,],
      testV273[273,],
      testV371[371,],
      testV481[481,],
      testV528[528,],
      testV556[556,],
      testV601[601,])

propagationbn3[c(206,280,273,371,481,528,556,601),]

class(testV206[206,1])
testframe[206,]
test[c(284,285),]
propagationbn3[c(206,280,273,371,481,528,556,601),]
system.time(Propagation(fitted,284,285,0.5,1), gcFirst = TRUE)

#################################################################################
# test B
#################################################################################
example <- hc_edges_loglik_10d_1400_1600i$networks
dfRMS <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d), rms = TRUE)
fitted <- bn.fit(example, df)

length(nodes(fitted))
test <- Propagation(BN= fitted,startnode= 284,endnode = 285,eventvalue = 1,evidencevalue = 0.5)

testV206b <- Propagation(BN= fitted,startnode= 206,endnode = 206,eventvalue = 1,evidencevalue = 0.5)
testV280b <- Propagation(BN= fitted,startnode= 280,endnode = 280,eventvalue = 1,evidencevalue = 0.5)
testV273b <- Propagation(BN= fitted,startnode= 273,endnode = 273,eventvalue = 1,evidencevalue = 0.5)
testV371b <- Propagation(BN= fitted,startnode= 371,endnode = 371,eventvalue = 1,evidencevalue = 0.5)
testV481b <- Propagation(BN= fitted,startnode= 481,endnode = 481,eventvalue = 1,evidencevalue = 0.5)
testV556b <- Propagation(BN= fitted,startnode= 556,endnode = 556,eventvalue = 1,evidencevalue = 0.5)
testV528b <- Propagation(BN= fitted,startnode= 528,endnode = 528,eventvalue = 1,evidencevalue = 0.5)
testV601b <- Propagation(BN= fitted,startnode= 601,endnode = 601,eventvalue = 1,evidencevalue = 0.5)

testframeb <- data.frame(nodes = vector(mode = "double", length = 648), 
                         with = vector(mode = "numeric", length = 648),
                         withcomplement = vector(mode = "numeric", length = 648),
                         without = vector(mode = "numeric", length = 648))
testframeb[206,] <- testV206b[206,]
testframeb[280,] <- testV280b[280,]
testframeb[273,] <- testV273b[273,]
testframeb[371,] <- testV371b[371,]
testframeb[481,] <- testV481b[481,]
testframeb[528,] <- testV528b[528,]
testframeb[556,] <- testV556b[556,]
testframeb[601,] <- testV601b[601,]

testframeb

rbind(testV206b[206,],
      testV280b[280,],
      testV273b[273,],
      testV371b[371,],
      testV481b[481,],
      testV528b[528,],
      testV556b[556,],
      testV601b[601,])

testframe
propagationbn2[c(206,280,273,371,481,528,556,601),]


#################################################################################
# Comando propagation
# data without Root Mean Square. 
#################################################################################
df <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
example1 <- hc_edges_loglik_10d_200_400i$networks
example2 <- hc_edges_loglik_10d_600_800i$networks
example3 <- hc_edges_loglik_10d_1400_1600i$networks
bn1 <- bn.fit(example1, df)
bn2 <- bn.fit(example2, df)
bn3 <- bn.fit(example3, df)

testbn2c <- Propagation(BN = bn2, startnode = 647, endnode = length(nodes(bn2)), evidencevalue = 0.5)

propagationbn1 <- Propagation(BN = bn1,startnode = 1, endnode = length(nodes(bn1)), evidencevalue = 0.5)
propagationbn2 <- Propagation(BN = bn2,startnode = 1, endnode = length(nodes(bn2)), evidencevalue = 0.5)
propagationbn3 <- Propagation(BN = bn3,startnode = 1, endnode = length(nodes(bn3)), evidencevalue = 0.5)

save(df,example1,example2,example3,bn1,bn2,bn3,propagationbn1,propagationbn2,propagationbn3, 
     file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propationbns.rda")


testbn2 <- Propagation(bn2,640,length(nodes(bn2)),0.5)
testbn2b <-  Propagation(bn2,258,265,0.5)
testbn2c <- Propagation(BN = bn2, startnode = 647, endnode = length(nodes(bn2)),eventvalue =1 , evidencevalue = 0.5)

system.time(Propagation(bn2,640,length(nodes(bn2)),0.5))
testbn2[640:length(nodes(bn2)),]
testbn2b[258:265,]
testbn2c[647:648,]

#################################################################################
# Comando propagation
# data withRoot Mean Square. 
#################################################################################
dfRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
example1 <- hc_edges_loglik_10d_200_400i$networks
example2 <- hc_edges_loglik_10d_600_800i$networks
example3 <- hc_edges_loglik_10d_1400_1600i$networks
bn1 <- bn.fit(example1, as.data.frame(dfRMS))
bn2 <- bn.fit(example2, as.data.frame(dfRMS))
bn3 <- bn.fit(example3, as.data.frame(dfRMS))

testbn2c <- Propagation(BN = bn2, startnode = 647, endnode = length(nodes(bn2)), evidencevalue = 1)

propagationbn1 <- Propagation(BN = bn1,startnode = 1, endnode = length(nodes(bn1)), evidencevalue = 0.5)
propagationbn2 <- Propagation(BN = bn2,startnode = 1, endnode = length(nodes(bn2)), evidencevalue = 0.5)
propagationbn3 <- Propagation(BN = bn3,startnode = 1, endnode = length(nodes(bn3)), evidencevalue = 0.5)

save(df,example1,example2,example3,bn1,bn2,bn3,propagationbn1,propagationbn2,propagationbn3, 
     file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/propationbns.rda")


testbn2 <- Propagation(bn2,640,length(nodes(bn2)),0.5)
testbn2b <-  Propagation(bn2,258,265,0.5)
testbn2c <- Propagation(BN = bn2, startnode = 647, endnode = length(nodes(bn2)), evidencevalue = 0.5)

system.time(Propagation(bn2,640,length(nodes(bn2)),0.5))
testbn2[640:length(nodes(bn2)),]
testbn2b[258:265,]
testbn2c[647:648,]


#################################################################################
# Visualize Bayesianpropagation with plotClimatology
#################################################################################
str(propagationbn1)
propagationbn1$without
propagationbn1$with
matwithout <- matrix(propagationbn1$without, nrow = 1)
matwith <- matrix(propagationbn1$with, nrow = 1)
matdif <- matwith - matwithout
propbn1_tas_10d_without <- tas_ncep_10d
propbn1_tas_10d_with <- tas_ncep_10d
propbn1_tas_10d_dif <- tas_ncep_10d
propbn1_tas_10d_without$Data <- mat2Dto3Darray(matwithout, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
propbn1_tas_10d_with$Data <- mat2Dto3Darray(matwith, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
propbn1_tas_10d_dif$Data <- mat2Dto3Darray(matdif, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
attr(propbn1_tas_10d_without$Data, "climatology:fun") <- "without evidence"
attr(propbn1_tas_10d_with$Data, "climatology:fun") <- "with evidence"
attr(propbn1_tas_10d_dif$Data, "climatology:fun") <- "with - without evidence"
plotClimatology(propbn1_tas_10d_without,backdrop.theme = "coastline")
plotClimatology(propbn1_tas_10d_with,backdrop.theme = "coastline")
plotClimatology(propbn1_tas_10d_dif,backdrop.theme = "coastline")

str(propagationbn2)
propagationbn2$without
propagationbn2$with
matwithout2 <- matrix(propagationbn2$without, nrow = 1)
matwith2 <- matrix(propagationbn2$with, nrow = 1)
matdif2 <- matwith2 - matwithout2
propbn2_tas_10d_without <- tas_ncep_10d
propbn2_tas_10d_with <- tas_ncep_10d
propbn2_tas_10d_dif <- tas_ncep_10d
propbn2_tas_10d_without$Data <- mat2Dto3Darray(matwithout2, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
propbn2_tas_10d_with$Data <- mat2Dto3Darray(matwith2, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
propbn2_tas_10d_dif$Data <- mat2Dto3Darray(matdif2, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
attr(propbn2_tas_10d_without$Data, "climatology:fun") <- "without evidence"
attr(propbn2_tas_10d_with$Data, "climatology:fun") <- "with evidence"
attr(propbn2_tas_10d_dif$Data, "climatology:fun") <- "with - without evidence"
plotClimatology(propbn2_tas_10d_without,backdrop.theme = "coastline")
plotClimatology(propbn2_tas_10d_with,backdrop.theme = "coastline")
plotClimatology(propbn2_tas_10d_dif,backdrop.theme = "coastline")

str(propagationbn3)
propagationbn3$without
propagationbn3$with

matwithout3 <- matrix(propagationbn3$without, nrow = 1)
matwith3 <- matrix(propagationbn3$with, nrow = 1)
matdif3 <- matwith3 - matwithout3

propbn3_tas_10d_without <- tas_ncep_10d
propbn3_tas_10d_with <- tas_ncep_10d
propbn3_tas_10d_dif <- tas_ncep_10d
propbn3_tas_10d_without$Data <- mat2Dto3Darray(matwithout3, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
propbn3_tas_10d_with$Data <- mat2Dto3Darray(matwith3, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
propbn3_tas_10d_dif$Data <- mat2Dto3Darray(matdif3, x = attr(data.dag, "Xcoords"), y = attr(data.dag, "Ycoords"))
attr(propbn3_tas_10d_without$Data, "climatology:fun") <- "without evidence"
attr(propbn3_tas_10d_with$Data, "climatology:fun") <- "with evidence"
attr(propbn3_tas_10d_dif$Data, "climatology:fun") <- "with - without evidence"

plotClimatology(propbn3_tas_10d_without,backdrop.theme = "coastline")
plotClimatology(propbn3_tas_10d_with, backdrop.theme = "coastline")
plotClimatology(propbn3_tas_10d_dif,backdrop.theme = "coastline")


##############################################################################################
# Tryout three variables 
# Conditional independence ?
# P(x,y | z) = P(x | y,z) ?
##############################################################################################

cpquery(fitted3, event = (V286 >= 1) & (V458 >= 1), evidence = (V81 >= 1))
cpquery(fitted3, event = (V286 >= 1), evidence = ((V458 >= 1) & (V81 >= 1)))

cpquery(fitted3, event = (V286 >= 1) & (V315 >= 1), evidence = (V81 >= 1))
cpquery(fitted3, event = (V286 >= 1), evidence = ((V315 >= 1) & (V81 >= 1)))

cpquery(fitted3, event = (V286 <= -1) & (V458 >= 1), evidence = (V81 >= 1))
cpquery(fitted3, event = (V286 <= -1), evidence = ((V458 >= 1) & (V81 >= 1)))

cpquery(fitted3, event = (V286 <= -1) & (V315 >= 1), evidence = (V81 >= 1))
cpquery(fitted3, event = (V286 <= -1), evidence = ((V315 >= 1) & (V81 >= 1)))

cpquery(fitted3, event = (V286 <= -1), evidence = (V205>= 1))
cpquery(fitted3, event = (V286 >= 1), evidence = (V205>= 1))
cpquery(fitted3, event = (V205 >= 1), evidence = (V286 <= -1)) #BIG
cpquery(fitted3, event = (V205 <= -1), evidence = (V286 <= -1))
cpquery(fitted3, event = (V205 <= -1), evidence = (V286 <= -1))
cpquery(fitted3, event = (V205 >= 1), evidence = (V286 <= -1) & (V81 >=1))
cpquery(fitted3, event = (V205 >= 1) & (V286 <= -1) , evidence = (V81 >=1))
cpquery(fitted3, event = (V205 >= 1) , evidence = (V81 >=1))
cpquery(fitted3, event = (V205 >= 1) , evidence = TRUE)
cpquery(fitted3, event = (V205 >= 1) , evidence = (V81 >=1))

dsep(fitted3, 'V49', 'V42', c('V81','V99','V63','V45','V27'))

cpquery(fitted3, event = (V423>=1), evidence = (V81>= 1))
cpquery(fitted3, event = (V423>=1), evidence = ((V81>= 1) & (V63 <= 0.5) & (V63 >= 0)))
