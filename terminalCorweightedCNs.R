############################################################################################
# Creation of Correlation weighted Correlation networks (distance /weights / largeweights)
############################################################################################
rm(list = ls())
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/propagationFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")

load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
gridGraphsCN <- gridGraphs
rm(gridGraphs)
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(gridGraphsCN) <- as.character(numberofedgesCN)
listcn <- GraphsCN[1:100]
narcs <- numberofedgesCN[1:100]

listcnstrengths <- lapply(listcn, str_vs_dis_2, perm = NULL, data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
listcnstrengths <- listcnstrenghts
save(listcnstrengths, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNstrengths.rda")
listcnweights <- lapply(listcnstrengths, igraph.weights, type = "cn", fromdist = 0, perm = NULL)
save(listcnweights, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphsCNcorweights.rda")
listcnweights[[5]]


####################################################################################################
# Creation of BIC -  weighted Bayesian networks for all permutations (distance / weights / largeweights)
####################################################################################################
rm(list = ls())
library(igraph)
library(bnlearn)
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations20.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations20.rda")
################################################
# Load BNs
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm2sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm3sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm4sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm5sort.rda")
permsorts <- list(perm1sort,perm2sort,perm3sort,perm4sort,perm5sort)

for( i in 3:5){
permstrenghts <- lapply(permsorts[[i]], bn_to_igraph.strengths, perm = permutations[[i]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
permdists <- lapply(permstrenghts, igraph.distances, perm = permutations[[i]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 2000)
assign(paste0("perm",i,"weights"), permweights)
save(list = paste0("perm",i,"weights"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",i,"weights.rda"))
}

################################################
# Load BNs 6 tot en met 10
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm6sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm7sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm8sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm9sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm10sort.rda")
permsorts_6_10 <- list(perm6sort,perm7sort,perm8sort,perm9sort,perm10sort)
# first posibility to calculate
for( i in 1:5){
  permnumber <- i + 5
  permstrenghts <- lapply(permsorts_6_10[[i]], bn_to_igraph.strengths, perm = permutations20[[i]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
  permdists <- lapply(permstrenghts, igraph.distances, perm = permutations20[[i]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
  permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 2000)
  assign(paste0("perm",permnumber,"weights"), permweights)
  save(list = paste0("perm",permnumber,"weights"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",permnumber,"weights.rda"))
}
# or same but other written:
for( i in 6:10){
  permnumber <- i 
  listnumber <- i -5
  perm20number <- i-5
  permstrenghts <- lapply(permsorts_6_10[[listnumber]], bn_to_igraph.strengths, perm = permutations20[[perm20number]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
  permdists <- lapply(permstrenghts, igraph.distances, perm = permutations20[[perm20number]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
  permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 2000)
  assign(paste0("perm",permnumber,"weights"), permweights)
  save(list = paste0("perm",permnumber,"weights"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",permnumber,"weights.rda"))
}
################################################
# Load BNs 11 tot en met 25
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm11sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm12sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm13sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm14sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm15sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm16sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm17sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm18sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm19sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm20sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm21sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm22sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm23sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm24sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm25sort.rda")

permsorts_11_25 <- list(perm11sort,perm12sort,perm13sort,perm14sort,perm15sort,
                        perm16sort,perm17sort,perm18sort,perm19sort,perm20sort,
                        perm21sort,perm22sort,perm23sort,perm24sort,perm25sort)
for( i in 11:25){
  permnumber <- i
  listnumber <- i-10
  perm20number <- i-5
  permstrenghts <- lapply(permsorts_11_25[[listnumber]], bn_to_igraph.strengths, perm = permutations20[[perm20number]],data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
  permdists <- lapply(permstrenghts, igraph.distances, perm = permutations20[[perm20number]] ,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
  permweights <- lapply(permdists, igraph.weights, type = "bn", fromdist = 2000)
  assign(paste0("perm",permnumber,"weights"), permweights)
  save(list = paste0("perm",permnumber,"weights"), file = paste0("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm",permnumber,"weights.rda"))
}