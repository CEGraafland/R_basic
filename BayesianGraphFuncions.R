# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Makes graph with iamb of data with treshold hyp.alpha for correlations
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
graphmaker <- function(data, data.prop = 1, hyp.alpha = 0.05){
  
  if(!is.null(names(data))){nodenames <- names(data)} 
  
  prop <- data.prop
  nobs <- base::nrow(data)
  
  train <- as.data.frame(matrix(data = NA, nrow = prop*nobs, ncol = ncol(data)))
  names(train) <- nodenames
  for(i in 1:ncol(train)){
    train[1:(prop*nobs),i] <- data[1:(prop*nobs),i]
  }
  
  train.pdag <- iamb(train, test = "cor", alpha = hyp.alpha)
  train.dag <- cextend(train.pdag)
  
  return(train.dag)
  
}

# Compares the loglikelikhoods of the train and testdata
compare <- function(graph, data, train.prop = 0.5){
  if(!is.null(names(data))){nodenames <- names(data)
  } else {nodenames <- nodes(graph)}
  
  prop <- train.prop
  nobs <- base::nrow(data)
  
  train <- as.data.frame(matrix(data = NA, nrow = prop*nobs, ncol = length(nodenames)))
  names(train) <-nodenames
  
  for(i in names(train)){
    train[1:(prop*nobs),i] <- data[1:(prop*nobs),i]
  }
  
  testdata <- as.data.frame(matrix(data = NA, nrow = ((1-prop)*nobs), ncol = length(nodenames)))
  names(testdata)<- nodenames
  
  for(i in names(testdata)){
    testdata[1:((1-prop)*nobs),i] <- data[(prop*nobs+1):(nobs),i]
  }
  
  train.gbn <- bn.fit(graph, train)
  
  liktrain <- logLik(train.gbn, train, nodenames)
  liktestdata <- logLik(train.gbn, testdata, nodenames)
  
  #scoretrain <- score(graph, train)
  #scoretest <- score(train.gbn, testdata)
  
  print(paste0("Loglikelihood train is: ",liktrain)) 
  print(paste0("Loglikelihood test is: ",liktestdata))
  
  return(c(liktrain,liktestdata))
  
}



# alpha_loglik_pc2 <- function(data, alpha){
#  degreestat <- list(C = data$correlation, n = nrow(data$data_coords))
#  eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data$data_coords), alpha = alpha)
#  amata <- as(eqclass,"amat")
  
#  if(isValidGraph(amata, type = "pdag")){
#    dag <- pcalg::pdag2dag(eqclass@graph)
#    bn.dag <- as.bn(dag$graph)
#    newdata <- as.data.frame(data$data_coords)
#    names(newdata) <- nodes(bn.dag)
#    gbn <- bn.fit(bn.dag, newdata)
#    likgbn <- logLik(gbn, newdata)
#  }
#  else likgbn <- 0
#  
#  return(likgbn)
#}


#mypc2 <- function(data, alpha){
#  datastat <- list(C = data$correlation, n = nrow(data$data_coords))
#  eqclas<- pc(datastat, indepTest = gaussCItest, p = ncol(data$data_coords), alpha = alpha)
#  amata <- as(eqclas,"amat")
#  vcpdag <- isValidGraph(amata, type = "cpdag")
#  vpdag <- isValidGraph(amata, type = "pdag")
  
#  bneqclas <- as.bn(eqclas@graph)
#  nund <- length(undirected.arcs(bneqclas))
#  nedges <- narcs(bneqclas)
#  return(data.frame(nund, nedges, vcpdag, vpdag))
#}
################## cbind imp dataframe! NOT CHECKED YET + /2 not checked yet
alpha_loglik_pc3_simple <- function(data, correlation = NULL, method = c("spearman","kendall","pearson"), alpha){
  if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  degreestat <- list(C = correlation, n = nrow(data))
  eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
  amata <- as(eqclass,"amat")
  bn.eqclass <- as.bn(eqclass@graph)
  nedgeseqclass <- narcs(bn.eqclass)
  nundeqclass <- nrow(undirected.arcs(bn.eqclass))/2
  
  if(isValidGraph(amata, type = "pdag")){
    dag <- pcalg::pdag2dag(eqclass@graph)
    bn.dag <- as.bn(dag$graph)
    nedgesdag <- narcs(bn.dag)
    
    newdata <- as.data.frame(data)
    names(newdata) <- nodes(bn.dag)
    gbn <- bn.fit(bn.dag, newdata)
    likgbn <- logLik(gbn, newdata)
  }
  else {likgbn <- NA
  nedgesdag <- NA}
  
  
  return(cbind(alpha,
                    nundeqclass,
                    nedgeseqclass,
                    nedgesdag,
                    likgbn))
}

sapply(myalphastot, alpha_loglik_pc3_simple, data = gaussian.test, correlation = NULL, method = "spearman")

data10d <- TimeCoordsAnom_from_Grid(tas_ncep_10d)
corrspear10d <- cor(data10d, method = "spearman")
loglik10d_19_17 <- sapply(myalphas9, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_17_15 <- sapply(myalphas8, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")
loglik10d_15_13 <- sapply(myalphas7, alpha_loglik_pc3_simple, data = data10d, correlation = corrspear10d, method ="spearman")


alpha_loglik_pc3 <- function(data, alpha){
 
  degreestat <- list(C = data$correlation, n = nrow(data$data_coords))
  eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data$data_coords), alpha = alpha)
  amata <- as(eqclass,"amat")
  bn.eqclass <- as.bn(eqclass@graph)
  nedgeseqclass <- narcs(bn.eqclass)
  nundeqclass <- nrow(undirected.arcs(bn.eqclass))/2
  
  if(isValidGraph(amata, type = "pdag")){
    dag <- pcalg::pdag2dag(eqclass@graph)
    bn.dag <- as.bn(dag$graph)
    nedgesdag <- narcs(bn.dag)
    
    newdata <- as.data.frame(data$data_coords)
    names(newdata) <- nodes(bn.dag)
    gbn <- bn.fit(bn.dag, newdata)
    likgbn <- logLik(gbn, newdata)
  }
  else {likgbn <- NA
        nedgesdag <- NA}
  
  
  return(data.frame(alpha,
           nundeqclass,
           nedgeseqclass,
           nedgesdag,
           likgbn))
}

# Makes eqclass with PC algorithm
# Convert eqclass in DAG IF possible
# returns quantities as loglik of resulting dag
# returns resulting eqclass and dag
# !!! Right now works only for extendable alphas!!!!
alpha_loglik_pc3_withdag <- function(data, correlation = NULL, method = "spearman", alpha){
  if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  degreestat <- list(C = correlation, n = nrow(data))
  eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data), alpha = alpha)
  amata <- as(eqclass,"amat")
  bn.eqclass <- as.bn(eqclass@graph,check.cycles = FALSE)
  nedgeseqclass <- narcs(bn.eqclass)
  nundeqclass <- nrow(undirected.arcs(bn.eqclass))/2
  
  if(isValidGraph(amata, type = "pdag")){
    dag <- pcalg::pdag2dag(eqclass@graph)
    bn.dag <- as.bn(dag$graph)
    nedgesdag <- narcs(bn.dag)
    
    newdata <- as.data.frame(data)
    names(newdata) <- nodes(bn.dag)
    gbn <- bn.fit(bn.dag, newdata)
    likgbn <- logLik(gbn, newdata)
  }
  else {likgbn <- NA
  nedgesdag <- NA}
  
  
  return(list(networkdata = data.frame(alpha,
                    nundeqclass,
                    nedgeseqclass,
                    nedgesdag,
                    likgbn), networks = list(eqclass = bn.eqclass, dag = bn.dag)))
}#alpha_loglik_pc3_withdag


alpha_loglik_pc3_withdag(degree10d, alpha = 4e-108)

myalphas0 <- as.vector(seq(0.0001, 0.001, 0.0001))
myalphas1 <- as.vector(seq(0.001, 0.01, 0.001))
myalphas2 <- as.vector(seq(0.01, 0.1, 0.01))
myalphastot <- as.vector(seq(0.001,0.1,0.001))
myalphas3 <- as.vector(seq(0.00001,0.001,0.00001))
myalphas4 <- as.vector(seq(1e-7, 1e-5, 1e-7))
myalphas4b <- as.vector(seq(1e-9, 1e-7, 1e-9))
myalphas5 <- as.vector(seq(1e-11,1e-9,1e-11))
myalphas6 <- as.vector(seq(1e-13,1e-11,1e-13))
myalphas7 <- as.vector(seq(1e-15,1e-13,1e-15))
myalphas8 <- as.vector(seq(1e-17,1e-15,1e-17))
myalphas9 <- as.vector(seq(1e-19,1e-17,1e-19))

mapply(alpha_loglik_pc3, list(degree30d), myalphas1)
loglik9045d <- mapply(alpha_loglik_pc3, list(degree9045d), myalphastot)
loglik30  <- mapply(alpha_loglik_pc3, list(degree30d), myalphastot) #hamitad
loglik20  <- mapply(alpha_loglik_pc3, list(degree20d), myalphastot) #sale nada
loglik20b <- mapply(alpha_loglik_pc3, list(degree20d), myalphas3) # hamidad
loglik20c <- mapply(alpha_loglik_pc3, list(degree20d), myalphas4) # toda
loglik20d <- mapply(alpha_loglik_pc3, list(degree20d), myalphas4b) # todas
loglik20e <- mapply(alpha_loglik_pc3, list(degree20d), myalphas5)
loglik10  <- mapply(alpha_loglik_pc3, list(degree10d), myalphas3) # nada
loglik10b <- mapply(alpha_loglik_pc3, list(degree10d), myalphas4) # nada
loglik10c <- mapply(alpha_loglik_pc3, list(degree10d), myalphas5) # nada
loglik10d <- mapply(alpha_loglik_pc3, list(degree10d), myalphas6) # Solo uno
loglik10e <- mapply(alpha_loglik_pc3, list(degree10d), myalphas7) # Todas
loglik10f <- mapply(alpha_loglik_pc3, list(degree10d), myalphas8) # Todas
loglik10g <- mapply(alpha_loglik_pc3, list(degree10d), myalphas9) # Todas

loglik10d

degree10_alpha19_13 <-cbind(loglik10g,loglik10f, loglik10e)
degree10_alpha17_13 <-cbind(loglik10f, loglik10e)
degree20_alpha11_3 <- cbind(loglik20e, loglik20d,loglik20c, loglik20b)
degree20_alpha11_3[,1]

testpc3degree10 <- alpha_loglik_pc3(degree10d, alpha = 1e-13) #yes
testpc3degree10a <- alpha_loglik_pc3(degree10d, alpha = 2e-13) #no
testpc3degree10a

testgaussianpc3 <- alpha_loglik_pc3()

#plot nedges vs loglik for the pc algorithm 
# pc3obj = object from function alpha_loglik_pc3
plot.nedge.vs.loglik_pc <- function(pc3obj) {
  xeje <- unlist(pc3obj['nedgeseqclass',])
  yeje <- unlist(pc3obj['likgbn',])
  
  miss <- !is.na(yeje)
  
  xeje[which(miss)]
  graph <- plot(xeje[which(miss)], yeje[which(miss)], 
                main = "loglikelihood vs number of edges",
                sub = paste("PC algorithm for ", deparse(substitute(pc3obj)), sep=" "), 
                xlab = "number of edges in DAG", 
                ylab = "loglik()" )
  
}#plot.nedge.vs.loglik_pc



plot.nedge.vs.loglik_pc(loglik9045d)
plot.nedge.vs.loglik_pc(loglik20b)
plot.nedge.vs.loglik_pc(degree20_alpha11_3)
plot.nedge.vs.loglik_pc(loglik30) #niet opgeslagen nog doen
plot.nedge.vs.loglik_pc(loglik10e)
plot.nedge.vs.loglik_pc(loglik10f)
plot.nedge.vs.loglik_pc(degree10_alpha17_13)
plot.nedge.vs.loglik_pc(degree10_alpha19_13)

# Combi plot HC AND PC
plot.nedge.vs.loglik_combi <- function(pc3obj,it_hc_obj,xlim = NULL) {
  xejepc <- unlist(pc3obj['nedgeseqclass',])
  yejepc <- unlist(pc3obj['likgbn',])
  
  miss <- !is.na(yejepc)
  
  xejehc <- unlist(it_hc_obj['edges'])
  yejehc <- unlist(it_hc_obj['logliks'])
  
  plot(xejehc, yejehc, 
       main = "loglikelihood vs number of edges",
       sub = paste("HC and PC algorithm for ", deparse(substitute(it_hc_obj)), sep=" "), 
       xlab = "number of edges in DAG", 
       ylab = "log P(d|(DAG,P))",
       xlim = xlim)
  points(xejepc[which(miss)], yejepc[which(miss)])
  
}#plot.nedge.vs.loglik_combi

pcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/pcalg/pc10d_a19_13_hc10d_400_600i.pdf"
pdf(pcplotname)
plot.nedge.vs.loglik_combi(degree10_alpha19_13,
                           hc_edges_loglik_10d_400_600i$networkdata)
dev.off()


plot.nedge.vs.loglik_CNvsBN <- function(pc3obj,it_hc_obj,xlim = NULL) {
  xejepc <- unlist(pc3obj['nedgeseqclass',])
  yejepc <- unlist(pc3obj['likgbn',])
  
  miss <- !is.na(yejepc)
  
  xejehc <- unlist(it_hc_obj['edges'])
  yejehc <- unlist(it_hc_obj['logliks'])
  
  plot(xejehc, yejehc, 
       main = paste0("Information vs number of edges\nHC and PC algorithm until ", deparse(substitute(it_hc_obj)),"\nComplex networks until ",xlim[2]), 
       xlab = "number of edges in complex graph G and Bayesian graph DAG", 
       ylab = "log P(d|(DAG,P))",
       xlim = xlim)
  points(xejepc[which(miss)], yejepc[which(miss)])
  
}#plot.nedge.vs.loglik_CNvsBN

###opslaan -----------
#opslaan plot 10d met alphas 10^-19 tot en met 10 ^-13
pcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/pcalg/pc10d_a19_13.pdf"
pdf(pcplotname)
plot.nedge.vs.loglik_pc(degree10_alpha19_13)
dev.off()

#maken netwerkplaatjes voor PC 10d alpha = c(3.6e-18,1e-15,1e-13). VERBETEREN! NU elke keer uitrekenen
for(a in c(3.6e-18,1e-15,1e-13)){
  pc10da <- alpha_loglik_pc3_withdag(degree10d,a) #pc10da3.6_18$networkdata$nedgesdag #class(pc10da3.6_18$networks$dag)

  hcplotname <- paste0("/Users/lisettegraafland/Documents/R_practice/plots/pcalg/network10d",as.numeric(pc10da$networkdata["alpha"]),".pdf")  
  pdf(hcplotname)
  plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d),pc10da$networks$dag)
  title(paste0("number of edges = ",pc10da$networkdata$nedgesdag,"\nalpha = ",pc10da$networkdata$alpha), line = -4)
  dev.off()
}




#verwijdering zware grafieken.
gridcomplexgraphs <- list(a,b,c,d,e,f)
rm(a,b,c,d,e,f)
testbayesiangraphs <- list(test.cpdag, test.cpdagNEL, test.dag.igraph,test.iamb,test.dag, test.dagNEL, struDegree20, struDegree20un, struDegree30a,struDegree30mmhc, struDegree10,struDegree20test,struDegree30hc, struDegree9045, stru.hc.Degree9045igraph, stru.hc.Degree20, stru.fast.iamb.Degree20, stru.gs.Degree20, stru.hc.Degree9045, stru.hc.Degree9045igraph,iGraphDegree20.dir, iGraphDegree20.simp, igraphDegree9045, iGraphDegree20, iGraphDegree20.pdag, iGraphDegree20.sp)
save(test.cpdag, 
     test.cpdagNEL, 
     test.dag.igraph,
     test.iamb,
     test.dag, 
     test.dagNEL, 
     struDegree20, 
     struDegree20un, 
     struDegree30a,
     struDegree30mmhc, 
     struDegree10,
     struDegree20test,
     struDegree30hc, 
     struDegree9045,
     stru.hc.Degree9045igraph, 
     stru.hc.Degree20, 
     stru.fast.iamb.Degree20, 
     stru.gs.Degree20,
     stru.hc.Degree9045, 
     stru.hc.Degree9045igraph,
     iGraphDegree20.dir, 
     iGraphDegree20.simp, 
     igraphDegree9045, 
     iGraphDegree20, 
     iGraphDegree20.pdag, 
     iGraphDegree20.sp,
     file="/Users/lisettegraafland/Documents/R_practice/Data/data_testbayesiangraphs.rda")
load("/Users/lisettegraafland/Documents/R_practice/Data/data_testbayesiangraphs.rda")
save(a,b,c,d,e,f,file="/Users/lisettegraafland/Documents/R_practice/Data/data_gridcomplexgraphs.rda")

rm(test.cpdag, 
     test.cpdagNEL, 
     test.dag.igraph,
     test.iamb,
     test.dag, 
     test.dagNEL, 
     struDegree20, 
     struDegree20un, 
     struDegree30a,
     struDegree30mmhc, 
     struDegree10,
     struDegree20test,
     struDegree30hc, 
     struDegree9045,
     stru.hc.Degree9045igraph, 
     stru.hc.Degree20, 
     stru.fast.iamb.Degree20, 
     stru.gs.Degree20,
     stru.hc.Degree9045, 
     iGraphDegree20.dir, 
     iGraphDegree20.simp, 
     igraphDegree9045, 
     iGraphDegree20, 
     iGraphDegree20.pdag, 
     iGraphDegree20.sp)
warnings()
