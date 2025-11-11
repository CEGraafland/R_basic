##################################################################################################
# Iteration function hill climbing 
# Include parameters
##################################################################################################
# Applies hill climbing algorithm step by step.
# entries:  data, number of iterations of de algorithm, 
#           network from which next best will be produced.
# Outcome:  nedges in intermediate networks, scores of intermediate networks, final network

iteration_loglik_hc <- function(data, iterations, start = NULL){
  logliks = list()
  scores = list()
  edges = list()
  networks = NULL
  
  step = hc(data, max.iter = 1, start = start)
  
  edges[1] <- narcs(step)
  scores[1] <- score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = hc(data, max.iter = 1, start = step)
    nedges.stepi <- narcs(stepi)
    edges[i+1] <- nedges.stepi
    score.stepi <- score(stepi, data)
    scores[i+1] <- score.stepi
    loglik.stepi <- logLik(stepi, data)
    logliks[i+1] <- loglik.stepi
    
    if (!unlist(scores[i+1])>unlist(scores[i])) {
      print(paste0("is broken at", i+1))
      break
    }
    else{step <-stepi
    i <- i+1}
  }
  
  networks <- stepi
  return(list(networkdata = data.frame(edges = unlist(edges),logliks = unlist(logliks)), networks = networks))
}
###########################################################################################
# Iteration with number of parameters. 
###########################################################################################
iteration_hc_old <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- NULL
  
  step = hc(data, max.iter = 1, start = start)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  
  i<-1
  while (i < iterations){
    stepi = hc(data, max.iter = 1, start = step)
    nedges.stepi <- narcs(stepi)
    edges[i+1] <- nedges.stepi
    score.stepi <- bnlearn::score(stepi, data)
    scores[i+1] <- score.stepi
    loglik.stepi <- logLik(stepi, data)
    logliks[i+1] <- loglik.stepi
    params.stepi <- nparams(stepi, data)
    params[i+1] <- params.stepi
    tests.stepi <- ntests(stepi)
    tests[i+1] <- tests.stepi
    
    if (!scores[i+1]>scores[i]) {
      print(paste0("is broken at", i+1))
      break
    }
    else{step <-stepi
    i <- i+1}
  }
  
  networks <- stepi
  output <- list(networkdata = data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests), networks = networks)
  return(output)
}

###############################################################################################
# iets verzinnen voor list append!
# IDEE: modificar de ratios
###############################################################################################
iteration_tabu <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  step = tabu(data, score = "bic-g", max.iter = 1, start = start)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = tabu(data, score = "bic-g", max.iter = 1, start = step)
    nedges.stepi <- narcs(stepi)
    edges[i+1] <- nedges.stepi
    score.stepi <- bnlearn::score(stepi, data)
    scores[i+1] <- score.stepi
    loglik.stepi <- logLik(stepi, data)
    logliks[i+1] <- loglik.stepi
    params.stepi <- nparams(stepi, data)
    params[i+1] <- params.stepi
    tests.stepi <- ntests(stepi)
    tests[i+1] <- tests.stepi
    
    ###########################
    if (samplesize/params[i] >= 0.5& samplesize/params[i + 1]< 0.5) {
    networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
    attr(stepi, "n/par") <- (samplesize/params[i + 1])
    networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    }
    
    
    if (!scores[i+1] > scores[i]) {
      message("is broken at", i+1)
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
      break
    } else {
      step <-stepi
      i <- i+1
    }
  }
  
  if (scores[i] > scores[i-1]) {
  networks.names[i] <- paste0("tabu_10d_",i,"i")
  attr(stepi, "n/par") <- (samplesize/params[i])
  networks[[i]] <- stepi
  }
  
  ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
  networks <- networks[ind0]
  names(networks) <- networks.names[ind0]

  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

iteration_hc <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  step = hc(data, score = "bic-g", max.iter = 1, start = start)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = hc(data, score = "bic-g", max.iter = 1, start = step)
    nedges.stepi <- narcs(stepi)
    edges[i+1] <- nedges.stepi
    score.stepi <- bnlearn::score(stepi, data)
    scores[i+1] <- score.stepi
    loglik.stepi <- logLik(stepi, data)
    logliks[i+1] <- loglik.stepi
    params.stepi <- nparams(stepi, data)
    params[i+1] <- params.stepi
    tests.stepi <- ntests(stepi)
    tests[i+1] <- tests.stepi
    
    ###########################
    if (samplesize/params[i] >= 0.5& samplesize/params[i + 1]< 0.5) {
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    }
    
    
    if (!scores[i+1] > scores[i]) {
      message("is broken at", i+1)
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
      break
    } else {
      step <-stepi
      i <- i+1
    }
  }
  
  if (scores[i] > scores[i-1]) {
    networks.names[i] <- paste0("tabu_10d_",i,"i")
    attr(stepi, "n/par") <- (samplesize/params[i])
    networks[[i]] <- stepi
  }
  
  ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
  networks <- networks[ind0]
  names(networks) <- networks.names[ind0]
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

iteration_mmhc <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  step = mmhc(data, maximize.args = list(max.iter = 1))
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = mmhc(data, restrict.args = list(test = "bic-gt"), maximize.args = list(max.iter = i+1))
    nedges.stepi <- narcs(stepi)
    edges[i+1] <- nedges.stepi
    score.stepi <- bnlearn::score(stepi, data)
    scores[i+1] <- score.stepi
    loglik.stepi <- logLik(stepi, data)
    logliks[i+1] <- loglik.stepi
    params.stepi <- nparams(stepi, data)
    params[i+1] <- params.stepi
    tests.stepi <- ntests(stepi)
    tests[i+1] <- tests.stepi
    
    ###########################
    if (samplesize/params[i] >= 0.5& samplesize/params[i + 1]< 0.5) {
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    }
    
    if (!scores[i+1] > scores[i]) {
      message("is broken at", i+1)
      networks.names[i] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i+1])
      networks[[i+1]] <- stepi
      break
    } else {
      step <-stepi
      i <- i+1
    }
  }
  
  
  
  
  if (scores[i] > scores[i-1]) {
    networks.names[i] <- paste0("mmhc_10d_",i,"i")
    attr(stepi, "n/par") <- (samplesize/params[i])
    networks[[i]] <- stepi
  }
  
  ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
  networks <- networks[ind0]
  names(networks) <- networks.names[ind0]
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

iteration_mmhc_split <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  begin = mmpc(data, test = "bic-gt")
  
  fulligraph <- make_full_graph(ncol(data))
  fullbnlearn <- as.bn(igraph.to.graphNEL(fulligraph))
  nodes(fullbnlearn) <- names(data)
  
  igraph1 <- igraph.from.graphNEL(as.graphNEL(fullbnlearn))
  igraph2 <- igraph.from.graphNEL(as.graphNEL(begin))
  
  blackigraph <- difference(igraph1, igraph2)
  blackbngraph <- as.bn(igraph.to.graphNEL(blackigraph))
  
  step <- hc(data, score = "bic-g", blacklist = arcs(blackbngraph), max.iter = 1)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi <- hc(data, score = "bic-g", blacklist = arcs(blackbngraph), max.iter = 1, start = step)
    nedges.stepi <- narcs(stepi)
    edges[i+1] <- nedges.stepi
    score.stepi <- bnlearn::score(stepi, data)
    scores[i+1] <- score.stepi
    loglik.stepi <- logLik(stepi, data)
    logliks[i+1] <- loglik.stepi
    params.stepi <- nparams(stepi, data)
    params[i+1] <- params.stepi
    tests.stepi <- ntests(stepi)
    tests[i+1] <- tests.stepi
    
    ###########################
    if (samplesize/params[i] >= 0.5& samplesize/params[i + 1]< 0.5) {
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    }
    
    if (!scores[i+1] > scores[i]) {
      message("is broken at", i+1)
      networks.names[i] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i+1])
      networks[[i+1]] <- stepi
      break
    } else {
      step <-stepi
      i <- i+1
    }
  }
  
  
  
  
  if (scores[i] > scores[i-1]) {
    networks.names[i] <- paste0("mmhc_10d_",i,"i")
    attr(stepi, "n/par") <- (samplesize/params[i])
    networks[[i]] <- stepi
  }
  
  ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
  networks <- networks[ind0]
  names(networks) <- networks.names[ind0]
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata, begin = begin)
  return(output)
}


alpha_loglik_pc3_withdag <- function(data, correlation = NULL, method = "spearman", alpha){
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
  
  
  return(list(networkdata = data.frame(alpha,
                                       nundeqclass,
                                       nedgeseqclass,
                                       nedgesdag,
                                       likgbn), networks = list(eqclass = bn.eqclass, dag = bn.dag)))
}#alpha_loglik_pc3_withdag


iteration_pc <- function(data, correlation = NULL, method = "spearman", alpha){
  if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  degreestat <- list(C = correlation, n = nrow(data))
  
  names <- c()
  for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  alphas <- c()
  nedgeseqclass <- c()
  nundeqclass <- c()
  edges <- c()
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  i <- 1
  while (i <= length(alpha)){
    eqclass <- pc(degreestat, indepTest = gaussCItest, p = ncol(data), alpha = alpha[i])
    tests[i] <- sum(eqclass@n.edgetests)
    amata <- as(eqclass,"amat")
    bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)
    nedgeseqclass[i] <- narcs(bn.eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(bn.eqclass))/2
  
    if(isValidGraph(amata, type = "pdag")){
      dag <- pcalg::pdag2dag(eqclass@graph)
      bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(bn.dag)
    
      # data maken variable names!
      newdata <- as.data.frame(data)
      names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(bn.dag, newdata)
      
      alphas[i] <- alpha[i]
      logliks[i] <- logLik(gbn, newdata)
      params[i] <- nparams(gbn)
      scores[i] <- score(bn.dag, newdata)
    } else {alpha[i] <- NA
      logliks[i] <- NA
      edges[i] <- NA
      params[i] <- NA
      scores[i] <- NA}
    
    
  
    ###########################
    if (!i == 1){
      if(!is.na(params[i-1]) & !is.na(params[i])) {
        if (samplesize/params[i-1] >= 0.5 & samplesize/params[i] <= 0.5) { 
          if ((samplesize/params[i-1] - 0.5) < (0.5 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.5")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.4 & samplesize/params[i] <= 0.4) { 
          if ((samplesize/params[i-1] - 0.4) < (0.4 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.4")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.3 & samplesize/params[i] <= 0.3) { 
          if ((samplesize/params[i-1] - 0.3) < (0.3 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.3")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.2 & samplesize/params[i] <= 0.2) { 
          if ((samplesize/params[i-1] - 0.2) < (0.2 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.2")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.1 & samplesize/params[i] <= 0.1) { 
          if ((samplesize/params[i-1] - 0.1) < (0.1 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.1")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.05 & samplesize/params[i] <= 0.05) { 
          if ((samplesize/params[i-1] - 0.05) < (0.05 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.05")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        }
      }
    }
    
    gbnmin1 <- gbn
    i <- i +1
    
  }
  if(!length(networks) == 0){
  ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
  networks <- networks[ind0]
  names(networks) <- networks.names[ind0]}
  
  
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

data <- gaussian.test
alpha <- c(1,5,10,50)

iteration_pcstable <- function(data, alpha){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  alphas <- c()
  nedgeseqclass <- c()
  nundeqclass <- c()
  edges <- c()
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  i <- 1
  while (i <= length(alpha)){
    eqclass <- pc.stable(data, test = "bic-gt", B = alpha[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)

    if(acyclic(eqclass, directed = TRUE)){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      alphas[i] <- alpha[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
    } else {alpha[i] <- NA
    logliks[i] <- NA
    edges[i] <- NA
    params[i] <- NA
    scores[i] <- NA}
    
    
    
    ###########################
    if (!i == 1){
      if(!is.na(params[i-1]) & !is.na(params[i])) {
        if (samplesize/params[i-1] >= 0.5 & samplesize/params[i] <= 0.5) { 
          if ((samplesize/params[i-1] - 0.5) < (0.5 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.5")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.4 & samplesize/params[i] <= 0.4) { 
          if ((samplesize/params[i-1] - 0.4) < (0.4 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.4")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.3 & samplesize/params[i] <= 0.3) { 
          if ((samplesize/params[i-1] - 0.3) < (0.3 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.3")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.2 & samplesize/params[i] <= 0.2) { 
          if ((samplesize/params[i-1] - 0.2) < (0.2 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.2")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.1 & samplesize/params[i] <= 0.1) { 
          if ((samplesize/params[i-1] - 0.1) < (0.1 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.1")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.05 & samplesize/params[i] <= 0.05) { 
          if ((samplesize/params[i-1] - 0.05) < (0.05 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pc_10d_0.05")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "alpha") <- (alpha[l])
          networks[[l]] <- add
        }
      }
    }
    
    gbnmin1 <- gbn
    i <- i +1
    
  }
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

output

##################################################################################################
# Visualize loglik vs nedges.
# Parameters ??
##################################################################################################
# Plots number of edges of dag with respect to loglik dag
# Entries:  outcome $networkdata of 
#           iteration_loglik_hc: number of edges and logliks.
# Outcome:  Plot

plot.nedge.vs.loglik_hc <- function(it_hc_obj) {
  xeje <- unlist(it_hc_obj['edges'])
  yeje <- unlist(it_hc_obj['logliks'])
  
  graph <- plot(xeje, yeje, 
                main = "loglikelihood vs number of edges",
                sub = paste("HC algorithm for ", deparse(substitute(it_hc_obj)), sep=" "), 
                xlab = "number of edges in DAG", 
                ylab = "loglik()" )
}
