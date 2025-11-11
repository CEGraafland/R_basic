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

iteration_tabu <- function(data, iterations, gamma = 0, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  step = tabu(data, score = "bic-g", k = gamma, max.iter = 1, start = start)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = tabu(data, score = "bic-g", k = gamma, max.iter = 1, start = step)
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
    } else if (samplesize/params[i] >=0.1641 & samplesize/params[i + 1]< 0.1641){
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

################################################################
# start will give other amount of ntests
# if using gamma = 0 USE BNLEARN
# if using !gamma = 0 USE modified_4.3
################################################################
gamma_tabu <- function(data, iterations, gamma = 0, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  
  if(!gamma == 0){
    step = tabu(data, score = "bic-g", k = gamma, max.iter = iterations, start = start)
  }
  if(gamma == 0){
    step = tabu(data, score = "bic-g", max.iter = iterations, start = start)
    #scores[1] = logLik(step, data)-log(nrow(data))/2*nparams(step, data = data)
  }
  
  networks.names<- paste0("tabu_10d_g",gamma)
  
  # test<-tabu(gaussian.test, score = "bic-g", k = 5, max.iter = 1000, start = NULL)
  # test$learning
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  networks[[1]] <- step
  names(networks) <- networks.names
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gamma = c(gamma))
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}



# test1 <- iteration_tabu(data = datapermutations[[2]], iterations = 20, gamma = 25)
# test2<- iteration_tabu(data = datapermutations[[2]], iterations = 20, gamma = 20)
# test3 <- iteration_tabu(data = datapermutations[[2]], iterations = 20, gamma = 15)
# test4 <- iteration_tabu(data = datapermutations[[2]], iterations = 20, gamma = 10)
# test5 <- iteration_tabu(data = datapermutations[[2]], iterations = 20, gamma = 5)
# test1$networkdata
# test2$networkdata
# test3$networkdata
# test4$networkdata
# test5$networkdata
# all.equal(test1$networkdata$logliks, test3$networkdata$logliks, test5$networkdata$logliks)
# test6 <- iteration_tabu(data = datapermutations[[2]], iterations = 20, gamma = 40)
# all.equal(test5$networkdata$logliks[1:3], test6$networkdata$logliks)
# test6$networkdata

iteration_tabu_sel <- function(data, iterations, select, gamma, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  step = tabu(data, score = "bic-g", max.iter = 1, k = gamma, start = start)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = tabu(data, score = "bic-g", max.iter = 1, k = gamma, start = step)
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
    if (edges[i+1] >= select[1] & edges[i]< select[1]) {
      networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    }
  
    if (length(scores) >= 10) {    
      if (!scores[i+1] > scores[i-8]) {
        message("is broken at", i+1)
        networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
        attr(stepi, "n/par") <- (samplesize/params[i + 1])
        networks[[i+1]] <- stepi
        break
      } else {
        step <-stepi
        i <- i+1}
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

##############################################3
# Gamma == 0 is only usable with package bnlearn normal
# !Gamma == 0 is only usable with package bnlearn_4.3_modified
# if you use the start function the ntests will not be realiable
#################################################
gamma_hc <- function(data, iterations, gamma = 0, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  if(!gamma == 0){
    step = hc(data, score = "bic-g", k = gamma, max.iter = iterations, start = start)
  }
  if(gamma == 0){
    step = hc(data, score = "bic-g", max.iter = iterations, start = start)
    #scores[1] = logLik(step, data)-log(nrow(data))/2*nparams(step, data = data)
    }
  networks.names<- paste0("hc_10d_g",gamma)
  
  scores[1] <- bnlearn::score(step, data)
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  logliks[1] <- logLik(step, data)
  
  networks[[1]] <- step
  names(networks) <- networks.names
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gamma = c(gamma))
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

#################################################
# The number of tests is not reliable because of the loose of memory in the iterations.
###################################################

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
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1641 & samplesize/params[i + 1]< 0.1641){
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    }
    
    
    if (!scores[i+1] > scores[i]) {
      message("is broken at", i+1)
      networks.names[i+1] <- paste0("hc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
      break
    } else {
      step <-stepi
      i <- i+1
    }
  }
  
  if (scores[i] > scores[i-1]) {
    networks.names[i] <- paste0("hc_10d_",i,"i")
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

###################################
# The number of tests output is not realible because of the iterations.
##################################
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
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
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
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
      attr(stepi, "n/par") <- (samplesize/params[i + 1])
      networks[[i+1]] <- stepi
    } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05){
      networks.names[i+1] <- paste0("mmhc_10d_",i+1,"i")
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

gamma_mmhc_split <- function(data, iterations, gamma, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  begin = mmpc(data, test = "bic-gt", B = gamma)
  
  fulligraph <- make_full_graph(ncol(data))
  fullbnlearn <- as.bn(igraph.to.graphNEL(fulligraph))
  nodes(fullbnlearn) <- names(data)
  
  igraph1 <- igraph.from.graphNEL(as.graphNEL(fullbnlearn))
  igraph2 <- igraph.from.graphNEL(as.graphNEL(begin))
  
  blackigraph <- difference(igraph1, igraph2)
  blackbngraph <- as.bn(igraph.to.graphNEL(blackigraph))
  
  step <- hc(data, score = "bic-g", k = gamma, blacklist = arcs(blackbngraph), max.iter = iterations)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  

  networks.names <- paste0("mmhc_10d_g",gamma)
  
  networks[[1]] <- step
  names(networks) <- networks.names
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gamma = c(gamma))
  output <- list(networks = networks, networkdata = networkdata, begin = begin)
  return(output)
}

iteration_mmhc_split_simple <- function(data, iterations, begin){
  
  fulligraph <- make_full_graph(ncol(data))
  fullbnlearn <- as.bn(igraph.to.graphNEL(fulligraph))
  nodes(fullbnlearn) <- names(data)
  
  igraph1 <- igraph.from.graphNEL(as.graphNEL(fullbnlearn))
  igraph2 <- igraph.from.graphNEL(as.graphNEL(begin))
  
  blackigraph <- difference(igraph1, igraph2)
  blackbngraph <- as.bn(igraph.to.graphNEL(blackigraph))
  
  step <- hc(data, score = "bic-g", blacklist = arcs(blackbngraph), max.iter = iterations)
  
  return(step)
}

iteration_h2pc_split <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  begin <- hpc(data, test = "bic-gt")
  
  fulligraph <- make_full_graph(ncol(data))
  fullbnlearn <- as.bn(igraph.to.graphNEL(fulligraph))
  nodes(fullbnlearn) <- names(data)
  
  igraph1 <- igraph.from.graphNEL(as.graphNEL(fullbnlearn))
  igraph2 <- igraph.from.graphNEL(as.graphNEL(begin))
  
  blackigraph <- difference(igraph1, igraph2)
  blackbngraph <- as.bn(igraph.to.graphNEL(blackigraph))
  
  step <- hc(data, score = "bic-g", blacklist = arcs(blackbngraph), max.iter = iterations)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  
  networks.names <- "h2pc_10d_g0"
  
  networks[[1]] <- step
  names(networks) <- networks.names
  
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata, begin = begin)
  return(output)
}

iteration_h2pc_split_part1 <- function(data){
  begin <- hpc(data, test = "bic-gt")
  return(begin)}
  
iteration_h2pc_split_part2 <- function(data, iterations, start = NULL){
  
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  fulligraph <- make_full_graph(ncol(data))
  fullbnlearn <- as.bn(igraph.to.graphNEL(fulligraph))
  nodes(fullbnlearn) <- names(data)
  
  igraph1 <- igraph.from.graphNEL(as.graphNEL(fullbnlearn))
  igraph2 <- igraph.from.graphNEL(as.graphNEL(start))
  
  blackigraph <- difference(igraph1, igraph2)
  blackbngraph <- as.bn(igraph.to.graphNEL(blackigraph))
  
  step <- hc(data, score = "bic-g", blacklist = arcs(blackbngraph), max.iter = iterations)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  
  networks.names <- "h2pc_10d_g0"
  
  networks[[1]] <- step
  names(networks) <- networks.names
  
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata, begin = start)
  return(output)
}

gamma_h2pc_split <- function(data, iterations, gamma, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  begin = hpc(data, test = "bic-gt", B = gamma)
  
  fulligraph <- make_full_graph(ncol(data))
  fullbnlearn <- as.bn(igraph.to.graphNEL(fulligraph))
  nodes(fullbnlearn) <- names(data)
  
  igraph1 <- igraph.from.graphNEL(as.graphNEL(fullbnlearn))
  igraph2 <- igraph.from.graphNEL(as.graphNEL(begin))
  
  blackigraph <- difference(igraph1, igraph2)
  blackbngraph <- as.bn(igraph.to.graphNEL(blackigraph))
  
  step <- hc(data, score = "bic-g", k = gamma, blacklist = arcs(blackbngraph), max.iter = iterations)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  
  networks.names <- paste0("h2pc_10d_g",gamma)
  
  networks[[1]] <- step
  names(networks) <- networks.names
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gamma = c(gamma))
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

iteration_pcstable <- function(data, gammas){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  gammalist <- c()
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
  while (i <= length(gammas)){
    eqclass <- pc.stable(data, test = "bic-gt", B = gammas[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    amat <- amat(eqclass)
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)

    if(isValidGraph(amat, type = "pdag")){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      gammalist[i] <- gammas[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
    } else {gammalist[i] <- gammas[i]
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
          networks.names[l] <- paste0("pcstable_10d_0.5")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.4 & samplesize/params[i] <= 0.4) { 
          if ((samplesize/params[i-1] - 0.4) < (0.4 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.4")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.3 & samplesize/params[i] <= 0.3) { 
          if ((samplesize/params[i-1] - 0.3) < (0.3 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.3")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.2 & samplesize/params[i] <= 0.2) { 
          if ((samplesize/params[i-1] - 0.2) < (0.2 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.2")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.19 & samplesize/params[i] <= 0.19) { 
          if ((samplesize/params[i-1] - 0.19) < (0.19- samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.19")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.18 & samplesize/params[i] <= 0.18) { 
          if ((samplesize/params[i-1] - 0.18) < (0.18 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.18")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.17 & samplesize/params[i] <= 0.17) { 
          if ((samplesize/params[i-1] - 0.17) < (0.17 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.17")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.16 & samplesize/params[i] <= 0.16) { 
          if ((samplesize/params[i-1] - 0.16) < (0.16 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.16")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.15 & samplesize/params[i] <= 0.15) { 
          if ((samplesize/params[i-1] - 0.15) < (0.15 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.15")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.1 & samplesize/params[i] <= 0.1) { 
          if ((samplesize/params[i-1] - 0.1) < (0.1 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.1")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.05 & samplesize/params[i] <= 0.05) { 
          if ((samplesize/params[i-1] - 0.05) < (0.05 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.05")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        }
        gbnmin1 <- gbn
      }
    }
    
    #gbnmin1 <- gbn
    i <- i +1
    
  }
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gammas = gammalist)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

iteration_pcstable_2 <- function(data, gammas){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  gammalist <- c()
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
  while (i <= length(gammas)){
    eqclass <- pc.stable(data, test = "bic-gt", B = gammas[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    amat <- amat(eqclass)
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)
    
    if(pdag2dag(getGraph(amat))$success == TRUE){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      gammalist[i] <- gammas[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
    } else {gammalist[i] <- gammas[i]
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
          networks.names[l] <- paste0("pcstable_10d_0.5")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.4 & samplesize/params[i] <= 0.4) { 
          if ((samplesize/params[i-1] - 0.4) < (0.4 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.4")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.3 & samplesize/params[i] <= 0.3) { 
          if ((samplesize/params[i-1] - 0.3) < (0.3 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.3")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.2 & samplesize/params[i] <= 0.2) { 
          if ((samplesize/params[i-1] - 0.2) < (0.2 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.2")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.19 & samplesize/params[i] <= 0.19) { 
          if ((samplesize/params[i-1] - 0.19) < (0.19- samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.19")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.18 & samplesize/params[i] <= 0.18) { 
          if ((samplesize/params[i-1] - 0.18) < (0.18 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.18")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.17 & samplesize/params[i] <= 0.17) { 
          if ((samplesize/params[i-1] - 0.17) < (0.17 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.17")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.16 & samplesize/params[i] <= 0.16) { 
          if ((samplesize/params[i-1] - 0.16) < (0.16 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.16")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.15 & samplesize/params[i] <= 0.15) { 
          if ((samplesize/params[i-1] - 0.15) < (0.15 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.15")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.1 & samplesize/params[i] <= 0.1) { 
          if ((samplesize/params[i-1] - 0.1) < (0.1 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.1")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.05 & samplesize/params[i] <= 0.05) { 
          if ((samplesize/params[i-1] - 0.05) < (0.05 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.05")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        }
        gbnmin1 <- gbn
      }
    }
    
    #gbnmin1 <- gbn
    i <- i +1
    
  }
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gammas = gammalist)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

# data <- datapermutations[[1]]
# gammas <- c(25)
iteration_pcstable_2_withnet <- function(data, gammas){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  gammalist <- c()
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
  while (i <= length(gammas)){
    eqclass <- pc.stable(data, test = "bic-gt", B = gammas[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    amat <- amat(eqclass)
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)
    
    if(pdag2dag(getGraph(amat))$success == TRUE){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      gammalist[i] <- gammas[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
      networks[[i]] <- dag
      networks.names[i] <- paste0("pcstable_10d_g",gammas[i])
    } else {gammalist[i] <- gammas[i]
    logliks[i] <- NA
    edges[i] <- NA
    params[i] <- NA
    scores[i] <- NA
    networks[[i]] <- eqclass
    networks.names[i] <- paste0("pcstable_10d_g",gammas[i])
    }
    
    i <- i +1
  }
  
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gammas = gammalist)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}




iteration_GS <- function(data, gammas){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  gammalist <- c()
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
  while (i <= length(gammas)){
    eqclass <- gs(data, test = "bic-gt", B = gammas[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    amat <- amat(eqclass)
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)
    
    if(isValidGraph(amat, type = "pdag")){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      gammalist[i] <- gammas[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
    } else {gammalist[i] <- gammas[i]
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
          networks.names[l] <- paste0("pcstable_10d_0.5")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.4 & samplesize/params[i] <= 0.4) { 
          if ((samplesize/params[i-1] - 0.4) < (0.4 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.4")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.3 & samplesize/params[i] <= 0.3) { 
          if ((samplesize/params[i-1] - 0.3) < (0.3 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.3")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.2 & samplesize/params[i] <= 0.2) { 
          if ((samplesize/params[i-1] - 0.2) < (0.2 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.2")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.19 & samplesize/params[i] <= 0.19) { 
          if ((samplesize/params[i-1] - 0.19) < (0.19- samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.19")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.18 & samplesize/params[i] <= 0.18) { 
          if ((samplesize/params[i-1] - 0.18) < (0.18 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.18")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.17 & samplesize/params[i] <= 0.17) { 
          if ((samplesize/params[i-1] - 0.17) < (0.17 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.17")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.16 & samplesize/params[i] <= 0.16) { 
          if ((samplesize/params[i-1] - 0.16) < (0.16 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.16")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.15 & samplesize/params[i] <= 0.15) { 
          if ((samplesize/params[i-1] - 0.15) < (0.15 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.15")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.1 & samplesize/params[i] <= 0.1) { 
          if ((samplesize/params[i-1] - 0.1) < (0.1 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.1")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.05 & samplesize/params[i] <= 0.05) { 
          if ((samplesize/params[i-1] - 0.05) < (0.05 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.05")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        }
        gbnmin1 <- gbn
      }
    }
    
    #gbnmin1 <- gbn
    i <- i +1
    
  }
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gammas = gammalist)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

iteration_GS_2 <- function(data, gammas){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  gammalist <- c()
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
  while (i <= length(gammas)){
    eqclass <- gs(data, test = "bic-gt", B = gammas[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    amat <- amat(eqclass)
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)
    
    if(pdag2dag(getGraph(amat))$success == TRUE){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      gammalist[i] <- gammas[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
    } else {gammalist[i] <- gammas[i]
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
          networks.names[l] <- paste0("pcstable_10d_0.5")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.4 & samplesize/params[i] <= 0.4) { 
          if ((samplesize/params[i-1] - 0.4) < (0.4 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.4")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.3 & samplesize/params[i] <= 0.3) { 
          if ((samplesize/params[i-1] - 0.3) < (0.3 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.3")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.2 & samplesize/params[i] <= 0.2) { 
          if ((samplesize/params[i-1] - 0.2) < (0.2 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.2")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.19 & samplesize/params[i] <= 0.19) { 
          if ((samplesize/params[i-1] - 0.19) < (0.19- samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.19")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.18 & samplesize/params[i] <= 0.18) { 
          if ((samplesize/params[i-1] - 0.18) < (0.18 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.18")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.17 & samplesize/params[i] <= 0.17) { 
          if ((samplesize/params[i-1] - 0.17) < (0.17 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.17")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.16 & samplesize/params[i] <= 0.16) { 
          if ((samplesize/params[i-1] - 0.16) < (0.16 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.16")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.15 & samplesize/params[i] <= 0.15) { 
          if ((samplesize/params[i-1] - 0.15) < (0.15 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.15")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.1 & samplesize/params[i] <= 0.1) { 
          if ((samplesize/params[i-1] - 0.1) < (0.1 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.1")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        } else if (samplesize/params[i-1] >= 0.05 & samplesize/params[i] <= 0.05) { 
          if ((samplesize/params[i-1] - 0.05) < (0.05 - samplesize/params[i])){
            l <- i-1; add <- gbnmin1} else {l <- i; add <- gbn}
          networks.names[l] <- paste0("pcstable_10d_0.05")
          attr(add, "n/par") <- (samplesize/params[l])
          attr(add, "gamma") <- (gammas[l])
          networks[[l]] <- add
        }
        gbnmin1 <- gbn
      }
    }
    
    #gbnmin1 <- gbn
    i <- i +1
    
  }
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gammas = gammalist)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

iteration_gs_2_withnet <- function(data, gammas){
  # if (is.null(correlation)) {correlation <- cor(data, method = method)} 
  # degreestat <- list(C = correlation, n = nrow(data))
  # 
  # names <- c()
  # for(i in 1:ncol(data)) names <- append(names,paste0("V",i))
  
  gammalist <- c()
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
  while (i <= length(gammas)){
    eqclass <- gs(data, test = "bic-gt", B = gammas[i])
    tests[i] <- ntests(eqclass)
    nedgeseqclass[i] <- narcs(eqclass)
    nundeqclass[i] <- nrow(undirected.arcs(eqclass))/2
    amat <- amat(eqclass)
    # amata <- as(eqclass,"amat")
    # bn.eqclass <- as.bn(eqclass@graph, check.cycles = FALSE)
    
    if(pdag2dag(getGraph(amat))$success == TRUE){
      dag <- cextend(eqclass)
      # bn.dag <- as.bn(dag$graph)
      edges[i]<- narcs(dag)
      
      # data maken variable names!
      # newdata <- as.data.frame(data)
      # names(newdata) <- nodes(bn.dag)
      gbn <- bn.fit(dag, data)
      
      gammalist[i] <- gammas[i]
      logliks[i] <- logLik(gbn, data)
      params[i] <- nparams(gbn)
      scores[i] <- score(dag, data)
      networks[[i]] <- dag
      networks.names[i] <- paste0("gs_10d_g",gammas[i])
    } else {gammalist[i] <- gammas[i]
    logliks[i] <- NA
    edges[i] <- NA
    params[i] <- NA
    scores[i] <- NA
    networks[[i]] <- eqclass
    networks.names[i] <- paste0("gs_10d_g",gammas[i])
    }
    
    i <- i +1
  }
  
  if(!length(networks) == 0){
    ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
    networks <- networks[ind0]
    names(networks) <- networks.names[ind0]}
  
  networkdata <- data.frame(edgeseqclass = nedgeseqclass, undeqclass = nundeqclass , edges = edges,logliks = logliks, scores = scores, params = params, tests = tests, gammas = gammalist)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}



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
