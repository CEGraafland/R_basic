###HC code change-------------

shuffle  <- function(m,n) {
  deck <- 1:m
  w <- numeric(n)
  for(i in 1:n){
    x <- sample(deck, m)
    w[i] <- 1*!(sum(deck==x)==0)
  }
  return(w)
}

r = list()
repeat{
  
  r[[length(r)+1]] = list(x,j)
}


# unified hill climbing implementation (both optimized and by spec).
hill.climbing2 = function(x, start, whitelist, blacklist, score, extra.args,
                         restart, perturb, max.iter, maxp, optimized, debug = FALSE) {
  #LIST
  r = list()
  # cache nodes' labels.
  nodes = names(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # set the iteration counter.
  iter = 1
  # check whether the score is score-equivalent.
  score.equivalence = is.score.equivalent(score, nodes, extra.args)
  # check whether the score is decomposable.
  score.decomposability = is.score.decomposable(score, nodes, extra.args)
  # allocate the cache matrix.
  cache = matrix(0, nrow = n.nodes, ncol = n.nodes)
  # nodes to be updated (all of them in the first iteration).
  updated = seq_len(n.nodes) - 1L
  # set the number of random restarts.
  restart.counter = restart
  
  # set the reference score.
  reference.score = per.node.score(network = start, score = score,
                                   targets = nodes, extra.args = extra.args, data = x)
  
  # convert the blacklist to an adjacency matrix for easy use.
  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)
  
  # convert the whitelist to an adjacency matrix for easy use.
  if (!is.null(whitelist))
    wlmat = arcs2amat(whitelist, nodes)
  else
    wlmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)
  
  if (debug) {
    
    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(start)
    cat("* current score:", sum(reference.score), "\n")
    cat("* whitelisted arcs are:\n")
    if (!is.null(whitelist)) print(whitelist)
    cat("* blacklisted arcs are:\n")
    if (!is.null(blacklist)) print(blacklist)
    
    # set the metadata of the network; othewise the debugging output is
    # confusing and not nearly as informative.
    start$learning$algo = "hc"
    start$learning$ntests = 0
    start$learning$test = score
    start$learning$args = extra.args
    start$learning$optimized = optimized
    
  }#THEN
  
  repeat {
    
    # build the adjacency matrix of the current network structure.
    amat = arcs2amat(start$arcs, nodes)
    # compute the number of parents of each node.
    nparents = colSums(amat)
    
    # set up the score cache (BEWARE: in place modification!).
    .Call("hc_cache_fill",
          nodes = nodes,
          data = x,
          network = start,
          score = score,
          extra = extra.args,
          reference = reference.score,
          equivalence = score.equivalence && optimized,
          decomposability = score.decomposability,
          updated = (if (optimized) updated else seq(length(nodes)) - 1L),
          env = environment(),
          amat = amat,
          cache = cache,
          blmat = blmat,
          debug = debug)
    
    # select which arcs should be tested for inclusion in the graph (hybrid
    # learning algorithms should hook the restrict phase here).
    to.be.added = arcs.to.be.added(amat = amat, nodes = nodes,
                                   blacklist = blmat, whitelist = NULL, nparents = nparents,
                                   maxp = maxp, arcs = FALSE)
    
    # get the best arc addition/removal/reversal.
    bestop = .Call("hc_opt_step",
                   amat = amat,
                   nodes = nodes,
                   added = to.be.added,
                   cache = cache,
                   reference = reference.score,
                   wlmat = wlmat,
                   blmat = blmat,
                   nparents = nparents,
                   maxp = maxp,
                   debug = debug)
    
    # the value FALSE is the canary value in bestop$op meaning "no operation
    # improved the network score"; break the loop.
    if (bestop$op == FALSE) {
      
      if ((restart.counter >= 0) && (restart > 0)) {
        
        if (restart.counter == restart) {
          
          # store away the learned network at the first restart.
          best.network = start
          best.network.score = sum(reference.score)
          
        }#THEN
        else {
          
          # if the network found by the algorithm is better, use that one as
          # starting network for the next random restart; use the old one
          # otherwise.
          if (best.network.score > sum(reference.score)) {
            
            start = best.network
            
            if (debug) {
              
              cat("----------------------------------------------------------------\n")
              cat("* the old network was better, discarding the current one.\n")
              cat("* now the current network is :\n")
              print(start)
              cat("* current score:", best.network.score, "\n")
              
            }#THEN
            
          }#THEN
          else {
            
            # store away the learned network at the first restart.
            best.network = start
            best.network.score = sum(reference.score)
            
          }#ELSE
          
        }#ELSE
        
        # this is the end; if the network found by the algorithm is the best one
        # it's now the 'end' one once more.
        if (restart.counter == 0) break
        
        # don't try to do anything if there are no more iterations left.
        if (iter >= max.iter) {
          
          if (debug)
            cat("@ stopping at iteration", max.iter, "ignoring random restart.\n")
          
          break
          
        }#THEN
        # increment the iteration counter.
        iter = iter + 1
        
        # decrement the number of random restart we still have to do.
        restart.counter = restart.counter - 1
        
        if (debug) {
          
          cat("----------------------------------------------------------------\n")
          cat("* doing a random restart,", restart.counter, "of", restart, "left.\n")
          
        }#THEN
        
        # perturb the network.
        start = perturb.backend(network = start, iter = perturb, nodes = nodes,
                                amat = amat, whitelist = whitelist, blacklist = blacklist,
                                maxp = maxp, debug = debug)
        
        # update the cached values of the end network.
        start$nodes = cache.structure(nodes, arcs = start$arcs)
        
        # update the scores of the nodes as needed.
        if (best.network.score > sum(reference.score)) {
          
          # all scores must be updated; this happens when both the network
          # resulted from the random restart is discarded and the old network
          # is perturbed.
          reference.score = per.node.score(network = start, score = score,
                                           targets = nodes, extra.args = extra.args, data = x)
          
          updated = which(nodes %in% nodes) - 1L
          
        }#THEN
        else {
          
          # the scores whose nodes' parents changed must be updated.
          reference.score[names(start$updates)] =
            per.node.score(network = start, score = score, targets = names(start$updates),
                           extra.args = extra.args, data = x)
          
          updated = which(nodes %in% names(start$updates)) - 1L
          
        }#THEN
        
        # nuke start$updates from orbit.
        start$updates = NULL
        
        next
        
      }#THEN
      else
        break
      
    }#THEN
    
    # update the network structure.
    start = arc.operations(start, from = bestop$from, to = bestop$to,
                           op = bestop$op, check.cycles = FALSE, update = TRUE,
                           debug = FALSE)
    
    # set the nodes whose cached score deltas are to be updated.
    if (bestop$op == "reverse")
      updated = which(nodes %in% c(bestop$from, bestop$to)) - 1L
    else
      updated = which(nodes %in% bestop$to) - 1L
    
    if (debug) {
      
      # update the test counter of the network; very useful to check how many
      # score comparison has been done up to now.
      start$learning$ntests = test.counter()
      
      cat("----------------------------------------------------------------\n")
      cat("* best operation was: ")
      if (bestop$op == "set")
        cat("adding", bestop$from, "->", bestop$to, ".\n")
      else if (bestop$op == "drop")
        cat("removing", bestop$from, "->", bestop$to, ".\n")
      else
        cat("reversing", bestop$from, "->", bestop$to, ".\n")
      cat("* current network is :\n")
      print(start)
      cat("* current score:", sum(reference.score), "\n")
      
    }#THEN
    
    # check the current iteration index against the max.iter parameter.
    if (iter >= max.iter) {
      
      if (debug)
        cat("@ stopping at iteration", max.iter, ".\n")
      
      break
      
    }#THEN
    else iter = iter + 1
    
  gbn <- bn.fit(start, x)
  likgbn <- logLik(gbn, x)  
  r[[length(r)+1]] = list(likgbn)  
  }#REPEAT
  
  return(list(start,
              r))
  
}#HILL.CLIMBING



# score-based learning algorithms.
greedy.search2 = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
                         score = "bic-g", heuristic = "hc2", expand, optimized = FALSE, debug = FALSE) {
  
  # check the data are there.
  #check.data(x)
  # check the algorithm.
  #check.learning.algorithm(heuristic, class = "score")
  # check the score label.
  #score = check.score(score, x)
  # check debug.
  #check.logical(debug)
  
  # check unused arguments in misc.args.
  misc.args = expand[names(expand) %in% method.extra.args[[heuristic]]]
  extra.args = expand[names(expand) %in% score.extra.args[[score]]]
  check.unused.args(expand, c(method.extra.args[[heuristic]], score.extra.args[[score]]))
  
  # expand and check the max.iter parameter (common to all algorithms).
  max.iter = check.max.iter(misc.args$max.iter)
  # expand and check the maxp parameter (common to all algorithms).
  maxp = check.maxp(misc.args$maxp, data = x)
  
  if (heuristic == "hc") {
    
    # expand and check the number of random restarts.
    restart = check.restart(misc.args$restart)
    # expand and check the magnitude of the perturbation when random restarts
    # are effectively used.
    perturb = ifelse((restart > 0), check.perturb(misc.args$perturb), 0)
    
  }#THEN
  else if (heuristic == "tabu") {
    
    # expand and check the arguments related to the tabu list.
    tabu = check.tabu(misc.args$tabu)
    max.tabu = check.max.tabu(misc.args$max.tabu, tabu)
    
  }#THEN
  
  # sanitize whitelist and blacklist, if any.
  whitelist = build.whitelist(whitelist, names(x))
  blacklist = build.blacklist(blacklist, whitelist, names(x))
  # if there is no preseeded network, use an empty one.
  if (is.null(start))
    start = empty.graph(nodes = names(x))
  else {
    
    # check start's class.
    check.bn(start)
    # check the preseeded network against the data set.
    check.bn.vs.data(start, x)
    # check the preseeded network against the maximum number of parents.
    nparents = sapply(start$nodes, function(x) length(x$parents))
    if (any(nparents > maxp))
      stop("nodes ", paste(names(which(nparents > maxp)), collapse = " "),
           " have more than 'maxp' parents.")
    
  }#ELSE
  
  # apply the whitelist to the preseeded network.
  if (!is.null(whitelist)) {
    
    for (i in 1:nrow(whitelist))
      start$arcs = set.arc.direction(whitelist[i, "from"],
                                     whitelist[i, "to"], start$arcs)
    
  }#THEN
  
  # apply the blacklist to the preseeded network.
  if (!is.null(blacklist)) {
    
    blacklisted = apply(start$arcs, 1, function(x){ is.blacklisted(blacklist, x) })
    start$arcs = start$arcs[!blacklisted, , drop = FALSE]
    
  }#THEN
  
  # be sure the graph structure is up to date.
  start$nodes = cache.structure(names(start$nodes), arcs = start$arcs)
  # no party if the graph is partially directed.
  if (is.pdag(start$arcs, names(start$nodes)))
    stop("the graph is only partially directed.")
  # check whether the graph is acyclic.
  if (!is.acyclic(arcs = start$arcs, nodes = names(start$nodes)))
    stop("the preseeded graph contains cycles.")
  
  # sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = start,
                                data = x, extra.args = extra.args)
  
  # reset the test counter.
  reset.test.counter()
  
  # call the right backend.
  if (heuristic == "hc") {
    
    res = hill.climbing2(x = x, start = start, whitelist = whitelist,
                        blacklist = blacklist, score = score, extra.args = extra.args,
                        restart = restart, perturb = perturb, max.iter = max.iter,
                        maxp = maxp, optimized = optimized, debug = debug)
    
  }#THEN
  else if (heuristic == "tabu"){
    
    res = tabu.search(x = x, start = start, whitelist = whitelist,
                      blacklist = blacklist, score = score, extra.args = extra.args,
                      max.iter = max.iter, optimized = optimized, tabu = tabu,
                      maxp = maxp, debug = debug)
    
  }#THEN
  
  # set the metadata of the network in one stroke.
  res$learning = list(whitelist = whitelist, blacklist = blacklist,
                      test = score, ntests = test.counter(),
                      algo = heuristic, args = extra.args, optimized = optimized)
  
  invisible(res)
  
}#GREEDY.SEARCH



# Hill Climbing greedy search frontend.
hc2 = function(x, start = NULL, whitelist = NULL, blacklist = NULL,
              score = NULL, ..., debug = FALSE, restart = 0, perturb = 1,
              max.iter = Inf, maxp = Inf, optimized = TRUE) {
  
  greedy.search2(x = x, start = start, whitelist = whitelist,
                blacklist = blacklist, score = score, heuristic = "hc2", debug = debug,
                expand = c(list(...), restart = restart, perturb = perturb,
                           max.iter = max.iter, maxp = maxp), optimized = optimized)
  
}#HC
hc2(gaussian.test)





### HC easy network store --------------------

scores = list()
edges = list()

step = hc(gaussian.test, max.iter = 1)
nedges.step <- narcs(step)
edges[1] = nedges.step
score.step <- score(step, gaussian.test)
scores[1] <- score.step
logLik.step <- logLik(step, gaussian.test)


step1 = hc(gaussian.test, max.iter = 1, start = step)
nedges.step1 <- narcs(step1)
edges[2] = nedges.step1
score.step1 <- score(step1, gaussian.test)
scores[2] <- score.step1
logLik.step1 <- logLik(step1, gaussian.test)


if(score.step1 > score.step){
  
}

step2 = hc(gaussian.test, max.iter = 1, start = step2)
nedges.step3 <- narcs(step3)
score.step3 <- score(step3, gaussian.test)
logLik.step3 <- logLik(step3, gaussian.test)



laststep = hc(gaussian.test)
nedges.laststep <- narcs(laststep)
logLik(laststep, gaussian.test)

score.step1 <- score(step1, gaussian.test)
score.step1
logLik(step1, gaussian.test)

laststep

hc_opt_step

#nieuw
logliks = list()
scores = list()
edges = list()
step = hc(gaussian.test, max.iter = 1)

edges[1] <- narcs(step)
scores[1] <- score(step, gaussian.test)
logliks[1] <- logLik(step, gaussian.test)

i<-1
while (i <80){
  stepi = hc(gaussian.test, max.iter = 1, start = step)
  nedges.stepi <- narcs(stepi)
  edges[i+1] <- nedges.stepi
  score.stepi <- score(stepi, gaussian.test)
  scores[i+1] <- score.stepi
  loglik.stepi <- logLik(stepi, gaussian.test)
  logliks[i+1] <- loglik.stepi

  if (!unlist(scores[i+1])>unlist(scores[i])) {
    print(paste0("is broken at", i+1))
    break}
  else{step <-stepi
      i <- i+1}
}
########################################################################
# Hillclimbing iterations
########################################################################
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

hc_edges_loglik_gaussian_5i <- iteration_loglik_hc(gaussian.test, 5)
hc_edges_loglik_gaussian_5i$networks
hc_edges_loglik_gaussian_7i <- iteration_loglik_hc(gaussian.test, 2, start = hc_edges_loglik_gaussian_5i$networks)
hc_edges_loglik_gaussian_7i
hc_edges_loglik_gaussian <- iteration_loglik_hc(gaussian.test, 100)
hc_edges_loglik_gaussian$networks

# Manual hill climbing 30d
data30d <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_30d))
hc_edges_loglik_30d_100i <- iteration_loglik_hc(data30d, 100)

# Manual hill climbing 20d
data20d <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_20d2))
data20d
hc_edges_loglik_20d_100i <- iteration_loglik_hc(data20d, 100, start = NULL)
hc_edges_loglik_20d_200i <- iteration_loglik_hc(data20d, 200)
hc_edges_loglik_20d_200_400i <- iteration_loglik_hc(data20d, 200, start = hc_edges_loglik_20d_200i$networks)
hc_edges_loglik_20d_400_600i <- iteration_loglik_hc(data20d, 200, start = hc_edges_loglik_20d_200_400i$networks)

# Manual hill climbing 10d
data10d <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
hc_edges_loglik_10d_100i <- iteration_loglik_hc(data10d, 100)
hc_edges_loglik_10d_200i <- iteration_loglik_hc(data10d, 200)
hc_edges_loglik_10d_700i <- iteration_loglik_hc(data10d, 700)
hc_edges_loglik_10d_200_400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_200i$networks)
hc_edges_loglik_10d_400_600i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_200_400i$networks)
hc_edges_loglik_10d_600_800i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_400_600i$networks)
hc_edges_loglik_10d_800_1000i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_600_800i$networks)
hc_edges_loglik_10d_1000_1200i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_800_1000i$networks)
hc_edges_loglik_10d_1200_1400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_1000_1200i$networks)
hc_edges_loglik_10d_1400_1600i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_1200_1400i$networks)
hc_edges_loglik_10d_1600_1800i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_1400_1600i$networks)
hc_edges_loglik_10d_1800_2000i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_1600_1800i$networks)
hc_edges_loglik_10d_2000_2200i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_1800_2000i$networks)
hc_edges_loglik_10d_2200_2400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_2000_2200i$networks)
hc_edges_loglik_10d_2400_2600i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_2200_2400i$networks)
hc_edges_loglik_10d_2600_2800i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_2400_2600i$networks)
hc_edges_loglik_10d_2800_3000i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_2600_2800i$networks)
hc_edges_loglik_10d_3000_3200i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_2800_3000i$networks)
hc_edges_loglik_10d_3200_3400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_3000_3200i$networks)
hc_edges_loglik_10d_3400_3600i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_3200_3400i$networks)
hc_edges_loglik_10d_3600_3800i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_3400_3600i$networks)
hc_edges_loglik_10d_3800_4000i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_3600_3800i$networks)
hc_edges_loglik_10d_4000_4200i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_3800_4000i$networks)
hc_edges_loglik_10d_4200_4400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_4000_4200i$networks)
hc_edges_loglik_10d_4400_4600i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_4200_4400i$networks)
hc_edges_loglik_10d_4600_4800i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_4400_4600i$networks)
hc_edges_loglik_10d_4800_5000i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_4600_4800i$networks)
hc_edges_loglik_10d_5000_5200i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_4800_5000i$networks)
hc_edges_loglik_10d_5200_5400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_5000_5200i$networks)
hc_edges_loglik_10d_5400_5600i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_5200_5400i$networks)
hc_edges_loglik_10d_5600_5800i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_5400_5600i$networks)
hc_edges_loglik_10d_5800_6000i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_5600_5800i$networks)
hc_edges_loglik_10d_6000_6200i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_5800_6000i$networks)
hc_edges_loglik_10d_6200_6400i <- iteration_loglik_hc(data10d, 200, start = hc_edges_loglik_10d_6000_6200i$networks)

# last one:
hc_edges_loglik_10d_8600_8695i <- iteration_loglik_hc(data10d, 95, start = hc_edges_loglik_10d_8400_8600i$networks)
hc_edges_loglik_10d_8600_8695i

# Calculate hill climbing with loop
# gayssian.test

 # i <-2
  # start <- iteration_loglik_hc(gaussian.test, 2)
  # while ( i < 10) {
    # j <- i +2
    # berekening <- iteration_loglik_hc(gaussian.test, 2, start = start$networks)
    # assign(paste0("gs_edges_loglik_10d_",i,"_",j,"i"),berekening)
    # start <- berekening
    # i <- i+2
  # }
  
# Automatically Hill climbing 10 d with loop

  i <- 6400
  start <- hc_edges_loglik_10d_6200_6400i
  while ( i < 6800) {
    j <- i+200
    berekening <- iteration_loglik_hc(data10d, 200, start = start$networks)
    assign(paste0("hc_edges_loglik_10d_",i,"_",j,"i"), berekening)
    start <- berekening
    i <- i+200
  }
  
  i <- 6800
  start <- hc_edges_loglik_10d_6600_6800i
  while ( i < 8000) {
    j <- i+200
    berekening <- iteration_loglik_hc(data10d, 200, start = start$networks)
    assign(paste0("hc_edges_loglik_10d_",i,"_",j,"i"), berekening)
    start <- berekening
    i <- i+200
  }
  
  ## Last sequence: is broken at 95 -> 94 the best. in 8600 - 8800 (2na laatste) (next two are broken at two)
  i <- 8000
  start <- hc_edges_loglik_10d_7800_8000i
  data10d <- as.data.frame(TimeCoordsAnom_from_Grid(tas_ncep_10d))
  while ( i < 9200) {
    j <- i+200
    berekening <- iteration_loglik_hc(data10d, 200, start = start$networks)
    assign(paste0("hc_edges_loglik_10d_",i,"_",j,"i"), berekening)
    start <- berekening
    i <- i+200
  }
  # variablenames
  # i<-2
  # variablenames[2]
  # variablenames[i-1]
  # for (i in 2:length(variablenames)){
    # assign(start, variablenames[i-1])
    # start
    # assign(variablenames[i],narcs(start$networks))
    # start <- i 
  # }
  
  # variablenames
  # edges_loglik_10d_5400_5600i
  # variablenames <- c(variablenames,1)
# variablenames
  
# Bind the networkdata
hc_edges_loglik_20d_400i <- rbind(hc_edges_loglik_10d_200i$networkdata, hc_edges_loglik_10d_200_400i$networkdata)
hc_edges_loglik_20d_600i <- rbind(hc_edges_loglik_20d_200i$networkdata, 
                                 hc_edges_loglik_20d_200_400i$networkdata,
                                 hc_edges_loglik_20d_400_600i$networkdata)

hc_edges_loglik_10d_1000i <- rbind(hc_edges_loglik_10d_200i$networkdata,
      hc_edges_loglik_10d_200_400i$networkdata,
      hc_edges_loglik_10d_400_600i$networkdata,
      hc_edges_loglik_10d_600_800i$networkdata, 
      hc_edges_loglik_10d_800_1000i$networkdata)

hc_edges_loglik_10d_1400i <- rbind(hc_edges_loglik_10d_200i$networkdata,
                                   hc_edges_loglik_10d_200_400i$networkdata,
                                   hc_edges_loglik_10d_400_600i$networkdata,
                                   hc_edges_loglik_10d_600_800i$networkdata, 
                                   hc_edges_loglik_10d_800_1000i$networkdata,
                                   hc_edges_loglik_10d_1000_1200i$networkdata,
                                   hc_edges_loglik_10d_1200_1400i$networkdata)

hc_edges_loglik_10d_1600i <- rbind(hc_edges_loglik_10d_1400i,hc_edges_loglik_10d_1400_1600i$networkdata)
hc_edges_loglik_10d_1800i <- rbind(hc_edges_loglik_10d_1600i,hc_edges_loglik_10d_1600_1800i$networkdata)
hc_edges_loglik_10d_2000i <- rbind(hc_edges_loglik_10d_1800i,hc_edges_loglik_10d_1800_2000i$networkdata)
hc_edges_loglik_10d_2200i <- rbind(hc_edges_loglik_10d_2000i,hc_edges_loglik_10d_2000_2200i$networkdata)
hc_edges_loglik_10d_2400i <- rbind(hc_edges_loglik_10d_2200i,hc_edges_loglik_10d_2200_2400i$networkdata)
hc_edges_loglik_10d_2600i <- rbind(hc_edges_loglik_10d_2400i,hc_edges_loglik_10d_2400_2600i$networkdata)

hc_edges_loglik_10d_600_2000i <- rbind(hc_edges_loglik_10d_600_800i$networkdata, 
                                       hc_edges_loglik_10d_800_1000i$networkdata,
                                       hc_edges_loglik_10d_1000_1200i$networkdata,
                                       hc_edges_loglik_10d_1200_1400i$networkdata, 
                                       hc_edges_loglik_10d_1400_1600i$networkdata,
                                       hc_edges_loglik_10d_1600_1800i$networkdata,
                                       hc_edges_loglik_10d_1800_2000i$networkdata)

hc_edges_loglik_10d_600_3000i <- rbind(hc_edges_loglik_10d_600_2000i,
                                       hc_edges_loglik_10d_2000_2200i$networkdata,
                                       hc_edges_loglik_10d_2200_2400i$networkdata,
                                       hc_edges_loglik_10d_2400_2600i$networkdata,
                                       hc_edges_loglik_10d_2600_2800i$networkdata,
                                       hc_edges_loglik_10d_2800_3000i$networkdata)

hc_edges_loglik_10d_3000i <- rbind(hc_edges_loglik_10d_2600i,
                                   hc_edges_loglik_10d_2600_2800i$networkdata,
                                   hc_edges_loglik_10d_2800_3000i$networkdata)

hc_edges_loglik_10d_3600i <- rbind(hc_edges_loglik_10d_3000i,
                                   hc_edges_loglik_10d_3000_3200i$networkdata,
                                   hc_edges_loglik_10d_3200_3400i$networkdata,
                                   hc_edges_loglik_10d_3400_3600i$networkdata)
hc_edges_loglik_10d_4600i <- rbind(hc_edges_loglik_10d_3600i,
                                   hc_edges_loglik_10d_3600_3800i$networkdata,
                                   hc_edges_loglik_10d_3800_4000i$networkdata,
                                   hc_edges_loglik_10d_4000_4200i$networkdata,
                                   hc_edges_loglik_10d_4200_4400i$networkdata,
                                   hc_edges_loglik_10d_4400_4600i$networkdata)
hc_edges_loglik_10d_5200i <- rbind(hc_edges_loglik_10d_4600i,
                                   hc_edges_loglik_10d_4600_4800i$networkdata,
                                   hc_edges_loglik_10d_4800_5000i$networkdata,
                                   hc_edges_loglik_10d_5000_5200i$networkdata)

hc_edges_loglik_10d_6400i <- rbind(hc_edges_loglik_10d_5200i,
                                   hc_edges_loglik_10d_5200_5400i$networkdata,
                                   hc_edges_loglik_10d_5400_5600i$networkdata,
                                   hc_edges_loglik_10d_5600_5800i$networkdata,
                                   hc_edges_loglik_10d_5800_6000i$networkdata,
                                   hc_edges_loglik_10d_6000_6200i$networkdata,
                                   hc_edges_loglik_10d_6200_6400i$networkdata)

# Bind with loop
i <- 6400
vector <- hc_edges_loglik_10d_6400i

while ( i < 8000) {
  j <- i+200
  vector <- rbind(vector, eval(parse(text = paste0("hc_edges_loglik_10d_",i,"_",j,"i$networkdata"))))
  vector
  i <- i+200
}
vector
hc_edges_loglik_10d_8000i <- vector
rm(vector)

i <- 8000
vector <- hc_edges_loglik_10d_8000i

while ( i < 8600) {
  j <- i+200
  vector <- rbind(vector, eval(parse(text = paste0("hc_edges_loglik_10d_",i,"_",j,"i$networkdata"))))
  vector
  i <- i+200
}
vector
hc_edges_loglik_10d_8600i <- vector
rm(vector)

hc_edges_loglik_10d_8695i <- rbind(hc_edges_loglik_10d_8600i,hc_edges_loglik_10d_8600_8695i$networkdata)

# make list with fres Hill climbing networks
i <- 200
newnetworks <- c()
while ( i < 9200) {
  j <- i+200
  newnetworks <- append(newnetworks, paste0("hc_edges_loglik_10d_",i,"_",j,"i"))
  newnetworks
  i <- i+200
}
newnetworks

###########################################################################
# Script old and new for iterations marco
###########################################################################
iteration_tabu_dev <- function(data, iterations, start = NULL){
  tests <- c()
  params<- c()
  logliks <- c()
  scores <- c()
  edges <- c()
  networks <- list()
  networks.names <- character()
  samplesize <- nrow(data)
  
  step = tabu(data, max.iter = 1, start = start)
  
  tests[1] <- ntests(step)
  params[1] <- nparams(step, data)
  edges[1] <- narcs(step)
  scores[1] <- bnlearn::score(step, data)
  logliks[1] <- logLik(step, data)
  
  i<-1
  while (i < iterations){
    stepi = tabu(data, max.iter = 1, start = step)
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
    # i <- 3
    
    ###########################
    if (samplesize/params[i] >= 400& samplesize/params[i + 1]< 400) {
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
      break
    } else {
      step <-stepi
      i <- i+1
    }
    #########################
    
    # if (samplesize/params[i] >= 400& samplesize/params[i + 1]< 400) {
    #   # save
    #   assign(paste0("tabu_10d_",i+1,"i"),stepi)
    #   networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
    #   attr(parse(text = paste0("tabu_10d_",i+1,"i")), "n/par") <- (samplesize/params[i + 1])
    # } else if (samplesize/params[i] >= 0.4 & samplesize/params[i + 1]< 0.4) {
    #   assign(paste0("tabu_10d_",i+1,"i"),stepi)
    #   networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
    # } else if (samplesize/params[i] >= 0.3 & samplesize/params[i + 1]< 0.3) {
    #   assign(paste0("tabu_10d_",i+1,"i"),stepi)
    #   networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
    # } else if (samplesize/params[i] >= 0.2 & samplesize/params[i + 1]< 0.2) {
    #   assign(paste0("tabu_10d_",i+1,"i"),stepi)
    #   networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
    # } else if (samplesize/params[i] >=0.1 & samplesize/params[i + 1]< 0.1) {
    #   assign(paste0("tabu_10d_",i+1,"i"),stepi)
    #   networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
    # } else if (samplesize/params[i] >= 0.05 & samplesize/params[i + 1] < 0.05) {
    #   assign(paste0("tabu_10d_",i+1,"i"),stepi)
    #   networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
    # }
    
    # if (!scores[i+1] > scores[i]) {
    #   print(paste0("is broken at", i+1))
    #   break
    # } else {
    #   step <-stepi
    #   i <- i+1
    # }
  }
  
  networks.names[i+1] <- paste0("tabu_10d_",i+1,"i")
  attr(stepi, "n/par") <- (samplesize/params[i + 1])
  networks[[i+1]] <- stepi
  
  ###################################
  ind0 <- which(unlist(lapply(networks, function(x) !is.null(x))))
  networks <- networks[ind0]
  names(networks) <- networks.names[ind0]
  #################################
  str(networks)
  
  
  # assign(paste0("tabu_10d_",i+1,"i"),stepi)
  # networks <- append(networks,paste0("tabu_10d_",i+1,"i"))
  
  networkdata <- data.frame(edges = edges,logliks = logliks, scores = scores, params = params, tests = tests)
  output <- list(networks = networks, networkdata = networkdata)
  return(output)
}

##################################################################
# PLOT
##################################################################

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
  
}#plot.nedge.vs.loglik_hc



plotgausian <- plot.nedge.vs.loglik_hc(hc_edges_loglik_gaussian)

plot30d100i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_30d_100i)

plot20d100i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_100i)
plot20d200i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_200i)
plot20d200_400i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_200_400i$networkdata)
plot20d600i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_600i)
plot20d400i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_400i)
plot20d400_600i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_400_600i$networkdata)

plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_9000_9200i$networkdata)
###opslaan --------------
#opslaan plot 20d met iteraties 
hcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/hc20d600i.pdf"
pdf(hcplotname)
plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_600i)
dev.off()

#opslaan plot 10d met 1000 iteraties
hcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/hc10d1000i.pdf"
pdf(hcplotname)
plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1000i)
dev.off()

#opslaan plot 10d met 1200 iteraties
hcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/hc10d1200i.pdf"
pdf(hcplotname)
plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1200i)
dev.off()

#difference plot 600 tot 3000i and tot 0 tot 600i and 0 tot 3000i
hcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/hc10d3000i.pdf"
pdf(hcplotname)
## divide the device into two rows and two columns
## allocate figure 1 all of row 1
## allocate figure 2 the intersection of column 2 and row 2
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
## show the regions that have been allocated to each plot
layout.show(3)
plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_3000i)
plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_700i)
plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_600_3000i)
dev.off()


plot10d100i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_100i)
plot10d700i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_700i)
plot10d600_800i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_600_800i$networkdata)
plot10d1400_1600i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1400_1600i$networkdata)
plot10d1600_1800i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1600_1800i$networkdata)
plot10d1800_2000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1800_2000i$networkdata)
plot10d2000_2200i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_2000_2200i$networkdata)
plot10d800_1000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_800_1000i$networkdata)
plot10d600_2000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_600_2000i)
plor10d600_3000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_600_3000i)
plot10d1000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1000i)
plot10d1200i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1200i)
plot10d1400i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1400i)
plot10d1600i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1600i)
plot10d1800i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_1800i)
plot10d2000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_2000i)
plot10d2200i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_2200i)
plot10d2400i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_2400i)
plot10d2600i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_2600i)
plot10d3000i <- plot.nedge.vs.loglik_hc(hc_edges_loglik_10d_3000i)

#opslaan netwerken 10d hill climbing
load(file = "/Users/lisettegraafland/Documents/R_practice/Data/data_hcnetworks10d.rda")

save(hc_edges_loglik_10d_200i,
      hc_edges_loglik_10d_200_400i,
      hc_edges_loglik_10d_400_600i,
      hc_edges_loglik_10d_600_800i,
      hc_edges_loglik_10d_800_1000i,
      hc_edges_loglik_10d_1000_1200i,
      hc_edges_loglik_10d_1200_1400i,
      hc_edges_loglik_10d_1400_1600i,
      hc_edges_loglik_10d_1600_1800i,
      hc_edges_loglik_10d_1800_2000i,
      hc_edges_loglik_10d_2000_2200i,
      hc_edges_loglik_10d_2200_2400i,
      hc_edges_loglik_10d_2400_2600i,
     hc_edges_loglik_10d_2600_2800i,
     hc_edges_loglik_10d_2800_3000i,
     hc_edges_loglik_10d_3000_3200i,
     hc_edges_loglik_10d_3200_3400i,
     hc_edges_loglik_10d_3400_3600i,
     hc_edges_loglik_10d_3600_3800i,
     hc_edges_loglik_10d_3800_4000i,
     hc_edges_loglik_10d_4000_4200i,
     hc_edges_loglik_10d_4200_4400i,
     hc_edges_loglik_10d_4400_4600i,
     hc_edges_loglik_10d_4600_4800i,
     hc_edges_loglik_10d_4800_5000i,
     hc_edges_loglik_10d_5000_5200i,
     hc_edges_loglik_10d_5200_5400i,
     hc_edges_loglik_10d_5400_5600i,
     hc_edges_loglik_10d_5600_5800i,
     hc_edges_loglik_10d_5800_6000i,
     hc_edges_loglik_10d_6000_6200i,
     hc_edges_loglik_10d_6200_6400i,
     hc_edges_loglik_10d_6400_6600i,
     hc_edges_loglik_10d_6600_6800i,
      hc_edges_loglik_10d_100i,
      hc_edges_loglik_10d_1000i,
      hc_edges_loglik_10d_1400i,
      hc_edges_loglik_10d_1600i,
      hc_edges_loglik_10d_1800i,
      hc_edges_loglik_10d_2000i,
      hc_edges_loglik_10d_2200i,
      hc_edges_loglik_10d_2400i,
      hc_edges_loglik_10d_2600i,
      hc_edges_loglik_10d_3000i,
      hc_edges_loglik_10d_3600i,
      hc_edges_loglik_10d_700i,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda")

save(hc_edges_loglik_10d_8600_8695i,
     list = newnetworks,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/networks_hcnetworks10d.rda")

save(hc_edges_loglik_10d_200_400i,
     hc_edges_loglik_10d_400_600i,
     hc_edges_loglik_10d_600_800i,
     hc_edges_loglik_10d_1000_1200i,
     hc_edges_loglik_10d_1400_1600i,
     hc_edges_loglik_10d_2400_2600i,
     hc_edges_loglik_10d_4800_5000i,
     hc_edges_loglik_10d_8600_8695i,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/utilnetworks_hcnetworks10d.rda")

# save(hc_edges_loglik_10d_200i,
#      hc_edges_loglik_10d_100i,
#      hc_edges_loglik_10d_1000i,
#      hc_edges_loglik_10d_1400i,
#      hc_edges_loglik_10d_1600i,
#      hc_edges_loglik_10d_1800i,
#      hc_edges_loglik_10d_2000i,
#      hc_edges_loglik_10d_2200i,
#      hc_edges_loglik_10d_2400i,
#      hc_edges_loglik_10d_2600i,
#      hc_edges_loglik_10d_3000i,
#      hc_edges_loglik_10d_3600i,
#      hc_edges_loglik_10d_700i,
#      hc_edges_loglik_10d_8000i,
#      file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda")
# View(load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/data_hcnetworks10d.rda"))

rm(hc_edges_loglik_10d_200i,
   hc_edges_loglik_10d_200_400i,
   hc_edges_loglik_10d_400_600i,
   hc_edges_loglik_10d_600_800i,
   hc_edges_loglik_10d_800_1000i,
   hc_edges_loglik_10d_1000_1200i,
   hc_edges_loglik_10d_1200_1400i,
   hc_edges_loglik_10d_1400_1600i,
   hc_edges_loglik_10d_1600_1800i,
   hc_edges_loglik_10d_1800_2000i,
   hc_edges_loglik_10d_2000_2200i,
   hc_edges_loglik_10d_2200_2400i,
   hc_edges_loglik_10d_2400_2600i,
   hc_edges_loglik_10d_2600_2800i,
   hc_edges_loglik_10d_2800_3000i,
   hc_edges_loglik_10d_3000_3200i,
   hc_edges_loglik_10d_3200_3400i,
   hc_edges_loglik_10d_3400_3600i,
   hc_edges_loglik_10d_3600_3800i,
   hc_edges_loglik_10d_3800_4000i,
   hc_edges_loglik_10d_4000_4200i,
   hc_edges_loglik_10d_4200_4400i,
   hc_edges_loglik_10d_4400_4600i,
   hc_edges_loglik_10d_4600_4800i,
   hc_edges_loglik_10d_4800_5000i,
   hc_edges_loglik_10d_5000_5200i,
   hc_edges_loglik_10d_5200_5400i,
   hc_edges_loglik_10d_5400_5600i,
   hc_edges_loglik_10d_5600_5800i,
   hc_edges_loglik_10d_5800_6000i,
   hc_edges_loglik_10d_6000_6200i,
   hc_edges_loglik_10d_6200_6400i,
   hc_edges_loglik_10d_6400_6600i,
   hc_edges_loglik_10d_6600_6800i,
   hc_edges_loglik_10d_6800_7000i,
   hc_edges_loglik_10d_7000_7200i,
   hc_edges_loglik_10d_7200_7400i,
   hc_edges_loglik_10d_7400_7600i,
   hc_edges_loglik_10d_7600_7800i,
   hc_edges_loglik_10d_7800_8000i,
   hc_edges_loglik_10d_8000i,
   hc_edges_loglik_10d_100i,
   hc_edges_loglik_10d_1000i,
   hc_edges_loglik_10d_1400i,
   hc_edges_loglik_10d_1600i,
   hc_edges_loglik_10d_1800i,
   hc_edges_loglik_10d_2000i,
   hc_edges_loglik_10d_2200i,
   hc_edges_loglik_10d_2400i,
   hc_edges_loglik_10d_2600i,
   hc_edges_loglik_10d_3000i,
   hc_edges_loglik_10d_3600i,
   hc_edges_loglik_10d_700i)

for (i in list(hc_edges_loglik_10d_200i,
           hc_edges_loglik_10d_200_400i,
           hc_edges_loglik_10d_400_600i,
           hc_edges_loglik_10d_600_800i,
           hc_edges_loglik_10d_800_1000i,
           hc_edges_loglik_10d_1000_1200i)){

  hcplotname <- paste0("/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/network10d",length(i$networks$arcs)/2,".pdf")  
  pdf(hcplotname)
  plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d),i$networks)
  title(paste0("number of edges = ",length(i$networks$arcs)/2), line = -4)
  dev.off()
}

# plot different networks that came from different iterations
# in one frame. (and save it)
hcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/networks10d200_1600i.pdf"
pdf(hcplotname)

par(mar=c(1,1,2,2))
par(mfrow = c(4, 2))
for (i in list(hc_edges_loglik_10d_200i,
               hc_edges_loglik_10d_200_400i,
               hc_edges_loglik_10d_400_600i,
               hc_edges_loglik_10d_600_800i,
               hc_edges_loglik_10d_800_1000i,
               hc_edges_loglik_10d_1000_1200i,
               hc_edges_loglik_10d_1200_1400i,
               hc_edges_loglik_10d_1400_1600i)){
  
  plot.Meteodag(TimeCoordsAnom_from_Grid(tas_ncep_10d),i$networks)
  title(paste0("number of edges = ",length(i$networks$arcs)/2))
}
dev.off()

#oefenen met plotten
hcplotname <- "/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/hc20d600i.pdf"
pdf(hcplotname)
plot.nedge.vs.loglik_hc(hc_edges_loglik_20d_600i)
dev.off()
toString(substitute(hc_edges_loglik_10d_1000_1200i))
sprintf
hc_edges_loglik_10d_1000_1200i$networks

paste0("/Users/lisettegraafland/Documents/R_practice/plots/hillclimbing/network10d",length(i$networks$arcs)/2,".pdf")

# oefenen
nparams(hc_edges_loglik_10d_1000_1200i$networks, TimeCoordsAnom_from_Grid(tas_ncep_10d))
