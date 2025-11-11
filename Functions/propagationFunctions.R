#####################################################################################
# Evidence propagation function
# necesario timeCoordsAnom_from_Grid_rms              (same as FINAL FUNCTION in "old)
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
  
  
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} 
  else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
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
####################################################################################
# Propagation from correlation matrix or covariance matrix (package condnorm)
####################################################################################

propagationCorr <- function(cormatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- cormatrix
  # nodesEvents <- 1:648
  # valueEvent <- 1
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  m <- nrow(cormatrix)
  names <- c()
  names
  for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V >", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
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
    
    
    probname <- paste0("P(V >", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  upper <- rep(Inf, (m-length(nodesEvidence)))
  upperw <- rep(Inf, m)
  
  for(i in 1:(m-length(nodesEvidence))) {
    lower2 <- rep(-Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    lower2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
                        dependent.ind = dependent, given.ind = nodesEvidence, X.given = valueEvidence,
                        check.sigma= FALSE)
  }
  
 for(i in 1:length(valueEvidence)){
  if(valueEvidence[i]>=0){
    with[nodesEvidence[i]] <- 1
  } else {with[nodesEvidence[i]] <- 0}
   }
  
  for(i in 1:m){
    lower2w <- rep(-Inf,m)
    lower2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), sigma = cormatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V >", valueEvent,")")
  df <- data.frame(names = names[nodesEvents], with = with, without = without)
  return(df)
}

# werkt niet (alleen maar corr = ipv sigm = voor evt berekenen van n > 1000)
propagationCorrspec <- function(cormatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- cormatrix
  # nodesEvents <- 1:648
  # valueEvent <- 1
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  m <- nrow(cormatrix)
  names <- c()
  names
  for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V >", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
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
    
    
    probname <- paste0("P(V >", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  upper <- rep(Inf, (m-length(nodesEvidence)))
  upperw <- rep(Inf, m)
  
  for(i in 1:(m-length(nodesEvidence))) {
    lower2 <- rep(-Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    lower2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower2, upper=upper, rep(0, m), sigma = cormatrix,
                        dependent.ind = dependent, given.ind = nodesEvidence, X.given = valueEvidence,
                        check.sigma= FALSE)
  }
  
  for(i in 1:length(valueEvidence)){
    if(valueEvidence[i]>=0){
      with[nodesEvidence[i]] <- 1
    } else {with[nodesEvidence[i]] <- 0}
  }
  
  for(i in 1:m){
    lower2w <- rep(-Inf,m)
    lower2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), corr = cormatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V >", valueEvent,")")
  df <- data.frame(names = names[nodesEvents], with = with, without = without)
  return(df)
}




propagationnegCorr <- function(cormatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- cormatrix
  # nodesEvents <- 1:648
  # valueEvent <- 1
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  m <- nrow(cormatrix)
  names <- c()
  names
  for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V <", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
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
    
    
    probname <- paste0("P(V <", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  lower <- rep(-Inf, (m-length(nodesEvidence)))
  lowerw <- rep(-Inf, m)
  
  # upper <- rep(Inf, (m-length(nodesEvidence)))
  # upperw <- rep(Inf, m)
  
  for(i in 1:(m-length(nodesEvidence))) {
    # lower2 <- rep(-Inf,m-length(nodesEvidence))
    upper2 <- rep(Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    # lower2[c(i)] <- valueEvent
    upper2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower, upper=upper2, rep(0, m), cormatrix,
                        dependent.ind = dependent, given.ind = nodesEvidence, X.given = valueEvidence,
                        check.sigma= FALSE)
  }
  for(i in 1:length(valueEvidence)){
    if(valueEvidence[i]<=0){
      with[nodesEvidence[i]] <- 1
    } else {with[nodesEvidence[i]] <- 0}
  }
  for(i in 1:m){
    # lower2w <- rep(-Inf,m)
    # lower2w[c(i)] <- valueEvent
    upper2w <- rep(Inf,m)
    upper2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lowerw, upper = upper2w, mean = rep(0, m), sigma = cormatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V <", valueEvent,")")
  df <- data.frame(names = names[nodesEvents], with = with, without = without)
  return(df)
}
propagationCov <- function(covmatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- cormatrix
  # nodesEvents <- 1:648
  # valueEvent <- 1
  # nodesEvidence <- c(81,280)
  # valueEvidence <- c(1,1)
  m <- nrow(covmatrix)
  names <- c()
  names
  for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(V", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V >", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
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
    
    
    probname <- paste0("P(V >", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  upper <- rep(Inf, (m-length(nodesEvidence)))
  upperw <- rep(Inf, m)
  
  for(i in 1:(m-length(nodesEvidence))) {
    lower2 <- rep(-Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    lower2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower2, upper=upper, mean = rep(0, m), covmatrix,
                        dependent.ind = dependent, given.ind = nodesEvidence, X.given = valueEvidence,
                        check.sigma= FALSE)
  }
  with[nodesEvidence] <- 1
  
  for(i in 1:m){
    lower2w <- rep(-Inf,m)
    lower2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), sigma = covmatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V >", valueEvent,")")
  df <- data.frame(names = names[nodesEvents], with = with, without = without)
  return(df)
}


#######################################################################################
# Genes
#######################################################################################
PropagationExactGenesPerm <- function(baysnet, nodesEvents, valueEvent, nodesEvidence, valueEvidence, perm){
  
  # baysnet = hc_gene_fit
  # nodesEvents = nodes(hc_gene_fit)
  # valueEvent = ">=1"
  # nodesEvidence = "fadR"
  # valueEvidence = c(2)
  # perm = NULL
  # 
  if (is.null(perm)) {nodesEventsRef <- nodesEvents} else {
    nodesEventsRef <- c()
    for (i in 1:length(nodesEvents)){
      nodesEventsRef[i] <- which(perm == nodesEvents[i])
    }
  }
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
  }
  
  
  if (length(nodesEvidence) > 1){
    proves <- c()
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0(nodesEvidence[j]," = ", valueEvidence[j])
    }
    
    text <- "list("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],",")
    }
    text <- paste0(text,proves[length(nodesEvidence)],")")
    str2 <- text
    
    
    probname <- paste0("P(", valueEvent,"|")
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
    str <- paste0("(", l, valueEvent, ")")
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
  attr(without, "probability") <- paste0("P(V", valueEvent,")")
  df <- data.frame(names = nodesEventsRef, with = with, without = without)
  return(df)
  
}

# df[which(df$names =="fadE"),,]
propagationCorrGenes <- function(cormatrix, nodesEvents, valueEvent, nodesEvidence, valueEvidence){
  
  # cormatrix <- covariation
  # nodesEvents <- colnames(covariation)
  # valueEvent <- 1
  # nodesEvidence <- "fadR"
  # valueEvidence <- c(1)
  # m <- nrow(cormatrix)
  
  
  ind.nodesEvidence <- which(nodesEvents == nodesEvidence)

  # names <- c()
  names <- nodesEvents
  # for(i in 1:length(nodesEvents)) names <- append(names,paste0("V",nodesEvents[i]))
  names
  
  dependent <- 1:m
  dependent <- dependent[-c(ind.nodesEvidence)]
  
  with <- numeric(length = length(nodesEvents))
  without <- numeric(length = length(nodesEvents))
  
  if (length(nodesEvidence) == 1) {
    str2 <- paste0("list(", nodesEvidence[1]," = ", valueEvidence[1], ")")
    probname <- paste0("P(V >", valueEvent,"|",nodesEvidence[1]," = ", valueEvidence[1], ")")
  }
  
  
  if (length(nodesEvidence) > 1){
    proves <- c()
    for (j in 1:length(nodesEvidence)){
      proves[j]<- paste0(nodesEvidence[j]," = ", valueEvidence[j])
    }
    
    text <- "list("
    for (j in 1:(length(nodesEvidence)-1)){
      text <- paste0(text,proves[j],",")
    }
    text <- paste0(text,proves[length(nodesEvidence)],")")
    str2 <- text
    
    
    probname <- paste0("P(V >", valueEvent,"|")
    for (j in 1:(length(nodesEvidence)-1)){
      probname <- paste0(probname,proves[j],",")
    }
    probname <- paste0(probname,proves[length(nodesEvidence)],")")
  }
  
  upper <- rep(Inf, (m-length(nodesEvidence)))
  upperw <- rep(Inf, m)
  
  i <- 1
  for(i in 1:(m-length(nodesEvidence))) {
    lower2 <- rep(-Inf,m-length(nodesEvidence))
    l <- dependent[i]
    l
    lower2[c(i)] <- valueEvent
    with[l] <- pcmvnorm(lower= lower2, upper=upper, rep(0, m), cormatrix,
                        dependent.ind = dependent, given.ind = ind.nodesEvidence, X.given = valueEvidence,
                        check.sigma= FALSE)
  }
  
  for(i in 1:length(valueEvidence)){
    if(valueEvidence[i]>=0){
      with[nodesEvidence[i]] <- 1
    } else {with[nodesEvidence[i]] <- 0}
  }
  
  for(i in 1:m){
    lower2w <- rep(-Inf,m)
    lower2w[c(i)] <- valueEvent
    without[i] <- pmvnorm(lower = lower2w, upper = upperw, mean = rep(0, m), sigma = cormatrix)
  }
  
  attr(with, "probability") <- probname
  attr(without, "probability") <- paste0("P(V >", valueEvent,")")
  df <- data.frame(names = nodesEvents, with = with, without = without)
  return(df)
}
