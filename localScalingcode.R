#localScaling

localScaling <- function(grid,
                         base = NULL,
                         ref = NULL,
                         clim.fun = list(FUN = "mean", na.rm = TRUE),
                         by.member = TRUE,
                         time.frame = c("none", "monthly", "daily"),
                         parallel = FALSE,
                         max.ncores = 16,
                         ncores = NULL) {
  time.frame <- match.arg(time.frame, choices = c("none", "monthly", "daily"))
  if (time.frame == "none") {
    message("[", Sys.time(), "] - Scaling ...")
    out <- .localScaling(grid, base, ref, clim.fun, by.member, parallel, max.ncores, ncores)
    message("[", Sys.time(), "] - Done")
  } else if (time.frame == "monthly") {
    message("[", Sys.time(), "] - Scaling by months ...")
    months <- getSeason(grid)
    aux.list <- lapply(1:length(months), function(x) {
      grid1 <- subsetGrid(grid, season = months[x])
      if (!is.null(base)) {
        base1 <- subsetGrid(base, season = months[x])
      } else {
        base1 <- NULL
      }
      if (!is.null(ref)) {
        ref1 <- subsetGrid(ref, season = months[x])
      } else {
        ref1 <- NULL
      }
      .localScaling(grid1, base1, ref1, clim.fun, by.member, parallel, max.ncores, ncores)
    })
    out <- do.call("bindGrid.time", aux.list)
    message("[", Sys.time(), "] - Done")
  } else if (time.frame == "daily") {
    doys.grid <- grid %>% getRefDates() %>% substr(6,10) 
    doys.grid <- gsub("02-29", "02-28", doys.grid)
    if (!is.null(base)) {
      doys.base <- base %>% getRefDates() %>% substr(6, 10)
      doys.base <- gsub("02-29", "02-28", doys.base)
    }
    if (!is.null(ref)) {
      doys.ref <- ref %>% getRefDates() %>% substr(6, 10)
      doys.ref <- gsub("02-29", "02-28", doys.ref)
    }
    message("[", Sys.time(), "] - Scaling by julian days ...")
    aux.list <- lapply(unique(doys.grid), function(x) {
      grid1 <- subsetDimension(grid, dimension = "time", indices = which(doys.grid == x))
      if (!is.null(base)) {
        base1 <- subsetDimension(base, dimension = "time", indices = which(doys.base == x))
      } else {
        base1 <- base
      }
      if (!is.null(ref)) {
        ref1 <- subsetDimension(ref, dimension = "time", indices = which(doys.ref == x))
      } else {
        ref1 <- ref
      }
      .localScaling(grid1, base1, ref1, clim.fun, by.member, parallel, max.ncores, ncores)
    })
    out <- do.call("bindGrid.time", aux.list)
    message("[", Sys.time(), "] - Done")
  }
  invisible(out)
}

.localScaling <- function(grid, base, ref, clim.fun, by.member, parallel, max.ncores, ncores) {
  grid <- redim(grid)
  if (is.null(base)) {
    base <- suppressMessages({
      climatology(grid, clim.fun, by.member, parallel, max.ncores, ncores)
    }) %>% redim()
  } else {
    checkSeason(grid, base)
    checkDim(grid, base, dimensions = c("lat", "lon"))
    base <- suppressMessages({
      climatology(base, clim.fun, by.member, parallel, max.ncores, ncores)
    }) %>% redim()
  }
  if (!is.null(ref)) {
    checkDim(grid, ref, dimensions = c("lat", "lon"))
    checkSeason(grid, ref)
    ref <- suppressMessages({
      climatology(ref, clim.fun, by.member, parallel, max.ncores,ncores)
    }) %>% redim()
  } else {
    ref <- list()
    ref[["Data"]] <- array(0, getShape(base))
    attr(ref[["Data"]], "dimensions") <- getDim(base)
  }    
  parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
  lapply_fun <- selectPar.pplyFun(parallel.pars, .pplyFUN = "lapply")
  if (parallel.pars$hasparallel) on.exit(parallel::stopCluster(parallel.pars$cl))
  clim <- grid[["Data"]]
  dimNames <- getDim(grid)
  ind.time <- grep("^time", dimNames)
  n.times <- getShape(grid, "time")
  Xc <- base[["Data"]]
  Xref <- ref[["Data"]]
  aux.list <- lapply_fun(1:n.times, function(x) {
    X <- asub(clim, idx = x, dims = ind.time, drop = FALSE)
    out <- if (dim(X)[1] != dim(Xc)[1]) {
      aux <- lapply(1:dim(X)[1], function(i) X[i, , , , drop = FALSE] - Xc + Xref)
      do.call("abind", c(aux, along = 1)) %>% unname()
    } else {
      X - Xc + Xref
    }
    return(out)
  })
  X <- Xc <- Xref <- base <- ref <- NULL
  grid[["Data"]] <- do.call("abind", c(aux.list, along = ind.time)) %>% unname()
  aux.list <- NULL
  attr(grid[["Data"]], "dimensions") <- dimNames
  return(grid)
}
