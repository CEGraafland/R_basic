##########################################################################
# Store perhaps obtained iterations tabu_2_eBIC_g
##########################################################################
# for(j in c(2)){
#   pattern <- paste0("tabu_",j,"_eBIC_g")
#   filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", full.names = T, pattern = pattern)
#   filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", pattern = pattern)
#   filestabu1names <- gsub(".rda", "", filestabu1names)
#   filestabu1names
# 
#   tabu2eBICiterations <- list(it.networks = list())
#   names <- c()
# 
#   for (i in 1:length(filestabu1)){
#     variablepos <- get(load(filestabu1[i]))
#     tabu2eBICiterations$it.networks[[i]] <- variablepos
#     names(tabu2eBICiterations$it.networks)
#     names[i] <- filestabu1names[i]
#   }
# 
#   names(tabu2eBICiterations$it.networks)<- names
# }
# 
# save(tabu2eBICiterations, file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu2eBICiterations.rda")
##########################################################################
# Load tabu eBIC
##########################################################################
rm(list = ls())

for(j in c(1,2,3,4,5)){
  pattern <- paste0("tabu_",j,"_eBIC_g")
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1$networks[[i]] <- variablepos$networks[[1]]
    names[i] <- names(variablepos$networks)
    data.frames[i,] <- as.matrix(variablepos$networkdata)
  }
  colnames(data.frames) <- names(variablepos$networkdata)
  data.frames <- as.data.frame(data.frames)
  variablelisttabu1$networkdata<- data.frames
  variablelisttabu1$networkdata
  
  names(variablelisttabu1$networks) <- names
  assign(paste0("tabu_",j,"_eBIC"),variablelisttabu1)
}

##########################################################################
# Load HC eBIC
##########################################################################

for(j in c(1,2,3,4,5)){
  pattern <- paste0("hc_",j,"_eBIC_g")
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(),networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))
  
  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1$networks[[i]] <- variablepos$networks[[1]]
    names[i] <- names(variablepos$networks)
    data.frames[i,] <- as.matrix(variablepos$networkdata)
  }
  colnames(data.frames) <- names(variablepos$networkdata)
  data.frames <- as.data.frame(data.frames)
  variablelisttabu1$networkdata<- data.frames
  variablelisttabu1$networkdata
  
  names(variablelisttabu1$networks) <- names
  assign(paste0("hc_",j,"_eBIC"),variablelisttabu1)
}


##########################################################################
# Load mmhc eBIC
##########################################################################

for(j in c(1,2,3,4,5)){
  pattern <- paste0("mmhc_",j,"_eBIC_g")
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc", full.names = T, pattern = pattern)
  filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc", pattern = pattern)
  filestabu1names <- gsub(".rda", "", filestabu1names)
  filestabu1names
  
  variablelisttabu1 <- list(networks = list(), begin = list(), networkdata = data.frame())
  names <- c()
  data.frames <- matrix(ncol = 6, nrow = length(filestabu1))

  for (i in 1:length(filestabu1)){
    variablepos <- get(load(filestabu1[i]))
    variablelisttabu1$networks[[i]] <- variablepos$networks[[1]]
    names[i] <- names(variablepos$networks)
    variablelisttabu1$begin[[i]] <- variablepos$begin
    data.frames[i,] <- as.matrix(variablepos$networkdata)
  }
  colnames(data.frames) <- names(variablepos$networkdata)
  data.frames <- as.data.frame(data.frames)
  variablelisttabu1$networkdata<- data.frames
  variablelisttabu1$networkdata
  
  names(variablelisttabu1$networks) <- names
  assign(paste0("mmhc_",j,"_eBIC"),variablelisttabu1)
}



testssum <- function(x){
  f1 <- function(x, n){ sum(head(x,n))}
  sums <- c()
  for (i in 1:length(x)){
    sums[i] <- f1(x,i)
  }
  return(sums)
}

# nets <- mmhc_2
# attributes(nets$networks$pc_10d_0.2)
# mmhc_2$networkdata$testsplus <- mmhc_2$networkdata$tests +mmhc_2$begin$learning$ntests
# check <- function(x){networkdata <- x$networkdata; return(networkdata)}
# check(pc_2)
##########################################################################
# Load PCstable eBIC
##########################################################################

for(j in c(1,2,3,4,5)){
  pattern <- paste0("pcstable2_",j)
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable", full.names = T, pattern = pattern)
  # filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable", pattern = pattern)
  # filestabu1names <- gsub(".rda", "", filestabu1names)
  # filestabu1names
  
  load(filestabu1)
}

##########################################################################
# Load gs eBIC
##########################################################################

for(j in c(1,2,3,4,5)){
  pattern <- paste0("gs_",j)
  filestabu1 <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs", full.names = T, pattern = pattern)
  # filestabu1names <- list.files("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/gs", pattern = pattern)
  # filestabu1names <- gsub(".rda", "", filestabu1names)
  # filestabu1names
  filestabu1
  load(filestabu1)
}

####################################################################################
# functions
####################################################################################

ejes_gamma_tests <- function(nets, begin = NULL){
  networkdata <- nets$networkdata
  network <- nets$networks
  if (!is.null(begin)){
    begin <- nets$begin
    begintests<- sapply(begin, ntests)}
 
  xeje <- (nets$networkdata$gamma)*2
  miss <- !is.na(xeje)
  if (length(networkdata) == 5){
    sumtests <- networkdata$tests
    yeje <- log10(sumtests)
  } else if (length(nets) == 3){
    sumtests <- networkdata$tests
    sumtests <- sumtests + begintests
    yeje <- log10(sumtests)
  } else {yeje <- log10(networkdata$tests)}
  return(list(xeje,yeje))
}


ejes_gamma_log <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- (nets$networkdata$gamma)*2
  miss <- !is.na(xeje)
  if (length(networkdata) == 5){
    logliks <- networkdata$logliks
    yeje <- logliks
  } else if (length(nets) == 3){
    logliks <- networkdata$logliks
    yeje <- logliks
  } else {yeje <- networkdata$logliks}
  return(list(xeje,yeje))
}

ejes_gamma_edges <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- (nets$networkdata$gamma)*2
  miss <- !is.na(xeje)
  edges <- networkdata$edges
  yeje <- edges
  return(list(xeje,yeje))
}

ejes_gamma_size <- function(nets){
  networkdata <- nets$networkdata
  network <- nets$networks
  xeje <- (nets$networkdata$gamma)*2
  miss <- !is.na(xeje)
  edges <- 360/networkdata$params
  yeje <- edges
  return(list(xeje,yeje))
}

# plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/testsandloglik.pdf")
# pdf(plotname, width = 10, height = 5) 
# par(mfrow = c(1,2))
gammas <- c(0.25,0.5,1,5,6,7,8,9,10,25)
ejes_gamma_tests_tabu <- function(nets){
  len <- length(nets$networkdata$tests)
  y <- log10(nets$networkdata$tests[len])
  return(y)}
ejes_gamma_loglik_tabu <- function(nets){
  len <- length(nets$networkdata$loglik)
  y <- (nets$networkdata$loglik[len])
  return(y)}

hc_1_eBIC_g0$networkdata
hc_2_eBIC_g0$networkdata
tabu_1_eBIC_g0$networkdata
mmhc_1_eBIC_g0$networkdata

log10(mmhc_1_eBIC_g0.1$networkdata$tests + mmhc_1_eBIC_g0.1$begin$learning$ntests)
log10(mmhc_2_eBIC_g0.1$networkdata$tests + mmhc_2_eBIC_g0.1$begin$learning$ntests)
log10(mmhc_3_eBIC_g0.1$networkdata$tests + mmhc_3_eBIC_g0.1$begin$learning$ntests)
log10(mmhc_4_eBIC_g0.1$networkdata$tests + mmhc_4_eBIC_g0.1$begin$learning$ntests)
log10(mmhc_5_eBIC_g0.1$networkdata$tests + mmhc_5_eBIC_g0.1$begin$learning$ntests)

log10(mmhc_1_eBIC_g0$networkdata$tests + mmhc_1_eBIC_g0$begin$learning$ntests)
log10(mmhc_2_eBIC_g0$networkdata$tests + mmhc_2_eBIC_g0$begin$learning$ntests)
log10(mmhc_3_eBIC_g0$networkdata$tests + mmhc_3_eBIC_g0$begin$learning$ntests)
log10(mmhc_4_eBIC_g0$networkdata$tests + mmhc_4_eBIC_g0$begin$learning$ntests)
log10(mmhc_5_eBIC_g0$networkdata$tests + mmhc_5_eBIC_g0$begin$learning$ntests)

log10(mmhc_1_eBIC_g0.25$networkdata$tests + mmhc_1_eBIC_g0.25$begin$learning$ntests)
log10(mmhc_2_eBIC_g0.25$networkdata$tests + mmhc_2_eBIC_g0.25$begin$learning$ntests)
log10(mmhc_3_eBIC_g0.25$networkdata$tests + mmhc_3_eBIC_g0.25$begin$learning$ntests)
log10(mmhc_4_eBIC_g0.25$networkdata$tests + mmhc_4_eBIC_g0.25$begin$learning$ntests)
log10(mmhc_5_eBIC_g0.25$networkdata$tests + mmhc_5_eBIC_g0.25$begin$learning$ntests)

log10(hc_1_eBIC_g0$networkdata$tests)
log10(hc_2_eBIC_g0$networkdata$tests)
log10(hc_3_eBIC_g0$networkdata$tests)
log10(hc_4_eBIC_g0$networkdata$tests)
log10(hc_5_eBIC_g0$networkdata$tests)

log10(tabu_1_eBIC_g0$networkdata$tests)
log10(tabu_2_eBIC_g0$networkdata$tests)
log10(tabu_3_eBIC_g0$networkdata$tests)
log10(tabu_4_eBIC_g0$networkdata$tests)
log10(tabu_5_eBIC_g0$networkdata$tests)

hc_1_eBIC_g0$networkdata$logliks
hc_2_eBIC_g0$networkdata$logliks
hc_3_eBIC_g0$networkdata$logliks
hc_4_eBIC_g0$networkdata$logliks
hc_5_eBIC_g0$networkdata$logliks

tabu_1_eBIC_g0$networkdata$logliks
tabu_2_eBIC_g0$networkdata$logliks
tabu_3_eBIC_g0$networkdata$logliks
tabu_4_eBIC_g0$networkdata$logliks
tabu_5_eBIC_g0$networkdata$logliks

hc_1_eBIC_g0$networkdata$edges
hc_2_eBIC_g0$networkdata$edges
hc_3_eBIC_g0$networkdata$edges
hc_4_eBIC_g0$networkdata$edges
hc_5_eBIC_g0$networkdata$edges

hc_1_eBIC_g0.1$networkdata$edges
hc_2_eBIC_g0.1$networkdata$edges
hc_3_eBIC_g0.1$networkdata$edges
hc_4_eBIC_g0.1$networkdata$edges
hc_5_eBIC_g0.1$networkdata$edges

tabu_1_eBIC_g0$networkdata$edges
tabu_2_eBIC_g0$networkdata$edges
tabu_3_eBIC_g0$networkdata$edges
tabu_4_eBIC_g0$networkdata$edges
tabu_5_eBIC_g0$networkdata$edges

tabu_1_eBIC_g0.1$networkdata$edges
tabu_2_eBIC_g0.1$networkdata$edges
tabu_3_eBIC_g0.1$networkdata$edges
tabu_4_eBIC_g0.1$networkdata$edges
tabu_5_eBIC_g0.1$networkdata$edges


tabu_1_eBIC_g0$networkdata$logliks
tabu_2_eBIC_g0$networkdata$logliks
tabu_3_eBIC_g0$networkdata$logliks
tabu_4_eBIC_g0$networkdata$logliks
tabu_5_eBIC_g0$networkdata$logliks


log10(hc_1_eBIC_g0$networkdata$tests)
mmhc_1_eBIC$begin[[1]]
log10(tabu_1_eBIC_g0.1$networkdata$tests)
a <- mmhc_1_eBIC_g0.1$networkdata$tests
b <- mmhc_1_eBIC_g0.1$begin$learning$ntests
log10(a + b)
c <- sum(mmhc_1$networkdata$tests[1:length(mmhc_1$networkdata$tests)])
d <- mmhc_1$begin$learning$ntests
log10(c+d)
log10(sum(tabu_1$networkdata$tests[1:7999]))
#####################################################################################
# high gamma plot tests
######################################################################################
# ejesmmhc <- lapply(list(mmhc_1, mmhc_2, mmhc_3, mmhc_4, mmhc_5),ejes)
# ejestabu <- lapply(list(tabu_2, tabu_3, tabu_4, tabu_5),ejes)
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/testsvsgammaQuad.pdf")
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/testsvsgammaQuad2.pdf")
pdf(plotname, height = 5, width = 15)
par(xpd = FALSE)
par(mfrow = c(1,3))
par(mar = c(2, 4, 2, 2))
dev.off()

tabuBIC <- list(tabu_2_eBIC_g0.25, tabu_2_eBIC_g0.5, tabu_2_eBIC_g1, tabu_2_eBIC_g5,tabu_2_eBIC_g6,tabu_2_eBIC_g7,tabu_2_eBIC_g8,tabu_2_eBIC_g9,tabu_2_eBIC_g10,tabu_2_eBIC_g25)
ejestabuBIC <- lapply(list(tabu_1_eBIC,tabu_2_eBIC,tabu_3_eBIC,tabu_4_eBIC,tabu_5_eBIC),ejes_gamma_tests)
ejeshcBIC <- lapply(list(hc_1_eBIC, hc_2_eBIC, hc_3_eBIC,hc_4_eBIC,hc_5_eBIC),ejes_gamma_tests)
ejesmmhcBIC <- lapply(list(mmhc_1_eBIC,mmhc_2_eBIC,mmhc_3_eBIC,mmhc_4_eBIC,mmhc_5_eBIC), ejes_gamma_tests,1)
ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4,pcstable2_5),ejes_gamma_tests)
ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejes_gamma_tests)


plot(x = ejestabuBIC[[5]][[1]], y = ejestabuBIC[[5]][[2]], type = "p",
     pch = 23,
     bg = rgb(255,0,0, alpha = 255, maxColorValue = 255),
     main = "(a) Speed",
     xlab = expression(paste(gamma)), 
     ylab = "log10(calls to the statistical criterion)",
     xlim = c(0,50),
     ylim = c(5,7.1),
     yaxt = "n")
axis(side=2, at=c(5,5.5,6,6.5,7,7.5,8), labels=c(5,5.5,6,6.5,7,7.5,8))


points(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))

points(ejeshcBIC[[5]][[1]], ejeshcBIC[[5]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[2]][[1]], ejeshcBIC[[2]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[3]][[1]], ejeshcBIC[[3]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[4]][[1]], ejeshcBIC[[4]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[1]][[1]], ejeshcBIC[[1]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))

points(ejesmmhcBIC[[5]][[1]], ejesmmhcBIC[[5]][[2]], pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[2]][[1]], ejesmmhcBIC[[2]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[3]][[1]], ejesmmhcBIC[[3]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[4]][[1]], ejesmmhcBIC[[4]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[1]][[1]], ejesmmhcBIC[[1]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
col2rgb("navyblue")
points(ejespcstable[[5]][[1]], ejespcstable[[5]][[2]], pch = 23, bg = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))

points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], pch = 23, bg = rgb(0,0,255, alpha = 255, maxColorValue = 255))
points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))

par(xpd=T)
rect(-1, 7.00, 51, 7.15, col="white", border="NA")
rect(-3, 7, 0.2, 7.05, col="white", border="white")
lines(x=c(-3,0.2), y=c(7, 7))
lines(x=c(-3,0.2), y=c(7.05, 7.05))

points(x = c(0), y = c(7.15), pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(x = c(0.2), y = c(7.1), pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(x = c(0.5), y = c(7.08), pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
#########################################################################################
# low gamma plot tests
#########################################################################################
# plot(x = ejestabuBIC[[1]][[1]], y = ejestabuBIC[[1]][[2]], type = "p",
#      pch = 23,
#      bg = rgb(255,0,0, alpha = 255, maxColorValue = 255),
#      main = "Speed",
#      xlab = expression(paste(gamma)), 
#      ylab = "log10(calls to the statistical criterion)",
#      xlim = c(0,1),
#      ylim = c(6,8.5))
# points(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(ejestabuBIC[[5]][[1]], ejestabuBIC[[5]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# 
# points(ejeshcBIC[[1]][[1]], ejeshcBIC[[1]][[2]], pch = 23, bg = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[2]][[1]], ejeshcBIC[[2]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[3]][[1]], ejeshcBIC[[3]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[4]][[1]], ejeshcBIC[[4]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[5]][[1]], ejeshcBIC[[5]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# 
# points(ejesmmhcBIC[[1]][[1]], ejesmmhcBIC[[1]][[2]], pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[2]][[1]], ejesmmhcBIC[[2]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[3]][[1]], ejesmmhcBIC[[3]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[4]][[1]], ejesmmhcBIC[[4]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[5]][[1]], ejesmmhcBIC[[5]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# col2rgb("navyblue")
# points(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], pch = 23, bg = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[5]][[1]], ejespcstable[[5]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# 
# points(ejesgs[[1]][[1]], ejesgs[[1]][[2]], pch = 23, bg = rgb(0,0,255, alpha = 255, maxColorValue = 255))
# points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
# points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
# points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
# points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
# 
# 
# 
# test <- mmhc(datapermutations[[1]],restrict.args = list(test = "bic-gt", B = 5), maximize.args = list(score = "bic-g", k = 5))
# test$learning$ntests
# mmhc_1_eBIC_g5$networkdata$tests + mmhc_1_eBIC_g5$begin$learning$ntests
# mmhc_1_eBIC_g0.1
# 55
############################################################################################
# high gamma plot logs 
############################################################################################

tabuBIC <- list(tabu_2_eBIC_g0.25, tabu_2_eBIC_g0.5, tabu_2_eBIC_g1, tabu_2_eBIC_g5,tabu_2_eBIC_g6,tabu_2_eBIC_g7,tabu_2_eBIC_g8,tabu_2_eBIC_g9,tabu_2_eBIC_g10,tabu_2_eBIC_g25)
ejestabuBIC <- lapply(list(tabu_1_eBIC,tabu_2_eBIC,tabu_3_eBIC,tabu_4_eBIC,tabu_5_eBIC),ejes_gamma_log)
ejeshcBIC <- lapply(list(hc_1_eBIC, hc_2_eBIC, hc_3_eBIC,hc_4_eBIC,hc_5_eBIC),ejes_gamma_log)
ejesmmhcBIC <- lapply(list(mmhc_1_eBIC,mmhc_2_eBIC, mmhc_3_eBIC,mmhc_4_eBIC,mmhc_5_eBIC), ejes_gamma_log)
ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4,pcstable2_5),ejes_gamma_log)
ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejes_gamma_log)

par(xpd=F)
plot(ejestabuBIC[[5]][[1]], ejestabuBIC[[5]][[2]], type = "p",
     #lwd = 1,
     pch = 23,
     bg = rgb(255,0,0, alpha = 255, maxColorValue = 255),
     main = "(b) Score (log-likelihood)",
     xlab = expression(paste(gamma)), 
     ylab = expression(paste("Log P(G|(x,",theta,")")),
     xlim = c(0,50),
     ylim = c(-350000,-140000)
)

points(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))

points(ejeshcBIC[[5]][[1]], ejeshcBIC[[5]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[2]][[1]], ejeshcBIC[[2]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[3]][[1]], ejeshcBIC[[3]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[4]][[1]], ejeshcBIC[[4]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[1]][[1]], ejeshcBIC[[1]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))

points(ejesmmhcBIC[[5]][[1]], ejesmmhcBIC[[5]][[2]], pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[2]][[1]], ejesmmhcBIC[[2]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[3]][[1]], ejesmmhcBIC[[3]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[4]][[1]], ejesmmhcBIC[[4]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[1]][[1]], ejesmmhcBIC[[1]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
col2rgb("navyblue")
points(ejespcstable[[5]][[1]], ejespcstable[[5]][[2]], pch = 23, bg = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))

points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], pch = 23, bg = rgb(0,0,255, alpha = 255, maxColorValue = 255))
points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))

# par(xpd=T)
# load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/utilnetworks_hcnetworks10d.rda")
# y <- hc_edges_loglik_10d_8600_8695i$networkdata$logliks[length(hc_edges_loglik_10d_8600_8695i$networkdata$logliks)]
# points(0,y, pch = 23, bg = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(0,y, col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
############################################################################################
# low gamma plot logs 
############################################################################################

# # ejesmmhc <- lapply(list(mmhc_1, mmhc_2, mmhc_3, mmhc_4, mmhc_5),ejes)
# # ejestabu <- lapply(list(tabu_2, tabu_3, tabu_4, tabu_5),ejes)
# ejestabuBIC <- lapply(list(tabu_1_eBIC,tabu_2_eBIC,tabu_3_eBIC,tabu_4_eBIC,tabu_5_eBIC),ejes_gamma_log)
# ejeshcBIC <- lapply(list(hc_1_eBIC, hc_2_eBIC, hc_3_eBIC,hc_4_eBIC,hc_5_eBIC),ejes_gamma_log)
# ejesmmhcBIC <- lapply(list(mmhc_1_eBIC,mmhc_3_eBIC,mmhc_4_eBIC,mmhc_5_eBIC), ejes_gamma_log)
# ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4),ejes_gamma_log)
# ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejes_gamma_log)
# 
# plot(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], type = "p",
#      #lwd = 1,
#      pch = 23,
#      bg = rgb(255,0,0, alpha = 255, maxColorValue = 255),
#      main = "logliklihood vs sample-parameter ratio",
#      xlab = expression(paste(gamma)), 
#      ylab = expression(paste("P(G|(x,",theta,")")),
#      xlim = c(0,1),
#      ylim = c(-230000,-150000)
# )
# 
# points(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# points(ejestabuBIC[[5]][[1]], ejestabuBIC[[5]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
# 
# points(ejeshcBIC[[1]][[1]], ejeshcBIC[[1]][[2]], pch = 23, bg = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[2]][[1]], ejeshcBIC[[2]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[3]][[1]], ejeshcBIC[[3]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[4]][[1]], ejeshcBIC[[4]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# points(ejeshcBIC[[5]][[1]], ejeshcBIC[[5]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
# 
# points(ejesmmhcBIC[[1]][[1]], ejesmmhcBIC[[1]][[2]], pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[2]][[1]], ejesmmhcBIC[[2]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[3]][[1]], ejesmmhcBIC[[3]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[4]][[1]], ejesmmhcBIC[[4]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# points(ejesmmhcBIC[[5]][[1]], ejesmmhcBIC[[5]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
# 
# col2rgb("navyblue")
# points(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], pch = 23, bg = rgb(0,0,128, alpha = 125, maxColorValue = 255))
# 
# points(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 100, maxColorValue = 255))
# points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
# points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
# points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 20, maxColorValue = 255))
# points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], pch = 23, bg = rgb(0,0,255, alpha = 255, maxColorValue = 255))

############################################################################################
# gamma plot edges
############################################################################################

# ejesmmhc <- lapply(list(mmhc_1, mmhc_2, mmhc_3, mmhc_4, mmhc_5),ejes)
# ejestabu <- lapply(list(tabu_2, tabu_3, tabu_4, tabu_5),ejes)
ejestabuBIC <- lapply(list(tabu_1_eBIC,tabu_2_eBIC,tabu_3_eBIC,tabu_4_eBIC,tabu_5_eBIC),ejes_gamma_edges)
ejeshcBIC <- lapply(list(hc_1_eBIC, hc_2_eBIC, hc_3_eBIC,hc_4_eBIC,hc_5_eBIC),ejes_gamma_edges)
ejesmmhcBIC <- lapply(list(mmhc_1_eBIC,mmhc_2_eBIC,mmhc_3_eBIC,mmhc_4_eBIC,mmhc_5_eBIC), ejes_gamma_edges)
ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4,pcstable2_5),ejes_gamma_edges)
ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejes_gamma_edges)

par(xpd=F)
plot(ejestabuBIC[[5]][[1]], ejestabuBIC[[5]][[2]], type = "p",
     #lwd = 1,
     pch = 23,
     bg = rgb(255,0,0, alpha = 255, maxColorValue = 255),
     main = "(c) Size",
     xlab = expression(paste(gamma)), 
     ylab = "|A| = number of arcs" ,
     # xlim = c(0,1),
     ylim = c(0,2600)
)

points(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))

points(ejeshcBIC[[5]][[1]], ejeshcBIC[[5]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[2]][[1]], ejeshcBIC[[2]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[3]][[1]], ejeshcBIC[[3]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[4]][[1]], ejeshcBIC[[4]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[1]][[1]], ejeshcBIC[[1]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))

points(ejesmmhcBIC[[5]][[1]], ejesmmhcBIC[[5]][[2]], pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[2]][[1]], ejesmmhcBIC[[2]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[3]][[1]], ejesmmhcBIC[[3]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[4]][[1]], ejesmmhcBIC[[4]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[1]][[1]], ejesmmhcBIC[[1]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
col2rgb("navyblue")
points(ejespcstable[[5]][[1]], ejespcstable[[5]][[2]], pch = 23, bg = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))

points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], pch = 23, bg = rgb(0,0,255, alpha = 255, maxColorValue = 255))
points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))

par(xpd=T)
rect(-1, 2500, 51, 2600, col="white", border="white")
rect(-3, 2500, 0.2, 2550, col="white", border="white")
lines(x=c(-3,0.2), y=c(2500, 2500))
lines(x=c(-3,0.2), y=c(2550, 2550))

points(x = c(0.2), y = c(2600),pch = 23,bg = rgb(255,0,0, alpha = 255, maxColorValue = 255))
points(x = c(0), y = c(2650),pch = 23,bg = rgb(255,0,0, alpha = 255, maxColorValue = 255))
points(x = c(0.2), y = c(2600), col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(x = c(0), y = c(2650), col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
############################################################################################
#  gamma plot size 
############################################################################################

# ejesmmhc <- lapply(list(mmhc_1, mmhc_2, mmhc_3, mmhc_4, mmhc_5),ejes)
# ejestabu <- lapply(list(tabu_2, tabu_3, tabu_4, tabu_5),ejes)
ejestabuBIC <- lapply(list(tabu_1_eBIC,tabu_2_eBIC,tabu_3_eBIC,tabu_4_eBIC,tabu_5_eBIC),ejes_gamma_size)
ejeshcBIC <- lapply(list(hc_1_eBIC, hc_2_eBIC, hc_3_eBIC,hc_4_eBIC,hc_5_eBIC),ejes_gamma_size)
ejesmmhcBIC <- lapply(list(mmhc_1_eBIC,mmhc_2_eBIC, mmhc_3_eBIC,mmhc_4_eBIC,mmhc_5_eBIC), ejes_gamma_size)
ejespcstable <- lapply(list(pcstable2_1, pcstable2_2, pcstable2_3,pcstable2_4,pcstable2_5),ejes_gamma_size)
ejesgs <- lapply(list(gs_1, gs_2, gs_3, gs_4, gs_5), ejes_gamma_size)

plot(ejestabuBIC[[5]][[1]], ejestabuBIC[[5]][[2]], type = "p",
     #lwd = 1,
     pch = 23,
     bg = rgb(255,0,0, alpha = 255, maxColorValue = 255),
     main = "Size",
     xlab = expression(paste(gamma)), 
     ylab = expression(paste("n/",theta))
     # xlim = c(0,1),
     # ylim = c(-230000,-150000)
)
points(ejestabuBIC[[2]][[1]], ejestabuBIC[[2]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[3]][[1]], ejestabuBIC[[3]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[4]][[1]], ejestabuBIC[[4]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))
points(ejestabuBIC[[1]][[1]], ejestabuBIC[[1]][[2]], col = rgb(255,0,0, alpha = 125, maxColorValue = 255))

points(ejeshcBIC[[5]][[1]], ejeshcBIC[[5]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[2]][[1]], ejeshcBIC[[2]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[3]][[1]], ejeshcBIC[[3]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[4]][[1]], ejeshcBIC[[4]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))
points(ejeshcBIC[[1]][[1]], ejeshcBIC[[1]][[2]], col = rgb(255,165,0, alpha = 255, maxColorValue = 255))

points(ejesmmhcBIC[[5]][[1]], ejesmmhcBIC[[5]][[2]], pch = 23, bg = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[2]][[1]], ejesmmhcBIC[[2]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[3]][[1]], ejesmmhcBIC[[3]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[4]][[1]], ejesmmhcBIC[[4]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
points(ejesmmhcBIC[[1]][[1]], ejesmmhcBIC[[1]][[2]], col = rgb(0,255,0, alpha = 255, maxColorValue = 255))
col2rgb("navyblue")
points(ejespcstable[[5]][[1]], ejespcstable[[5]][[2]], pch = 23, bg = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[2]][[1]], ejespcstable[[2]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[3]][[1]], ejespcstable[[3]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[4]][[1]], ejespcstable[[4]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))
points(ejespcstable[[1]][[1]], ejespcstable[[1]][[2]], col = rgb(0,0,128, alpha = 125, maxColorValue = 255))

points(ejesgs[[5]][[1]], ejesgs[[5]][[2]], pch = 23, bg = rgb(0,0,255, alpha = 255, maxColorValue = 255))
points(ejesgs[[2]][[1]], ejesgs[[2]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[3]][[1]], ejesgs[[3]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[4]][[1]], ejesgs[[4]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))
points(ejesgs[[1]][[1]], ejesgs[[1]][[2]], col = rgb(0,0,255, alpha = 125, maxColorValue = 255))

dev.off()
##############################################################################################
# plot mean
##############################################################################################
# method 1: nino indices
subset <- tas_ncep_10d$Data
time.coords.matrix <- array3Dto2Dmat(subset)
subtcm <- time.coords.matrix[ind1981_2010,]
meansubtcm <- colMeans(subtcm)
climnino <- quantity2clim(meansubtcm, ref.grid = tas_ncep_10d, what = "mean")
plotClimatology(climnino, backdrop.theme = "coastline", main = list("niño months",cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))
# method 2: warm tercile
subset <- tas_ncep_10d$Data
time.coords.matrix <- array3Dto2Dmat(subset)
subtcm <- time.coords.matrix[warmmonths,]
meansubtcm <- colMeans(subtcm)
climnino <- quantity2clim(meansubtcm, ref.grid = tas_ncep_10d, what = "mean")
plotClimatology(climnino, backdrop.theme = "coastline")

nosubtcm <- time.coords.matrix[no.ind1981_2010,]
nomeansubtcm <- colMeans(nosubtcm)
climnonino <- quantity2clim(nomeansubtcm, ref.grid = tas_ncep_10d, what = "mean")
plotClimatology(climnonino, backdrop.theme = "coastline", main = list("no niño months",cex = 0.5),colorkey = list(width = 0.6, lables = list(cex = 0.5)))

nino_tas_ncep_10d <- tas_ncep_10d
nino_tas_ncep_10d$Data <- tas_ncep_10d$Data[ind1981_2010,,]
nino_tas_ncep_10d$Dates <- tas_ncep_10d$Dates[ind1981_2010, outside = TRUE]
subsetGrid(tas_ncep_10d, season = ind1981_2010, outside = TRUE)
climnino <- climatology(nino_tas_ncep_10d)
clim <- climatology(tas_ncep_10d)
clim$Data
plotClimatology(clim, backdrop.theme = "coastline")
clim <- climatology()
########################################################################################
# Niño months other method 97-98
########################################################################################
grid <- tas_ncep_10d
seas <- getSeason(grid)
coords <- getCoordinates(grid)
x <- coords$x
y <- coords$y
ref.coords <- expand.grid(y, x)[2:1]
names(ref.coords) <- c("x", "y")
ref.dates <- getRefDates(grid)
seas.list <- lapply(1:length(seas), function(i) {
  subsetGrid(grid, season = seas[i]) %>% localScaling() %>% redim(drop = TRUE)
})
seas.list
grid <- NULL
aux <- do.call("bindGrid.time", seas.list) %>% redim(drop = TRUE)
aux <-subsetGrid(aux, years = 1998, season = 1)
climtot <- climatology(aux)




sub <- subsetGrid(tas_ncep_10d, years = 1998, season = 1)%>% localScaling() %>% redim(drop = TRUE)
climsub <- climatology(sub)
plotClimatology(climsub, backdrop.theme = "coastline")

totJan <- subsetGrid(tas_ncep_10d, season = 1)
climtotJan <- climatology(totJan)
plotClimatology(climtotJan, backdrop.theme = "coastline")

tot <- tas_ncep_10d
climtot <- climatology(tot)
medplot <- plotClimatology(climtot, backdrop.theme = "coastline", main = list("Mean monthly temperature 1981:2010", cex = 0.8),colorkey = list(width = 0.6, lables = list(cex = 0.6)))
medplot
dif <- quantity2clim(climsub$Data-climtotJan$Data, what = "difference",tas_ncep_10d)


library('RColorBrewer')
display.brewer.all()
col.blue <- rev(brewer.pal(8,"Blues"))
col.blue
col.red <- brewer.pal(8,"Reds")
col.b <- c(col.blue,col.red)
length(col.b)

difplot <- plotClimatology(dif,backdrop.theme = "coastline", set.max = 3,col.regions = col.b,  main = list("Anomaly January 1998", cex = 0.8),colorkey = list(width = 0.6, lables = list(cex = 0.6)), at = seq(-6,6,1))


plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/ninoyear.pdf")
pdf(plotname, width = 10, height = 4)

lay = cbind(c(1),c(2))
grid::gpar(mar = c(0.0, 0.0, 0.0, 0.0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0))
grid.arrange(grobs = list(medplot,difplot), layout_matrix = lay, gp =  gpar(mar = c(0.0, 0.0, 0.0, 0.0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0)))

medplot
difplot
dev.off()

pdf(plotname, width = 14, height = 8)
op <- par(no.readonly = TRUE)
lay <- rbind(c(1,0,0),
             c(2,3,0),
             c(4,5,6))
graphics::layout(as.matrix(lay))
par(mar = c(0.0, 0.0, 0.0, 0.0))
par(mgp = c(0, 0, 0))
par(oma = c(0, 0, 0, 0))
# for (i in 1:20) {
#   +     plot(diat.max[, i] ~ env[, "TP"], ann = F, cex = 0.6, cex.axis = 0.8,
#              +         tcl = -0.4, las = 1, pch = 19)
#   +     mtext(colnames(diat.max)[i], side = 3, line = 0.2, cex = 0.8)
#   + }
# mtext("Diatom distribution vs. log 10 Total Phosphorus", outer = TRUE,
#         +     side = 3, cex = 1.2, line = 1)


#######################################################################################
# long distance gamma plot
#######################################################################################
plotname <- paste0("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/Writing/paper/figures/comparegraphsdistances.pdf")
pdf(plotname, width = 14, height = 8)
op <- par(no.readonly = TRUE)
lay <- rbind(c(1,0,0),
             c(2,3,0),
             c(4,5,6))
graphics::layout(as.matrix(lay))
par(mar = c(0.0, 0, 1, 0.0))
par(mgp = c(0, 0, 0))
par(oma = c(0, 0, 0, 0))

#permutations
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pc_1_eBIC_g5.rda")
# best of pcstable
plot_long_distances(pc_1_eBIC_g5, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(0,0,255,alpha = 100, maxColorValue = 255), perm = permutations[[1]],title = FALSE)
title(mar = c(0,0,0,0), main = paste0("|E| = ",narcs(pc_1_eBIC_g5)))
# pcstablbest and best of mmhc
plot_long_distances(mmhc_1_eBIC_g5$networks$mmhc_10d_g5, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(0,255,0,alpha = 100, maxColorValue = 255) ,perm = permutations[[1]], title = FALSE)
title(mar = c(0,0,0,0), main = paste0("|E| = ",narcs(mmhc_1_eBIC_g5$networks$mmhc_10d_g5)))
plot_long_distances(mmhc_2$networks[[2]], data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(0,255,0,alpha = 100, maxColorValue = 255), perm = permutations[[2]],title = FALSE)
title(mar = c(0,0,0,0), main = paste0("|E| = ",narcs(mmhc_2$networks[[2]])))
# pcstablebest, mmhcbest, hcbien. 
plot_long_distances(tabu_1_eBIC_g5$networks$tabu_10d_g5, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(200,100,100, alpha = 100, maxColorValue = 255), perm = permutations[[1]],title = FALSE, remove = TRUE)
title(mar = c(0,0,0,0), main = paste0("|E| = ",narcs(tabu_1_eBIC_g5$networks$tabu_10d_g5)))
plot_long_distances(tabu_mmhc_simple2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(200,100,100, alpha = 50,maxColorValue = 255), perm = permutations[[2]], title = FALSE, remove = TRUE)
title(mar = c(0,0,0,0), main = paste0("|E| = ",narcs(tabu_mmhc_simple2)))
plot_long_distances(tabu_simple_1600_2, data.dag = TimeCoordsAnom_from_Grid_std(tas_ncep_10d),minimdist = 5000, smallcol = rgb(200,100,100, alpha = 50,maxColorValue = 255),perm = permutations[[2]], title = FALSE, remove = TRUE)
title(mar = c(0,0,0,0), main = paste0("|E| = ",narcs(tabu_simple_1600_2)))
par(op)
dev.off()


load(file ="/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/pcstable/pc_1_eBIC_g5.rda")
