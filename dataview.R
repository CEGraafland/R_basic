install.packages(pkgs = "/home/catharina/Documents/Installations/modified_4.3.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)


load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/HC/hc_2.rda")
str(hc_2)     
hc_2$networkdata
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_2.rda")
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/mmhc/mmhc_3.rda")
str(mmhc_2)
mmhc_2$networks$
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/PCalg/pc_2.rda")
attributes(pc_2$networks$pc_10d_0.2)
360/pc_2$networkdata$params
pc_2$networkdata
c(2:10 %o% 10^(-112:-2))[841:999]
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/Tabu/tabu_2.rda")
tabu_2$networks$tabu_10d_505i
hc_2$networks$tabu_10d_505i
all.equal(tabu_2$networks$tabu_10d_505i,hc_2$networks$tabu_10d_505i)
360/mmhc_2$networkdata$params

mmhc_2$begin
mmhc_2$netwok
360/mmhc_3$networkdata$params

test <- pc.stable(gaussian.test, test ="bic-gt", B = 30)
test
net
cextend(net)

library(transformeR)
library(magrittr)
library(pcalg)
library(igraph)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions.R")

data <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
stabletest1 <- pc.stable(data, test ="bic-gt", B = 100)
stabletest2 <- pc.stable(data, test ="bic-gt", B = 0.5)
stabletest3 <- pc.stable(data, test ="bic-gt", B = 0.75)
stabletest4
cextend(stabletest2) # no
acyclic(stabletest3, directed = FALSE)
cextend(stabletest3) # no
save(stabletest1, stabletest2, stabletest3,
     file = "/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/tests/stabletests.rda")
