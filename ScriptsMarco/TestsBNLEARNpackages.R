x <- hc_2$networkdata
x <- hc_1$networkdata
hc_2_eBIC_g0$networkdata
tail(x, 1)
hc_2$networks$tabu_10d_8000i$learning
data <- as.data.frame(TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = FALSE))
data[permutations[[5]]]
datapermutations[[]]

install.packages("bnlearn")
library(bnlearn)
test <- gamma_hc(data = data, iterations = 2, gamma = 0, start = hc_2$networks$tabu_10d_8000i)
logLik(hc_2$networks$tabu_10d_8000i,data)
test$networkdata
bnlearn::score(hc_2$networks$tabu_10d_8000i,data = data, type = "bic-g")

install.packages(pkgs = "/home/catharina/Documents/Installations/modified_bnlearn.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)


install.packages(pkgs = "/home/catharina/Documents/Installations/modified_4.3.tar.gz", repos = NULL,type ='source')
library(bnlearn, quietly = TRUE, warn.conflict = FALSE)
test2 <- gamma_hc(data = data, iterations = 2, gamma = 0, start = hc_2$networks$tabu_10d_8000i)

test2$networkdata
test$networkdata
all.equal(test,test2)


test3 <-  gamma_hc(data = data, iterations = 4, gamma = 0, start = hc_2$networks$tabu_10d_8000i)
test4 <-  gamma_hc(data = data[permutations[[5]]], iterations = 4, gamma = 0, start = hc_5$networks$hc_10d_8000i)
load("/run/user/1000/gvfs/sftp:host=ui01.macc.unican.es/vols/seal/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
test4$networkdata
score(data = data[permutations[[5]]],test4$networks$hc_10d_g0)
hc_5$networkdata
x5 <- hc_5$networkdata
tail(x5, 1)
hc_5$networkdata
gamma_hc(data = data[permutations[[5]]], iterations = 4, gamma = 25, start = NULL)
hc_5_eBIC_g25$networkdata
hc_1_eBIC_g0$networkdata

-140000-log(nrow(data))/2*8000

test5 <-hc(x = data[permutations[[5]]], max.iter = 4, start = hc_5$networks$hc_10d_8000i)
score(hc_5$networks$hc_10d_8000i, data = data[permutations[[5]]],type = "bic-g")
hc_5$networkdata[8000,]
score(test5, data = data[permutations[[5]]], type = "bic-g")
logLik
bnlearn::score(data[permutations[[5]]], x = test5, type = "bic-g")
a <- bnlearn::score(data = data[permutations[[5]]], x = hc_5$networks$hc_10d_8000i, type = "bic-g")
hc_2$networks$tabu_10d_8000i
score(a)
a


aold <- hc(x = gaussian.test, score = "bic-g")
bold <- hc(x = gaussian.test, score = "bic-g", k = 10)
scoreaold <- score(aold, data = gaussian.test, type = "bic-g")
scorebold <- score(bold, data = gaussian.test, type = "bic-g", k = 10)
  logaold <- logLik(aold, data = gaussian.test)
  logbold <- logLik(bold, data = gaussian.test)

  nparams(data = gaussian.test, aold)
  
  
  amed <- hc(x = gaussian.test, score = "bic-g")
  bmed <- hc(x = gaussian.test, score = "bic-g", k = 10)
  scoreamed <- score(amed, data = gaussian.test, type = "bic-g")
  scorebmed <- score(bmed, data = gaussian.test, type = "bic-g", k = 10)
  logamed <- logLik(amed, data = gaussian.test)
  logbmed <- logLik(bmed, data = gaussian.test)  
  
  nparams(data = gaussian.test, amed)
  amed
  
anew <- hc(x = gaussian.test, score = "bic-g")
bnew <- hc(x = gaussian.test, score = "bic-g", k = 10)
scoreanew<- score(anew, data = gaussian.test, type = "bic-g")
scorebnew <- score(bnew, data = gaussian.test, type = "bic-g", k = 10)
loganew <- logLik(anew, data = gaussian.test)
logbnew <- logLik(bnew, data = gaussian.test)
  
nparams(data = gaussian.test, anew)

anew
bnew
aold
bold

b <- gamma_hc(data = gaussian.test, iterations = 4, gamma = 0)

score(a$networks$hc_10d_g0, data = gaussian.test, type = "bic-g")
b

b <- log(nrow(data))/2*nparams(data = data[permutations[[5]]],hc_5$networks$hc_10d_8000i)
b
a -b
bnlearn::score(test5,data)
logLik(test5, data)
hc_5$networkdata

hc_5_eBIC_g25$networkdata
test6 <-  gamma_hc(data = data[permutations[[5]]], iterations = 2000, gamma = 50, start = NULL)
test6$networkdata
gamma = 25
packageDescription("bnlearn", lib.loc = NULL, fields = NULL,
                   drop = TRUE, encoding = "")
nparams()
hc_1$networkdata
tabu_1$networkdata
tabu_5$networkdata
tabu_2_eBIC_g0.1$networkdata
tabu_2_eBIC_g0$networkdata
hc_edges_loglik_10d_1000i
score(tabu_2_eBIC_g0$networks$tabu_10d_g0)
