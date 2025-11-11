x <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d)
hugeObj <- huge(x, method = "glasso")
hugeObj$loglik
hugeObj$df
hugeObj$sparsity
hugeSEL <- huge.select(hugeObj, criterion = "ebic", ebic.gamma =0.5)
hugeSEL
hugeSEL$loglik
hugeSEL$df
hugeSEL$refit
hugeSEL$opt.icov
hugeSEL$ebic.score
plot(hugeSEL)
plot(hugeObj)


# original data with one very large y value
x <- c(1, 2, 3, 4, 5)
y <- c(1000, 120, 110, 108, 104)
plot(x,y)

# truncated y-axis plot
y <- c(140, 120, 110, 108, 104)
plot(x,y, yaxt="n")
axis(side=2, at=c(100, 105, 110, 115, 120, 140), labels=c(100, 105, 110, 115, 120, 1000))
rect(0, 130, 1, 131, col="white", border="white")
par(xpd=T)
lines(x=c(0.7,1), y=c(130, 130))

lines(x=c(0.7,1), y=c(131, 131))

