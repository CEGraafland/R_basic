install.packages("psych")
install.packages("qgraph")
library("psych")
library("qgraph")
data(bfi)
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
dataRMS <- TimeCoordsAnom_from_Grid_rms(tas_ncep_10d, rms = TRUE)
# Compute correlations:
cor1 <- cor(dataRMS)
corpos1 <- make.positive.definite(cor1, tol = 0.06)
CorMat <- cor_auto(cp)

# Compute graph with tuning = 0 (BIC):
BICgraph <- EBICglasso(S = corpos1, n = 360, 0, threshold = TRUE)

# Compute graph with tuning = 0.5 (EBIC)
EBICgraph <- EBICglasso(CorMat, nrow(bfi), 0.5, threshold = TRUE)
EBICgraph
# Plot both:
layout(t(1:2))
BICgraph <- qgraph(BICgraph, layout = "spring", title = "BIC", details = TRUE)
EBICgraph <- qgraph(EBICgraph, layout = "spring", title = "EBIC")

# Compare centrality and clustering:
par(mfrow = c(2,1))

centralityPlot(list(EBIC = EBICgraph))
centralityPlot(list(BIC = BICgraph))
centralityPlot(list(BIC = BICgraph, EBIC = EBICgraph))
clusteringPlot(list(BIC = BICgraph))
# }