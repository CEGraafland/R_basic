##############################################################################
# Arcstrength by distances.
##############################################################################
rm(list = ls())
library(gridExtra)
library(gridGraphics)
library(grid)
library("bnlearn")
################################################
# For MAC:
# source("/Users/lisettegraafland/Desktop/mnt/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
################################################
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/BasicNetworkFunctions.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/HillclimbingFunctions2.R")
source("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/R/Functions/CN_ConstructionandMeasuresFunctions.R")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/tas_ncep_10d.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/permutations.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/datapermutations.rda")
################################################
# Load BNs
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm2sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm3sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm4sort.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm5sort.rda")
################################################
# Load CNs
################################################
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/gridGraphs.rda")
gridGraphsCN <- gridGraphs
rm(gridGraphs)
GraphsCN <- lapply(gridGraphsCN, function(x){x$graph})
edgelistsCN <- lapply(GraphsCN, E)
numberofedgesCN <- sapply(edgelistsCN, length)
names(gridGraphsCN) <- as.character(numberofedgesCN)
#############################################################################################
# Graph representations.
#############################################################################################
# perm1strenghts <- lapply(perm1sort, bn_to_igraph.strengths, perm = NULL,data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
# perm1dists <- lapply(perm1strenghts, igraph.distances, perm = NULL,data.igraph = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
# perm1weights <- lapply(perm1dists, igraph.weights, type = "bn", fromdist = 2000)
#(perm1weights, file = "/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1weights.rda")
load("/oceano/gmeteo/WORK/lisette/Trabajo/R_practice/Data/Struct_learn/hciterations/perm1weights.rda")
arcs1 <- lapply(perm1sort, arcs)
narcs1 <- sapply(arcs1, nrow)

edge.attributes(perm1weights[[86]])





listcn <- GraphsCN[35:70]
narcs <- numberofedgesCN[35:70]
narcs[[]]
listcnstrenghts <- lapply(listcn, str_vs_dis_2, perm = permutations[[1]], data.dag = TimeCoordsAnom_from_Grid_rms(tas_ncep_10d))
listcnweights <- lapply(listcnstrenghts, igraph.weights, type = "cn", fromdist = 0, perm = permutations[[1]])

# Choose some:
# bn:
order13 <- node.ordering(perm1sort[[13]])

order13perm <- c()
for (i in 1:ncol(datapermutations[[1]])){
  int <- which(colnames(datapermutations[[1]]) == order13[i])
order13perm[i] <- int
}
order13perm


  backorder13perm <- c()
  for (i in 1:ncol(datapermutations[[1]])){
    intback <- which(order13 == colnames(datapermutations[[1]][i]))
    backorder13perm[i] <- intback
  }
 
backorder13perm

narcs1[c(13,50,86)]
bn1 <- perm1weights[[13]]
bn2 <- perm1weights[[50]]
bn3 <- perm1weights[[86]]
# cn:
numberofedgesCN[c(69,42,35)]
narcs[c(35,8,1)]
cn1 <- listcnweights[[35]]
cn2 <- listcnweights[[8]]
cn3 <- listcnweights[[1]]

#representations
# layout as star
dev.off()
par(mfrow = c(2,2))

random <- sample(1:648,648)
parents <- backorder13perm
grid <- 1:648


plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/asStarBNCN111.pdf")
pdf(plotname)
par(mfrow = c(2,2))
dev.off()
plot.igraph(bn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            # vertex.label.cex =0.8,
            edge.arrow.width = 0.1,
            edge.arrow.size = 0.1,
            layout = layout_as_star(bn1, center = V(bn1)[81], order = random),
            main = paste0("bn: ",length(E(bn1)), " random"))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            arrow.width = 0.1,
            edge.arrow.size = 0.1,
            layout = layout_as_star(cn1, center = V(cn1)[81], order = random),
            main = paste0("cn: ",length(E(cn1)), " random"))

plot.igraph(bn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_as_star(bn1, center = V(bn1)[81], order = parents),
            main = paste0("bn: ",length(E(bn1)), " parents"))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_as_star(cn1, center = V(cn1)[81], order = parents),
            main = paste0("cn: ",length(E(cn1)), " parents"))

plot.igraph(bn1, 
            vertex.size = 0.1,
            # vertex.label.cex =0.5,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_as_star(bn1, center = V(bn1)[81], order = grid),
            main = paste0("bn: ",length(E(bn1)), " grid"))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_as_star(cn1, center = V(cn1)[81], order = grid),
            main = paste0("cn: ",length(E(cn1)), " grid"))

dev.off()
plot.igraph(bn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_as_star(bn1, center = V(bn1)[81] ),
            main = paste0(length(E(bn1))))

 plot.igraph(bn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_as_star(bn1, center = V(bn1)[81], order = backorder13perm),
            main = paste0(length(E(bn1))))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_star(cn1, center = V(cn1)[81]),
            main = paste0(length(E(cn1))))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_star(cn1, center = V(cn1)[81], order = backorder13perm),
            main = paste0(length(E(cn1))))

plot.igraph(bn2, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_star(bn2),
            edge.arrow.size = 0.1,
            main = paste0(length(E(bn2))))
plot.igraph(cn2, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_star(cn2),
            edge.arrow.size = 0.1,
            main = paste0(length(E(cn2))))


# layout as tree


plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(bn1,circular = TRUE),
            edge.arrow.size = 0.1,
            main = paste0("bn: ",length(E(bn1))))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(cn1,circular = TRUE),
            edge.arrow.size = 0.1,
            main = paste0("cn: ",length(E(cn1))))



plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(bn1,circular = FALSE,root = c(81,280),rootlevel = c(2,2)),
            edge.arrow.size = 0.1,
            main = paste0("bn: ",length(E(bn1))))

plot.igraph(cn1, 
            vertex.size = 0.1,
            # vertex.label = NA,
            layout = layout_as_tree(cn1,circular = FALSE,root = c(81,280),rootlevel = c(2,2)),
            edge.arrow.size = 0.1)

plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(bn2,circular = TRUE),
            edge.arrow.size = 0.1,
            main = paste0("bn: ",length(E(bn2))))

plot.igraph(cn2, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(cn2,circular = TRUE),
            edge.arrow.size = 0.1,
            main = paste0("cn: ",length(E(cn2))))

dev.off()

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/asTreeBNCN123.pdf")
pdf(plotname)
par(mfrow = c(3,2))


plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(bn1,circular = FALSE, mode = "all"),
            edge.arrow.size = 0.1,
            main = paste0("bn: ",length(E(bn1))))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(cn1,circular = FALSE),
            edge.arrow.size = 0.1,
            main = paste0("cn: ",length(E(cn1))))

plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(bn2,circular = FALSE),
            edge.arrow.size = 0.1,
            main = paste0("bn: ",length(E(bn2))))

plot.igraph(cn2, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(cn2,circular = FALSE),
            edge.arrow.size = 0.1,
            main = paste0("cn: ",length(E(cn2))))



plot.igraph(bn3,
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(bn3,circular = FALSE),
            edge.arrow.size = 0.1,
            main = paste0("bn: ",length(E(bn3))))

plot.igraph(cn3, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_as_tree(cn3,circular = FALSE),
            edge.arrow.size = 0.1,
            main = paste0("cn: ",length(E(cn3))))

dev.off()

# layout on grid
plot.igraph(permute.vertices(bn1, backorder13perm), 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_on_grid(permute.vertices(bn1, backorder13perm)),
            edge.arrow.size = 0.1)

plot.igraph(bn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_on_grid(bn1),
            edge.arrow.size = 0.1)

# layout in circle
# With node ordering bn: No borders: links are directed to neighbours or cross over. 
plot.igraph(bn1,
  #permute.vertices(bn1, backorder13perm), 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_in_circle(bn1, order = backorder13perm))

plot.igraph(bn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_in_circle(bn1))

plot.igraph(cn1, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_in_circle(cn1))

plot.igraph(bn2, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_in_circle(bn2),
            edge.arrow.size = 0.1)

plot.igraph(cn2, 
            vertex.size = 0.1,
            vertex.label = NA,
            layout = layout_in_circle(cn2),
            edge.arrow.size = 0.1)

# layout Fruchterman Raingold (in nicely < 1000 edges)
plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/asfrBNCN123.pdf")
pdf(plotname)

par(mfrow = c(3,2))
plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_fr(bn1, weights = E(bn1)$weights),
            main = paste0("bn: ",length(E(bn1))))

plot.igraph(cn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_fr(cn1, weights = E(cn1)$weights),
            main = paste0("cn: ",length(E(cn1))))

plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_fr(bn2, weights = E(bn2)$weights),
            main = paste0("bn: ",length(E(bn2))))

plot.igraph(cn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_fr(cn2, weights = E(cn2)$weights),
            main = paste0("cn: ",length(E(cn2))))

plot.igraph(bn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_fr(bn3, weights = E(bn3)$weights),
            main = paste0("cn: ",length(E(bn3))))

plot.igraph(cn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_fr(cn3, weights = E(cn3)$weights),
            main = paste0("cn: ",length(E(cn3))))
dev.off()
# layout with lgl: 
plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_lgl(bn1))
plot.igraph(cn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_lgl(cn1))
plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_lgl(bn2))
plot.igraph(cn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_lgl(cn2))


# layout with drl:
# with weights 

plotname <- paste0("/oceano/gmeteo/WORK/lisette/Trabajo/Resumen4/figures/asdrlBNCN123.pdf")
pdf(plotname)
par(mfrow = c(3,2))

plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(bn1, weights = E(bn1)$weights),
            main = paste0("bn: ",length(E(bn1))))
plot.igraph(cn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(cn1, weights = E(cn1)$weights),
            main = paste0("cn: ",length(E(cn1))))

plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_drl(bn2,  weights = E(bn2)$weights),
            main = paste0("bn: ",length(E(bn2))))

plot.igraph(cn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(cn2, weights = E(cn2)$weights),
            main = paste0("cn: ",length(E(cn2))))

plot.igraph(bn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_drl(bn3,  weights = E(bn3)$weights),
            main = paste0("cn: ",length(E(bn3))))

plot.igraph(cn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(cn3, weights = E(cn3)$weights),
            main = paste0("cn: ",length(E(cn3))))
dev.off()
# layout with drl:
# without weights 
plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(bn1, weights = E(bn1)$weight))
plot.igraph(cn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(cn1, weights = E(cn1)$weight))



plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(bn1, weights = E(bn2)$weight))
plot.igraph(cn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(cn1, weights = E(cn2)$weight))


plot.igraph(bn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(bn3, weights = E(bn3)$weight))
plot.igraph(cn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(cn3, weights = E(cn3)$weight))


# layout with kk with weights (don't function)

plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(bn1, weights = E(bn1)$weights),
            main = paste0("bn: ",length(E(bn1))))
plot.igraph(cn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(cn1, weights = E(cn1)$weights),
            main = paste0("cn: ",length(E(cn1))))

plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_kk(bn2,  weights = E(bn2)$weights),
            main = paste0("bn: ",length(E(bn2))))

plot.igraph(cn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(cn2, weights = E(cn2)$weights),
            main = paste0("cn: ",length(E(cn2))))

plot.igraph(bn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_kk(bn3,  weights = E(bn3)$weights),
            main = paste0("cn: ",length(E(bn3))))

plot.igraph(cn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(cn3, weights = E(cn3)$weights),
            main = paste0("cn: ",length(E(cn3))))

# layout with kk without wegith: Does function
plot.igraph(bn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(bn1, weights = E(bn1)$weight),
            main = paste0("bn: ",length(E(bn1))))
plot.igraph(cn1,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(cn1, weights = E(cn1)$weight),
            main = paste0("cn: ",length(E(cn1))))

plot.igraph(bn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_kk(bn2,  weights = E(bn2)$weight),
            main = paste0("bn: ",length(E(bn2))))

plot.igraph(cn2,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(cn2, weights = E(cn2)$weight),
            main = paste0("cn: ",length(E(cn2))))

plot.igraph(bn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.1,
            layout = layout_with_kk(bn3,  weights = E(bn3)$weight),
            main = paste0("cn: ",length(E(bn3))))

plot.igraph(cn3,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_kk(cn3, weights = E(cn3)$weight),
            main = paste0("cn: ",length(E(cn3))))



# rest
oefen <- perm1sort[[26]]
oefeni <- igraph.from.graphNEL(as.graphNEL(oefen))
nodenames <- nodes(oefen) 


plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_as_star(oefeni))

plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_fr(oefeni))


plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_as_tree(oefeni,circular = TRUE))

plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_as_tree(oefeni,circular = FALSE))

plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_in_circle(oefeni))

plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_lgl(oefeni))

layout1 <- layout_with_sugiyama(oefeni, layers = NULL, hgap = 1, vgap = 1, maxiter = 100, weights = NULL, attributes = c("default", "all", "none"))
layout1$extd_graph
layout1$layout
plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout1$layout)


plot.igraph(oefeni,
            vertex.size = 0.1,
            vertex.label = NA,
            edge.arrow.size = 0.3,
            layout = layout_with_drl(oefeni))


            