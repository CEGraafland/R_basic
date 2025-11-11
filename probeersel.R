install.packages("data.table")
library(data.table)

names <- c()
for(i in 1:648) {names <- append(names,paste0("V",i))}
limitsl <- names[1:18]
limitsl <- t(limitsl)
limitsr <- names[631:648]
limitsr <- t(limitsr)

edges <- as_edgelist(igraph)
black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
blacks <- rbind(black,black2)

m1 <- edges
m2 <- blacks
colnames(m2) <- paste0("V", seq(len=ncol(m2)))
DT1 <- data.table(m1)
DT1
DT2 <- data.table(cbind(m2, 0), key=paste0("V", seq(len=ncol(m2))))
DT2
setnames(DT2, c(head(names(DT2), -1L), "found"))
a <- DT2[DT1, list(found=ifelse(is.na(found), 0, 1))]
a <-as.matrix(a)
a

edges2 <- edges[-which(a[,1] == 1, arr.ind = TRUE),]
edges2






  names <- c()
  for(i in 1:648) {names <- append(names,paste0("V",i))}
  limitsl <- names[1:18]
  limitsl <- t(limitsl)
  limitsr <- names[631:648]
  limitsr <- t(limitsr)

  black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
  black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))

  indlimits <- c()

  for (i in 1:nrow(black)){
    int <- intersect(which(as_edgelist(igraph)[,1] == black2[i,1]),
                     which(as_edgelist(igraph)[,2] == black2[i,2]))
    indlimits[i] <- int
  }
  indsame
  black2[i,]


  for (i in 1:length(limitsl)){
    which(limits[,1] == limitsl[i,1] & limits[,2]== limitsl[i,2])
  }
duplicated(c(limits,black))
crossing(limitsl,limitsr)
black <- as.matrix(expand.grid(limitsl,limitsr, stringsAsFactors = FALSE))
black2 <- as.matrix(expand.grid(limitsr,limitsl, stringsAsFactors = FALSE))
combn(letters[4:1],2)

install.packages("gtable")
library(gtable)
library(gridExtra)
a <- gtable(grid::unit(1:3, c("cm")), grid::unit(5, "cm"))
a
gtable_show_layout(a)
# Add a grob:
rect <- rectGrob(gp = gpar(fill = "black"))
a <- gtable_add_grob(a, rect, 1, 1)
a
plot(a)
# gtables behave like matrices:
dim(a)
t(a)
plot(t(a))
# when subsetting, grobs are retained if their extents lie in the
# rows/columns that retained.
b <- gtable(grid::unit(c(2, 2, 2), "cm"), grid::unit(c(2, 2, 2), "cm"))
b <- gtable_add_grob(b, rect, 2, 2)
b[1, ]
b[, 1]
b[2, 2]
# gtable have row and column names
rownames(b) <- 1:3
rownames(b)[2] <- 200
colnames(b) <- letters[1:3]
dimnames(b
)
