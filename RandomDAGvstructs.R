f <- function(a,n){(a-2)/((n*(n-1)/2)-2)}
f(898,648)
prop <- function(a,n){a/(n*(n-1)/2)}
twopaths <- function(n,prop){choose(n,3)*prop^2}
triangles <- function(n,prop){choose(n,3)*prop^3}

twop <-twopaths(648,p)

clus <- function(n,a){triangles(n,prop(a,n))/(0.25*twopaths(n,prop(a,n)))}
clus(648,898)

a <- 7863
b <- 648
p <- prop(a,b)
triangles(b,p)/twopaths(b,p)*0.25


twopaths(b,p)*0.25
triangles(b,p)
(0.25*twopaths(b,p) - triangles(b,p))/twopaths(b,p)
