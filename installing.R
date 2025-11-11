#For installing rJava:
#https://github.com/snowflakedb/dplyr-snowflakedb/wiki/Configuring-R-rJava-RJDBC-on-Mac-OS-X
devtools::install_github(c("SantanderMetGroup/loadeR.java", "SantanderMetGroup/loadeR"))

library(loadeR)
Sys.getenv("JAVA_HOME")

# For installing climate4R things:
devtools::install_github('SantanderMetGroup/transformeR')

# For installing Bioconductor dependencies: Rgrapviz rgbl graph
source("https://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")