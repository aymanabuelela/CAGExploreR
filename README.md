CAGExploreR
===========

#To install:
#First install packages from Bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite(c("edgeR","biomaRt","GenomicFeatures"))

#then packages from CRAN
install.packages(c("R2HTML","data.table","rbamtools"))
install.packages("CAGExploreR")

#load package:
library(CAGExploreR)

#To view PDF guide after loading:
vignette("CAGExploreR")

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/edimont/cagexplorer/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
