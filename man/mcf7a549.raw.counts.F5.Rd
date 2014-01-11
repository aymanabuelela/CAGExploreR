\name{mcf7a549.raw.counts.F5}
\alias{mcf7a549.raw.counts.F5}
\docType{data}
\title{
Promoter and Gene expression from whole-cell CAGE-Seq applied to MCF7 breast and A549 lung cancer cell lines
}
\description{
CAGE was applied to a set of cell lines as part of the ENCODE project to determine promoter-level expression. This data can be reproduced by running \code{countTags} on the proper BAM and BAI input files. 
}
\usage{data(mcf7a549.raw.counts.F5)}
\format{
  The format is:
List of 2
 $ depth : Named num [1:4] 28855537 18239960 52947757 36259720
  ..- attr(*, "names")= chr [1:4] "a545.1" "a545.2" "mcf7.1" "mcf7.2"
 $ counts:Classes 'data.table' and 'data.frame':	91460 obs. of  6 variables:
  ..$ region: chr [1:91460] "chr10:100143322..100174982,-" "chr10:100174900..100174956,-" "chr10:100174957..100174982,-" "chr10:100216834..100995635,-" ...
  ..$ gene  : chr [1:91460] "PYROXD2@" "PYROXD2" "PYROXD2" "HPSE2@" ...
  ..$ a545.1: int [1:91460] 457 5 282 27 0 0 0 0 98 31 ...
  ..$ a545.2: int [1:91460] 330 4 200 33 0 0 0 0 91 27 ...
  ..$ mcf7.1: int [1:91460] 1154 1 925 185 0 0 0 0 246 10 ...
  ..$ mcf7.2: int [1:91460] 3530 21 648 54 0 0 0 0 127 17 ...
  ..- attr(*, "sorted")= chr "region"
  ..- attr(*, ".internal.selfref")=<externalptr> 
}
\details{
Raw data can be downloaded from the ENCODE website http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRikenCage/ . The files are:
\enumerate{
\item wgEncodeRikenCageA549CellPapAlnRep1.bam - A549 whole-cell RNA CAGE rep1
\item wgEncodeRikenCageA549CellPapAlnRep2.bam - A549 whole-cell RNA CAGE rep2
\item wgEncodeRikenCageMcf7CellPapAlnRep1.bam - MCF7 whole-cell RNA CAGE rep1
\item wgEncodeRikenCageMcf7CellPapAlnRep2.bam - MCF7 whole-cell RNA CAGE rep2
}
}
\source{
ENCODE Consortium (2012). An integrated encyclopedia of DNA elements in the human genome. Nature. 489: 57-74.
}
\references{
~~ possibly secondary sources and usages ~~
}

\keyword{datasets}
