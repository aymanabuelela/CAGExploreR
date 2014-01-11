\name{F5.hg19.promoters}
\alias{F5.hg19.promoters}
\alias{mpromdb.hg19.promoters}
\docType{data}
\title{
Promoter and Gene Definitions from FANTOM5 and MPromDB
}
\description{
These are collections of promoters from the FANTOM5 and MPromDB databases using the hg19 build. Collections of genes are included too for matching promoters. 
}
\usage{data(F5.hg19.promoters)}
\format{
  A data frame with 91460 observations on the following 5 variables.
  \describe{
    \item{\code{chr}}{a character vector. Chromosome.}
    \item{\code{strand}}{a character vector. Strand.}
    \item{\code{start}}{a numeric vector. Start position.}
    \item{\code{end}}{a numeric vector. End position.}
    \item{\code{gene}}{a character vector. Associated gene.}
  }
}
\details{
Gene regions are included from the ENSEMBL resource. Genes are identified with '@' following their name, e.g. SOX2@. Promoters associated to a gene do not have the '@' sign. Only multi-promoter genes are included. See references for details on how promoter regions are determined. Start and end positions are 0- and 1-based respectively. 
}
\source{
FANTOM5 Consortium (2013). A promoter level mammalian expression atlas. Under review at Nature.

Gupta, R., Bhattacharyya, A., Agosto-Perez, F. J., Wickramasinghe, P. and Davuluri, R. V. (2011). MPromDb
update 2010: An integrated resource for annotation and visualization of mammalian gene promoters and
ChIP-seq experimental data. Nucleic Acids Research, Vol. 39, D92-97.
}
\references{
~~ possibly secondary sources and usages ~~
}
\keyword{datasets}
