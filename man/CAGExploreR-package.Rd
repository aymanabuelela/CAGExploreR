\name{CAGExploreR-package}
\alias{CAGExploreR-package}
\alias{CAGExploreR}
\docType{package}
\title{
Detection, Quantification and Visualization of Differential Promoter Composition in CAGE-Seq Data
}
\description{
This package is a differential composition analysis tool for Next-Generation Sequencing technologies. It can take any number of different samples or conditions (at least 2) and compare transcription (promoter) composition between them. Both MPromDB and FANTOM5 promoter set definitions are included with the package. Any number of replicates can be used for each condition, even 1, although this is highly discouraged. 
}
\details{
\tabular{ll}{
Package: \tab CAGExploreR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-01-05\cr
License: \tab MIT\cr
}
A typical workflow begins with a set of CAGE-Seq BAM/BAI files and a set of promoter and gene region definitions generated using \code{definePromoters()} and \code{defineGenes()}. \code{countTags()} is used to quantify transcription levels at these promoter and gene regions by using the BAM files. A \code{DGEList} object is created from this. pool() is then called in order to combine replicates and generate a usable file for downstream analysis. \code{diffcomp()} calculates differential promoter composition metrics and \code{plotcomp()} visualizes the data for single genes. \code{html.report()} generates multiple such figures in a convenient HTML report. 
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)

Maintainer: Emmanuel Dimont
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
\code{\link{countTags}},\code{\link{definePromoters}},\code{\link{defineGenes}},\code{\link{pool}},\code{\link{diffcomp}},\code{\link{plotcomp}},\code{\link{html.report}},\code{\link{DGEList}},\code{\link{info2osc}},\code{\link{osc2info}},\code{\link{select}}
}
