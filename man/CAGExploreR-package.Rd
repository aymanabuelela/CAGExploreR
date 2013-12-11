\name{CAGExploreR-package}
\alias{CAGExploreR-package}
\alias{CAGExploreR}
\docType{package}
\title{
Detection, Quantification and Visualization of Differential Promoter Composition in CAGE-Seq Data
}
\description{
This package is designed as a differential composition tool for Next-Generation Sequencing technologies. It can take any number of different samples or conditions (at least 2) and compare transcription (promoter) composition between them. Both MPromDB and FANTOM5 promoter set definitions are included with the package. Any number of replicates can be used for each condition, even 1, although this is highly discouraged. 
}
\details{
\tabular{ll}{
Package: \tab CAGExploreR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-14\cr
License: \tab MIT\cr
}
A typical workflow begins with a set of CAGE-Seq BAM files and a set of promoter region definitions generated using definePromoters(). countTags() is used to quantify transcription levels at these promoter regions by using the BAM files. A DGEList object is created from this. pool() is then called in order to combine replicates and generate a usable file for downstream analysis. diffcomp() calculates differential promoter composition metrics and plotcomp() visualizes the data for single genes. html.report() generates multiple such figures in a convenient HTML report. 
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
\code{\link{countTags}},\code{\link{definePromoters}},\code{\link{pool}},\code{\link{diffcomp}},\code{\link{plotcomp}},\code{\link{html.report}},\code{\link{DGEList}}
}
\examples{
\dontrun{
}}
