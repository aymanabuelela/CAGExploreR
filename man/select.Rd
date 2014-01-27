\name{Select}
\alias{Select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract a substrings from a larger string that contains delimiters
}
\description{
A long character string may contain data that is concatenated together with the use of delimiters. \code{Select} is a convenient wrapper for the \code{substr}-locate-\code{unlist} routine. 
}
\usage{
Select(x, split, which)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
character. Vector of strings to be split. 
}
  \item{split}{
character. Single character used as the delimiter on which to split. 
}
  \item{which}{
integer. After splitting, this number indicates which piece to keep. 
}
}
\details{
A combination of calls to \code{Select} can be used to perform quite complicated string parsing. 
}
\value{
character vector of substring extracted. 
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}

\seealso{
\code{\link{osc2info}},\code{\link{info2osc}}
}
\examples{
\dontrun{

chr = paste("chr",1:22,sep="")
strand = sample(c("+","-"),22,replace=TRUE)
chr.strand = paste(chr,strand,sep=".")

chr.extracted = Select(chr.strand,".",1)
strand.extracted = Select(chr.strand,".",2)


}}
