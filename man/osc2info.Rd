\name{osc2info}
\alias{osc2info}
\alias{info2osc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between RIKEN OSC genomic region format and list with extracted data
}
\description{
The RIKEN OSC format for genomic regions is in the form 'chr:start..end,strand'. \code{osc2info} converts a vector of such regions into a list with separate entries for chr, start, end and strand. \code{info2osc} does the opposite, converts a list (or data.frame) with chr, start, end and strand data to compact OSC form. 
}
\usage{
osc2info(x)
info2osc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Either a data.frame with chr, start, end and strand columns (for \code{info2osc}) or a character vector of RIKEN OSC format regions (for \code{osc2info}). 
}
}
\details{
When using \code{osc2info}, a list is generated, however this can be converted to a data.frame by applying the \code{as.data.frame()} function to it. 
}
\value{
A list with chromosome, strand, start and end columns when using \code{osc2info}.
A character vector of RIKEN OSC format gene regions when using \code{info2osc}.
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}

\seealso{
\code{\link{Select}}
}
\examples{
\dontrun{
data(mcf7a549.raw.counts.F5)

annotations = osc2info(mcf7a549.raw.counts.F5$counts$region)
data.frame.format = head(as.data.frame(annotations))
osc.format = info2osc(data.frame.format)

)
}}
