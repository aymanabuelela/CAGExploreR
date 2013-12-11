\name{html.report}
\alias{html.report}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates an HTML report with a list of promoter composition plots for convenience
}
\description{
html.report is a wrapper for plotcomp() that repeatedly calls the plotting function and saves the output in an HTML file. 
}
\usage{
html.report(data.pooled, k, fig.dir = "Figures", report.name = "Switch Report")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.pooled}{
DGEList from pool() function call. 
}
  \item{k}{
character. Gene name(s) that are present in the pooled DGEList object. Usually a vector of gene names that need to be in the HTML report. 
}
  \item{fig.dir}{
character. Name of the directory that will be created to hold the image files generated. 
}
  \item{report.name}{
character. Name of the HTML report file. 
}
}
\details{
Gene names can be specified either manually one by one, or by referring to a certain number of the rownames of a diffcomp() function call. This way only the top K significant genes are displayed in the report. 
}
\value{
Nothing is returned. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}
\note{
Just like plotcomp(), this function can take a long time to run. Each figure could take a few seconds to generate, so budget time accordingly. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{countTags}},\code{\link{definePromoters}},\code{\link{pool}},\code{\link{diffcomp}},\code{\link{plotcomp}},\code{\link{html.report}},\code{\link{DGEList}}
}
\examples{
\dontrun{
#Load a count table
data(mcf7a549.raw.counts.F5)
#Create a DGEList object
my.data = DGEList(
	counts = mcf7a549.raw.counts.F5$counts[,5:16],
	lib.size = mcf7a549.raw.counts.F5$depth,
	group = sapply(strsplit(colnames(mcf7a549.raw.counts.F5$counts[,5:16]),".",fixed=TRUE)
		,function(x) x[1]),
	genes = data.frame(
	chr = mcf7a549.raw.counts.F5$counts$chr,
	strand = mcf7a549.raw.counts.F5$counts$strand,
	start = mcf7a549.raw.counts.F5$counts$start,
	end = mcf7a549.raw.counts.F5$counts$end,
	gene = mcf7a549.raw.counts.F5$counts$gene
),
	remove.zeros = FALSE)

#OPTIONALLY prior to pooling:
#Normalize samples
my.data = calcNormFactors(my.data)
#Estimate negative binomial dispersion (common or tagwise)
my.data = estimateCommonDisp(my.data, verbose = T)
my.data = estimateTagwiseDisp(my.data)

#Finally pool the DGEList object
data.pooled = pool(my.data)

#Make plot
plotcomp(data.pooled,"SYT1")

results <- diffcomp(data.pooled)
html.report(data.pooled, rownames(results)[1:3])
}}
