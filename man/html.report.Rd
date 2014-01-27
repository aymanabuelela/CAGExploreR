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
html.report(x, k, anot_list, fig.dir = "Figures", report.name = "Switch Report")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
DGEList from pool() function call. 
}
  \item{k}{
character. Gene name(s) that are present in the pooled DGEList object. Usually a vector of gene names that need to be in the HTML report. 
}
  \item{fig.dir}{
character. Name of the directory that will be created to hold the image files generated. 
}
\item{anot_list}{
list. Annotation list as output from \code{prepEnsemblAnot} 
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
 ~put references to the literature/web site here ~
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
data(mcf7a549.raw.counts.F5)
annotations = osc2info(mcf7a549.raw.counts.F5$counts$region)
my.ids <- c("a549.1", "a549.2", "mcf7.1", "mcf7.2")

my.data = DGEList(
	counts = mcf7a549.raw.counts.F5$counts[,-c(1,2),with=FALSE], #remove columns 1&2
	lib.size = mcf7a549.raw.counts.F5$depth,
	group = my.ids, #will not pool replicates
	#group = c("a549","a549","mcf7","mcf7"), #pools replicates
	#group = Select(my.ids,".",1), #pools replicates (alternative)
	genes = data.frame(
	chr = annotations$chr,
	strand = annotations$strand,
	start = annotations$start,
	end = annotations$end,
	gene = mcf7a549.raw.counts.F5$counts$gene
),
	remove.zeros = FALSE
)

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
