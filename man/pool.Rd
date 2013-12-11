\name{pool}
\alias{pool}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pools tag counts across replicates
}
\description{
If replicates are present in CAGE data following a countTags() function call and the creation of a DGEList object, this function will pool tag counts across replicates by summing. Replicates are determined by common sample name prefix and different numeric suffix. More conveniently however, they are determined from the DGEList object "samples" element. 
}
\usage{
pool(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
DGEList object. This is typically the output from countTags() function call followed by a DGEList object creation, or a user-specified count DGEList object.
}
}
\details{
Pooling adds to additional elements to a DGEList object, specifically: separates gene regions and promoter regions into separate elements, computes effective library sizes and creates pooled count tables for both gene regions and promoter regions to be used in downstream analysis. 
}
\value{
An augmented DGEList object is returned.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}

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

}}
