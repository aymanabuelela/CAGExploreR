\name{pool}
\alias{pool}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pools tag counts across replicates
}
\description{
If replicates are present in CAGE data following a \code{countTags()} function call and the creation of a DGEList object, this function will pool tag counts across replicates by summing. Replicates are determined by common sample name prefix and different numeric suffix. More conveniently however, they are determined from the DGEList object "samples" element. 
}
\usage{
pool(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
DGEList object. This is typically the output from \code{countTags()} function call followed by a DGEList object creation, or a user-specified count DGEList object.
}
}
\details{
Pooling adds to additional elements to a DGEList object, specifically: separates gene regions and promoter regions into separate elements, computes effective library sizes and creates pooled count tables for both gene regions and promoter regions to be used in downstream analysis. This step is required even if no replicates are to be pooled!
}
\value{
An augmented DGEList object is returned.
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}

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
	#group = select(my.ids,".",1), #pools replicates (alternative)
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

}}
