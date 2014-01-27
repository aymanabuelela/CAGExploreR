\name{diffcomp}
\alias{diffcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Differential Promoter Composition (DPC). Gene-level and Promoter-level statistics are available.
}
\description{
diffcomp generates a sorted data frame with genes that have the strongest differential promoter composition across samples at the very top. DPC is quantified by a proportional entropy reduction score (Theil's U). pool() must be called before calling this function. 
}
\usage{
diffcomp(x,detailed=F,top="all",gene=NA,B=1000,seed=1,mc.adjust="fdr",text=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{DGEList from pool() function call. }
\item{detailed}{boolean. Whether to display brief gene-based DPC statistics (default) or detailed promoter-level switching event statistics}
\item{top}{integer or "all". If detailed=TRUE and gene=NA, number of top genes whose promoter-level statistics be displayed. By default all genes are displayed}
\item{gene}{character. If detailed=TRUE, vector of genes to include in promoter-level statistics display}
\item{B}{integer. Number of iterations to use for Monte Carlo p-value estimation. Larger values take longer to compute but are more accurate}
\item{seed}{integer. Seed for random number generator in Monte Carlo p-value generation. This can be used for sensitivity analysis}
\item{mc.adjust}{one of: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none" method for adjusting p-values for multiple comparisons}
\item{text}{logical. If detailed=TRUE, whether to display which promoters are being compared. By default, this will be displayed. Can be turned off to save memory on large data frames.}
}
\details{
Entropy reduction varies from 0 to 1, with an entropy reduction of 1 corresponding to e.g. condition A where promoter P1 is exclusively used out of a set of possible promoters, while in condition B, promoter P2 is used exclusively, i.e. a perfect switch in promoter usage between conditions. Entropy reduction of 0 corresponds to equal relative transcription from promoters between 2 or more conditions. Gene heterogeneity measures are entropy based and also range from 0 to 1. Values of 0 correspond to no change while vaues of 1 correspond to maximal change between conditions. These measures give an idea of the relative changes in gene expression v.s. changes in promoter composition. For promoter-level statistics, the odds ratio is the main effect measure. NOTE: odds ratios are calculated after 0.5 is added to all counts to avoid dividing by zero. 
}
\value{
if detailed=FALSE (default), a data.frame is produced with gene-level statistics. Each row corresponds to a gene, and genes are sorted by decreasing entropy reduction. The following columns are included:
\item{Entropy.Reduction}{0 to 1. measures promoter switching effect size (see details above). If NaN, this means that one of the samples had 0 expression and so entropy reduction could not be calculated}
\item{pvalue}{0 to 1. Monte Carlo based P-value for promoter switching, testing if Entropy Reduction is 0 under the null hypothesis}
\item{qvalue}{0 to 1. Multiplicity adjusted p-value for promoter switching}
\item{GeneHetero}{0 to 1. Gene expression heterogeneity. (1 - normalized entropy) of gene expression across samples. Higher values show larger changes across samples}
\item{Coverage}{numeric. Proportion of CAGE-Seq tags mapping to defined promoter regions relative to all tags mapping to the entire gene region. Values close to 1 signify good promoter set definitions. Values over 1 point to the possibility of defined promoter regions that overlap one another. Values close to 0 signify poor coverage by promoters, possibility of large level of background noise or novel undefined promoter. NA if gene boundaries are unknown and hence, coverage cannot be computed.}
\item{dominant.promoter.switch}{character. displays which promoter(s) switch from being dominant in one condition to another}
\item{RepAgree}{integer. Specified how well replicates cluster together within their respective samples. A value of 1 means that replicates cluster perfectly together, larger values mean that replicates do not. This number is the maximum number of discordant sample labels that cluster together after cutting the dendrogram of the promoter composition vectors hierarchical clustering results.}

If detailed=TRUE, a more detailed data.frame is produced with promoter-level statistics. By default, every pairwise comparison of promoters within genes for all genes will be displayed, but as the number of samples simultaneously compared increases, this value can be decreased, since the number of rows in the output grows very rapidly and can quickly lead to memory overload problems.
\item{comparison}{character. if text=TRUE, Shows pair of promoters and pair of samples being compared}
\item{gene}{character. Gene name}
\item{theilU}{numeric. Entropy reduction score}
\item{OR}{numeric. Odds ratio measuring switching effect size. Value of 1 means no switching}
\item{log2OR}{numeric. Log-base-2 odds ratio measuring switching effect size. Value of 0 means no switching. Results within genes are sorted on this value}
\item{pvalue}{numeric. P-value associated with the odds ratio}
\item{qvalue}{numeric. Multiplicity corrected p-value. Note: if top is used as an argument, then this is based only on the p-values in the current table only, and not on all results for all genes.}
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

#Obtain DPC results
results = diffcomp(data.pooled)
results = diffcomp(data.pooled,detailed=TRUE)
}}
