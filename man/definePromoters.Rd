\name{definePromoters}
\alias{definePromoters}
\alias{defineGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate gene and promoter regions usable by CAGExploreR
}
\description{
\code{definePromoters} and \code{defineGenes} convert built-in FANTOM5, MPromDB or user-supplied promoter and gene regions into a usable objects by CAGExploreR for downstream analysis, specifically for quantifying the number of CAGE tags mapping to these regions. This is one of the first steps that needs to be performed before any other analysis can be performed. All promoters must be associated with a gene. At least one gene in the list must have 2 or more promoters. Analysis cannot be performed on mono-promoter genes or promoters not assigned to a gene and these will be dropped. Currently, only non-overlapping promoters and non-overlapping gene regions are supported. 
}
\usage{
definePromoters(x)
defineGenes(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
R object. This is the name gene-promoter specification file.This should be a name of an object already loaded into the workspace. 
}
}
\details{
Promoter-gene region definition files should have the following columns: 1. chromosome (e.g. chr1, chr20, chrX, chrM, etc.)(character), 2. strand (i.e. +,-)(character), 3. start position (integer), 4. end position (integer), 5. HGNC gene symbol (character). Column names must not end with a number. An input file can contain both promoter and gene region definitions. To differentiate between the two: gene regions should have the gene specified with an @ sign at the end of the gene name (e.g. SOX2@). 
}
\value{
A data.table is returned with 5 entries:
  \item{chr}{character. Chromosome}
  \item{strand}{character. Strand}
  \item{position}{integer. Genomic position (either start or end)}
  \item{gene}{character. Name of Gene}  
  \item{prom}{character. concatenation of chr, strand and position in OSC (RIKEN) format}
This data.table is not meant for direct use, but as an input for downstream analysis. 
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}
\note{
Raw BAM and BAI files from ENCODE can be found at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRikenCage/
}

\seealso{
\code{\link{countTags}},\code{\link{pool}},\code{\link{diffcomp}},\code{\link{plotcomp}},\code{\link{html.report}},\code{\link{DGEList}}
}
\examples{
\dontrun{
#This example uses the following BAM and BAI files 
#in the local folder from the ENCODE project website
#wgEncodeRikenCageA549CellPapAlnRep1.bam
#wgEncodeRikenCageA549CellPapAlnRep2.bam   
#wgEncodeRikenCageMcf7CellPapAlnRep1.bam 
#wgEncodeRikenCageMcf7CellPapAlnRep2.bam

files <- dir()[grep("bam", dir())]
my.bai.files <- dir()[grep("bam.bai", dir())]
my.bam.files <- setdiff(files, my.bai.files)
my.ids <- c("a549.1", "a549.2", "mcf7.1", "mcf7.2")

data(F5.hg19.promoters)
head(F5.hg19.promoters)
my.promoters = definePromoters(F5.hg19.promoters)
my.genes = defineGenes(F5.hg19.promoters)

mcf7a549.raw.counts.F5 = countTags(my.bam.files, my.bai.files, my.ids, my.genes, my.promoters)
# This takes 2-10 minutes per BAM file
}}
