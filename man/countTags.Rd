\name{countTags}
\alias{countTags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantify transcription from defined genomic regions in CAGE data
}
\description{
countTags quantifies the level of transcription from defined regions (usually promoters) by counting the number of transcription start sites (TSSs) that map to pre-specified genomic regions. Because CAGE-Seq captures the 5' end of mRNA and generates a short tag after converting to cDNA, the first base position of such a read is taken as the TSS. The function needs as an input a BAM alignment file together with its index BAI file and a promoter region definition file from a definePromoters() function call.  
}
\usage{
countTags(bamfiles, idx, ids, GENES, PROMS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamfiles}{
character. Vector of BAM alignment file names in the current working directory or file paths to them. 
}
  \item{idx}{
character. Vector of BAI index file names in the current working directory or file paths to them. 
}
  \item{ids}{
character. Vector of short identifiers for the different samples compared. The names should have common prefixes for samples from the same experimental condition and replicates should be differentiated by a different number at the end of the name. If the name ends with a number (e.g. cell lines such as MCF7, A549, etc.) then a dot should separate this prefix from the replicate number at the end. See examples below. 
}
  \item{GENES}{
R object from \code{defineGenes()} call. 
}
  \item{PROMS}{
R object from \code{definePromoters()} call. 
}
}
\details{
Each BAM file should have a corresponding BAI index file specified. At least 2 BAM files are needed that correspond to 2 different conditions being studied (single-replicate setting). However, additional replicates are highly recommended to account for biological and technical variation. NOTE: if only 1 replicate is used for any of the conditions being compared, estimated levels of transcription can be unstable and unreliable for generalization. ids are used for downstream analysis steps in order to decide how replicates are pooled together. Samples with the same id prefix names and different suffix replicate numbers are to be pooled together at a later stage in the pipeline. These pooling schemes can be modified later on. 
}
\value{
A list is returned with 2 entries:
	\item{depth}{numeric. Vector with the total number of tags mapped to the genome in each library}
	\item{counts}{data table. Contains annotation columns from the promoter region definition file used, together with columns for each BAM file specified with counts for each promoter region in those columns.} 
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}
\note{
Counting tags can be slow! This process takes approximately 10 minutes per genome-wide BAM file on an Intel Core2 Duo 2GHz laptop with 4GB RAM running Windows7. Raw BAM and BAI files from ENCODE can be found at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRikenCage/
}

\seealso{
\code{\link{countTags}},\code{\link{definePromoters}},\code{\link{pool}},\code{\link{diffcomp}},\code{\link{plotcomp}},\code{\link{html.report}},\code{\link{DGEList}}
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
