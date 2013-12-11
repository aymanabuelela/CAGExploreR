\name{definePromoters}
\alias{definePromoters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate gene promoter regions usable by CAGExploreR
}
\description{
definePromoters converts built-in FANTOM5, MPromDB or user-supplied gene promoter regions into a usable object by CAGExploreR for downstream analysis, specifically for quantifying the number of CAGE tags mapping to these regions. This is one of the first steps that needs to be performed before any other analysis can be performed. All promoters must be associated with a gene. At least one gene in the list must have 2 or more promoters. Analysis cannot be performed on mono-promoter genes or promoters not assigned to a gene and these will be dropped. 
}
\usage{
definePromoters(file = NA, sep = "\t", header = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
character. This is the name (or full file path) of the user-supplied gene promoter specification file. If file=NA (default), this should be a name of an object already loaded into the workspace. 
}
  \item{sep}{
character. If a user-supplied file is used, this specifies the column delimiter (tab-delimited by default). Not used if file=NA. 
}
  \item{header}{
logical. If a user-supplied file is used, this specifies whether the first row in the file corresponds to column headers. Not used if file=NA. 
}
}
\details{
Promoter region definition files should have the following columns: 1. chromosome (e.g. chr20, chrX, chrM, etc.), 2. strand (i.e. +,-), 3. start position (integer), 4. end position (integer), 5. HGNC gene symbol (character). Column names must not end with a number. 
}
\value{
A list is returned with 5 entries:
  \item{which.pos}{RangesList object, with a chromosome for each entry in the list. Promoter regions on + strand are specified here. Each list entry is an IRanges object with 3 columns: start, end and region width}
  \item{which.neg}{RangesList object, with a chromosome for each entry in the list. Promoter regions on - strand are specified here. Each list entry is an IRanges object with 3 columns: start, end and region width}
  \item{CHROMOSOMES}{character. Vector of chromosome identifiers present in promoter definition file}
  \item{n.promoters}{integer. Number of promoter regions defined in input file}  
  \item{original}{data.frame. A copy of the original file submitted to definePromoters}
This list is not meant for direct use, but as an input for downstream analysis. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{countTags}},\code{\link{definePromoters}},\code{\link{pool}},\code{\link{diffcomp}},\code{\link{plotcomp}},\code{\link{html.report}},\code{\link{DGEList}}
}
\examples{
\dontrun{
if(FALSE){
This example uses the following BAM and BAI files in the local folder from the ENCODE project website at: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRikenCage/
wgEncodeRikenCageA549CellPapAlnRep1.bam
wgEncodeRikenCageA549CellPapAlnRep2.bam   
wgEncodeRikenCageA549CytosolPapAlnRep3.bam
wgEncodeRikenCageA549CytosolPapAlnRep4.bam
wgEncodeRikenCageA549NucleusPapAlnRep3.bam
wgEncodeRikenCageA549NucleusPapAlnRep4.bam
wgEncodeRikenCageMcf7CellPapAlnRep1.bam 
wgEncodeRikenCageMcf7CellPapAlnRep2.bam   
wgEncodeRikenCageMcf7CytosolPapAlnRep3.bam
wgEncodeRikenCageMcf7CytosolPapAlnRep4.bam
wgEncodeRikenCageMcf7NucleusPapAlnRep3.bam
wgEncodeRikenCageMcf7NucleusPapAlnRep4.bam}

data(F5.hg19.promoters)
my.prom.defs = definePromoters(F5.hg19.promoters)

files = dir()[grep("bam",dir())]
my.bai.files = dir()[grep("bam.bai",dir())]
my.bam.files = setdiff(files, my.bai.files)
my.ids = c("a545.1","a545.2","a545.3","a545.4","a545.5","a545.6","mcf7.1","mcf7.2","mcf7.3","mcf7.4","mcf7.5","mcf7.6")

#data = countTags(my.bam.files, my.bai.files, my.ids, my.prom.defs) #2 mins per BAM file
}}
