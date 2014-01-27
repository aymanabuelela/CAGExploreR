\name{getEnsemblAnot}
\alias{getEnsemblAnot}
\alias{prepEnsemblAnot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download and Prepare Human Gene Annotation Data from ENSEMBL Biomart
}
\description{
These functions are convenience wrappers for downloading relevant human gene annotation data from the ENSEMBL biomart in data.table format for use in CAGExploreR. 
}
\usage{
getEnsemblAnot()
prepEnsemblAnot(Ensembl_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ensembl_list}{
list. Output list from getEnsemblAnot. 
}
}
\details{
These functions need to be run at the beginning of an R session if plotting is required using CAGExploreR. A call to \code{getEnsemblAnot} will download gene, transcript, exon, and cds coordinate data from ENSEMBL via the package \code{biomaRt}. An output list of length 2 is generated. This is then passed to \code{prepEnsemblAnot} which returns a more detailed list, which stores data.table format tables that CAGExploreR can use when calling \code{plotcomp} for plotting. To avoid having to download data from an online biomart each time, the output from \code{prepEnsemblAnot} can be saved as an .RData file and loaded each time. 
}
\value{
  \item{getEnsemblAnot}{
list of length 2. The first element is a data.table with gene id, transcript id, exon rank, start, stop, chromosome and strand information. The second element is a data.table of ENSEMBL to HGNC and Entrezgene ID mappings. 
}
  \item{prepEnsemblAnot}{
returns a list, 5 data.tables are stored: (a) ENSEMBL, Entrez, HGNC ID gene mappings, (b) exons by gene, (c) cds by gene, (d) exons by transcript by gene, (e) exons by transcript.  
}
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Emmanuel Dimont - Hide Laboratory for Computational Biology. Department of Biostatistics. Harvard School of Public Health. (edimont@mail.harvard.edu)
}

\seealso{
\code{\link{plotcomp}}
}

\examples{
\dontrun{

Ensembl_list = getEnsemblAnot() #this can take 10 minutes
anot.tabs = prepEnsemblAnot(Ensembl_list)
save(anot.tabs,file="anot.tabs.RData") #save for future use
load("anot.tabs.RData") #faster to load than to download fresh each time

}}
