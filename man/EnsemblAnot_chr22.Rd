\name{EnsemblAnot_chr22}
\alias{EnsemblAnot_chr22}
\alias{EnsemblAnot}
\docType{data}
\title{
ENSEMBL hg19 gene definitions in data.table format (chromosome 22 only)
}
\description{
The data contains 4 human genomic coordinate (start and stop) datasets from ENSEMBL based on the hg19 build on chr22: (a) exons by gene, (b) CDS by gene, (c) exons by transcript, (d) exons by transcript by gene. The package also contains one dictionary dataset that converts gene names between ENSEMBL, Entrez and HGNC symbols. All datasets are in data table format (from package data.table) for faster access time. These data are used by R package CAGExploreR for making genomic plots.
}
\usage{data(EnsemblAnot_chr22)}
\format{a list of data.tables that contain various useful combinations of the following variables:

  \describe{
    \item{\code{GENEID}}{character. ENSEMBL gene name}
    \item{\code{EXONNAME}}{character. ENSEMBL exon name}
    \item{\code{EXONSTART}}{integer. Exon start position}
    \item{\code{EXONEND}}{integer. Exon end position}
    \item{\code{CDSID}}{integer. ENSEMBL CDS ID}
    \item{\code{CDSSTART}}{integer. CDS start position}
    \item{\code{CDSEND}}{integer. CDS end position}
    \item{\code{TXNAME}}{character. ENSEMBL transcript name}
    \item{\code{TXSTART}}{integer. Transcript start position}
    \item{\code{TXEND}}{integer. Transcript end position}
    \item{\code{ensembl_gene_id}}{character. ENSEMBL gene name}
    \item{\code{entrezgene}}{character. Entrez gene name}
    \item{\code{hgnc_symbol}}{character. HGNC gene name}
  }
}
\details{
Gene regions are included from the ENSEMBL resource, obtained using package biomaRt. 
}
\source{
http://www.ensembl.org/
}
\references{
~~ possibly secondary sources and usages ~~
}

\keyword{datasets}
\seealso{
\code{\link{plotcomp}},\code{\link{getEnsemblAnot}},\code{\link{prepEnsemblAnot}}
}