prepEnsemblAnot = function(Ensembl_list)
{
EXONNAME=EXONSTART=EXONEND=GENEID=CDSID=CDSSTART=CDSEND=TXNAME=TXSTART=TXEND=NULL
	humanEnsembl = Ensembl_list[[1]]
	hg19.ensembl.entrez.hgnc = Ensembl_list[[2]]
	
	#assign(envir,new.env(),envir=emptyenv())	
	hg19.exons.by.gene = humanEnsembl[,list(EXONNAME,EXONSTART,EXONEND),by=GENEID]
	hg19.cds.by.gene = humanEnsembl[,list(CDSID,CDSSTART,CDSEND),by=GENEID]
	hg19.exons.by.transcript.by.gene = humanEnsembl[,list(TXNAME,TXSTART,TXEND,EXONSTART,EXONEND),by=GENEID]
	hg19.exons.by.transcript = humanEnsembl[,list(EXONSTART,EXONEND),by=TXNAME]
	setkey(hg19.exons.by.transcript,TXNAME)
	
	to.return = list(
		hg19.ensembl.entrez.hgnc
		,hg19.exons.by.gene
		,hg19.cds.by.gene
		,hg19.exons.by.transcript.by.gene
		,hg19.exons.by.transcript
	)
	names(to.return) = c("hg19.ensembl.entrez.hgnc"
		,"hg19.exons.by.gene"
		,"hg19.cds.by.gene"
		,"hg19.exons.by.transcript.by.gene"
		,"hg19.exons.by.transcript"
	)
	return(to.return)
}
