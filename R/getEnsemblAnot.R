getEnsemblAnot = function()
{
CDSNAME=CDSCHROM=EXONCHROM=CDSSTRAND=EXONSTRAND=hgnc_symbol=NULL

	humanEnsembl.db = makeTranscriptDbFromBiomart(biomart="ensembl"
		,dataset="hsapiens_gene_ensembl")
	
	suppressWarnings({
	humanEnsembl = select(humanEnsembl.db
		,keys(humanEnsembl.db,"GENEID")
		,keytype="GENEID"
		,columns=columns(humanEnsembl.db))
	})

	humanEnsembl = data.table(humanEnsembl)
	humanEnsembl = humanEnsembl[,CDSNAME:=NULL]
	humanEnsembl = humanEnsembl[,CDSCHROM:=NULL]
	humanEnsembl = humanEnsembl[,EXONCHROM:=NULL]
	humanEnsembl = humanEnsembl[,CDSSTRAND:=NULL]
	humanEnsembl = humanEnsembl[,EXONSTRAND:=NULL]
	setnames(humanEnsembl,old="TXCHROM",new="CHROM")
	setnames(humanEnsembl,old="TXSTRAND",new="STRAND")
	setkeyv(humanEnsembl,c("GENEID","TXNAME","EXONRANK"))
	
	ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
	hg19.ensembl.entrez.hgnc = data.table(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),mart=ensembl))
	hg19.ensembl.entrez.hgnc[hg19.ensembl.entrez.hgnc$hgnc_symbol=="","hgnc_symbol"] = as.character(NA)
	setkey(hg19.ensembl.entrez.hgnc,hgnc_symbol) 
	
	return(list(humanEnsembl,hg19.ensembl.entrez.hgnc))
}