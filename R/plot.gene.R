plot.gene = function(gene,my.prom.defs,counts,lib.counts)
{
EXONSTART=EXONEND=CDSSTART=CDSEND=TXNAME=V1=.N=N=NULL
#data(hg19.ensembl.entrez.hgnc,package="DTGenAnnot",envir=environment())
#data(hg19.exons.by.gene,package="DTGenAnnot",envir=environment())
#data(hg19.cds.by.gene,package="DTGenAnnot",envir=environment())
#data(hg19.exons.by.transcript.by.gene,package="DTGenAnnot",envir=environment())
#data(hg19.exons.by.transcript,package="DTGenAnnot",envir=environment())

	my.gene = convert(gene,from="HGNC.symbol",to="Ensembl.Gene.ID",using=hg19.ensembl.entrez.hgnc)
	proms = my.prom.defs
	exons = hg19.exons.by.gene[my.gene][,list(EXONSTART,EXONEND)]
	cds = hg19.cds.by.gene[my.gene][,list(CDSSTART,CDSEND)]
	tx.names = hg19.exons.by.transcript.by.gene[my.gene][,unique(TXNAME)]

	txs = hg19.exons.by.transcript[tx.names]
	gene.thickness = 0.05
	leftmost = min(c(exons$EXONSTART,exons$EXONEND))
	rightmost = max(c(exons$EXONSTART,exons$EXONEND))
	window.width = rightmost - leftmost
	gene.ylevel = 1;prom.ylevel = 0
	tx.ylevel = seq(prom.ylevel+gene.thickness*2,gene.ylevel-gene.thickness*2,length.out=length(tx.names)) 
	tx.thickness = min(diff(tx.ylevel)[1] / 2 , gene.thickness)
	
	par(mar=c(5.1,4.1+5,4.1,2.1))
	par(mgp=c(3,0.5,0))
	plot(x=c(leftmost,rightmost),y=c(1,1),ylim=c(0-gene.thickness,1+gene.thickness),type="n",xlab="",ylab="",yaxt="n",xaxt="n")
	title(paste0(gene," on ",proms$chr[1],",",proms$strand[1]),line=2)
	axis(1,at=axTicks(1),labels=format(axTicks(1),big.mark=",",scientific=FALSE))
	axis(3,at=axTicks(1),labels=format(axTicks(1),big.mark=",",scientific=FALSE))
	segments(x0=leftmost,x1=rightmost,y0=gene.ylevel,y1=gene.ylevel, col="black")		
	rect(xleft=exons$EXONSTART,xright=exons$EXONEND,ybottom=gene.ylevel-gene.thickness/2,ytop=gene.ylevel+gene.thickness/2, col = "black",lwd=0,border="black")
	rect(xleft=cds$CDSSTART,xright=cds$CDSEND,ybottom=gene.ylevel-gene.thickness/2,ytop=gene.ylevel+gene.thickness/2, col = "red",lwd=0,border="red")
	mtext(gene,side=2,line=1,at=1,las=2)
	rect(xleft=proms$start,xright=proms$end,ybottom=prom.ylevel-gene.thickness/2,ytop=prom.ylevel+gene.thickness/2, col = "black",lwd=0,border="black")
	mtext("Promoters",side=2,line=1,at=0,las=2)

	leftmost = txs[,min(EXONSTART,EXONEND),by=TXNAME][,V1]
	rightmost = txs[,max(EXONSTART,EXONEND),by=TXNAME][,V1]
	segments(x0=leftmost,x1=rightmost,y0=tx.ylevel,y1=tx.ylevel, col="black")
	rect(xleft=txs[,EXONSTART,by=TXNAME][,EXONSTART],xright=txs[,EXONEND,by=TXNAME][,EXONEND],ybottom=rep(tx.ylevel-tx.thickness/2,txs[,.N,by=TXNAME][,N]),ytop=rep(tx.ylevel+tx.thickness/2,txs[,.N,by=TXNAME][,N]), col = "blue",lwd=0,border="blue")
	mtext(tx.names,side=2,line=1,at=tx.ylevel,las=2,cex=1/(length(tx.names))^(1/3))
	my.colors = rainbow(nrow(proms),alpha=0.3)
	rect(xleft=proms$start-window.width/100,xright=proms$end+window.width/100,ybottom=prom.ylevel-gene.thickness/2,ytop=(rowSums(t(t(counts)/lib.counts))/sum(t(t(counts)/lib.counts))), col = my.colors,lwd=0,border=NA)
}
