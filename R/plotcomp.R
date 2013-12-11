plotcomp <-
function(x,gene)
{
	gene = gene[1]
	if(is.na(match(gene,as.character(x$pooled.sub.counts$gene))[1])) {cat(paste0("Sorry there's no gene ",gene,"\n","Did you mean: ",paste(grep(gene,unique(as.character(x$pooled.sub.counts$gene)),value=T,ignore.case=T),collapse=","),"?"));return(1)}
	PROMS = plot.switch(x$pooled.sub.counts[x$pooled.sub.counts$gene==gene,-match(c("gene","dispersion"),colnames(x$pooled.sub.counts))]
	,lib.counts=x$pooled.samples[,"effective.lib.size"]
	,coverage=x$pooled.super.counts[,-match("dispersion",colnames(x$pooled.super.counts))]
	,GENE=gene)
	
	return(PROMS)
}
