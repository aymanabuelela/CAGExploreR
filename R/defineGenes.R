defineGenes = function(x){
	gene.regions = grep("@",x$gene)
	genes = x[gene.regions,]
	return(convert2pos(genes))
}