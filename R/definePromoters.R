definePromoters = function(x){
	gene.regions = grep("@",x$gene)
	proms = x[-gene.regions,]
	return(convert2pos(proms))
}