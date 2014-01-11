convert2pos = function(x){
	proms = x
	prom.names = info2osc(proms)
	proms$prom = prom.names
	proms = proms[rep(1:nrow(proms),rep(2,nrow(proms))),]
	proms$start[seq(2,nrow(proms),2)] = proms$end[seq(1,nrow(proms)-1,2)]
	proms$end = NULL
	colnames(proms)[match("start",colnames(proms))] = "position"
	proms = data.table(proms)
	setkeyv(proms,c("chr","strand","position"))
	temp = data.frame(proms)
	temp = split(temp,list(temp$chr,temp$strand))
	to.remove = lapply(temp,function(x) which(diff(x$position)<2))
	to.remove = lapply(to.remove,function(x) c(x,x+1))
	temp = mapply(function(x,y) {
		if(length(y)>=1) {return(x[-y,])}
		if(length(y)<1) {return(x)}
		},temp,to.remove,SIMPLIFY=FALSE)
	tryout = lapply(temp,function(x) {
		emptpos = integer(nrow(x)+1)
		emptpos[1] = 0
		emptpos[seq(2,length(emptpos)-1,2)] = x$position[seq(1,nrow(x)-1,2)]-1
		emptpos[seq(3,length(emptpos),2)] = x$position[seq(2,nrow(x),2)]+1
		data.frame(chr=x$chr[1],strand=x$strand[1],position=emptpos,gene="intergenic",prom="intergenic")
	})
	tryout = rbindlist(tryout)
	PROMS = data.table(rbind(proms,tryout))
	setkeyv(PROMS,c("chr","strand","position"))
	return(PROMS)
}