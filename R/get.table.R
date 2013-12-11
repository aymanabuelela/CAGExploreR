get.table = function(gene,data,theil,text=T)
{
	data = data[order(-rowSums(data)),]
	mydata = t(t(as.matrix(data))/colSums(as.matrix(data))) #percentages
	mydata.raw = as.matrix(data)
	mydata.raw = mydata.raw[,order(mydata[1,])]
	mydata = mydata[,order(mydata[1,])]
	all.combs.samples = unlist(apply(combn(ncol(mydata),2),2,list),recursive=F)
	all.combs.proms = unlist(apply(combn(nrow(mydata),2),2,list),recursive=F)
	all2by2  = mapply(function(x,y) mydata.raw[x,y],
		rep(all.combs.proms,length(all.combs.samples)),
		all.combs.samples[rep(1:length(all.combs.samples),rep(length(all.combs.proms),
		length(all.combs.samples)))]
	,SIMPLIFY=F)
	temp = do.call(rbind,all2by2) + 0.5
	odd.numbers = seq(1,nrow(temp)-1,2)
	even.numbers = odd.numbers+1
	odds.ratios = as.numeric(temp[odd.numbers,1]*temp[even.numbers,2]/(temp[odd.numbers,2]*temp[even.numbers,1]))
	log.odds.ratios = log(odds.ratios,base=2)
	log.odds.se = as.numeric(sqrt(1/temp[odd.numbers,1]+1/temp[odd.numbers,2]+1/temp[even.numbers,1]+1/temp[even.numbers,2]))
	log.odds.z = log.odds.ratios/log.odds.se
	log.odds.pvalue = 2*pnorm(abs(log.odds.z),lower.tail=F)
		
	if(text){
	prom1 = unlist(lapply(all2by2,function(x) rownames(x)[1]))
	prom2 = unlist(lapply(all2by2,function(x) rownames(x)[2]))
	lib1 = unlist(lapply(all2by2,function(x) colnames(x)[1]))
	lib2 = unlist(lapply(all2by2,function(x) colnames(x)[2]))
	comparison = paste(paste(prom1,prom2,sep="|"),paste(lib1,lib2,sep="|"),sep="#")
	to.return = data.frame(comparison,gene,theilU=theil,OR=odds.ratios,log2OR=log.odds.ratios,pvalue=log.odds.pvalue)
	to.return = to.return[order(abs(to.return$log2OR),decreasing=T),]
	return(to.return)	
	}

	to.return = data.frame(gene,theilU=theil,OR=odds.ratios,log2OR=log.odds.ratios,pvalue=log.odds.pvalue)
	to.return = to.return[order(abs(to.return$log2OR),decreasing=T),]
	return(to.return)	
}

