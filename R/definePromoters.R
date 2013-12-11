definePromoters <-
function(file=NA,sep="\t",header=T)
{
	if(!is.character(file)) prom.file = file
	if(is.character(file)) prom.file = read.table(file,sep=sep,header=header)
	n.promoters = nrow(prom.file)
	prom.file.pos = prom.file[prom.file$strand=="+",]
	prom.file.neg = prom.file[prom.file$strand=="-",]

	CHROMOSOMES = unique(prom.file$chr)
	PROMOTERS = list();length(PROMOTERS) = 2;names(PROMOTERS) = c("+","-")
	PROMOTERS[[1]] = list();length(PROMOTERS[[1]]) = length(CHROMOSOMES);names(PROMOTERS[[1]]) = CHROMOSOMES
	PROMOTERS[[2]] = list();length(PROMOTERS[[2]]) = length(CHROMOSOMES);names(PROMOTERS[[2]]) = CHROMOSOMES

	for(i in 1:length(CHROMOSOMES))
	{
		PROMOTERS[["+"]][[i]] = as.matrix(prom.file.pos[prom.file.pos$chr==CHROMOSOMES[i],c("start","end")])
		PROMOTERS[["-"]][[i]] = as.matrix(prom.file.neg[prom.file.neg$chr==CHROMOSOMES[i],c("start","end")])
	}
	IR.pos = lapply(PROMOTERS[["+"]],function(x) IRanges(start=x[,1],end=x[,2]))
	IR.neg = lapply(PROMOTERS[["-"]],function(x) IRanges(start=x[,1],end=x[,2]))
	
	which.pos = RangesList(IR.pos)
	which.neg = RangesList(IR.neg)
	return(list(which.pos=which.pos,which.neg=which.neg,CHROMOSOMES=CHROMOSOMES,n.promoters=n.promoters,original=prom.file))
}
