countTags = function(bamfiles,idx,ids,GENES,PROMS)
{
prom=gene=width=cigar=revstrand=strand=.BY=chr=J=N=.N=NULL


all.regions = data.table(prom=unique(c(unique(GENES$prom),unique(PROMS$prom))))
all.regions = all.regions[-which(all.regions$prom=="intergenic")]
setkey(all.regions,"prom")
genes.n.promoters = data.table(rbind(GENES,PROMS))
setkey(genes.n.promoters,"prom")
table.start = genes.n.promoters[all.regions,mult="first"]
table.start = table.start[,list(prom,gene)]
setkey(table.start,"prom")

counter = 1
for(bam in bamfiles)
{
	reader<-bamReader(bam)
	load.index(reader,idx[counter])
	rdf<-getRefData(reader)
	
	for(i in 1:nrow(rdf))
	{
	coords<-as.integer(c(rdf$ID[i],0,rdf[rdf$ID==rdf$ID[i],"LN"]))
	range<-bamRange(reader,coords)
	temp = data.table(as.data.frame(range)[,c("position","revstrand","cigar","flag")])
	setkeyv(temp,"revstrand","position")
	temp[,width := sapply(strsplit(temp[,cigar],"M|I|D"),function(x) sum(as.integer(x)))]
	temp[,position := (position + revstrand*(width-1))]
	temp[,strand := ifelse(.BY,"-","+"), by = revstrand]
	temp[,c("revstrand","flag","width","cigar") := NULL]
	temp[,chr := rdf$SN[i] ]
	setkeyv(temp,c("chr","strand","position"))
	gene.counts = GENES[J(temp),roll=TRUE,allow.cartesian=TRUE][,.N,by=list(chr,strand,gene,prom)]
	prom.counts = PROMS[J(temp),roll=TRUE,allow.cartesian=TRUE][,.N,by=list(chr,strand,gene,prom)]
	rm(range,temp);gc()
	if(i==1) results = rbindlist(list(gene.counts,prom.counts))
	if(i>1) results = rbindlist(list(results,gene.counts,prom.counts))
	rm(gene.counts,prom.counts);gc()
	}
	results = results[,list(prom,N)]
	setnames(results,"N",ids[counter])
	#setkeyv(results,c("chr","strand","gene"))
	setkey(results,"prom")
	depth = sum(results[,2,with=FALSE])
	results = results[table.start]
	suppressWarnings({for (i in seq_along(results)) set(results, i=which(is.na(results[[i]])), j=i, value=0)})
	results = results[,c(1,3,2),with=FALSE]
	bamClose(reader)
	if(counter==1) {to.return = list()
		to.return$depth = numeric()
		to.return$depth[counter] = depth
		to.return$counts = results}
	if(counter>1) {to.return$depth[counter] = depth
		to.return$counts = cbind(to.return$counts,results[,3,with=FALSE])}
	rm(results);gc()
	counter = counter + 1
}	
names(to.return$depth) = ids
setnames(to.return$counts,"prom","region")
return(to.return)
}
