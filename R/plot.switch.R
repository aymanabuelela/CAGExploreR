plot.switch <-
function(x,lib.counts,coverage,GENE)
{
x = x[order(-rowSums(x)),]
mydata = t(t(as.matrix(x))/colSums(as.matrix(x))) #percentages
mydata[is.na(mydata)] = 0
mydata.raw = as.matrix(x)
dendro = as.dendrogram(hclust(dist(t(mydata))))
mydata.raw = mydata.raw[,order.dendrogram(dendro)]
my.lib.counts = lib.counts[order.dendrogram(dendro)]
mydata = mydata[,order.dendrogram(dendro)]

layout(matrix(c(1,2,3,4,4,4),2,3,byrow=T),widths=c(1,2,2),heights=c(1,1))

par(omd=c(0.01,0.99,0,0.95),mar=c(3,0.5,2,1),mgp=c(1.5,0.5,0))
#legend("center",inset=c(0,0),cex=1,fill=rainbow(nrow(mydata)),bty="n",legend=rownames(mydata))
plot(dendro,horiz=T,yaxt="n",leaflab="none")
par(mar=c(3,3,2,0))
b=barplot(mydata,las=2,horiz=T,col=rainbow(nrow(mydata)),main="Promoters",xlab="proportion",space=0)

ypos <- apply(mydata, 2, cumsum)
ypos <- ypos - mydata/2
ypos <- t(ypos)

text(ypos,b,t(mydata.raw),cex=ifelse(round(t(mydata)*100)>5,1,0.001),
	font=ifelse(t(apply(round(t(mydata)*100),1,function(x) x==max(x))),2,1),
	col=ifelse(t(apply(round(t(mydata)*100),1,function(x) x==max(x))),"white","black"))

par(mar=c(3,3,2,0))
if(!GENE %in% coverage$gene) barplot(colSums(mydata.raw)/my.lib.counts*1e6,las=2,horiz=T,names.arg="",main="Gene",xlab="tpm",space=0)
if(GENE %in% coverage$gene) {
coverage.temp = as.numeric(coverage[coverage$gene==GENE,-match("gene",colnames(coverage))])
barplot(coverage.temp[order.dendrogram(dendro)]/my.lib.counts*1e6,las=2,horiz=T,names.arg="",main="Gene",xlab="tpm",space=0)
}
title(paste(GENE),outer=T)
prom.data = data.frame(CAGExploreR:::osc2info(rownames(x)))

plot.gene(gene=GENE,my.prom.defs=prom.data,counts=x,lib.counts=lib.counts)

return(rownames(mydata))
}
