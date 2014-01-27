repagree = function(x,ncond){

samples=cluster=V1=NULL

	x = x[order(-rowSums(x)),]
	mydata = t(t(as.matrix(x))/colSums(as.matrix(x))) #percentages
	mydata[is.na(mydata)] = 0
	test = cutree(hclust(dist(t(mydata))),ncond)
	test2 = data.table(cluster=test,samples=Select(names(test),".",1))
	return(test2[,length(unique(samples)),by=cluster][,max(V1)])
}
