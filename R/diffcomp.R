diffcomp <-
function(x,detailed=F,top="all",gene=NA,B=1000,seed=1,mc.adjust="fdr",text=TRUE) 
{
ncond = length(unique(Select(as.character(x$samples$group),".",1)))
lib.counts = x$pooled.samples[,"effective.lib.size"]
n.lib = length(lib.counts)

all.data.w.phi = split(x$pooled.sub.counts,x$pooled.sub.counts$gene)
all.data.w.phi = lapply(all.data.w.phi,function(x) as.matrix(x[,-match(c("gene"),colnames(x))]))
all.data = lapply(all.data.w.phi,function(x) as.matrix(x[,-match("dispersion",colnames(x))]))

if(!is.null(x$pooled.super.counts)) coverage.data = x$pooled.super.counts[,-match("dispersion",colnames(x$pooled.super.counts))]

if(!detailed){
all.data.prop = lapply(all.data,function(x) x/sum(x))
all.data.row.prop = lapply(all.data,function(x) rowSums(x)/sum(x))
all.data.col.prop = lapply(all.data,function(x) colSums(x)/sum(x))

theilU.num = mapply(function(r,p,c) {
	temp = (log(diag(1/r)%*%p%*%diag(1/c)))*p
	temp[is.nan(temp)]=0
	temp=-sum(temp)
	return(temp)},all.data.row.prop,all.data.prop,all.data.col.prop)
theilU.den = unlist(lapply(all.data.row.prop,function(c) {
	temp = entropy(c)
	temp[is.nan(temp)]=0
	return(sum(temp))
}))
theilU = theilU.num/theilU.den
theilU[is.nan(theilU)] = 0 #all zeroes, or only one promoter with nonzero counts
theilU[theilU<0] = 0 #one condition expressed, the other isnt
theilU[theilU==Inf] = 0 #only one promoter expressed

u.pvalues = mapply(function(x,u) {
	dat.only = x[,-grep("dispersion",colnames(x))] 
	rs = rowSums(dat.only)
	rs[rs<1] = 1
	mr = max(rs)
	norm.x = dat.only * mr / rs
	to.ret = gen.U(C=(ncol(x)-1),P=nrow(x),B=B,mu=rowMeans(norm.x),cutoff=u,seed=seed,phi=x[,"dispersion"])
	return(to.ret)
	#gen.U(C=(ncol(x)-1),P=nrow(x),B=B,mu=rowMeans(x[,-grep("dispersion",colnames(x))]),cutoff=u,seed=seed,phi=x[,"dispersion"])
	#gen.U(C=(ncol(x)-1),P=nrow(x),B=B,mu=as.numeric(t(x[,-grep("dispersion",colnames(x))])),cutoff=u,seed=seed,phi=x[,"dispersion"])
	#gen.U(C=(ncol(x)-1),P=nrow(x),B=B,mu=rowSums(x[,-grep("dispersion",colnames(x))]),cutoff=u,seed=seed,phi=x[,"dispersion"])
	#gen.U(C=(ncol(x)-1),P=nrow(x),B=B,mu=rowMeans(x[,-grep("phi",colnames(x))]),cutoff=u,seed=seed,phi=rep(0.8,nrow(x)))
},all.data.w.phi,theilU)

u.qvalues = p.adjust(u.pvalues,mc.adjust)

geneHetero = unlist(mapply(function(X,y,z) 
{	if(!y) gene.prop = colSums(X)/lib.counts
	if(y) gene.prop = t(t(coverage.data[match(z,coverage.data$gene),-match("gene",colnames(coverage.data))])/lib.counts)
	temp=gene.prop/sum(gene.prop)
	return(1 - sum(entropy(temp))/log(1/n.lib))
},all.data,as.list(names(all.data)%in%x$pooled.super.counts$gene),names(all.data),SIMPLIFY=F))#0: gene expressed equally across all samples, 1: DEG

Coverage = unlist(mapply(function(X,y,z) 
{	if(!y) return(NA)
	gene.total = as.numeric(coverage.data[match(z,coverage.data$gene),-match("gene",colnames(coverage.data))])
	return(mean(colSums(X)/gene.total,na.rm=T))
},all.data,as.list(names(all.data)%in%x$pooled.super.counts$gene),names(all.data),SIMPLIFY=F))

qvalue.text = mc.adjust
if(mc.adjust=="none") qvalue.text = "pvalue"

Dominating = unlist(lapply(all.data,function(x) dominating(x)))
RepAgree = unlist(lapply(all.data,function(x) repagree(x,ncond)))

significance.pars = data.frame(entropy.Reduction=theilU,pvalue=u.pvalues,U.FDR=u.qvalues,geneHetero=geneHetero,coverage=Coverage,dominant.promoter.switch=Dominating,RepAgree=RepAgree)
colnames(significance.pars)[3] = qvalue.text
significance.pars = significance.pars[order(-significance.pars[,1],significance.pars[,2]),]
return(significance.pars)
} #end if(!detailed)

if(detailed) {

all.data.prop = lapply(all.data,function(x) x/sum(x))
all.data.row.prop = lapply(all.data,function(x) rowSums(x)/sum(x))
all.data.col.prop = lapply(all.data,function(x) colSums(x)/sum(x))

theilU.num = mapply(function(r,p,c) {
temp = (log(diag(1/r)%*%p%*%diag(1/c)))*p
temp[is.nan(temp)]=0
temp=-sum(temp)
return(temp)},all.data.row.prop,all.data.prop,all.data.col.prop)
theilU.den = unlist(lapply(all.data.row.prop,function(c) {
temp = entropy(c)
temp[is.nan(temp)]=0
return(sum(temp))
}))
theilU = theilU.num/theilU.den
theilU[is.nan(theilU)] = 0
theilU[theilU<0] = 0 #one condition expressed, the other isnt
theilU[theilU==Inf] = 0 #only one promoter expressed

if(top=="all") top=length(all.data)
if(is.na(gene[1])) {results.list = mapply(function(x,y,z,v) {get.table(gene=x,data=z,theil=y,text=v)}
				,x=as.list(names(theilU[1:top])),y=as.list(theilU[1:top]),z=all.data[1:top],v=text,SIMPLIFY=F)}
if(!is.na(gene[1])) {results.list = mapply(function(x,y,z,v) {get.table(gene=x,data=z,theil=y,text=v)}
				,x=as.list(gene),y=as.list(theilU[match(gene,names(theilU))]),z=all.data[match(gene,names(all.data))],v=text,SIMPLIFY=F)}

results.list = as.data.frame(rbindlist(results.list))
results.list$local.qvalue = p.adjust(results.list$pvalue,mc.adjust)
qvalue.text = mc.adjust
if(mc.adjust=="none") qvalue.text = "pvalue"
colnames(results.list)[ncol(results.list)] = qvalue.text
results.list = results.list[order(-results.list$theilU,results.list$gene,-abs(results.list$log2OR)),]
rownames(results.list) = NULL
return(results.list)
}

}
