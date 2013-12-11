pool = function(x)
{
#currently requires super regions
#currently assumes replicates exist even if they don't, should work!

	x$samples$effective.lib.size = round(x$samples$lib.size * x$samples$norm.factors)
	if(is.null(x$common.dispersion)) x$common.dispersion = 0
	if(is.null(x$tagwise.dispersion)) x$tagwise.dispersion = rep(x$common.dispersion,nrow(x$counts))

	if(length(grep("@",x$genes$gene))>0) super.regions.exist = TRUE
	if(length(grep("@",x$genes$gene))==0) super.regions.exist = FALSE

	if(super.regions.exist) {
		super.regions = grep("@",x$genes$gene)
		x$sub.counts = x$counts[-super.regions,]
		rownames(x$sub.counts) = NULL
		x$sub.genes = x$genes[-super.regions,]
		x$sub.genes$gene = as.character(x$sub.genes$gene)
		rownames(x$sub.genes) = NULL
		if(!is.null(x$tagwise.dispersion)) {
			x$sub.tagwise.dispersion = x$tagwise.dispersion[-super.regions]
			x$super.tagwise.dispersion = x$tagwise.dispersion[super.regions]
		}
		x$super.counts = x$counts[super.regions,]
		rownames(x$super.counts) = NULL
		x$super.genes = x$genes[super.regions,]
		x$super.genes$gene = sapply(strsplit(as.character(x$super.genes$gene),"@"),function(x) x[1])
		rownames(x$super.genes) = NULL

		if(max(table(x$samples$group))>1) replicates.exist = TRUE
		if(max(table(x$samples$group))<=1) replicates.exist = FALSE

		#if(replicates.exist) {
		if(TRUE) {
			conditions = as.character(unique(x$samples$group))
			temp = list();length(temp) = length(conditions)
			for(i in 1:length(temp)) {
				temp[[i]] = rowSums(x$sub.counts[,which(x$samples$group %in% conditions[i]),drop=F])
			}
			x$pooled.sub.counts = do.call(cbind,temp)
			colnames(x$pooled.sub.counts) = conditions
			x$pooled.sub.counts = data.frame(x$pooled.sub.counts,gene=x$sub.genes$gene,stringsAsFactors=FALSE,dispersion=x$sub.tagwise.dispersion)
			rownames(x$pooled.sub.counts) = info2osc(x$sub.genes)
			x$sub.genes = NULL
			x$sub.tagwise.dispersion = NULL

			temp = list();length(temp) = length(conditions)
			for(i in 1:length(temp)) {
				temp[[i]] = rowSums(x$super.counts[,which(x$samples$group %in% conditions[i]),drop=F])
			}
			x$pooled.super.counts = do.call(cbind,temp)
			colnames(x$pooled.super.counts) = conditions
			x$pooled.super.counts = data.frame(x$pooled.super.counts,gene=x$super.genes$gene,stringsAsFactors=FALSE,dispersion=x$super.tagwise.dispersion)
			rownames(x$pooled.super.counts) = info2osc(x$super.genes)
			x$super.genes = NULL
			x$super.tagwise.dispersion = NULL

			temp = split(x$samples[,-match(c("group","norm.factors"),colnames(x$samples))],as.character(x$samples$group))
			x$pooled.samples = do.call(rbind,lapply(temp,function(x) colSums(x)))
		}
	}

	if(!super.regions.exist) {
		if(max(table(x$samples$group))>1) replicates.exist = TRUE
		if(max(table(x$samples$group))<=1) replicates.exist = FALSE

		#if(replicates.exist) {
		if(TRUE) {
			conditions = as.character(unique(x$samples$group))
			temp = list();length(temp) = length(conditions)
			for(i in 1:length(temp)) {
				temp[[i]] = rowSums(x$counts[,which(x$samples$group %in% conditions[i]),drop=F])
			}
			x$pooled.counts = do.call(cbind,temp)
			colnames(x$pooled.counts) = conditions
			x$pooled.counts = data.frame(x$pooled.counts,gene=as.character(x$genes$gene),stringsAsFactors=FALSE,dispersion=x$tagwise.dispersion)
			rownames(x$pooled.counts) = info2osc(x$genes)
			x$genes = NULL
			x$tagwise.dispersion = NULL

			temp = split(x$samples[,-match(c("group","norm.factors"),colnames(x$samples))],as.character(x$samples$group))
			x$pooled.samples = do.call(rbind,lapply(temp,function(x) colSums(x)))
		}
	}
	
	return(x)
}