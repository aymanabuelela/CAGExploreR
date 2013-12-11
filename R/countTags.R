countTags <-
function(bam.files,bai.files,ids,prom.defs)
{
	results = matrix(nrow=prom.defs[["n.promoters"]],ncol=length(bam.files))
	which.pos = prom.defs[["which.pos"]]
	which.neg = prom.defs[["which.neg"]]
	CHROMOSOMES = prom.defs[["CHROMOSOMES"]]
	depth = numeric(length(ids));names(depth) = ids
	for(i in 1:length(bam.files))
	{
		counts = list();length(counts) = 2;names(counts) = c("+","-")
		counts[["+"]] = counts[["-"]] = list()
		counter.pos = counter.neg = 1
		depth[i] = 0
		for(chr in 1:length(CHROMOSOMES))
		{
			no.chr.pos = no.chr.neg = F
	
			param.pos = ScanBamParam(flag = scanBamFlag(isMinusStrand=F,isUnmappedQuery=F,isNotPrimaryRead=F,isNotPassingQualityControls=F),which=which.pos[chr],what=c("pos"))
			param.neg = ScanBamParam(flag = scanBamFlag(isMinusStrand=T,isUnmappedQuery=F,isNotPrimaryRead=F,isNotPassingQualityControls=F),which=which.neg[chr],what=c("pos","qwidth"))
			tss.pos = tryCatch({toDF(scanBam(bam.files[i],index=bai.files[i], param=param.pos),param.pos)[,1]},error=function(e) no.chr.pos<<-T)
			test.neg = tryCatch({toDF(scanBam(bam.files[i],index=bai.files[i], param=param.neg),param.neg)},error=function(e) no.chr.neg<<-T)

			if(no.chr.pos & length(which.pos[chr][[1]])>0) 
			{
				counts[["+"]][[counter.pos]] = data.frame(chr=as.character(CHROMOSOMES[chr]),strand="+",start=as.data.frame(which.pos[[chr]])$start,end=as.data.frame(which.pos[[chr]])$end,count=0)
				counter.pos = counter.pos+1
			}
			if(no.chr.neg & length(which.neg[chr][[1]])>0) 
			{
				counts[["-"]][[counter.neg]] = data.frame(chr=as.character(CHROMOSOMES[chr]),strand="-",start=as.data.frame(which.neg[[chr]])$start,end=as.data.frame(which.neg[[chr]])$end,count=0)
				counter.neg = counter.neg+1
			}
			if(no.chr.pos & no.chr.neg) next
	
			tss.neg = as.integer(test.neg$pos + test.neg$qwidth - 1)

			tss.pos.IR = IRanges(start=tss.pos,width=1)
			tss.neg.IR = IRanges(start=tss.neg,width=1)
			depth[i] = depth[i] + length(tss.pos.IR) + length(tss.neg.IR)			

			n.blocks.pos = length(tss.pos.IR)/5e5
			counts.pos = matrix(nrow=length(which.pos[[chr]]),ncol=ceiling(n.blocks.pos))
			for(block in 1:ceiling(n.blocks.pos)) 
			{
				A = 1+5e5*(block-1)
				if(block<ceiling(n.blocks.pos)) B = 5e5*block
				if(block==ceiling(n.blocks.pos)) B = length(tss.pos.IR)
				counts.pos[,block] = countOverlaps(which.pos[[chr]],tss.pos.IR[A:B,])
			}
			counts.pos = as.integer(rowSums(counts.pos))
			#counts.pos = countOverlaps(which.pos[[chr]],tss.pos.IR)
			
			n.blocks.neg = length(tss.neg.IR)/5e5
			counts.neg = matrix(nrow=length(which.neg[[chr]]),ncol=ceiling(n.blocks.neg))
			for(block in 1:ceiling(n.blocks.neg)) 
			{
				A = 1+5e5*(block-1)
				if(block<ceiling(n.blocks.neg)) B = 5e5*block
				if(block==ceiling(n.blocks.neg)) B = length(tss.neg.IR)
				counts.neg[,block] = countOverlaps(which.neg[[chr]],tss.neg.IR[A:B,])
			}
			counts.neg = as.integer(rowSums(counts.neg))			
			#counts.neg = countOverlaps(which.neg[[chr]],tss.neg.IR)

			if(length(which.pos[[chr]])>0) counts[["+"]][[counter.pos]] = data.frame(chr=as.character(CHROMOSOMES[chr]),strand="+",start=as.data.frame(which.pos[[chr]])$start,end=as.data.frame(which.pos[[chr]])$end,count=counts.pos)
			if(length(which.neg[[chr]])>0) counts[["-"]][[counter.neg]] = data.frame(chr=as.character(CHROMOSOMES[chr]),strand="-",start=as.data.frame(which.neg[[chr]])$start,end=as.data.frame(which.neg[[chr]])$end,count=counts.neg)

			if(length(which.pos[[chr]])>0) counter.pos = counter.pos + 1
			if(length(which.neg[[chr]])>0) counter.neg = counter.neg + 1
		}

		#tag.counts = do.call(rbind,unlist(counts,recursive=F))
		tag.counts = rbindlist(unlist(counts,recursive=F))
		rownames(tag.counts) = NULL
		tag.counts = as.data.frame(tag.counts[order(tag.counts$chr,tag.counts$strand,tag.counts$start),])
		results[,i] = as.integer(tag.counts$count)
	}
	results = data.frame(tag.counts[,-grep("count",colnames(tag.counts))],results)
	colnames(results) = c("chr","strand","start","end",ids)
	temp = prom.defs[["original"]]
	#temp$strand = relevel(temp$strand,"+")
	temp = temp[order(temp$chr,temp$strand,temp$start,temp$end),]
	gene = temp[,-unlist(GREP(c("chr","strand","start","end"),colnames(temp)))]
	results = cbind(results,gene)
	return(list(depth=depth,counts=results))
	#save(tag.counts,file="tag.counts.RData")
}
