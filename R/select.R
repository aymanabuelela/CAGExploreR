Select = function(x,split,which)
{
	sapply(strsplit(x,split,fixed=TRUE),function(x) x[which])
}
