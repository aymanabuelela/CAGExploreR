convert = function(what,from,to,using=NA)
{
	if(key(using)[1]!=from) setkeyv(using,from)
	return(as.data.frame(using[what][,to,with=FALSE])[,1])
}