info2osc <-
function(x)
{
paste(paste(paste(x$chr,x$start,sep=":"),x$end,sep=".."),x$strand,sep=",")
}
