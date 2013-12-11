GREP <- function(pattern,x,...)
{
	temp = Vectorize("grep","pattern",SIMPLIFY=F)
	return(temp(pattern,x,...))
}