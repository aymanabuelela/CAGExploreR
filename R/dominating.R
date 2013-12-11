dominating = function(df)
{
	res = apply(df,2,function(x) (x==max(x) & x>0 ))
	temp = apply(res,1,any) & !apply(res,1,prod)
	return(paste(names(temp)[temp],collapse="|"))
}
