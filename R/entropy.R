entropy <-
function(x)
{
	ent = x*log(x)
	ent[is.nan(ent)] = 0
	return(ent)
}
