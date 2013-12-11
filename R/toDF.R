toDF <-
function(x,param)
{
x = unname(x)
elts = setNames(bamWhat(param),bamWhat(param))
lst = lapply(elts,function(elt) my.unlist(lapply(x,"[[",elt)))
as.data.frame(lst)
}
