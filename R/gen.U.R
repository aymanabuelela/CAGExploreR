gen.U <- function(C,P,B=1000,mu,cutoff,pvalue=T,seed=1,phi)
{
	lambda = mu[rep(1:P,C)] #uses promoter means
	phi = phi[rep(1:P,C)]
	size=1/phi
	set.seed(seed)
	if(sum(phi)>(1e-12)) Z = matrix(rnbinom(n=C*P*B,mu = lambda,size=size),nrow=C*P)
	if(sum(phi)<(1e-12)) Z = matrix(rpois(C*P*B,lambda = lambda),nrow=C*P)
	S = colSums(Z)
	Pij = t(t(Z)/S)
	c = matrix(NA,nrow=C,ncol=B)
	for(k in 1:C) c[k,] = t(t(colSums(Z[seq(1+P*(k-1),P*k),]))/S)	
	r = matrix(NA,nrow=P,ncol=B)
	for(k in 1:P) r[k,] = t(t(colSums(Z[seq(k,C*P-P+k,P),]))/S)
	
	RR = r[rep(1:P,C),]
	CC = c[rep(1:C,rep(P,C)),]

	temp = Pij*log((RR*CC)/Pij)
	temp[is.nan(temp)] = 0
	U.num = colSums(temp)

	#temp2 = c*log(c)
	temp2 = r*log(r)
	temp2[is.nan(temp2)] = 0
	U.den = colSums(temp2)
	U = U.num/U.den
	U[is.nan(U)] = 0
	U[U<0] = 0
	U[U==Inf] = 0
	if(!pvalue) return(U)

	#if(cutoff==1) return(sum(U==1)/length(U))

	if(sum(U)<(1e-12)) return(1) #sometimes always 0 U
	xbar = mean(U,na.rm=T)
	vbar = var(U,na.rm=T)
	alpha = xbar*((xbar*(1-xbar))/vbar - 1)
	beta = (1-xbar)*((xbar*(1-xbar))/vbar - 1)
	if(alpha<0 | beta<0) return(xbar) #this happens if U is 0 or 1 only
	epsilon = 1e6
	if(cutoff==1) return(pbeta(q=(cutoff-epsilon),shape1=alpha,shape2=beta,lower.tail=F))
	return(pbeta(q=cutoff,shape1=alpha,shape2=beta,lower.tail=F))
}

