# Mortality rate function:
mx.fun     = function(x,th, idm=3, log=FALSE){
	mx    = th[,3] + exp(th[,4] + th[,5]*x)
	if(idm==3) mx = mx + exp(th[,1]-th[,2]*x)
	if(log) mx = log(mx)
	return(mx)
}
