# pdf of ages at death:
fx.fun         = function(x,th, idm=3, log=FALSE){
	fx       = mx.fun(x, th, idm=idm, log=TRUE) + Sx.fun(x, th, idm=idm, log=TRUE)
	if(!log) fx = exp(fx)
	return(fx)
}
