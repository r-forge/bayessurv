# Survival probability function
Sx.fun         = function(x,th, idm=3, log=FALSE){
	Sx    =  x * (-th[,3]) + exp(th[,4])/th[,5] * (1-exp(th[,5]*x))
	if(idm==3) Sx = Sx + (exp(th[,1])/th[,2] * (exp(-th[,2]*x)-1))
	if(!log) Sx = exp(Sx)
	return(Sx)
} 
