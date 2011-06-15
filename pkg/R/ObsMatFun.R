ObsMatFun <-
function(f, l, Tm){
	Fm    = Tm - f; Fm[Fm>=0] =  1; Fm[Fm<0] = 0
	Lm    = Tm - l; Lm[Lm<=0] = -1; Lm[Lm>0] = 0
	return(Fm * (-Lm))	
}

