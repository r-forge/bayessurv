c.low <-
function(th, idm=3){
	if(idm==1) cl = 0
	if(idm==2) cl = ifelse(th[5] > 0, -exp(th[4]), 0)
	if(idm==3){
		x.minf = (th[1]+log(th[2]) - th[4]-log(th[5]))/(th[2] + th[5])
		cl     = -exp(th[1]-th[2]*(x.minf)) - exp(th[4]+th[5]*(x.minf))
	}
	return(cl)
}

