BayesSurvPlots <-
function(outBS){
	require(RColorBrewer)
	if(.Platform$OS.type=="unix") devtype=quartz else devtype=windows
	for(i in dev.list()) dev.off()

	# PLOTS:
	# Trace plots for parameters:
	parname     = c("theta","pi")
	np          = length(parname)
	nz          = ncol(outBS$Z)
	nth         = sum(outBS$modm[outBS$idm,])
	npi         = ceiling((ncol(outBS$pi)+3)/2)
	postname    = expression(paste(p,"(",theta," | ",X,")"), paste(p,"(",X[0]," | ",X[1],",",theta,",",pi,",",")"), paste(p,"(",X[0],",",theta,",",pi,"| ... ",")"))
	for(p in 1:np){
		x        = outBS[[parname[p]]]
		if(!is.element(as.character(p+1), dev.list())) devtype(width=8, height=6)
		dev.set(as.character(p+1))
		par(mfrow=c(c(nth, npi)[p], nz), mar=c(3,3,2,1))
		for(i in 1:ncol(x)){
			if(ncol(x)==1 ) cname = "pi" else cname = colnames(x)
			if(length(dim(x))<3){
				plot(x[,i], type='l', xlab="Iteration", ylab="", main=cname[i], frame.plot=FALSE)
			} else {
				yl    = range(x[,i,])
				plot(c(1,nrow(x)), yl, col=NA, xlab="Iteration", ylab="", main=cname[i], frame.plot=FALSE)
				for(j in 1:dim(x)[3]) lines(x[,i,j], type='l', col=brewer.pal(12, "Paired")[j])
			}
		}
		if(p==np){			
			x        = outBS$post
			for(i in 1:ncol(x)){
				if(length(dim(x))<3){
					plot(x[,i], type='l', xlab="Iteration", main=postname[i])
				} else {
					yl       = range(x[,i,])
					plot(c(1,nrow(x)), yl, col=NA, xlab="Iteration", main=postname[i],  frame.plot=FALSE)
					for(j in 1:dim(x)[3]) lines(x[,i,j], type='l', col=brewer.pal(12, "Paired")[j])
				} 
			}
		}		
	}
	


	# PLOTS:
	# Resulting survival and mortality:
	if(!is.element("4", dev.list())) devtype(width=6, height=6)
	dev.set(4)
	xv         = seq(0, max(outBS$xqsum), 0.1)
	Bord       = brewer.pal(12, "Paired")[round(seq(1,12, length=ncol(outBS$Z)))]
	Cols       = adjustcolor(Bord, alpha.f=0.5)
	mx         = outBS$mxsum
	Sx         = outBS$Sxsum
	ylmx       = c(0, ifelse(nz>1, max(unlist(mx)), max(mx)))
	
	par(mfrow=c(2,1))
	plot(range(xv), c(0,1), col=NA, xlab="", ylab=expression(S(x)), main="Survival probability", frame.plot=FALSE)
	if(nz>1){
		for(i in 1:nz) polygon(c(xv, rev(xv)), c(Sx[[i]][2,], rev(Sx[[i]][3,])), col=Cols[i], border=Bord[i])
		legend('topright', colnames(outBS$Z), pch=15, pt.cex=2, col=Cols, bty='n')
	} else {
		polygon(c(xv, rev(xv)), c(Sx[2,], rev(Sx[3,])), col=Cols[1], border=Bord[1])
	}

	plot(range(xv), ylmx, col=NA, xlab="Age (x)", ylab=expression(mu(x)), main="Mortality rate", frame.plot=FALSE)
	if(nz>1){
		for(i in 1:nz) polygon(c(xv, rev(xv)), c(mx[[i]][2,], rev(mx[[i]][3,])), col=Cols[i], border=Bord[i])	
	} else {
		polygon(c(xv, rev(xv)), c(mx[2,], rev(mx[3,])), col=Cols[1], border=Bord[1])
	}

}

