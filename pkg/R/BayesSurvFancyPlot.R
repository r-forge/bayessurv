BayesSurvFancyPlot    = function(outBS){

	thing      = outBS$thint
	nth        = 5
	nz         = ncol(outBS$Z)
	zname      = colnames(outBS$Z)
	nthn       = sum(outBS$modm[outBS$idm,])
	pname      = paste(rep(colnames(outBS$modm)[outBS$modm[outBS$idm,]==1],each=nz), "[",rep(zname, nthn),"]", sep="")
	parname    = colnames(outBS$modm)[outBS$modm[outBS$idm,]==1]
	expar     = expression(alpha[1],beta[1],c,alpha[2],beta[2])[outBS$modm[outBS$idm,]==1]
	Bord       = brewer.pal(12, "Paired")[round(seq(1,12, length=nz))]
	Cols       = adjustcolor(Bord, alpha.f=0.5)
	if(length(dim(outBS$theta))>2){
		nsim   = dim(outBS$theta)[3]
		parr   = matrix(0,0,nthn*nz); colnames(parr) = pname
		for(i in 1:nsim) parr  = rbind(parr, outBS$theta[,,i])
	} else {
		parr   = outBS$theta[outBS$thint,]
	}


	# Summary Survival and mortality functions:
	xm          = rep(NA, nz)
	for(i in 1:nz){
		xm[i]  = ceiling(max(outBS$xqsum[1,outBS$Z[,i]==1])*1.1)
	}
	S.x         = function(th) Sx.fun(xv, matrix(th,1,nth), idm=outBS$idm)
	m.x         = function(th) mx.fun(xv, matrix(th,1,nth), idm=outBS$idm)
	Sx          = list()
	mx          = list()
	for(i in 1:nz){
		xv       = seq(0,xm[i],0.1)
		pmat     = parr[,paste(parname,"[",zname[i],"]", sep="")]
		if(outBS$idm==1) pmat = cbind(0,0,0,pmat) else if(outBS$idm==2) pmat = cbind(0,0,pmat)
		Sx[[zname[i]]] = apply(apply(pmat,1,S.x),1, quantile, c(0.5,0.025,0.975))
		mx[[zname[i]]] = apply(apply(pmat,1,m.x),1, quantile, c(0.5,0.025,0.975))
	}
	ylmx       = round(c(0, ifelse(nz>1, max(unlist(mx)), max(mx))))
	mxv        = ceiling(max(xm)/5)*5
	
	layout(matrix(c(rep(1:nthn, each=2), rep(rep(c(nthn+1, nthn+2),each=nthn), 2)), nthn*2, 3), widths = rep(2, 3), heights=rep(1, nthn))
	par(mar=c(3,3,0.5,0.5))
	for(i in 1:nthn){
		xz     = list()
		dez    = list()
		ylz    = rep(NA, nz)
		xlz    = matrix(0,nz,2)
		for(j in 1:nz){
			xz[[zname[j]]]  = parr[,paste(parname[i],"[",zname[j],"]",sep="")]
			dez[[zname[j]]] = density(xz[[zname[j]]])
			ylz[j]          = max(dez[[zname[j]]]$y)
			xlz[j,]         = range(dez[[j]]$x)
		}
		xr     = range(xlz)
		xl     = c(floor(xr[1]*10)/10, ceiling(xr[2]*10)/10)
		xd     = ceiling(diff(xl)*10)/10
		plot(dez[[1]], xlab="", ylab="", xlim=xl, ylim=c(0, max(ylz)), lwd=3, axes=FALSE, main="", col=NA)
		for(j in 1:nz) polygon(c(dez[[j]]$x, dez[[j]]$x[1]), c(dez[[j]]$y, dez[[j]]$y[1]), col=Cols[j], border=Bord[j], lwd=1.5)
		axis(1, at=seq(xl[1], xl[2], length=5), line=0.5, labels=NA, tcl=0.4)
		axis(1, at=seq(xl[1], xl[2], length=3), lwd=NA)
		mtext(expar[i], 2, line=0, at=max(ylz)*0.8, las=2, cex=1.25)
	}

	par(mar=c(4, 7, 0.5, 0.5))
	# Plot survival probability:
	plot(c(0,mxv), range(0,1), col=NA, axes=FALSE, xlab="", ylab="")
	for(i in 1:nz){
		xv       = seq(0,xm[i],0.1)
		polygon(c(xv, rev(xv)), c(Sx[[i]][2,], rev(Sx[[i]][3,])), col=Cols[i], border=Bord[i])
		lines(xv, Sx[[i]][1,], col=Bord[i], lty=3)
	} 
	legend('topright', zname, pch=15, pt.cex=3, cex=1.5, col=Cols, bty='n')

	axis(1, seq(0,mxv, 5), labels=NA, tcl=0.4, line=0.5)
	axis(2, seq(0, 1, 0.2), tcl=0.4, las=2, cex.axis=1.2)
	mtext(expression(paste("Survival prob. ", italic(S(x)))), 2, line=3.5, cex=1.25)

	# Plot mortality rates:
	plot(c(0,mxv), ylmx, col=NA, axes=FALSE, xlab="", ylab="")
	for(i in 1:nz){
		xv       = seq(0,xm[i],0.1)
		polygon(c(xv, rev(xv)), c(mx[[i]][2,], rev(mx[[i]][3,])), col=Cols[i], border=Bord[i])
		lines(xv, mx[[i]][1,], col=Bord[i], lty=3)
	} 

	axis(1, seq(0,mxv, 5), labels=NA, tcl=0.4, line=0.5)
	axis(1, at=seq(0,mxv, 5), lwd=NA, cex.axis=1.2)
	axis(2, seq(0, ylmx[2], length=5), tcl=0.4, las=2, cex.axis=1.2)
	mtext(expression(paste("Mortality rate ", italic(mu(x)))), 2, line=3.5, cex=1.25)
	mtext(expression(paste("Age ", italic(x), " (years)")), 1, cex=1.25, line=3)

}