plot.SM <-
function(outBS){
	
	if(length(which(outBS$full.run)) != length(outBS$full.run)){
		stop("MCMC runs on BaSuMo object did not finish.\n Survival and mortality plots cannot be constructed, verify model input and run again\n", call.=FALSE)
	}

	requires(RColorBrewer)
	
	# PLOTS:
	# Resulting survival and mortality:
#	if(!is.element("4", dev.list())) devtype(width=6, height=6)
#	dev.set(4)
	xv         = seq(0, max(outBS$xqsum), 0.1)
	Bord       = brewer.pal(12, "Paired")[round(seq(1,12, length=ncol(outBS$Z)))]
	Cols       = adjustcolor(Bord, alpha.f=0.5)
	mx         = outBS$mxsum
	Sx         = outBS$Sxsum
	ylmx       = c(0, ifelse(nz>1, max(unlist(mx)), max(mx)))
	
	par(mfrow=c(2,1), mar=c(4,4,3,2))
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

