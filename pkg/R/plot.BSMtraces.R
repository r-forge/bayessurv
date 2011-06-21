plot.BSMtraces <-
function(outBS, tracename = "theta"){
	require(RColorBrewer)
	if(.Platform$OS.type=="unix") devtype=quartz else devtype=windows
	plotname     = c("theta","pi", "post")
	if(!is.element(tracename, plotname)) stop(paste("Wrong 'tracename' argument. Valid arguments are:", paste(paste("'",plotname,"'", sep=""), collapse=", "),".\n"), call.=FALSE)

	# PLOTS:
	# Trace plots for parameters:
	p           = which(plotname == tracename)
	nz          = ncol(outBS$Z)
	ydim        = c(sum(outBS$modm[outBS$idm,]), ceiling((ncol(outBS$pi))/2), 2)
	xdim        = c(nz, 2, 2)

	x           = outBS[[tracename]]
	Main        = list(theta = colnames(outBS$theta), pi = colnames(outBS$pi), post = expression(paste(p,"(",theta," | ",X,")"), paste(p,"(",X[0]," | ",X[1],",",theta,",",pi,",",")"), paste(p,"(",X[0],",",theta,",",pi,"| ... ",")")))
	if(ncol(x)==1){
		Main$pi = "pi"
		ydim[2] = 1
		xdim[2] = 1
	} 

	par(mfrow=c(ydim[p], xdim[p]), mar=c(3,3,2,1))
	for(i in 1:ncol(x)){
		
		if(length(dim(x))<3){
			plot(x[,i], type='l', xlab="Iteration", ylab="", main=Main[[p]][i], frame.plot=FALSE)
		} else {
			yl    = range(x[,i,], na.rm=TRUE)
			plot(c(1,nrow(x)), yl, col=NA, xlab="Iteration", ylab="", main=Main[[p]][i], frame.plot=FALSE)
			for(j in 1:dim(x)[3]) lines(x[,i,j], type='l', col=brewer.pal(12, "Paired")[j])
		}
	}
}

