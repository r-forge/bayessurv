BayesSurvPlots <-
function(outBS){
	require(RColorBrewer)
	if(.Platform$OS.type=="unix") devtype=quartz else devtype=windows
for(i in dev.list()) dev.off()
# PLOTS:
# Trace plots for parameters:
if(!is.element("2", dev.list())) devtype(width=8, height=6)
dev.set(2)
par(mfrow=c(sum(outBS$modm[outBS$idm,]), ncol(outBS$Z)), mar=c(3,2,2,1))
for(i in 1:ncol(outBS$theta)) plot(outBS$theta[,i], type='l', xlab="Iteration", ylab="", main=colnames(outBS$theta)[i])

if(!is.element("3", dev.list())) devtype(width=8, height=6)
dev.set(3)
par(mfrow=c(ceiling(ncol(outBS$pi)/2+1),2), mar=c(3,2,2,1))
for(i in 1:ncol(outBS$pi)) plot(outBS$pi[,i], type='l', xlab="", ylab="", main=paste("pi", i))
for(i in 1:ncol(outBS$post)) plot(outBS$post[,i], type='l', xlab="Iteration", ylab=expression(p(theta|X), p(X[0]|X[1],theta,pi))[i])


# PLOTS:
# Resulting survival and mortality:
if(!is.element("4", dev.list())) devtype(width=6, height=6)
dev.set(4)
xv         = seq(0, max(outBS$xqsum), 0.1)
Bord       = brewer.pal(12, "Paired")[round(seq(1,12, length=ncol(outBS$Z)))]
Cols       = adjustcolor(Bord, alpha.f=0.5)
par(mfrow=c(2,1))
plot(range(xv), c(0,1), col=NA, xlab="Age", ylab=expression(S(x)), main="Survival probability")
for(i in 1:ncol(outBS$Z)) polygon(c(xv, rev(xv)), c(outBS$Sxsum[[i]][2,], rev(outBS$Sxsum[[i]][3,])), col=Cols[i], border=Bord[i])
legend('topright', colnames(outBS$Z), pch=15, pt.cex=2, col=Cols, bty='n')

plot(range(xv), c(0,max(unlist(outBS$mxsum))), col=NA, xlab="Age", ylab=expression(mu(x)), main="Mortality rate")
for(i in 1:ncol(outBS$Z)) polygon(c(xv, rev(xv)), c(outBS$mxsum[[i]][2,], rev(outBS$mxsum[[i]][3,])), col=Cols[i], border=Bord[i])

}

