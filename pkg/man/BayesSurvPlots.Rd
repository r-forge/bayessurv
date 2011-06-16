\name{BayesSurvPlots}
\alias{BayesSurvPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to produce summary plots from a BayesSurv analysis}
\description{

}
\usage{
BayesSurvPlots(outBS)
}
\arguments{
  \item{outBS}{A BayesSurv object ouput from \code{\link{BayesSurv}} or \code{\link{BayesSurvParall}}
}
}
\details{
This function plots parameter traces from the MCMC analysis, the survival and mortality trajectories and ...
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Colchero and Clark 2011 Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. Journal of Animal Ecology. 
Colchero, Jones and Rebke. BayesSurv - a package to estimate mortality parameters from incomplete mark-recapture data with covariates. Unpublished Manuscript.
}
\author{
Fernando Colchero \email{colchero@mpg.demogr.de}
}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BayesSurv}}, \code{\link{BayesSurvParall}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ methods }
\keyword{ hplot }
