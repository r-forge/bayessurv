\name{TestIniPars}
\alias{TestIniPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TestIniPars(Data, ini.pars = NULL, ststart, stend, model = "SI", jumps = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{ini.pars}{
%%     ~~Describe \code{ini.pars} here~~
}
  \item{ststart}{
%%     ~~Describe \code{ststart} here~~
}
  \item{stend}{
%%     ~~Describe \code{stend} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{jumps}{
%%     ~~Describe \code{jumps} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Data, ini.pars=NULL, ststart, stend, model="SI", jumps=NULL){

	# Data formating:
	Ti          = ststart
	Tf          = stend
	st          = Ti:Tf
	nt          = length(st)
	idnames     = Data[,1]
	n           = nrow(Data)
	bd          = as.matrix(Data[,2:3])
	Y           = as.matrix(Data[,1:nt+3]); colnames(Y) = st
	if(ncol(Data)>nt+3){
		Z  = as.matrix(Data[,(nt+4):ncol(Data)])
	} else {
		Z  = matrix(1, n, 1)
	}
	nth         = 5
	nz          = ncol(Z)

	# Model Matrix and boundary values for parameters:
	nth         = 5
	modm        = matrix(1, 3, nth)
	modm[1,1:3] = 0
	modm[2,1:2] = 0
	dimnames(modm) = list(c("GO", "GM", "SI"), c("alpha1", "beta1","c","alpha2","beta2"))
	pname       = paste(rep(colnames(modm),each=nz), "[",rep(colnames(Z), nth),"]", sep="")
	idm         = which(rownames(modm)==model)
	idth        = which(modm[rep(idm, nz),]==1)
	nthm        = sum(modm[idm,])
	th.low      = matrix(-Inf, nrow(modm),nth, dimnames=dimnames(modm))
	th.low["SI",c("beta1","beta2")] = 0
	low         = matrix(th.low[idm,],nz, nth, byrow=TRUE)
	dimnames(low) = list(colnames(Z), colnames(modm))

	if(is.null(jumps)){
		thj       = t(t(modm[rep(idm, nz),])*c(0.005, 0.005, 0.02, 0.0075, 0.001))
	} else {
		ljumps    = length(jumps)
		if(ljumps/nthm != round(ljumps/nthm)){
			stop(paste("\nLength of jumps not multiple of", nthm,"\n")) 
		} else if(ljumps > nth*nz) {
			stop(paste("\nlength of jumps larger than", nthm * nz,"\n"))		}
		thj       = modm[rep(idm, nz), ]
		thj[,modm[idm,]==1] = matrix(jumps,nz, sum(modm[idm,]), byrow=TRUE)
	}

	if(is.null(ini.pars)){
		thg       = t(t(modm[rep(idm, nz),])*c(-1, 0.001, 0, -1, 0.001))
	} else {
		lini.pars = length(ini.pars)
		if(lini.pars/nthm != round(lini.pars/nthm)){
			stop(paste("\nLength of starting parameters not multiple of", nthm,"\n")) 
		} else if(lini.pars > nth*nz) {
			stop(paste("\nLength of starting parameters larger than", nthm * nz,"\n"))		}
		thg       = modm[rep(idm, nz), ]
		thg[,modm[idm,]==1] = matrix(ini.pars,nz, sum(modm[idm,]), byrow=TRUE)
	}	


	#The code:
	# Extract times of birth and death:
	bi          = bd[,1]
	di          = bd[,2]

	# Define study duration:
	Dx          = (st[2]-st[1])
	Tm          = matrix(st, n, nt, byrow=TRUE)

	# Calculate first and last time observed:
	ytemp       = t(t(Y) * st)
	li          = c(apply(ytemp,1,max))
	ytemp[ytemp==0] = 10000
	fi          = c(apply(ytemp,1,min))
	fi[fi==10000] = 0
	rm("ytemp")

	# Calculate number of times detected:
	oi          = Y \%*\% rep(1, nt)

	# FUNCTIONS:
	# Survival and mortality:
	m.g         = function(x,th) exp(th[,1]-th[,2]*x) * modm[idm,1]*modm[idm,2] + th[,3] + exp(th[,4] + th[,5]*x)
	S.g         = function(x, th){
		Sg          =  x * (-th[,3]) + exp(th[,4])/th[,5] * (1-exp(th[,5]*x))
		if(idm==3) Sg = Sg + (exp(th[,1])/th[,2] * (exp(-th[,2]*x)-1))
		return(exp(Sg))
  }
	f.g         = function(x,th) m.g(x,th) * S.g(x,th)
	S.x         = function(th) S.g(xv, matrix(th,1,nth))
	m.x         = function(th) m.g(xv, matrix(th,1,nth))

	# Lower bounds for parameter c:
	c.low       = function(th){
		if(idm==1) cl = 0
		if(idm==2) cl = ifelse(th[5] > 0, -exp(th[4]), 0)
		if(idm==3){
			x.minf = (th[1]+log(th[2]) - th[4]-log(th[5]))/(th[2] + th[5])
			cl     = -exp(th[1]-th[2]*(x.minf)) - exp(th[4]+th[5]*(x.minf))
		}
	return(cl)
	}

	# STARTING VALUES:
	# Survival parameters
	Zthg        = Z \%*\% thg

	# Times of birth and death:
	bi0         = which(bi==0)
	bg          = bi
	bg[bi==0 & fi>0] = fi[bi==0 & fi>0] - 1
	bg[bi==0 & fi==0 & di>0] = di[bi==0 & fi==0 & di>0] - 1

	di0         = which(di==0)
	dg          = di
	dg[di==0 & li>0]  = li[di==0 & li>0] + 1
	dg[di==0 & li==0] = bi[di==0 & li==0] + 1
	dg[dg<Ti] = Ti+1

	xg          = dg - bg

	idtrg       = which(bg<Ti)
	
	# Calculate likelihood:
	px          = sum(fx.fun(xg + 0.5*Dx, Zthg, idm=idm, log=TRUE)) - sum(Sx.fun(Ti-bg[idtrg] + 0.5*Dx, Zthg[idtrg,], idm=idm, log=TRUE))
	
	if(is.na(px) | px==-Inf) ModifyPars = TRUE else ModifyPars = FALSE
	
	# Calculate likelihood with range for inipars:
	thnl        = matrix(qtnorm(rep(c(0.995, 0.005, 0.995,0.995, 0.995),each=nz), mean=c(thg), sd=c(thj), lower=c(low)), nz, nth) * modm[rep(idm,nz), ]
	
	Zthnl       = Z \%*\% thnl
	pxl         = sum(fx.fun(xg + 0.5*Dx, Zthnl, idm=idm, log=TRUE)) - sum(Sx.fun(Ti-bg[idtrg] + 0.5*Dx, Zthnl[idtrg,], idm=idm, log=TRUE))

	thnu        = matrix(qtnorm(rep(c(0.005, 0.995, 0.005, 0.005, 0.005),each=nz), mean=c(thg), sd=c(thj), lower=c(low)), nz, nth)* modm[rep(idm,nz), ]
	Zthnu       = Z \%*\% thnu
	pxu         = sum(fx.fun(xg + 0.5*Dx, Zthnu, idm=idm, log=TRUE)) - sum(Sx.fun(Ti-bg[idtrg] + 0.5*Dx, Zthnu[idtrg,], idm=idm, log=TRUE))

	if(is.na(pxl) | pxl==-Inf | is.na(pxu) | pxu==Inf) ModifyJumps = TRUE else ModifyJumps = FALSE

	return(list(Pars = ModifyPars, Jumps = ModifyJumps, values=c(px, pxl, pxu)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
