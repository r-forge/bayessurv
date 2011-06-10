\name{BayesSurv}
\alias{BayesSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BayesSurv(Data, ststart, stend, model = "SI", niter = 50000, burnin = 5001, thinning = 50, rptp = ststart, jumps = NULL, ini.pars = NULL, priors = NULL, lifetable = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
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
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{thinning}{
%%     ~~Describe \code{thinning} here~~
}
  \item{rptp}{
%%     ~~Describe \code{rptp} here~~
}
  \item{jumps}{
%%     ~~Describe \code{jumps} here~~
}
  \item{ini.pars}{
%%     ~~Describe \code{ini.pars} here~~
}
  \item{priors}{
%%     ~~Describe \code{priors} here~~
}
  \item{lifetable}{
%%     ~~Describe \code{lifetable} here~~
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
function(Data, ststart, stend, model="SI", niter=50000, burnin=5001, thinning=50, 
  rptp = ststart, jumps=NULL, ini.pars=NULL, priors=NULL, lifetable = TRUE){

	require(msm)
	# Basic error checking:
	tempcheck  = DataCheck(Data, ststart, stend, autofix = FALSE, silent=TRUE)
	if(tempcheck[[1]] == FALSE){ stop("You have an error in Dataframe 'Data',\nplease use function 'DataCheck'") }

    #Check that niter, burnin, and thinning are compatible.

	# Data formating:
	Ti         = ststart
	Tf         = stend
	st         = Ti:Tf
	nt         = length(st)
	idnames    = Data[,1]
	n          = nrow(Data)
	bd         = as.matrix(Data[,2:3])
	Y          = as.matrix(Data[,1:nt+3]); colnames(Y) = st
	if(ncol(Data)>nt+3){
		Z  = as.matrix(Data[,(nt+4):ncol(Data)])
	} else {
		Z  = matrix(1, n, 1)
	}
	nth       = 5
	nz        = ncol(Z)

# Model Matrix and boundary values for parameters:
nth       = 5
modm      = matrix(1, 3, nth)
modm[1,1:3] = 0
modm[2,1:2] = 0
dimnames(modm) = list(c("GO", "GM", "SI"), 
                 c("alpha1", "beta1","c","alpha2","beta2"))
pname     = paste(rep(colnames(modm),each=nz), 
            "[",rep(colnames(Z), nth),"]", sep="")
idm       = which(rownames(modm)==model)
idth      = which(modm[rep(idm, nz),]==1)
nthm      = sum(modm[idm,])
th.low    = matrix(-Inf, nrow(modm),nth, dimnames=dimnames(modm))
th.low["SI",c("beta1","beta2")] = 0
low       = matrix(th.low[idm,],nz, nth, byrow=TRUE)
dimnames(low) = list(colnames(Z), colnames(modm))

	# Model variables:
	diffrec    = rptp 
	ng         = niter
	bng        = burnin
	thint      = thinning

	if(is.null(jumps)){
		thj  = c(0.005, 0.005, 0.02, 0.0075, 0.001)
	} else {
		ljumps   = length(jumps)
		if(ljumps/nthm != round(ljumps/nthm)){
			stop(paste("length of jumps not multiple of", nthm)) 
		} else if(ljumps > nth*nz) {
			stop(paste("length of jumps larger than", nthm * nz))		}
		thj  = jumps
	}

	if(is.null(ini.pars)){
		thg  = c(-1, 0.001, 0, -1, 0.001)
	} else {
		thg  = ini.pars
	}

	if(is.null(priors)){
		thp  = c(-5,0.1,-1,0.001,0.005)
	} else {
		thp  = priors
	}


#The code:
# Extract times of birth and death:
bi        = bd[,1]
di        = bd[,2]

# Define study duration:
Dx        = (st[2]-st[1])
Tm        = matrix(st, n, nt, byrow=TRUE)

# Calculate first and last time observed:
ytemp     = t(t(Y) * st)
li        = c(apply(ytemp,1,max))
ytemp[ytemp==0] = 10000
fi        = c(apply(ytemp,1,min))
fi[fi==10000] = 0
rm("ytemp")

# Calculate number of times detected:
oi        = Y \%*\% rep(1, nt)


# FUNCTIONS:
# Survival and mortality:
m.g     = function(x,th) exp(th[,1]-th[,2]*x) * modm[idm,1]*modm[idm,2] + 
          th[,3] + exp(th[,4] + th[,5]*x)
S.g     = function(x, th){
  Sg      = exp((exp(th[,1])/th[,2] * 
           (exp(-th[,2]*x)-1)) * modm[idm,1]*modm[idm,2] + 
           x * (-th[,3]) +
           exp(th[,4])/th[,5] * (1-exp(th[,5]*x)))
  return(Sg)
  }
f.g       = function(x,th) m.g(x,th) * S.g(x,th)
S.x       = function(th) S.g(xv, matrix(th,1,nth))
m.x       = function(th) m.g(xv, matrix(th,1,nth))


# Lower bounds for parameter c:
c.low     = function(th){
  if(idm==1) cl = 0
  if(idm==2) cl = ifelse(th[5] > 0, -exp(th[4]), 0)
  if(idm==3){
    x.minf = (th[1]+log(th[2]) - th[4]-log(th[5]))/(th[2] + th[5])
    cl     = -exp(th[1]-th[2]*(x.minf)) - exp(th[4]+th[5]*(x.minf))
  }
  return(cl)
  }

# Observation matrices:
ObsMatFun = function(f, l){
  Fm    = Tm - f; Fm[Fm>=0] =  1; Fm[Fm<0] = 0
  Lm    = Tm - l; Lm[Lm<=0] = -1; Lm[Lm>0] = 0
  return(Fm * (-Lm))	
  }

# DEFINE PRIORS:
# Priors for survival parameters:
thp       = matrix(thp*modm[idm,], nz, nth, byrow=TRUE)
dimnames(thp) = dimnames(low)
Zthp      = Z \%*\% thp
thv       = 0.5

# Prior for age distribution:
dxx       = 0.001
xx        = seq(0,100,dxx)
zz        = cbind(1,matrix(0, length(xx), nz-1))
Ex        = sum(xx*f.g(xx,zz \%*\% thp)*dxx)
v.x       = function(x) S.g(x,Zthp)/Ex

# Prior parameters for detection probability:
idpi      = findInterval(st, diffrec); names(idpi) = st
npi       = length(unique(idpi))
rho1      = 0.1
rho2      = 0.1

# DEFINE JUMP VARIANCES:
thj       = matrix(thj*modm[idm,], nz, nth, byrow=TRUE)

# STARTING VALUES:
# Survival parameters
thg       = matrix(thg*modm[idm,], nz, nth, byrow=TRUE) 
dimnames(thg) = dimnames(low)
Zthg       = Z \%*\% thg

# Recapture probability:
pig       = rep(0.5, npi)
Pig       = pig[idpi]

# Times of birth and death:
bi0       = which(bi==0)
bg        = bi
bg[bi==0 & fi>0] = fi[bi==0 & fi>0] - 1
bg[bi==0 & fi==0 & di>0] = di[bi==0 & fi==0 & di>0] - 1

di0       = which(di==0)
dg        = di
dg[di==0 & li>0]  = li[di==0 & li>0] + 1
dg[di==0 & li==0] = bi[di==0 & li==0] + 1
dg[dg<Ti] = Ti+1

xg        = dg - bg

# FULL OBSERVATION MATRIX:
Fg        = c(apply(cbind(Ti, bg+1), 1, max))
Lg        = c(apply(cbind(Tf, dg-1), 1, min))
Og        = ObsMatFun(Fg, Lg)
fii       = fi; fii[bi>0 & bi>=Ti] = bi[bi>0 & bi>=Ti]+1; fii[bi>0 & bi<Ti] = Ti
lii       = li; lii[di>0 & di<=Tf] = di[di>0 & di<=Tf]-1; lii[di>0 & di>Tf] = Tf
lfi       = ObsMatFun(fii, lii)

# OUTPUT TABLES:
thing     = seq(bng, ng, by=thint)
thgibbs   = matrix(NA,ng,length(idth))
colnames(thgibbs) = pname[idth]
pigibbs   = matrix(NA, ng, npi)
bgibbs    = matrix(NA,0,n)
dgibbs    = bgibbs
postm     = matrix(NA, ng, 2)

# RUN GIBBS SAMPLER:
if(.Platform$OS.type=="unix") devtype=quartz else devtype=windows
devtype(width=2, height=0.5)
par(mar=rep(0,4))

#Start     = Sys.time()
for(g in 1:ng){

if(g==1){cat("MCMC is running...\n")}
	
  # 1.- SAMPLING:
  # a) Sample survival parameters:
  thn         = thg * 0
  thn[idth]   = c(rtnorm(length(idth), thg[idth], thj[idth], lower=low[idth]))
  if(idm>1){
    low[,3]     = apply(thn, 1, c.low)
    idcl        = which(thn[,3] < low[,3])
    if(length(idcl)>0){
      for(cc in idcl) thn[cc,3]   = c(rtnorm(1, thg[cc,3], thj[cc,3], lower=low[cc,3]))
    }
  }
  
  Zthn        = Z \%*\% thn
  idtrg       = which(bg<Ti)
	
  p.thg       = sum(log(f.g(xg + 0.5*Dx, Zthg))) - 
                sum(log(S.g(Ti-bg[idtrg] + 0.5*Dx, Zthg[idtrg,]))) + 
                sum(dtnorm(c(thg), c(thp), thv, lower=c(low), log=TRUE)) 

  p.thn       = sum(log(f.g(xg + 0.5*Dx, Zthn))) - 
                sum(log(S.g(Ti-bg[idtrg] + 0.5*Dx, Zthn[idtrg,]))) + 
                sum(dtnorm(c(thn), c(thp), thv, lower=c(low), log=TRUE)) 

  r           = exp(p.thn-p.thg)
  z           = runif(1,0,1)
  
  if(is.na(r) & g==1){ 
  	stop("PDF equal to 0 for some ages\nModify starting survival parameters")
  } else if(is.na(r) & g > 1){
  	stop("Proposed pdf equal to 0 for some ages\nModify jumps for survival parameters")  
  } else {
  	if(r>z){
  	 thg   = thn
  	 Zthg  = Zthn
  	 p.thg = p.thn
  }
  }

  # b) Sample times of birth and death:
  bn          = bg 
  bn[bi0]     = bg[bi0] + sample(-1:1, length(bi0), replace=TRUE) 
  bn[bi0]     = apply(cbind(bn[bi0],fi[bi0]-1),1,min)

  idtrn       = which(bn<Ti)
  dn          = dg 
  dn[di0]     = dg[di0] + sample(-1:1, length(di0), replace=TRUE) 
  dn[di0]     = apply(cbind(dn[di0],bn[di0],li[di0]+1),1,max) 

  xn          = dn - bn
  
  Fn          = c(apply(cbind(Ti, bn+1), 1, max))
  Ln          = c(apply(cbind(Tf, dn-1), 1, min))
  On          = ObsMatFun(Fn, Ln)
    
  p.bdg       = log(f.g(xg + 0.5*Dx, Zthg)) + 
                (Og - lfi) \%*\% log(1-Pig) +
                log(v.x(xg + 0.5*Dx))

  p.bdn       = log(f.g(xn + 0.5*Dx, Zthg)) + 
                (On - lfi) \%*\% log(1-Pig) +
                log(v.x(xn + 0.5*Dx))
	
  r           = exp(p.bdn-p.bdg)
  if(length(which(is.na(r)))>0){ 
  	stop("Proposed pdf equal to 0 for some ages\nReduce jump sizes for survival parameters")  
  } else {
  z           = runif(n, 0, 1)
  bg[r>z]     = bn[r>z]
  dg[r>z]     = dn[r>z]
  xg[r>z]     = xn[r>z]
  p.bdg[r>z]  = p.bdn[r>z]
  Og[r>z,]    = On[r>z,]
  }
  # c) Sample recapture probability:
  rho1g       = rho1 + t(t(Y)\%*\% rep(1,n))
  rho2g       = rho2 + t(t(Og - Y)\%*\% rep(1,n))
  Rho1        = tapply(rho1g, idpi, sum)
  Rho2        = tapply(rho2g, idpi, sum)
  pig         = rbeta(npi, Rho1, Rho2)
  if(1 \%in\% pig){
  	pig[pig==1] = 1-1e-5
  	warning("Some recapture probabilities are equal to 1\nThey have been constraint to be fractionally less than 1 for computational reasons")
  	} 
  Pig         = pig[idpi]
		
  # 2.- STORE RESULTS:
  # Parameters and latent states:
  thgibbs[g,] = thg[idth]
  pigibbs[g,] = pig
  if(g \%in\% thing){
    bgibbs     = rbind(bgibbs,bg)
    dgibbs     = rbind(dgibbs,dg)
  }

  # Conditional posteriors:
  postm[g,]   = c(p.thg, sum(p.bdg))
  
  # Progress plot:
  if(g \%in\% round(seq(1,ng,length=100))){
	par(mar=rep(0,4))
	plot(c(0,ng*1.1), c(0,1), axes=FALSE, col=NA, xlab="", ylab="")
	polygon(c(0,ng,ng,0), c(0.35,0.35,0.65,0.65), col=NA, border='dark red')
	polygon(c(0,g,g,0), c(0.35,0.35,0.65,0.65), col='dark red', border='dark red')
	text(ng/2, 0.85, "MCMC progress", cex=0.9)
	text(g, 0.15, paste(round(g/ng*100), "\%", sep=""), cex=0.8)
  	}

  }

# Results:
# RESULTS SUMMARY:
# Mean, standard error and 95\% credible 
# intervals for survival parameters:
thq       = cbind(apply(thgibbs[thing,],2,mean), 
            apply(thgibbs[thing,],2, sd), 
            t(apply(thgibbs[thing,], 2, 
            quantile, c(0.025, 0.975))))
colnames(thq) = c("Mean", "se", "2.5\%", "97.5\%")

# Mean, standard error and 95\% credible 
# intervals for recapture parameters:
if(npi>1){
piq       = cbind(apply(pigibbs[thing,],2,mean), 
            apply(pigibbs[thing,],2, sd), 
            t(apply(pigibbs[thing,], 2, 
            quantile, c(0.025, 0.975))))
colnames(piq) = colnames(thq)
  } else {
piq       = c(mean(pigibbs[thing]), 
            sd(pigibbs[thing]), 
            quantile(pigibbs[thing], c(0.025, 0.975)))
names(piq) = colnames(thq)	
  }

# Median and upper and lower 95\% credible 
#intervals for latent states (i.e. ages at death):
xq        = apply(dgibbs-bgibbs, 2, 
            quantile, c(0.5, 0.025, 0.975))
bq        = apply(bgibbs, 2, 
            quantile, c(0.5, 0.025, 0.975))

# Median and 95\% predictive intervals for survival and mortality:
pmat      = matrix(0, length(thing), nth * nz); colnames(pmat) = pname
pmat[,colnames(thgibbs)] = thgibbs[thing,]
xv        = seq(0, max(xq), 0.1)
Sxq       = lapply(colnames(Z), function(zz) apply(apply(pmat[,paste(colnames(modm), 
            "[",zz,"]", sep="")],1,S.x),1, quantile, c(0.5,0.025,0.975)))
mxq       = lapply(colnames(Z), function(zz) apply(apply(pmat[,paste(colnames(modm), 
            "[",zz,"]", sep="")],1,m.x),1, quantile, c(0.5,0.025,0.975)))


#Return a list object
output = list(bd = bd,Y = Y,Z = Z, post=postm, g=g, ng=ng, bng = bng, thint = thint, theta=thgibbs, pi = pigibbs, bis = bgibbs, dis = dgibbs, thsum = thq, pisum=piq, xqsum=xq, Sxsum = Sxq, mxsum = mxq, modm=modm, idm=idm)

#Add the lifetable to the output list if required.
#To Do: Calculate ax from the model outputs to use as defaults

if(lifetable) {output$LT = CohortLT(xq[1,bq[1,]>=Ti],ax=0.5,n=1)}

return(output)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
