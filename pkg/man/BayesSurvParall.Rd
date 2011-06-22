\name{BayesSurvParall}
\alias{BayesSurvParall}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to conduct a BayesSurv analyis with the option of running in parallel.}

\description{
A function to perform Bayesian estimation of latent stats (birth and death year) using an undelying mortality model.
Further details required!!
}

\usage{
function(Data, ststart, stend, nsim=5, parallel=FALSE, ncpus=2, ini.pars.mat=NULL, model="SI", niter=50000, burnin=5001, thinning=50, rptp = ststart, jumps=NULL, priors=NULL, lifetable = FALSE)
}


\arguments{
  \item{Data}{
A \code{data.frame} to be used as an input data file for BayesSurv. The first column is the individual's ID, the second and third columns are birth and death years respectively. Columns 4-(n-1) represent the observation window of n years. This is followed (optionally) by columns for covariates, which must be of 'factor' or 'numeric' class.
}
  \item{ststart}{
The first year of the observation matrix.
}
  \item{stend}{
The last year of the observation matrix.
}
  \item{nsim}{
Definition required.
}
  \item{parallel}{
A logical argument indicating whether to run the analyses in parallel (TRUE) or serial (FALSE).
}
  \item{ncpus}{
The number of CPUs to use in the analysis
}
  \item{ini.pars.mat}{
A matrix of initial parameters. Number of columns should be appropriate to the number of terms in the model (2 for Gompertz, 3 for Gompertz-Makeham, 5 for Siler).
}
  \item{model}{
The underlying mortality model used. "GO"= Gompertz, "GM" = Gompertz-Makeham and "SI" = Siler.
}
  \item{niter}{
The number of iterations for the MCMC process.
}
  \item{burnin}{
Number of iterations for the burn in.
}
  \item{thinning}{
The thinning parameter.
}
  \item{rptp}{
RPTP is a vector defining the recapture probability transition points. These are points (years) where the recapture probability is thought to change.  The default setting is for the recapture probability to be constant throughout the study, so the RPTP is simply defined as a single element vector of the first year of the observation period (e.g. c(1985)).
If one wanted to set a different recapture probability for years after (say) 1990, the RPTP could be defined as c(1985, 1990) etc..
}
  \item{jumps}{
A vector defining the size of jumps for each parameter. The jump specifications depend on the model selected. for Siler are specified in the following order:   \code{c(alpha1, beta1, c, alpha2, beta2)}
The jumps for Gompertz-Makeham are specified as \code{c(c,alpha,beta)}, and the jumps for Gompertz are specified as \code{c(alpha,beta)}. 
Any additional elements in the vector are ignored (e.g. if a Gompertz model is specified, but a vector with 5 elements is provided, only the first two will be used (as alpha and beta jumps respectively)).
}
  \item{priors}{
A vector defining the priors for each parameter. The vector takes the same format as described for \code{jumps}.
}
  \item{lifetable}{
A logical argument specifying whether a cohort-based life table should be produced or not. The life table produced has the same format as produced by \code{CohortLT}.
}
}

\details{
Further details will be added here.
}

\value{
The function returns a "BayesSurv" list object with the following elements:
(Copied from BayesSurv - therefore needs changing)

\item{bd }{Description}
\item{Y }{Description}
\item{Z }{Description}
\item{postm }{Description}
\item{g }{Description}
\item{ng }{Description}
\item{bng }{Description}
\item{thint }{Description}
\item{theta }{Description}
\item{pi }{Description}
\item{bis }{Description}
\item{dis }{Description}
\item{thsum }{Description}
\item{pisum }{Description}
\item{xqsum }{Description}
\item{Sxsum }{Description}
\item{mxsum }{Description}
\item{modm }{Description}
\item{idm }{Description}
\item{jumps }{Description}
\item{ini.pars }{Description}
\item{priors }{Description}
\item{LT }{Description}
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
\code{\link{BayesSurv}}, \code{\link{BayesSurvPlots}}, \code{\link{BayesSurvDIC}}, \code{\link{DataCheck}}, \code{\link{CohortLT}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Data, ststart, stend, nsim=5, parallel=FALSE, ncpus=2, ini.pars.mat=NULL, model="SI", niter=50000, burnin=5001, thinning=50, rptp = ststart, jumps=NULL, priors=NULL, lifetable = FALSE){

	# Packages:
	require(msm)

	# Basic error checking:
	tempcheck   = DataCheck(Data, ststart, stend, autofix = FALSE, silent=TRUE)
	if(tempcheck[[1]] == FALSE){ stop("You have an error in Dataframe 'Data',\nplease use function 'DataCheck'") }

	# Model Matrix and boundary values for parameters:
	nth         = 5
	modm        = matrix(1, 3, nth)
	modm[1,1:3] = 0
	modm[2,1:2] = 0
	dimnames(modm) = list(c("GO", "GM", "SI"), c("alpha1", "beta1","c","alpha2","beta2"))
	idm         = which(rownames(modm)==model)
	nthm        = sum(modm[idm,])
	th.low      = matrix(-Inf, nrow(modm),nth, dimnames=dimnames(modm))
	th.low["SI",c("beta1","beta2")] = 0

	
	# Construct initial parameters matrix:
	if(is.null(ini.pars.mat)){
		ini.pars.mat    = matrix(0, nsim, nth)
		for(i in 1:nsim){
			test.pj       = TRUE
			while(test.pj){
				ini.pars.mat[i,] = rtnorm(nth, c(-1, 0.001, 0, -1, 0.001), 0.25, lower=th.low[idm,])
				if(idm>1){
					clow     = c.low(ini.pars.mat[i,])
					if(ini.pars.mat[i,3]<clow) ini.pars.mat[i,3] = rtnorm(1, 0, 0.5, lower=clow)
				}
				test.pars    = TestIniPars(Data=Data,ini.pars=ini.pars.mat[i,modm[idm,]==1], ststart=ststart,stend=stend,model=model)
				if(test.pars$P | test.pars$J) test.pj = TRUE else test.pj = FALSE
			}
		}
		ini.pars.mat  = ini.pars.mat[,modm[idm,]==1]
	} else {
		if(ncol(ini.pars.mat)!=nthm){
			stop(paste("\nWrong dimemnsions in 'ini.pars.mat'.\nFor model", model, "initial parameters matrix should have", nthm, "columns.\n"))
		} else {
			test.pars  = unlist(apply(ini.pars.mat,1, function(pars) TestIniPars(Data=Data,ini.pars=pars, ststart=ststart,stend=stend,model=model)$Pars))
			idchp      = which(test.pars)
			if(length(idchp)>0) stop(paste("\nInitial parameters for simulations", paste(idchp), collapse=", "), "produce pdf values equal to 0 for some individuals.\n")
		}
	}

	sets = c()
	for(i in 1:nsim){sets = c(sets,paste("set",i,collapse=""))}
	rownames(ini.pars.mat) = sets
	print("Set of initial parameters:")
	print(ini.pars.mat)
	
	# Parallel function:
	paralBS   = function(sim){
		if(parallel) for(ii in 1:(sim*2)){}
		outBS   = BayesSurv(Data=Data, ststart=ststart, stend=stend, model=model, niter=niter, burnin=burnin, thinning=thinning, rptp=rptp, jumps=jumps, ini.pars=ini.pars.mat[sim,], priors=priors, lifetable=lifetable, datacheck=FALSE)
		return(outBS)
	}

	# Run multiple BayesSurv simulations:
	if(parallel){
		availpkg     = available.packages()
		if(!is.element("snowfall", availpkg)){
			warning("\nPackage 'snowfall' is not installed.\nSimulations will not be ran in parallel (computing time will be longer...)\n")
			cat("\nMultiple simulations started.\n")
			OutBS   = lapply(1:nsim, paralBS)
			names(OutBS) = paste("Sim.", 1:nsim)
			cat("\nMultiple simulations finished.\n")
		} else {
			require(snowfall)
			sfInit(parallel=TRUE, cpus=ncpus)
			sfExport("Data", "ststart", "stend", "model", "niter", "burnin", "thinning", "rptp", "jumps", "ini.pars.mat", "priors", "BayesSurv", "DataCheck", "BayesSurvDIC")
			sfLibrary(msm)
			cat("\nMultiple simulations started.\n")
			OutBS = sfClusterApplyLB(1:nsim, paralBS)
			sfStop()
			names(OutBS) = paste("Sim.", 1:nsim)
			cat("\nMultiple simulations finished.\n")
		}
	} else {
		cat("\nMultiple simulations started.\n")
		OutBS   = lapply(1:nsim, paralBS)
		names(OutBS) = paste("Sim.", 1:nsim)
		cat("\nMultiple simulations finished.\n")
	}

	# Convergence diagnostics (potential scale reduction):
	tnth       = ncol(OutBS[[1]]$theta)
	thint      = seq(burnin, niter, thinning)
	nthin      = length(thint)
	parr       = array(NA, dim=c(nthin, tnth, nsim))
	for(i in 1:nsim) parr[,,i]  = OutBS[[i]]$theta[thint,]
	Means      = t(apply(parr, c(2,3), mean))
	Vars       = t(apply(parr, c(2,3), var))
	meanall    = apply(Means,2,mean)
	B          = nthin/(nsim-1)*apply(t((t(Means)-meanall)^2),2,sum)
	W          = 1/nsim*apply(Vars,2,sum)
	Varpl      = (nthin-1)/nthin * W + 1/nthin*B
	Rhat       = sqrt(Varpl/W)
	conv       = cbind(B,W,Varpl,Rhat)
	rownames(conv) = colnames(OutBS[[1]]$theta)

	# Final output:
	OutBS$Ini.Pars = ini.pars.mat
	OutBS$Converge = conv
	return(OutBS)
  }
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
