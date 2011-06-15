BayesSurvParall <-
function(Data, ststart, stend, nsim=5, parallel=FALSE, ncpus=2, ini.pars.mat=NULL, model="SI", niter=50000, burnin=5001, thinning=50, rptp = ststart, jumps=NULL, priors=NULL, lifetable = FALSE){

	# Packages:
	require(msm)

    #Check that niter, burnin, and thinning are compatible.
    if(burnin>niter) stop("Burnin larger than niter.")
    if(thinning>niter) stop("Thinning larger than niter.")

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
			sfExport("Data", "ststart", "stend", "model", "niter", "burnin", "thinning", "rptp", "jumps", "ini.pars.mat", "priors", "BayesSurv", "DataCheck", "BayesSurvDIC", "fx.fun", "mx.fun", "Sx.fun", "ObsMatFun", "c.low")
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

