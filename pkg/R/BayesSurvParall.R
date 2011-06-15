BayesSurvParall <-
function(Data, ststart, stend, nsim=5, parallel=FALSE, ncpus=2, ini.pars.mat=NULL, model="SI", niter=50000, burnin=5001, thinning=50, rptp = ststart, jumps=NULL, priors=NULL, lifetable = FALSE){

	# Packages:
	require(msm)

    #Check that niter, burnin, and thinning are compatible.
    if(burnin>niter) stop("\nObject 'burnin' larger than 'niter'.")
    if(thinning>niter) stop("\nObject 'thinning' larger than 'niter'.")

	# Check that nsim is larger than 1:
	if(nsim<=1) stop("\nObject nsim needs to be larger than 1.")

	# Basic error checking:
	tempcheck   = DataCheck(Data, ststart, stend, autofix = FALSE, silent=TRUE)
	if(tempcheck[[1]] == FALSE){ stop("You have an error in Dataframe 'Data',\nplease use function 'DataCheck'") }

	# Model Matrix and boundary values for parameters:
	nth         = 5
	modm        = matrix(c(0,0,1,0,0,1,0,rep(1,8)), 3, nth, dimnames=list(c("GO", "GM", "SI"), c("alpha1", "beta1","c","alpha2","beta2")))
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

	rownames(ini.pars.mat) = paste('set.', 1:nsim, sep="")
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
			sfExport("Data", "ststart", "stend", "model", "niter", "burnin", "thinning", "rptp", "jumps", "ini.pars.mat", "priors", "BayesSurv", "DataCheck", "BayesSurvDIC", "fx.fun", "mx.fun", "Sx.fun", "ObsMatFun", "c.low","parallel","lifetable")
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


	# Collect results:
	simname    = paste("Sim.", 1:nsim, sep="")
	tnth       = ncol(OutBS[[1]]$theta)
	tnpi       = ncol(OutBS[[1]]$pi)
	tnni       = ncol(OutBS[[1]]$bis)
	tnpo       = ncol(OutBS[[1]]$post)
	thing      = seq(burnin, niter, thinning)
	nthin      = length(thing)
	thmat      = array(NA, dim=c(nthin, tnth, nsim), dimnames=list(NULL, colnames(OutBS[[1]]$theta), simname)) 
	pimat      = array(NA, dim=c(nthin, tnpi, nsim), dimnames=list(NULL, colnames(OutBS[[1]]$pi), simname)) 
	bimat      = array(NA, dim=c(nthin, tnni, nsim), dimnames=list(NULL, NULL, simname))
	ximat      = array(NA, dim=c(nthin, tnni, nsim), dimnames=list(NULL, NULL, simname))
	pomat      = array(NA, dim=c(nthin, tnpo, nsim), dimnames=list(NULL, colnames(OutBS[[1]]$post), simname))
	DImat      = matrix(NA, nsim, 5, dimnames = list(simname, colnames(OutBS[[1]]$ModSel))) 
	for(i in 1:nsim){
		thmat[,,i]  = OutBS[[i]]$theta[thing,]
		pimat[,,i]  = OutBS[[i]]$pi[thing,]
		bimat[,,i]  = OutBS[[i]]$bis
		ximat[,,i]  = OutBS[[i]]$dis - OutBS[[i]]$bis
		pomat[,,i]  = OutBS[[i]]$post[thing,]
		DImat[i,]   = OutBS[[i]]$ModSel
	} 

	# RESULTS SUMMARY:
	# Mean, standard error and 95% credible 
	# intervals for survival parameters:
	thq       = cbind(apply(thmat,2,mean), 
	            apply(thmat,2, function(x) sd(c(x))), 
	            t(apply(thmat, 2, 
	            quantile, c(0.025, 0.975))))
	colnames(thq) = c("Mean", "se", "2.5%", "97.5%")

	# Mean, standard error and 95% credible 
	# intervals for recapture parameters:
	if(ncol(pimat)>1){
		piq       = cbind(apply(pimat,2,mean), 
		            apply(pimat,2, function(x) sd(c(x))), 
		            t(apply(pimat, 2, 
		            quantile, c(0.025, 0.975))))
		colnames(piq) = colnames(thq)
	} else {
		piq       = c(mean(pimat), 
		            sd(pimat), 
		            quantile(pimat, c(0.025, 0.975)))
		names(piq) = colnames(thq)	
	}

	# Median and upper and lower 95% credible 
	#intervals for latent states (i.e. ages at death):
	xq        = apply(ximat, 2, 
	            quantile, c(0.5, 0.025, 0.975))
	bq        = apply(bimat, 2, 
	            quantile, c(0.5, 0.025, 0.975))

	# Summary Survival and mortality functions:
	S.x         = function(th) Sx.fun(xv, matrix(th,1,nth), idm=idm)
	m.x         = function(th) mx.fun(xv, matrix(th,1,nth), idm=idm)

	# Median and 95% predictive intervals for survival and mortality:
	nz        = ncol(OutBS[[1]]$Z)
	pmat      = matrix(0, 0, nth * nz); colnames(pmat) = colnames(thmat)
	for(i in 1:nsim) pmat = rbind(pmat, thmat[,,i])
#	pmat      = matrix(0, length(thing), nth * nz); colnames(pmat) = pname
#	pmat[,colnames(thgibbs)] = thmat
	xv        = seq(0, max(xq), 0.1)
	Sxq       = lapply(colnames(OutBS[[1]]$Z), function(zz) apply(apply(pmat[,paste(colnames(modm), 
	            "[",zz,"]", sep="")],1,S.x),1, quantile, c(0.5,0.025,0.975)))
	mxq       = lapply(colnames(OutBS[[1]]$Z), function(zz) apply(apply(pmat[,paste(colnames(modm), 
	            "[",zz,"]", sep="")],1,m.x),1, quantile, c(0.5,0.025,0.975)))

	# Convergence diagnostics (potential scale reduction):
	Means      = t(apply(thmat, c(2,3), mean))
	Vars       = t(apply(thmat, c(2,3), var))
	meanall    = apply(Means,2,mean)
	B          = nthin/(nsim-1)*apply(t((t(Means)-meanall)^2),2,sum)
	W          = 1/nsim*apply(Vars,2,sum)
	Varpl      = (nthin-1)/nthin * W + 1/nthin*B
	Rhat       = sqrt(Varpl/W)
	conv       = cbind(B,W,Varpl,Rhat)
	rownames(conv) = colnames(OutBS[[1]]$theta)


	#Return a list object
	output = list(bd = OutBS[[1]]$bd,Y = OutBS[[1]]$Y,Z = OutBS[[1]]$Z, post=pomat, ng=niter, bng = burnin, thint = thing, theta=thmat, pi = pimat, bis = bimat, xis = ximat, thsum = thq, pisum=piq, xqsum=xq, Sxsum = Sxq, mxsum = mxq, modm=modm, idm=idm, jumps=jumps, ini.pars=ini.pars.mat, priors=priors, convergence=conv)

	return(output)
}

