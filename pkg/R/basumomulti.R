basumomulti <-
function(Data, ststart, stend, nsim=5, parallel=FALSE, ncpus=2, ini.pars.mat=NULL, model="SI", niter=50000, burnin=5001, thinning=50, rptp = ststart, jumps=NULL, priors=NULL, lifetable = FALSE){

	# Packages:
	require(msm)

	# Data error checking:
	tempcheck   = DataCheck(Data, ststart, stend, silent=TRUE)
	if(tempcheck[[1]] == FALSE) stop("You have an error in Dataframe 'Data',\nplease use function 'DataCheck'\n", call.=FALSE)

    #Check that niter, burnin, and thinning are compatible.
    if(burnin>niter) stop("\nObject 'burnin' larger than 'niter'.")
    if(thinning>niter) stop("\nObject 'thinning' larger than 'niter'.")

	# Check that nsim is larger than 1:
	if(nsim<=1) stop("\nObject 'nsim' needs to be larger than 1\nFor a single run use function 'basumo'.")

	# Model Matrix and boundary values for parameters:
	nth         = 5
	modm        = matrix(c(0,0,1,0,0,1,0,rep(1,8)), 3, nth, dimnames=list(c("GO", "GM", "SI"), c("alpha1", "beta1", "c", "alpha2", "beta2")))
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
					clow     = c.low(ini.pars.mat[i,], idm=idm)
					if(ini.pars.mat[i,3]<clow) ini.pars.mat[i,3] = rtnorm(1, 0, 0.5, lower=clow)
				}
				test.pars    = TestIniPars(Data=Data,ini.pars=ini.pars.mat[i,modm[idm,]==1], ststart=ststart,stend=stend,model=model)
				if(test.pars$P | test.pars$J) test.pj = TRUE else test.pj = FALSE
			}
		}
		ini.pars.mat  = ini.pars.mat[,modm[idm,]==1]
	} else {
		if(ncol(ini.pars.mat)!=nthm){
			stop(paste("\nWrong dimemnsions in 'ini.pars.mat'.\nFor model", model, "initial parameters matrix should have", nthm, "columns.\n"), call.=FALSE)
		} else {
			test.pars  = unlist(apply(ini.pars.mat,1, function(pars) TestIniPars(Data=Data,ini.pars=pars, ststart=ststart,stend=stend,model=model)$Pars))
			idchp      = which(test.pars)
			if(length(idchp)>0) stop(paste("\nInitial parameters for simulation(s)", paste(idchp, collapse=", "), "produce pdf values equal to 0 for some individuals.\n"), call.=FALSE)
		}
	}

	# Parallel function:
	paralBSM  = function(sim){
		if(parallel) for(ii in 1:(sim*2)){}
		outbsm  = basumo(Data=Data, ststart=ststart, stend=stend, model=model, niter=niter, burnin=burnin, thinning=thinning, rptp=rptp, jumps=jumps, ini.pars=ini.pars.mat[sim,], priors=priors, lifetable=lifetable, datacheck=FALSE)
		return(outbsm)
	}

	# Run multiple BayesSurv simulations:
	Start      = Sys.time()
	cat("\nMultiple simulations started.\n")
	if(parallel){
		availpkg     = available.packages()
		if(!is.element("snowfall", availpkg)){
			warning("\nPackage 'snowfall' is not installed.\nSimulations will not be ran in parallel (computing time will be longer...)\n")
			outBSM   = lapply(1:nsim, paralBSM)
			names(outBSM) = paste("Sim.", 1:nsim)
		} else {
			require(snowfall)
			sfInit(parallel=TRUE, cpus=ncpus);
			sfExport("Data", "ststart", "stend", "model", "niter", "burnin", "thinning", "rptp", "jumps", "ini.pars.mat", "priors", "basumo", "DataCheck", "fx.fun", "mx.fun", "Sx.fun", "ObsMatFun", "c.low","parallel","lifetable","TestIniPars")
			sfLibrary(msm)
			outBSM = sfClusterApplyLB(1:nsim, paralBSM)
			sfStop()
			names(outBSM) = paste("Sim.", 1:nsim)
		}
	} else {
		outBSM   = lapply(1:nsim, paralBSM)
		names(outBSM) = paste("Sim.", 1:nsim)
	}

	# Report running time:
	End        = Sys.time()
	cat(paste("\nMultiple MCMC computing time: ", round(as.numeric(julian(End)-julian(Start))*24*60, 2), " minutes\n", sep=""))

	# Report if simo simulations failed:
	full.runs  = rep(0,nsim)
	last.steps = full.runs
	for(i in 1:nsim){
		full.runs[i]  = ifelse(outBSM[[i]]$full.run, 1, 0)
		last.steps[i] = outBSM[[i]]$last.step
	} 
	id.failed  = which(full.runs==0)
	id.ran     = which(full.runs==1)
	all.ran    = FALSE
	if(length(id.failed)>0 & length(id.failed)<nsim){
		cat("\nOne or more simulations failed\nConvergence diagnostics could not be calculates\n")
	} else if(length(id.failed)==nsim){
		cat("\nAll simulations failed\nConvergence diagnostics could not be calculates\n")
	} else {
		all.ran    = TRUE
		cat("\nMultiple simulations finished.\n")
	}
	

	# Collect results:
	nz         = ncol(outBSM[[1]]$Z)
	simname    = paste("Sim.", (1:nsim), sep="")
	tnth       = sum(modm[idm,])*nz
	tnpi       = length(rptp)
	tnni       = nrow(Data)
	tnpo       = 3
	thing      = seq(burnin, niter, thinning)
	nthin      = length(thing)
	
	thmat      = array(NA, dim=c(nthin, tnth, nsim), dimnames=list(NULL, colnames(outBSM[[1]]$theta), simname)) 
	pimat      = array(NA, dim=c(nthin, tnpi, nsim), dimnames=list(NULL, colnames(outBSM[[1]]$pi), simname)) 
	bimat      = array(NA, dim=c(nthin, tnni, nsim), dimnames=list(NULL, NULL, simname))
	ximat      = array(NA, dim=c(nthin, tnni, nsim), dimnames=list(NULL, NULL, simname))
	pomat      = array(NA, dim=c(nthin, tnpo, nsim), dimnames=list(NULL, colnames(outBSM[[1]]$post), simname))
	DImat      = matrix(NA, nsim, 5, dimnames = list(simname, colnames(outBSM[[1]]$ModSel))) 
	if(lifetable) ltmat = list()
	
	for(i in 1:nsim){
		thmat[,,i]  = outBSM[[i]]$theta[thing,]
		pimat[,,i]  = outBSM[[i]]$pi[thing,]
		bimat[,,i]  = outBSM[[i]]$bis
		ximat[,,i]  = outBSM[[i]]$xis
		pomat[,,i]  = outBSM[[i]]$post[thing,]
		if(full.runs[i]==0){
			DImat[i,]   = rep(NA, 5)
			if(lifetable) ltmat[[i]] = NA
		} else {
			DImat[i,]   = outBSM[[i]]$ModSel
			if(lifetable) ltmat[[i]] = outBSM[[i]]$LT
		}
	} 

	# ADD DIAGNOSTIC OF TOTAL RUN PER SIMULATION
	if(all.ran){
		# RESULTS SUMMARY:
		# Mean, standard error and 95% credible 
		# intervals for survival parameters:
		thq       = cbind(apply(thmat[,,id.ran],2,mean), 
		            apply(thmat[,,id.ran],2, function(x) sd(c(x))), 
		            t(apply(thmat[,,id.ran], 2, 
		            quantile, c(0.025, 0.975))))
		colnames(thq) = c("Mean", "se", "2.5%", "97.5%")

		# Mean, standard error and 95% credible 
		# intervals for recapture parameters:
		if(ncol(pimat)>1){
			piq       = cbind(apply(pimat[,,id.ran],2,mean), 
			            apply(pimat[,,id.ran],2, function(x) sd(c(x))), 
			            t(apply(pimat[,,id.ran], 2, 
			            quantile, c(0.025, 0.975))))
			colnames(piq) = colnames(thq)
		} else {
			piq       = c(mean(pimat[,,id.ran]), 
			            sd(pimat[,,id.ran]), 
			            quantile(pimat[,,id.ran], c(0.025, 0.975)))
			names(piq) = colnames(thq)	
		}

		# Median and upper and lower 95% credible 
		#intervals for latent states (i.e. ages at death):
		xq        = apply(ximat[,,id.ran], 2, 
		            quantile, c(0.5, 0.025, 0.975))
		bq        = apply(bimat[,,id.ran], 2, 
		            quantile, c(0.5, 0.025, 0.975))

		# Summary Survival and mortality functions:
		S.x         = function(th) Sx.fun(xv, matrix(th,1,nth), idm=idm)
		m.x         = function(th) mx.fun(xv, matrix(th,1,nth), idm=idm)

		# Median and 95% predictive intervals for survival and mortality:
		pname       = paste(rep(colnames(modm),each=nz), "[",rep(colnames(outBSM[[1]]$Z), nth),"]", sep="")
		pmat        = matrix(0, 0, nth * nz); colnames(pmat) = pname
		for(i in id.ran){
			npmat   = matrix(0, nthin, nth * nz)
			colnames(npmat) = pname
			npmat[,colnames(thmat)] = thmat[,,i]
			pmat    = rbind(pmat, npmat)
		} 
		xv          = seq(0, max(xq), 0.1)
		Sxq         = lapply(colnames(outBSM[[1]]$Z), function(zz) apply(apply(pmat[,paste(colnames(modm), 
		            "[",zz,"]", sep="")],1,S.x),1, quantile, c(0.5,0.025,0.975)))
		mxq         = lapply(colnames(outBSM[[1]]$Z), function(zz) apply(apply(pmat[,paste(colnames(modm), 
		            "[",zz,"]", sep="")],1,m.x),1, quantile, c(0.5,0.025,0.975)))

		# Convergence diagnostics (potential scale reduction):
		Means       = t(apply(thmat, c(2,3), mean))
		Vars        = t(apply(thmat, c(2,3), var))
		meanall     = apply(Means,2,mean)
		B           = nthin/(nsim-1)*apply(t((t(Means)-meanall)^2),2,sum)
		W           = 1/nsim*apply(Vars,2,sum)
		Varpl       = (nthin-1)/nthin * W + 1/nthin*B
		Rhat        = sqrt(Varpl/W)
		conv        = cbind(B,W,Varpl,Rhat)
		rownames(conv) = colnames(outBSM[[1]]$theta)

		# assess convergence:
		idnconv     = which(conv[,'Rhat']< 0.95 | conv[,'Rhat']>1.2)
		if(length(idnconv)>0){
			warning("Convergence not reached for some survival parameters", call.=FALSE)
		} else {
			cat("Survival parameters converged appropriately.")
		} 
		if(lifetable) LT = ltmat
	} else {
		thq       = NULL
		piq       = NULL
		xq        = NULL
		bq        = NULL
		Sxq       = NULL
		mxq       = NULL
		if(lifetable) LT = NULL
		ModSel    = NULL
		conv      = NULL
	}

	#Return a list object
	output = list(bd = outBSM[[1]]$bd,Y = outBSM[[1]]$Y,Z = outBSM[[1]]$Z, last.step = last.steps, niter=niter, burnin = burnin, thinning = thinning, theta=thmat, pi = pimat, bis = bimat, xis = ximat, post=pomat, thsum = thq, pisum=piq, xqsum=xq, Sxsum = Sxq, mxsum = mxq, modm=modm, idm=idm, jumps=jumps, ini.pars=ini.pars.mat, priors=priors, full.run=full.runs, ModSel = ModSel)
	if(lifetable) output$LT = ltmat
	output$convergence = conv
	return(output)
}

