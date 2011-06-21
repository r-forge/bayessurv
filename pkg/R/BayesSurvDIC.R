BayesSurvDIC    = function(outBS){
	
# Calculate DIC
thin	    = seq(outBS$bng,outBS$ng,outBS$thint)
modepost   = outBS$post[thin,3]
L          = length(thin)
Dm         = -2*modepost
Dmode      = -2*modepost[which(modepost==max(modepost))[1]]
Dave       = mean(Dm)
pD         = Dave - Dmode
k          = ncol(outBS$pi)+sum(outBS$modm[outBS$idm,]*ncol(outBS$Z))
DIC        = 2*Dave - Dmode
ModSel     = c(Dave, Dmode, pD, k, DIC)
names(ModSel) = c("D.ave", "D.mode", "pD", "k", "DIC")

	
#Return ModSel
return(ModSel)
	
}