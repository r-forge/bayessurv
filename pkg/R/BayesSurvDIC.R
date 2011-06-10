BayesSurvDIC <-
function(outBS){
	
# Calculate DIC
ModSel     = matrix(0, 5, 1)
dimnames(ModSel) = list(c("D.ave", "D.mode", "pD", "k", "DIC"))

thin	   = seq(outBS$bng,outBS$ng,outBS$thint)
L          = length(seq(outBS$bng,outBS$ng,outBS$thint))
Dm         = rep(0,L)

Dm   	   = -2*outBS$post[thin]
DMode	   = -2*modepost[which(modepost==max(modepost))[1]] #???

Dave       = mean(Dm)
pD         = Dave - DMode
k          = length(outBS$pi)+length(outBS$theta)
DIC        = 2*Dave - DMode
ModSel     = c(Dave, Dmode, pD, k, DIC)
	
#Return ModSel
return(ModSel)
	
}

