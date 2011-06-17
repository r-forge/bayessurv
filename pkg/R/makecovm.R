makecovm <-
function(Data,covs) {

#covs =c("Assig","rand") ; Data = newdata

nd = Data[,which(names(Data)%in%covs)]
Cl = as.vector(sapply(nd,class))

Z=NULL

for (x in 1:length(covs)){

if(Cl[x] == "numeric"){Z = cbind(Z,nd[,x])}

if(Cl[x] == "factor"){

    vL = length(nd[,x])
    nLevs = length(levels(nd[,x]))
    a = (vL * nLevs) - nLevs
    sv = seq(0, a, nLevs) + as.numeric(nd[,x])
    MV = rep(0, vL * nLevs)
    MV[sv] = 1
    M = matrix(MV, ncol = nLevs, nrow = vL, byrow = TRUE)
    colnames(M) = levels(nd[,x])
    
    Z = cbind(Z,M)

}
}

return(Z)

}




