makecovm <-
function(Data,covs) {

#covs =c("Sex") ; Data = newdata

nd = data.frame(Data[,which(names(Data)%in%covs)])
if(ncol(nd)==1) names(nd)=covs

Cl = as.vector(sapply(nd,class))

temp = c(Cl == "factor"|Cl == "numeric")
if(sum(temp)<length(Cl)) stop("Data types for the covariates should be 'factor' or 'numeric'")

factorCovs  = which(Cl=="factor")
numericCovs = which(Cl=="numeric")

if(length(factorCovs) >  1)  FactCovs = as.factor(apply(nd[,factorCovs],1,paste,collapse="."))
if(length(factorCovs) == 1) FactCovs = nd[,factorCovs]
if(length(factorCovs) == 0) FactCovs = NULL

if(length(numericCovs) != 0) NumCovs = nd[,which(Cl=="numeric")]
if(length(numericCovs) == 0) NumCovs = NULL

if(is.null(FactCovs)) {tempCovs = data.frame(NumCovs)} 
if(is.null(NumCovs)) {tempCovs = data.frame(FactCovs)}
if(!is.null(NumCovs)&!is.null(FactCovs)) {tempCovs = data.frame(FactCovs,NumCovs)}

Z=NULL

for (x in 1:ncol(tempCovs)){


if(class(tempCovs[,x]) == "numeric"){
    Z = cbind(Z,tempCovs[,x])
    }

if(class(tempCovs[,x]) == "factor"){

    vL = length(tempCovs[,x])
    nLevs = length(levels(tempCovs[,x]))
    a = (vL * nLevs) - nLevs
    sv = seq(0, a, nLevs) + as.numeric(tempCovs[,x])
    MV = rep(0, vL * nLevs)
    MV[sv] = 1
    M = matrix(MV, ncol = nLevs, nrow = vL, byrow = TRUE)
    colnames(M) = levels(tempCovs[,x])
    
    Z = cbind(Z,M)

}


}

#Assign column names to the numeric variables (if there are any)
if(length(numericCovs)!=0){
if(is.null(colnames(Z))) {colnames(Z) = rep("",length(numericCovs))}
colnames(Z)[colnames(Z)==""] = names(nd)[Cl=="numeric"]
}

#head(Z)

return(Z)

}




