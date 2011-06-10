MakeCovariateMatrix <-
function(x) {
    if (class(x) != "factor") 
        stop("x is not a factor")
    
    vL = length(x)
    nLevs = length(levels(x))
    a = (vL * nLevs) - nLevs
    sv = seq(0, a, nLevs) + as.numeric(x)
    MV = rep(0, vL * nLevs)
    MV[sv] = 1
    M = matrix(MV, ncol = nLevs, nrow = vL, byrow = TRUE)
    colnames(M) = levels(x)
    
    return(M)
}

