DataCheck <-
function(Data, ststart, stend, autofix = FALSE, silent=TRUE) {

#ToDo: add code to specify variations on what to do with the autofix

    Ti         = ststart
	Tf         = stend
	st         = Ti:Tf
	nt         = length(st)
	idnames    = Data[,1]
	n          = nrow(Data)
	bd         = as.matrix(Data[,2:3])
	Y          = as.matrix(Data[,1:nt+3]); colnames(Y) = st
	Tm         = matrix(st, n, nt, byrow=TRUE)
	
	if(ncol(Data)>nt+3){
		Z  = as.matrix(Data[,(nt+4):ncol(Data)])
	} else {
		Z  = matrix(1, n, 1)
	}

  
# 0. Death before observations start
    type0 = which(bd[,2] < Ti & bd[,2]!=0)
    if (length(type0) != 0) {
        cat("The following rows deaths occur before observations start:\n")
        print(type0)
        
        #Actions - remove those rows from bd, Y and Z
        if (autofix == TRUE) {
            bd = bd[-type0, ]
            Y = Y[-type0, ]
            idnames=idnames[-type0]
            Z = Z[-type0, ]
            n = nrow(Y)
        cat("These records have been removed from the Dataframe\n")
        }
    }
    
# 1. No birth/death AND no obervations
    type1 = which(rowSums(bd) + rowSums(Y) == 0)
    if (length(type1) != 0) {
        cat("The following rows have no Data (unknown birth, unknown death, and no observations):\n")
        print(type1)
        
        #Actions - remove those rows from bd, Y and Z
        if (autofix == TRUE) {
            bd = bd[-type1, ]
            Y = Y[-type1, ]
            idnames=idnames[-type1]
            Z = Z[-type1, ]
            n = nrow(Y)
        cat("These records have been removed from the Dataframe\n")
        }
    }
    
# 2. Birth after death 
#Need to check that this works - looks like it has changed.
    bd2 = bd
    bd2 = cbind(bd2, 1:n)
#    bd2 = subset(bd2, bd2[, 1] != 0 & bd2[, 2] != 0)
#    type2 = bd2[,3][which(bd2[, 1] > bd2[, 2])]
    type2 = which(bd[,1] > bd[,2] & bd[, 1] != 0 & bd[, 2] != 0)    
    if (length(type2) != 0) {
        cat("The following rows have birth dates that are later than their death dates:\n")
        print(type2)
        
        #Actions - remove the death record?
        if (autofix == TRUE) {
        bd[type2,2] = 0; cat("The death records have been replaced with 0.\n\n")#remove death record
    #   bd[type2,1] = 0; cat("The birth records have been replaced with 0\n") #remove birth record
    #   bd[type2,1:2] = 0; cat("The birth and death records have been replaced with 0\n") #remove birth and death record
        }
    }
    
# 3. Observations after death
    # Calculate first and last time observed: 
    # (ANOTHER ISSUE, DEATH YEARS HAPPENING BEFORE THE STUDY STARTED...)
    # I think I fixed this above (type0 error)
    st = Ti:Tf
    ytemp = t(t(Y) * st)
    lastObs = c(apply(ytemp, 1, max))
    tempDeath = bd[,2]; tempDeath[tempDeath==0] = Inf
    type3 = which(lastObs>tempDeath & tempDeath>=Ti)
    rm(tempDeath)
    
        if (length(type3) != 0) {
        cat("The following rows have observations that occur after the year of death:\n")
        print(type3)
        
        #Actions - remove spurious post-death observations
        if (autofix == TRUE) {
        	 Ymd  = ((Tm - bd[,2]) * Y)[type3,]
        	 Ymd[Ymd>0]  = 0
        	 Ymd[Ymd<0]   = 1
        	 Y[type3,]  = Ymd
        	 
#            dyr = bd[type3, 2]
#            dyr = dyr - Ti + 1
            
#            for (i in 1:length(type3)) {
#                Y[type3[i], (dyr[i]):ncol(Y)] = rep(0, length((dyr[i]):ncol(Y)))
#            }
       cat("Observations that post-date year of death have been removed.\n\n")
       }
    }

# 4. Observations before birth
    ytemp[ytemp == 0] = Inf
    firstObs = c(apply(ytemp, 1, min))
    type4 = which(firstObs < bd[, 1])
    
    if (length(type4) != 0) {
        cat("The following rows have observations that occur before the year of birth:\n")
        print(type4)
        
        #Actions - remove spurious pre-birth observations
        if (autofix == TRUE) {
        	 Ymd  = ((Tm - bd[,1]) * Y)[type4,]
        	 Ymd[Ymd>0]  = 1
        	 Ymd[Ymd<0]   = 0
        	 Y[type4,]  = Ymd
            
#            byr = bd[type4, 1]
#            byr = byr - Ti
            
#            for (i in 1:length(type4)) {
#                Y[type4[i], (1:byr[i])] = rep(0, length(1:byr[i]))
#            }
              cat("Observations that pre-date year of birth have been removed.\n\n")

       }
    }
    
# 5. Year of birth should be a zero in recapture matrix Y
    idb   = which(bd[,1]>0 & bd[,1]>=Ti & bd[,1] <=Tf)
    bcol  = bd[idb,1] - Ti + 1
    bpos  = (bcol-1)*n + idb
    type5 = which(Y[bpos]==1)

    if (length(type5) != 0) {
        cat("The following rows have a one in the recapture matrix in the birth year:\n")
        print(type5)
        
        #Actions - put a zero.
        if (autofix == TRUE) Y[bpos] = 0
   }
    
# 6. Year of death should be a zero in recapture matrix Y
    idd   = which(bd[,2]>0 & bd[,2]>=Ti)
    dcol  = bd[idd,2] - Ti + 1
    dpos  = (dcol-1)*n + idd
    type6 = which(Y[dpos]==1)
    if (length(type6) != 0) {
        cat("The following rows have a one in the recapture matrix in the death year:\n")
        print(type6)
        
        #Actions - put a zero.
        if (autofix == TRUE) Y[dpos] = 0
    }   

    
#All OK
    if (length(c(type0, type1, type2, type3, type4, type5, type6)) == 
        0) {
       cat("No problems were detected with the Data.\n\n")
    }
    ok = length(c(type0, type1, type2, type3, type4, type5, type6)) == 
        0
    
    
    if(!silent){
    cat(paste("*DataSummary*\nNumber of individuals             =",  nrow(bd), "\n"))
    cat(paste("Earliest observation in recapture matrix         =", min(ytemp), "\n"))
    cat(paste("Latest observation in recapture matrix           =", max(ytemp[ytemp != Inf]), "\n"))
    cat(paste("Total number of observations in recapture matrix =", sum(Y), "\n\n"))
    
    cat(paste("Earliest recorded birth year     =", min(bd[, 1][bd[, 1] != 0]), "\n"))
    cat(paste("Latest recorded birth year       =", max(bd[, 1][bd[, 1] != 0]), "\n"))
    cat(paste("Earliest recorded death year     =", min(bd[, 2][bd[, 2] != 0]), "\n"))
    cat(paste("Latest recorded death year       =", max(bd[, 2][bd[, 2] != 0]), "\n"))
    cat(paste("Number with known birth year           =", sum(bd[, 1] != 0), "\n"))
    cat(paste("Number with known death year           =", sum(bd[, 2] != 0), "\n"))
    cat(paste("Number with known birth AND death year =", sum(bd[, 2] != 0 & bd[, 1] != 0), "\n"))
	}

    if (ncol(Data)>nt+3) {
        return(list(ok=ok, newData=data.frame(idnames,bd,Y,Z), type0=type0, type1=type1, type2=type2, type3=type3, type4=type4, 
            type5=type5))
    } else {
        return(list(ok=ok, newData=data.frame(idnames,bd,Y), type0=type0, type1=type1, type2=type2, type3=type3, type4=type4, 
            type5=type5))
    }
}

