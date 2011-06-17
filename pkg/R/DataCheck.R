DataCheck <-
function(Data, ststart, stend, autofix = rep(0,7), silent=TRUE) {

#ToDo: add code to specify variations on what to do with the autofix
#Autofix = a vector that defines whether, and how, to fix a problematic dataset.
#It is a vector with 7 elements (one for each type of error), and the numbers 
#specify how to deal with the error.
#The default is c(0,0,0,0,0,0,0) - to do nothing.
#Type 1: 0 = do nothing; 1 = remove from dataframe
#Type 2: 0 = do nothing; 1 = remove from dataframe
#Type 3: 0 = do nothing; 1 = replace death records with 0; 2 = replace birth records with 0; 3 = replace both birth and death records with 0
#Type 4: 0 = do nothing; 1 = remove spurious post-death observations
#Type 5: 0 = do nothing; 1 = remove observations that pre-date year of birth 
#Type 6: 0 = do nothing; 1 = replace birth year element of observation matrix with 0
#Type 7: 0 = do nothing; 1 = replace death year element of observation matrix with 0

#if(autofix == TRUE) stop("Autofix specification should be a numerical vector of length 7.")

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

  
# 1. Death before observations start
    type1 = which(bd[,2] < Ti & bd[,2]!=0)
    if (length(type1) != 0) {
        cat("The following rows deaths occur before observations start:\n")
        print(type1)
        
        #Actions - remove those rows from bd, Y and Z
        if (autofix[1] == 1) {
            bd = bd[-type1, ]
            Y = Y[-type1, ]
            idnames=idnames[-type1]
            Z = Z[-type1, ]
            n = nrow(Y)
        cat("These records have been removed from the Dataframe\n")
        }
    }
    
# 2. No birth/death AND no obervations
    type2 = which(rowSums(bd) + rowSums(Y) == 0)
    if (length(type2) != 0) {
        cat("The following rows have no Data (unknown birth, unknown death, and no observations):\n")
        print(type2)
        
        #Actions - remove those rows from bd, Y and Z
        if (autofix[2] == 1) {
            bd = bd[-type2, ]
            Y = Y[-type2, ]
            idnames=idnames[-type2]
            Z = Z[-type2, ]
            n = nrow(Y)
        cat("These records have been removed from the Dataframe\n")
        }
    }
    
# 3. Birth after death 
#Need to check that this works - looks like it has changed.
    bd2 = bd
    bd2 = cbind(bd2, 1:n)
#    bd2 = subset(bd2, bd2[, 1] != 0 & bd2[, 2] != 0)
#    type3 = bd2[,3][which(bd2[, 1] > bd2[, 2])]
    type3 = which(bd[,1] > bd[,2] & bd[, 1] != 0 & bd[, 2] != 0)    
    if (length(type3) != 0) {
        cat("The following rows have birth dates that are later than their death dates:\n")
        print(type3)
        
        #Actions - remove the death, birth, both records?
        if (autofix[3] == 1) {bd[type3,2] = 0; cat("The death records have been replaced with 0.\n\n")}
        else if (autofix[3] == 2) {bd[type3,1] = 0; cat("The birth records have been replaced with 0\n")}
        else if (autofix[3] == 3) {bd[type3,1:2] = 0; cat("The birth and death records have been replaced with 0\n")}
        }
    
# 4. Observations after death
    # Calculate first and last time observed: 
    st = Ti:Tf
    ytemp = t(t(Y) * st)
    lastObs = c(apply(ytemp, 1, max))
    tempDeath = bd[,2]; tempDeath[tempDeath==0] = Inf
    type4 = which(lastObs>tempDeath & tempDeath>=Ti)
    rm(tempDeath)
    
        if (length(type4) != 0) {
        cat("The following rows have observations that occur after the year of death:\n")
        print(type4)
        
        #Actions - remove spurious post-death observations
        if (autofix[4] == 1) {
        	 Ymd  = ((Tm - bd[,2]) * Y)[type4,]
        	 Ymd[Ymd>0]  = 0
        	 Ymd[Ymd<0]   = 1
        	 Y[type4,]  = Ymd
        	 
#            dyr = bd[type4, 2]
#            dyr = dyr - Ti + 1
            
#            for (i in 1:length(type4)) {
#                Y[type4[i], (dyr[i]):ncol(Y)] = rep(0, length((dyr[i]):ncol(Y)))
#            }
       cat("Observations that post-date year of death have been removed.\n\n")
       }
    }

# 5. Observations before birth
    ytemp[ytemp == 0] = Inf
    firstObs = c(apply(ytemp, 1, min))
    type5 = which(firstObs < bd[, 1])
    
    if (length(type5) != 0) {
        cat("The following rows have observations that occur before the year of birth:\n")
        print(type5)
        
        #Actions - remove spurious pre-birth observations
        if (autofix[5] == 1) {
        	 Ymd  = ((Tm - bd[,1]) * Y)[type5,]
        	 Ymd[Ymd>0]  = 1
        	 Ymd[Ymd<0]   = 0
        	 Y[type5,]  = Ymd
            
#            byr = bd[type5, 1]
#            byr = byr - Ti
            
#            for (i in 1:length(type5)) {
#                Y[type5[i], (1:byr[i])] = rep(0, length(1:byr[i]))
#            }
              cat("Observations that pre-date year of birth have been removed.\n\n")

       }
    }
    
# 6. Year of birth should be a zero in recapture matrix Y
    idb   = which(bd[,1]>0 & bd[,1]>=Ti & bd[,1] <=Tf)
    bcol  = bd[idb,1] - Ti + 1
    bpos  = (bcol-1)*n + idb
    type6 = which(Y[bpos]==1)

    if (length(type6) != 0) {
        cat("The following rows have a one in the recapture matrix in the birth year:\n")
        print(type6)
        
        #Actions - put a zero.
        if (autofix[6] == 1) Y[bpos] = 0
   }
    
# 7. Year of death should be a zero in recapture matrix Y
    idd   = which(bd[,2]>0 & bd[,2]>=Ti)
    dcol  = bd[idd,2] - Ti + 1
    dpos  = (dcol-1)*n + idd
    type7 = which(Y[dpos]==1)
    if (length(type7) != 0) {
        cat("The following rows have a one in the recapture matrix in the death year:\n")
        print(type7)
        
        #Actions - put a zero.
        if (autofix[7] == 1) Y[dpos] = 0
    }   

    
#All OK
    if (length(c(type1, type2, type3, type4, type5, type6, type7)) == 
        0) {
       cat("No problems were detected with the Data.\n\n")
    }
    ok = length(c(type1, type2, type3, type4, type5, type6, type7)) == 
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
        return(list(ok=ok, newData=data.frame(idnames,bd,Y,Z), type1=type1, type2=type2, type3=type3, type4=type4, type5=type5, 
            type6=type6, type7=type7))
    } else {
        return(list(ok=ok, newData=data.frame(idnames,bd,Y), type1=type1, type2=type2, type3=type3, type4=type4, type5=type5, 
            type6=type6, type7=type7))
    }
}

