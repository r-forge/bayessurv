# This function coverts census data into a capture history.
# Input is 2 vectors: ID of the animal, and the Date it has been observed.
# Date can be entered as an R Posix (?) Date object, or alternatively you can specify 
# If d is already a date format dformat is ignored
# the format of the date can be anything, as long as it is consistent and specified e.g. datetype = "ddmmyyyy", "dd.mm.yyyy", "mm.dd.yyyy" etc. 
# dd = day, mm= month, yyyy =year. 


CensusToCaptHist=function(ID,d,dformat="dd/mm/yyyy"){

# Check data
if(class(ID) != "character") ID = as.character(ID)


if(class(d) == "Date") {yr = as.numeric(format(d, format = "%Y")) }



if(class(d) != "Date"){
    # Deduce the d format, if it is not already a d object
    st = strsplit(dformat, split = NULL)
    y  = which(st[[1]]=="y")

    # Extract the year
    yr = as.numeric(substr(d,y[1],y[length(y)]))
    }


# Add a dummy individual seen in all years, 
# to cope with years where there are no recaptures
dyr = min(yr):max(yr)
ID  = c(ID,rep("XdummyX",length(dyr)))
yr  = c(yr,dyr)

mat        = as.matrix(table(ID,yr))
mat[mat>0] = 1

#Remove the dummy row
mat        = mat[-which(rownames(mat)=="XdummyX"),]
return(mat)
}


