\name{DataCheck}
\alias{DataCheck}
\title{
A function to check the input data file for a BayesSurv analysis.
}
\description{
This function performs some basic error checking on the input files for a BayesSurv analysis. A number of checks are performed including (1) individuals dying before the observation window starts; (2) birth date recorded as being after death date; (3) no observations of any kind (phantom individuals); (4) observations before birth; (5) observations after death; (6) years of birth and deaths must appear as 0 in the observation matrix.
}
\usage{
DataCheck(Data, ststart, stend, autofix = FALSE, silent = TRUE)
}
\arguments{
  \item{Data}{
A data.frame to be used as an input data file for BayesSurv. The first column is the individual's ID, the second and third columns are birth and death years respectively. Columns 4-(n-1) represent the observation window of n years. This is followed (optionally) by a column for a covariate.
}
  \item{ststart}{
The start year of the observation window.
}
  \item{stend}{
The end year of the observation window.}
  \item{autofix}{A numerical vector with a length of 7 that defines whether, and how, to fix a problematic dataset.
Each of the 7 elements of the vector refers to one of the 7 possible error types 1-7. #specify how to deal with the error.
The default is c(0,0,0,0,0,0,0) - to do nothing. Users can specify how each error is dealt with as follows.

Type 1 (death before the observation window starts): 0 = do nothing; 1 = remove from dataframe
Type 2 (no birth/death AND no obervations): 0 = do nothing; 1 = remove from dataframe
Type 3 (births recorded after death): 0 = do nothing; 1 = replace death records with 0; 2 = replace birth records with 0; 3 = replace both birth and death records with 0
Type 4 (observations after death): 0 = do nothing; 1 = remove spurious post-death observations
Type 5 (observations before birth): 0 = do nothing; 1 = remove observations that pre-date year of birth 
Type 6 (year of birth is not a zero in the recapture matrix): 0 = do nothing; 1 = replace birth year element of observation matrix with 0
Type 7 (year of death is not a zero in the recapture matrix): 0 = do nothing; 1 = replace death year element of observation matrix with 0

This should be used with extreme caution. We recommend going back to the individual-based data and fixing each error "by hand".}
  \item{silent}{A logical argument indicating whether to print a detailed report to the screen or not.}
}
\details{
Further details of the autofix argument:
**To be completed.**

}

\value{
This function returns a list object with the following properties.
  \item{ok }{A logical indicator that indicates if the data are free of errors or not. i.e. TRUE = the data are ok.}
  \item{newData }{A corrected data frame.}
  \item{type1 }{A vector of row numbers in the original data frame where there are deaths before the observation window starts.}
  \item{type2 }{A vector of row numbers in the original data frame where there are no birth/death AND no obervations.}
  \item{type3 }{A vector of row numbers in the original data frame where there are births recorded after death.}
  \item{type4 }{A vector of row numbers in the original data frame where there are observations after death.}
  \item{type5 }{A vector of row numbers in the original data frame where there are observations before birth.}
  \item{type6 }{A vector of row numbers in the original data frame where the year of birth is not a zero in the recapture matrix.}
  \item{type7 }{A vector of row numbers in the original data frame where the year of death is not a zero in the recapture matrix.}
}
%\references{}

\author{
Owen Jones \email{jones@mpg.demogr.de}
Maren Rebke \email{rebke@mpg.demogr.de}
}

\note{We recommend using the \code{autofix} option only with extreme caution.
We recommend going back to the individual-based data and fixing each error "by hand".}


\seealso{\code{\link{BayesSurv}}, \code{\link{BayesSurvParall}}}
%\examples{Need to add a working example}

\keyword{ misc }
