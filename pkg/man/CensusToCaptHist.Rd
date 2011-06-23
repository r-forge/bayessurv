\name{CensusToCaptHist}

\alias{CensusToCaptHist}

\title{A function to generate a capture history from simple census data.}

\description{This function converts observation records of individuals into an observation matrix with a column for each year, and a row for each individual.}

\usage{CensusToCaptHist(ID, d, dformat = "dd/mm/yyyy")}

\arguments{
  \item{ID}{
A vector of unique identifier of individual animals.
}
  \item{d}{
A vector of dates, which can be of class "Date" or a character vector of another format.}
  \item{dformat}{
A description of the format of \code{d}. This is only used if \code{d} is not a \code{Date} object.
The default format is \code{dd/mm/yyyy}. The only part of the expression used by the function is \code{yyyy}. Examples include \code{yyyy}, \code{mm\yyyy}, \code{mmyyyy}, \code{yyyy-mm} etc.
}
}

%\details{}

\value{
Returns a matrix where each row represents an individual, and where each column represents a year. Observations are coded as 1, while 0 indicates that the animal was not seen during the year.}

\references{
Owen Jones \email{jones@mpg.demogr.de}
}

\author{
%%  ~~who you are~~
}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BayesSurv}}, \code{\link{DataCheck}} 
}

\examples{
#Generate a random dataset of IDs and observation dates
myCensus = data.frame(ID = paste("A",sample(20,100,replace=TRUE),sep=""), obs.date = as.Date(sample(35000:40000,100,replace=TRUE), origin="1900-01-01"))

head(myCenus)
Y = CensusToCaptHist(myCensus$ID,myCensus$obs.date)
Y
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ utilities }
