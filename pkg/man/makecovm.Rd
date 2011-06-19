\name{makecovm}
\alias{makecovm}

\title{
A function to produce a covariate matrix for BayesSurv}

\description{
This function produces a covariate matrix for a BayesSurv analysis from a given dataset with covariates.
The covariates must be 'numeric' or 'factor' class.
}

\usage{
makecovm(Data, covs)
}

\arguments{
  \item{Data}{
A data frame with covariates}
  \item{covs}{
A character vector with the names of the desired covariates (e.g. c("Sex","Location","Weight"))
}
}

%\details{}

\value{
Returns a matrix to be used as the covariate matrix for BayesSurv
}

%\references{ }

\author{
Owen Jones \email{jones@mpg.demogr.de}
}

\note{
This is really intended as an internal function to be called by \code{\link{BayesSurv}}}

\seealso{
\code{\link{BayesSurv}},\code{\link{BayesSurvParall}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}