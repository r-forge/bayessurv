\name{MakeCovariateMatrix}
\alias{MakeCovariateMatrix}

\title{A function to produce a covariate matrix for BayesSurv analysis.
}
\description{
This function produces an n-column covariate matrix from a factor-Class vector with n levels.}
\usage{
MakeCovariateMatrix(x)
}
\arguments{
  \item{x}{A factor-Class vector.}
}
\details{This function produces an n-column covariate matrix from a factor-Class vector with n levels.}

\value{Returns a matrix with a number of rows equal to the number of individuals, and a column for each level in the input vector.}

%\references{ }

\author{Owen Jones \email{jones@mpg.demogr.de}}

%\note{ }

\seealso{
  \code{\link{BayesSurv}},   \code{\link{BayesSurvParall}},   \code{\link{DataCheck}}
  }

\examples{
sex <- as.factor(c("M","F","M","M","M","F","F","F","M","M"))
MakeCovariateMatrix(sex)

}
\keyword{ methods }
