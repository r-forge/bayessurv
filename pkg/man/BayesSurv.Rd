\name{BayesSurv}
\alias{BayesSurv}

\title{A function to perform Bayesian estimation of latent stats (birth and death year) using an undelying mortality model.}
\description{ }
\usage{
BayesSurv(Data, ststart, stend, model = "SI", niter = 50000, burnin = 5001, thinning = 50, rptp = ststart, jumps = NULL, ini.pars = NULL, priors = NULL, lifetable = TRUE, datacheck = TRUE)
}

\arguments{
  \item{Data}{
A \code{data.frame} to be used as an input data file for BayesSurv. The first column is the individual's ID, the second and third columns are birth and death years respectively. Columns 4-(n-1) represent the observation window of n years. This is followed (optionally) by a column for a covariate.
}
  \item{ststart}{
The first year of the observation matrix.
}
  \item{stend}{
The last year of the observation matrix.
}
  \item{model}{
The underlying mortality model used. "GO"= Gompertz, "GM" = Gompertz-Makeham and "SI" = Siler.}
  \item{niter}{
The number of iterations for the MCMC process.
}
  \item{burnin}{
Number of iterations for the burn in.
}
  \item{thinning}{
The thinning parameter.
}
  \item{rptp}{
%%     ~~Describe \code{rptp} here~~
}
  \item{jumps}{
A vector defining the size of jumps for each parameter. The jump specifications depend on the model selected. for Siler are specified in the following order:   c(alpha1, beta1, c, alpha2, beta2)
The jumps for Gompertz-Makeham are specified as c(c,alpha,beta), and the jumps for Gompertz are specified as c(alpha,beta). Any additional elements in the vector are ignored (e.g. if a Gompertz model is specified, but a vector with 5 elements is provided, only the first two will be used (as alpha and beta jumps respectively).
}
  \item{ini.pars}{
A vector defining the initial parameter values for each parameter. The vector takes the same format as described for \code{jumps}.}
  \item{priors}{
A vector defining the priors for each parameter. The vector takes the same format as described for \code{jumps}.}
}
  \item{lifetable}{
A logical argument specifying whether a cohort-based life table should be produced or not. The life table produced has the same format as produced by \code{CohortLT}.
}
  \item{datacheck}{
A logical argument specifying whether the input data should be checked for errors or not.}
}
\details{

Further details will be added here.

}
\value{
The function returns a "BayesSurv" list object with the following elements:

\item{bd }{Description of 'comp1'}
\item{Y }{Description of 'comp1'}
\item{Z }{Description of 'comp1'}
\item{postm }{Description of 'comp1'}
\item{g }{Description of 'comp1'}
\item{ng }{Description of 'comp1'}
\item{bng }{Description of 'comp1'}
\item{thint }{Description of 'comp1'}
\item{theta }{Description of 'comp1'}
\item{pi }{Description of 'comp1'}
\item{bis }{Description of 'comp1'}
\item{dis }{Description of 'comp1'}
\item{thsum }{Description of 'comp1'}
\item{pisum }{Description of 'comp1'}
\item{xqsum }{Description of 'comp1'}
\item{Sxsum }{Description of 'comp1'}
\item{mxsum }{Description of 'comp1'}
\item{modm }{Description of 'comp1'}
\item{idm }{Description of 'comp1'}
\item{jumps }{Description of 'comp1'}
\item{ini.pars }{Description of 'comp1'}
\item{priors }{Description of 'comp1'}
\item{LT }{Description of 'comp1'}

}
\references{
Colchero and Clark 2011 Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. Journal of Animal Ecology. 
Colchero, Jones and Rebke. BayesSurv - a package to estimate mortality parameters from incomplete mark-recapture data with covariates. Unpublished Manuscript.
}
\author{
Fernando Colchero \email{colchero@mpg.demogr.de}
}
\note{}


\seealso{
\code{\link{BayesSurvParall}}, \code{\link{BayesSurvPlots}}, \code{\link{BayesSurvDIC}}, \code{\link{DataCheck}}, \code{\link{CohortLT}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

\keyword{ methods }

