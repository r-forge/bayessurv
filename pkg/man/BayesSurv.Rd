\name{BayesSurv}
\alias{BayesSurv}

\title{A function to perform Bayesian estimation of latent stats (birth and death year) using an undelying mortality model.
}

\description{
A function to perform Bayesian estimation of latent stats (birth and death year) using an undelying mortality model.
Further details required!!
}

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
The underlying mortality model used. "GO"= Gompertz, "GM" = Gompertz-Makeham and "SI" = Siler.
}

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
RPTP is a vector defining the recapture probability transition points. These are points (years) where the recapture probability is thought to change.  The default setting is for the recapture probability to be constant throughout the study, so the RPTP is simply defined as a single element vector of the first year of the observation period (e.g. c(1985)).
If one wanted to set a different recapture probability for years after (say) 1990, the RPTP could be defined as \code{c(1985, 1990)} etc.
}
  \item{jumps}{
A vector defining the size of jumps for each parameter. The jump specifications depend on the model selected. for Siler are specified in the following order:   \code{c(alpha1, beta1, c, alpha2, beta2)}.
The jumps for Gompertz-Makeham are specified as \code{c(c,alpha,beta)}, and the jumps for Gompertz are specified as c(alpha,beta). 
Any additional elements in the vector are ignored (e.g. if a Gompertz model is specified, but a vector with 5 elements is provided, only the first two will be used (as alpha and beta jumps respectively)).
}
  \item{ini.pars}{
A vector defining the initial parameter values for each parameter. The vector takes the same format as described for \code{jumps}.
}

  \item{priors}{
A vector defining the priors for each parameter. The vector takes the same format as described for \code{jumps}.
}

  \item{lifetable}{
A logical argument specifying whether a cohort-based life table should be produced or not. The life table produced has the same format as produced by \code{CohortLT}.
}

  \item{datacheck}{
A logical argument specifying whether the input data should be checked for errors or not.
}
}

\details{
Further details will be added here.
}

\value{
The function returns a "BayesSurv" list object with the following elements:

\item{bd }{Description}
\item{Y }{Description}
\item{Z }{Description}
\item{postm }{Description}
\item{g }{Description}
\item{ng }{Description}
\item{bng }{Description}
\item{thint }{Description}
\item{theta }{Description}
\item{pi }{Description}
\item{bis }{Description}
\item{dis }{Description}
\item{thsum }{Description}
\item{pisum }{Description}
\item{xqsum }{Description}
\item{Sxsum }{Description}
\item{mxsum }{Description}
\item{modm }{Description}
\item{idm }{Description}
\item{jumps }{Description}
\item{ini.pars }{Description}
\item{priors }{Description}
\item{LT }{Description}

}

\references{
Colchero and Clark 2011 Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. Journal of Animal Ecology. 
Colchero, Jones and Rebke. BayesSurv - a package to estimate mortality parameters from incomplete mark-recapture data with covariates. Unpublished Manuscript.
}

\author{
Fernando Colchero \email{colchero@mpg.demogr.de}
}

\seealso{
\code{\link{BayesSurvParall}}, \code{\link{BayesSurvPlots}}, \code{\link{BayesSurvDIC}}, \code{\link{DataCheck}}, \code{\link{CohortLT}}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

\keyword{ methods }

