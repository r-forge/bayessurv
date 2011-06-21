\name{basumomulti}

\alias{basumomulti}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to conduct a basumo analyis with the option of running in parallel.}


\description{
A function to perform Bayesian estimation of latent states (birth and death year) using an undelying mortality model.
Further details required!!}

\usage{
basumomulti(Data, ststart, stend, nsim = 5, parallel = FALSE, ncpus = 2, ini.pars.mat = NULL, model = "SI", niter = 50000, burnin = 5001, thinning = 50, rptp = ststart, jumps = NULL, priors = NULL, lifetable = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
A \code{data.frame} to be used as an input data file for \code{basumo}. The first column is the individual's ID, the second and third columns are birth and death years respectively. Columns 4-(n-1) represent the observation window of n years. This is followed (optionally) by columns for covariates, which must be of 'factor' or 'numeric' class.
}
  \item{ststart}{
The first year of the observation matrix.
}
  \item{stend}{
The last year of the observation matrix.
}
  \item{nsim}{
Definition required.
}
  \item{parallel}{
A logical argument indicating whether to run the analyses in parallel (TRUE) or serial (FALSE).
}
  \item{ncpus}{
The number of CPUs to use in the analysis
}
  \item{ini.pars.mat}{
A matrix of initial parameters. Number of columns should be appropriate to the number of terms in the model (2 for Gompertz, 3 for Gompertz-Makeham, 5 for Siler).
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
If one wanted to set a different recapture probability for years after (say) 1990, the RPTP could be defined as c(1985, 1990) etc..
}
  \item{jumps}{
A vector defining the size of jumps for each parameter. The jump specifications depend on the model selected. for Siler are specified in the following order:   \code{c(alpha1, beta1, c, alpha2, beta2)}
The jumps for Gompertz-Makeham are specified as \code{c(c,alpha,beta)}, and the jumps for Gompertz are specified as \code{c(alpha,beta)}. 
Any additional elements in the vector are ignored (e.g. if a Gompertz model is specified, but a vector with 5 elements is provided, only the first two will be used (as alpha and beta jumps respectively)).
}
  \item{priors}{
A vector defining the priors for each parameter. The vector takes the same format as described for \code{jumps}.
}
  \item{lifetable}{
A logical argument specifying whether a cohort-based life table should be produced or not. The life table produced has the same format as produced by \code{CohortLT}.
}
}

\details{
Further details will be added here.
}

\value{
The function returns a "basumo" list object with the following elements:
(Copied from basumo - probably needs changing)

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

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{basumo}}, \code{\link{DataCheck}}, \code{\link{CohortLT}}
}

\examples{}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
