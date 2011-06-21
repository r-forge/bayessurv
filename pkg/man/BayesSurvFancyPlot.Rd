\name{BayesSurvFancyPlot}
\alias{BayesSurvFancyPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to produce pretty plots of survival and mortality trajectories accompanied by parameter estimate probability distributions.}

\description{
A function to produce pretty plots of survival and mortality trajectories accompanied by parameter estimate probability distributions.
}

\usage{
BayesSurvFancyPlot(outBS)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outBS}{
A BayesSurv object (produced with the function \code{BayesSurv}).
}
}
\details{}

\value{A plot is produced.}

\references{
Colchero and Clark 2011 Bayesian inference on age-specific survival from capture-recapture data for censored and truncated data. Journal of Animal Ecology. 

Colchero, Jones and Rebke. BayesSurv - a package to estimate mortality parameters from incomplete mark-recapture data with covariates. Unpublished Manuscript.
}

\author{
Fernando Colchero \email{colchero@mpg.demogr.de}
}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{BayesSurv}}}

\examples{}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
