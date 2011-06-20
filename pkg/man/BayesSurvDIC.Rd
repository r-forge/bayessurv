\name{BayesSurvDIC}
\alias{BayesSurvDIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to calculate the Deviance Information Criterion (DIC) from a BayesSurv object.}
\description{
This function calculates the Deviance Information Criterion (DIC) from a BayesSurv object. The DIC is a measure of model fit and can be used to compare different nested models, where a lower DIC reflects a higher predictive power of the model.
}
\usage{
BayesSurvDIC(outBS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outBS}{The output object from \code{\link{BayesSurv}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
\item{D.ave }{Posterior mean deviance}
\item{D.mode }{Deviance at the parameter estimates}
\item{pD }{Effective number of parameters in the model}
\item{k }{Number of parameters in the model}
\item{DIC }{Deviance Information Criterion}
}
\references{
Spiegelhalter D.J., Best N.G. and Carlin B.P. (1998) Bayesian deviance, the effective number o parameters, and the comparison of arbitrarily complex models. Technical Report, MRC Biostatistics Unit, Cambridge, UK.  
}
\author{
Fernando Colchero \email{colchero@mpg.demogr.de}
Owen Jones \email{jones@mpg.demogr.de}
}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BayesSurv}}, \code{\link{BayesSurvParall}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
