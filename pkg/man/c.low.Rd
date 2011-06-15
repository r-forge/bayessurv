\name{c.low}
\alias{c.low}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
c.low(th, idm = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{th}{
%%     ~~Describe \code{th} here~~
}
  \item{idm}{
%%     ~~Describe \code{idm} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(th, idm=3){
	if(idm==1) cl = 0
	if(idm==2) cl = ifelse(th[5] > 0, -exp(th[4]), 0)
	if(idm==3){
		x.minf = (th[1]+log(th[2]) - th[4]-log(th[5]))/(th[2] + th[5])
		cl     = -exp(th[1]-th[2]*(x.minf)) - exp(th[4]+th[5]*(x.minf))
	}
	return(cl)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
