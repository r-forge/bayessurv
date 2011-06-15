\name{BayesSurvDIC}
\alias{BayesSurvDIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BayesSurvDIC(outBS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outBS}{
%%     ~~Describe \code{outBS} here~~
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
function(outBS){
	
# Calculate DIC
thin	    = seq(outBS$bng,outBS$ng,outBS$thint)
modepost   = outBS$post[thin,3]
L          = length(thin)
Dm         = -2*modepost
Dmode      = -2*modepost[which(modepost==max(modepost))[1]]
Dave       = mean(Dm)
pD         = Dave - Dmode
k          = ncol(outBS$pi)+sum(outBS$modm[outBS$idm,]*ncol(outBS$Z))
DIC        = 2*Dave - Dmode
ModSel     = c(Dave, Dmode, pD, k, DIC)
names(ModSel) = c("D.ave", "D.mode", "pD", "k", "DIC")

	
#Return ModSel
return(ModSel)
	
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
