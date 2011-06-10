\name{CohortLT}
\alias{CohortLT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CohortLT(DeathAges, ax = 0.5, n = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DeathAges}{
%%     ~~Describe \code{DeathAges} here~~
}
  \item{ax}{
%%     ~~Describe \code{ax} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
function(DeathAges, ax = 0.5, n = 1) {
    
    # DeathAges should be a numeric integer vector (e.g. years)
    #
    # ax is a vector describing the average proportion of year lived by an individual
    # dying during the interval. It is only necessary to put the first few values -
    # further values are extrapolated to the correct dimensions.
    # e.g. if ax is set as c(0.3, 0.5) and there are 5 rows in the life tables, ax
    # will become c(0.3, 0.5, 0.5, 0.5, 0.5).
    # If we have precise ages at death it is possible to calculate this directly from the data.
    # I could add this in later on if necessary.
    #
    # n defines the interval width. The default is 1 (yr), but you can set it to be anything.
    # If you want small intervals for early ages, and larger intervals for older ages you can put
    # a sequence such as c(1,1,5) to indicate that the first 2 intervals are 1 year long, then
    # the rest of the intervals are 5 years long.
    
    max.age = max(DeathAges, na.rm = TRUE)
    ageseq = c(0, cumsum(n))
    
    if (max(ageseq) > max.age) {
        wa = which(ageseq < max.age)
        ageseq = ageseq[c(wa, max(wa) + 1)]
    }
    else {
        fr = ageseq[length(ageseq)] + n[length(n)]
        to = ceiling((max.age + 1)/n[length(n)]) * n[length(n)]
        b = n[length(n)]
        ageseq = c(ageseq, seq(fr, to, b))
    }
    
    int = cbind(ageseq[1:(length(ageseq) - 1)], ageseq[2:(length(ageseq))])
    colnames(int) = c("StartAge", "EndAge")
    nx = int[, 2] - int[, 1]
    
    if (length(ax) > nrow(int)) {
        ax = ax[1:nrow(int)]
    }
    if (length(ax) < nrow(int)) {
        ax = c(ax, rep(ax[length(ax)], nrow(int) - length(ax)))
    }
    
    Nx = NULL
    Nx.plus.n = NULL
    
    for (x in 1:nrow(int)) {
        Nx[x] = length(DeathAges[DeathAges >= int[x, 1]])
        Nx.plus.n[x] = length(DeathAges[DeathAges >= int[x, 2]])
    }
    
    Dx = Nx - Nx.plus.n
    dx = Dx
    lx = Nx
    qx = dx/lx
    qx[length(qx)] = 1
    px = 1 - qx
    Lx = ((lx - dx) * nx) + (ax * dx)
    mx = dx/Lx
    Tx = rev(cumsum(rev(Lx)))
    ex = Tx/lx
    
#    return(data.frame(cbind(int, Nx, Dx, mx, ax, qx, px, lx, 
#        dx, Lx, Tx, ex)))
    return(data.frame(cbind(int, lx, dx, mx, ax, qx, px,  
        Lx, Tx, ex)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
