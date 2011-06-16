\name{CohortLT}
\alias{CohortLT}

\title{A function for calculating standard cohort-based life tables from age at death data.
}

\description{
This function calculates a standard cohort-based life table from death ages.
}

\usage{CohortLT(DeathAges, ax = 0.5, n = 1)}

\arguments{
  \item{DeathAges}{A vector of ages at death. These can be precise (e.g. 5.4 yrs) or rounded to the integer year (e.g. 5 yrs).}
  \item{ax}{A vector describing the averaged individual-years lived in the interval by those dying in the interval.   It is only necessary to put the first few values - further values are extrapolated to the correct dimensions. e.g. if ax is set as c(0.3, 0.5) and there are 5 rows in the life tables, ax will become c(0.3, 0.5, 0.5, 0.5, 0.5).}
  \item{n}{n defines the interval width. The default is 1 (yr), but it can be set to be anything numerical value(s). For example, if you want small intervals for early ages, and larger intervals for older ages, you can define n as a sequence such as c(1,1,5) to indicate that the first 2 intervals are 1 year long, then the rest of the intervals are 5 years long.}
}

\details{
Must add further details here.
}

\value{
The function returns a data.frame in the form of a life table with the following column headings.
Start and end of the interval, where start of the interval is defined as x, lx (number left alive at age x), dx (number dying between ages x and x+1), mx (death rate in the cohort between ages x and x+n), ax (average person-years lived in the interval by those dying in the interval), qx (probability of dying between ages x and x+1), px (probability of surviving between ages x and x+1), Lx (person years lived between ages x and x+1), Tx (person-years lived above age x), ex (expectation of life at age x).
}

\references{
Preston et al. 2001. Demography: Measuring and modeling population processes. Blackwell Publishers.
}

\author{Owen Jones \email{jones@mpg.demogr.de}}
%\note{}

%\seealso{ }

\examples{
CohortLT(rpois(100,3))
}

\keyword{ misc }
\keyword{ datagen }