\name{bStep2}
\alias{bStep2}
\title{Update Observation Numbers in Step 2}
\description{
Derives the next set of observation numbers for forsearch in linear mixed effects models
}
\usage{
bStep2(fixed, nOuter, mnf, mstart, nobs, yobs, fbg, n.f, s.o, ras, b.d, verbose)
}
\arguments{
  \item{fixed}{Fixed parameter formula of lm function}
  \item{nOuter}{Number of outer subgroups}
  \item{mnf}{Maximum number of observations in an outer subgroup}
  \item{mstart}{Number of observations in each outer subgroup}
  \item{nobs}{Number of observations in entire database}
  \item{yobs}{Column number of response variable}
  \item{fbg}{List of observation numbers by outer subgroup}
  \item{n.f}{Vector of number of observation in each outer subgroup}
  \item{s.o}{Original observation numbers prior to renumbering in each outer subgroup}
  \item{ras}{List of observation numbers in each outer subgroup}
  \item{b.d}{Indicator of place in code to begin diagnostic printouts}
  \item{verbose}{TRUE causes printing of function ID at beginning and end of run}
}
\details{Support function, usually not called independently}
\value{List of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
