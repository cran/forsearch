\name{bStep2}
\alias{bStep2}
\title{Update Observation Numbers in Step 2}
\description{
Derives the next set of observation numbers for forsearch in linear mixed effects models
}
\usage{
bStep2(fixed2, fulldata, random2, yf, nOuter, mstart, nobs, yobs, fbg, n.f, s.o, ras, b.d, 
          ufixdat, LLL, verbose)
}
\arguments{
  \item{fixed2}{Fixed parameter formula}
  \item{fulldata}{Complete data set}
  \item{random2}{Random parameter formula}
  \item{yf}{Indicator of presence of internal factors in model}
  \item{nOuter}{Number of outer subgroups}
  \item{mstart}{Number of observations in each outer subgroup}
  \item{nobs}{Number of observations in entire database}
  \item{yobs}{Column number of response variable}
  \item{fbg}{List of observation numbers by outer subgroup}
  \item{n.f}{Vector of number of observation in each outer subgroup}
  \item{s.o}{Original observation numbers prior to renumbering in each outer subgroup}
  \item{ras}{List of observation numbers in each outer subgroup}
  \item{b.d}{Indicator of place in code to begin diagnostic printouts}
  \item{ufixdat}{Unique outer or outer and inner grouping codes}
  \item{LLL}{Empty list of observation number sets}
  \item{verbose}{TRUE causes printing of function ID at beginning and end of run}
}
\details{Support function, usually not called independently}
\value{List of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
