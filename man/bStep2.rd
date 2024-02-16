\name{bStep2}
\alias{bStep2}
\title{Update Observation Numbers in Step 2}
\description{
Derives the set of Step 2 observation numbers for forsearch in linear mixed 
          effects models
}
\usage{bStep2(f2, dfa2, randm2, ms, finalm, fbg, b.d, rnk2)
}
\arguments{
  \item{f2}{Fixed parameter formula}
  \item{dfa2}{Complete data set with factor subset identification codes}
  \item{randm2}{Random parameter formula}
  \item{ms}{Number of observations beginning Step 2}
  \item{finalm}{List of expanding subset observation numbers}
  \item{fbg}{List of observation numbers by factor subgroup}
  \item{b.d}{Indicator of place in code to begin diagnostic printouts}
  \item{rnk2}{Rank of linear regression with factor variables eliminated}
}
\details{Support function, usually not called independently}
\value{List of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
