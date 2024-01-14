\name{cStep2}
\alias{cStep2}
\title{Update Observation Set in Step 2}
\description{
Derives the set of observation numbers for step 2 for forsearch in Cox proportional 
     hazard models
}
\usage{cStep2(f.e, finalm, dfa2, ms, rnk2, ss, b.d)
}
\arguments{
  \item{f.e}{Right hand side of formula}
  \item{finalm}{List of rows in model at each stage}
  \item{dfa2}{Complete data frame with factor subset indicator codes}
  \item{ms}{Number of observations in first stage of Step 2}
  \item{rnk2}{Rank of linear analysis with factor variables removed}
  \item{ss}{NULL or vector of observation numbers manually entered into Step 1}
  \item{b.d}{Indicator of starting point for diagnostic listings}
}
\details{Support function, usually not called independently}
\value{Vector of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
