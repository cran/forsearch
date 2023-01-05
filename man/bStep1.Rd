\name{bStep1}
\alias{bStep1}
\title{Create Set of Observation Numbers in Step 1 for Linear Mixed Effects 
         `Model Analysis
}
\description{
Derives the first set of observation numbers for forsearch_lme in linear mixed 
         effects models
}
\usage{
bStep1(fixed, yf, mnf, nOuter, yobs, s.o, nopl, nobs, i.s, fbg, b.d, verbose=TRUE)
}
\arguments{
  \item{fixed}{Fixed parameter formula of lm function}
  \item{yf}{Indicator of presence of factors in the X matrix}
  \item{mnf}{Maximum number of observations in any outer subgroup}
  \item{nOuter}{Number of outer subgroups}
  \item{yobs}{Column number of response variable}
  \item{s.o}{Vector of original observation numbers}
  \item{nopl}{Number of observations to take for each subgroup. Default is 1}
  \item{nobs}{Total number of observations in database}
  \item{i.s}{Number of random sets of observations in Step 1}
  \item{fbg}{Observation numbers for each outer subgroup}
  \item{b.d}{Numeric indicator of position in code to begin diagnostic prints}
  \item{verbose}{TRUE causes printing of function ID at beginning and end of run}
}
\details{Support function, usually not called independently
}
\value{
Produces set of observation numbers for Step 1. Accounts for presence of factors
    in the dataset
}
\author{William R. Fairweather
}
\keyword{ manip }
