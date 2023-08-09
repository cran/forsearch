\name{aStep1}
\alias{aStep1}
\title{
Create Set of Observation Numbers in Step 1 for Linear Model Analysis
}
\description{
Derives the first set of observation numbers for forsearch in linear models
}
\usage{
aStep1(yesfactor, data, inner.rank, initial.sample, formula, ycol, nopl)
}
\arguments{
  \item{yesfactor}{Logical. TRUE if there are factors in the X matrix
}
  \item{data}{Data frame being analyzed by forward search. 
}
  \item{inner.rank}{Rank of X matrix of lm analysis on entire database
}
  \item{initial.sample}{Number of random samples from which to take set of
       initial observations
}
  \item{formula}{Fixed parameter formula of lm function
}
  \item{ycol}{Response column number
}
  \item{nopl}{Number of observations per level of combined factor variables
}
}
\note{Presence of Observation column has no effect on outcome
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
