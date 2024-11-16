\name{aStep1}
\alias{aStep1}
\title{Create Set of Observation Numbers in Step 1 for Linear Model 
Analysis
}
\description{
Derives the first set of observation numbers for forsearch in linear models
}
\usage{
aStep1(yesfactor, df1, df1.ls, inner.rank, initial.sample, formulaA, 
   nofactform, ycol, b.d)
}
\arguments{
  \item{yesfactor}{Logical. TRUE if there are factors in the X matrix
}
  \item{df1}{Data frame being analyzed by forward search. 
}
  \item{df1.ls}{List, each element of which is a factor subset of df1
}
  \item{inner.rank}{Rank of X matrix of lm analysis on entire database
}
  \item{initial.sample}{Number of random samples from which to take set of
       initial observations
}
  \item{formulaA}{Fixed parameter formula of lm function
}
  \item{nofactform}{2-sided formula excluding factor variables
  }
  \item{ycol}{Response column number
}
  \item{b.d}{Index of point to begin diagnostic listings}
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
