\name{bStep1}
\alias{bStep1}
\title{Create Set of Observation Numbers in Step 1 for Linear Mixed
     Effects Model Analysis
}
\description{
Derives the first set of observation numbers for forsearch in linear mixed
     effects models
}
\usage{
bStep1(yesfactor, df1, df1.ls, inner.rank, initial.sample, formula, randform, 
     ycol, nopl, b.d)
}
\arguments{
  \item{yesfactor}{Logical. TRUE if there are factors in the X matrix
}
  \item{df1}{Data frame being analyzed by forward search. 
}
  \item{df1.ls}{List, each element of which is a factor subset of df1
}
  \item{inner.rank}{Rank of X matrix of lme analysis on entire database
}
  \item{initial.sample}{Number of random samples from which to take set of
       initial observations
}
  \item{formula}{Two-sided fixed parameter formula of lme function
}
  \item{randform}{One-sided random effects formula
  }
  \item{ycol}{Response column number
}
  \item{nopl}{Number of observations per level of combined factor variables
}
  \item{b.d}{Index of point to begin diagnostic listings}
}
\note{Presence of Observation column has no effect on outcome
}
\details{Support function, usually not called independently
}
\value{
Produces set of observation numbers for Step 1. Accounts for presence of 
     factors in the dataset 
}
\author{William R. Fairweather
}
\keyword{ manip }
