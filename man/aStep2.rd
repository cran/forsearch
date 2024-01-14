\name{aStep2}
\alias{aStep2}
\title{Update Observation Set in Step 2
}
\description{
Derives the set of observation numbers for forsearch in Step 2 for linear models
}
\usage{
aStep2(yesfactor, form.A2, finalm, rimbs, dfa2, ycol, mstart, rnk, b.d)
}
\arguments{
  \item{yesfactor}{True or False for presence of factors
}
  \item{form.A2}{Formula for analysis of entire dataset
}
  \item{finalm}{See VALUE above. finalm argument is the same but only for Step 1 values
}
  \item{rimbs}{List, each element is a matrix of obs numbers and corresponding subset codes
}
  \item{dfa2}{Data frame being analyzed by forward search. Presence of Observation column has
       no effect on output
}
  \item{ycol}{Response column number, including 1 for Observation 
}
  \item{mstart}{Number of first subset to be defined in Step 2
}
  \item{rnk}{Rank of X matrix. For factors, this is rank with factors removed.
}
  \item{b.d}{Number at which to begin diagnostic listings 
}
}

\details{Support function, usually not called independently
}
\value{Vector of integers corresponding to observation numbers}
\author{William R. Fairweather
}
\keyword{ manip }
