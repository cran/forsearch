\name{picksome}
\alias{picksome}
\title{Structured Sets of Random Samples of Observations 
}
\description{
Restricts Step 1 of forward search procedures to ensure that every possible
   combination of levels is included for observations containing factors
}
\usage{
picksome(subsetlist, nobs, initial.sample, n.obs.per.level, rank)
}
\arguments{
  \item{subsetlist}{List, each element is a data frame of 2 columns with code 
     indicating the highest possible level of interaction to which each 
     observation can belong. Usually, output from variablelist function.
}
  \item{nobs}{Number of observations in data frame containing observations of 
     forward search.
}
  \item{initial.sample}{Number of randomized sets of observations in Step 1 of
      forward search.
}
  \item{n.obs.per.level}{Number of observations to pull from each level.
}
  \item{rank}{Rank of the X matrix of the analytical function to be used on 
     data frame.
}
}
\details{
Support function, usually not called independently. Argument n.obs.per.level is
   set by user in forsearch_xxx function call.
}
\value{
Matrix, each row of which identifies observations in each set of random sample  
   of observations.
}
\author{William R. Fairweather
}
\keyword{ manip }
