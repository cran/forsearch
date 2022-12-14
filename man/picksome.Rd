\name{picksome}
\alias{picksome}
\title{Structured Sets of Random Samples of Observations 
}
\description{
Restricts Step 1 of forward search procedures to ensure that every possible
   combination of levels is included for observations containing factors
}
\usage{
picksome(subsetlist, nobs, initial.sample, n.obs.per.level, rank, 
     verbose)
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
  \item{verbose}{If TRUE, indicates beginning and end of function run.
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
\examples{
V1 <- rep(c(3,7,9,11),times=4);
V2 <- V1+3;
V3 <- V1+90;
V1 <- as.factor(V1);
V2 <- as.factor(V2);
V3 <- as.factor(V3);
Observation <- 1:16;
testdf <- data.frame(Observation, V1,V2,V3);
%print(testdf);
testoutdf <- variablelist(testdf);
picksome(testoutdf,nobs=16,initial.sample=11,n.obs.per.level=1,rank=4,
   verbose=FALSE)
}
\keyword{ manip }
