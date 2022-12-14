\name{variablelist}
\alias{variablelist}
\title{Identify Level(s) to Which Each Factor Observation Belongs
}
\description{
For a data frame with factor variables V1, V2, V3, etc having levels n1,
   n2, n3, etc, lists the n1*n2*n3*... possible interaction levels and
   identifies which of the observations of the data frame belong in which
   of these interaction levels.
}
\usage{
variablelist(datadf, verbose = TRUE)
}
\arguments{
  \item{datadf}{Data frame of independent variables in analysis. First column
       of data frame is Observation number}
  \item{verbose}{If TRUE, indicates beginning and end of function run}
}
\details{
Support function, usually not called independently
}
\value{
List, each element is a data frame of 2 columns with code indicating the highest
    possible level of interaction to which each observation can belong
}
\author{
William R. Fairweather
}
\examples{
V1 <- rep(c(3,7,9,11),times=4)
V2 <- V1+3
V3 <- V1+90
V1 <- as.factor(V1)
V2 <- as.factor(V2)
V3 <- as.factor(V3)
Observation<-1:16
testdf <- data.frame(Observation, V1,V2,V3)
variablelist(testdf)
}
\keyword{  manip }
