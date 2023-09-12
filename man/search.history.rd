\name{search.history}
\alias{search.history}
\title{Create Tabular History Of Forward Search}
\description{
The forward search functions output a list of vectors, each of which indicates which observations
are in the model at each stage of the search. This function processes that list to create a more easily
understood matrix of the observation numbers that are newly entered into the model and any that
were temporarily removed from the model over the course of the search.
}
\usage{
search.history(list1, verbose = TRUE)
}
\arguments{
   \item{list1}{Name of a forsearch_xxx output file}
   \item{verbose}{If TRUE, indicates beginning and end of function}
}
\value{Printout of matrix showing evolution of observations to enter or leave 
   the model during the course of the forward search}
\author{William R. Fairweather
}
\examples{
info3 <- system.file("extdata", "crossdata.for1.R", package="forsearch");
info3 <- source(info3);
info3 <- info3[[1]];
search.history(list1=info3, verbose=TRUE)
}
\keyword{ manip }
