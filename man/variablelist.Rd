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
variablelist(datadf, prank)
}
\arguments{
  \item{datadf}{Data frame of independent variables in analysis. First column
       of data frame is Observation number}
  \item{prank}{Number of continuous variables among independent variables}
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
\keyword{  manip }
