\name{cStep2}
\alias{cStep2}
\title{Update Observation Set in Step 2}
\description{
Derives the set of observation numbers for step 2 for forsearch in Cox proportional 
     hazard models
}
\usage{cStep2(fe, finalm, rimbs, dfa2, onlyfactor=FALSE, ycol, cphties,mstart, rnk, b.d)
}
\arguments{
  \item{fe}{Right hand side of formula}
  \item{finalm}{List of rows in model at each stage}
  \item{rimbs}{List, each element is a complete matrix of obs numbers and 
       corresponding subset codes}
  \item{dfa2}{Complete data frame with factor subset indicator codes}
  \item{onlyfactor}{Logical. TRUE if there are no continuous independent variables}
  \item{ycol}{Response column number}
  \item{cphties}{Character designation of method of handling ties}
  \item{mstart}{Number of observations in first stage of Step 2}
  \item{rnk}{Rank of linear analysis with factor variables removed}
  \item{b.d}{Indicator of starting point for diagnostic listings}
}
\details{Support function; usually not called independently}
\value{Vector of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
