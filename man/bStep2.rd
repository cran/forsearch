\name{bStep2}
\alias{bStep2}
\title{Update Observation Numbers in Step 2}
\description{
Derives the set of Step 2 observation numbers for forsearch in linear mixed 
          effects models
}
\usage{bStep2(yf, f2, dfa2, randm2, onlyfactor = FALSE,ms, ycol, initn, inc, finalm, fbg, b.d)
}
\arguments{
  \item{yf}{Logical. Indicates presence of factor variables}
  \item{f2}{Fixed parameter formula}
  \item{dfa2}{Complete data set with factor subset identification codes}
  \item{randm2}{Random parameter formula}
  \item{onlyfactor}{TRUE if there are no continuous independent variables in 
         the model}
  \item{ms}{Number of observations beginning Step 2}
  \item{ycol}{Column number of response variable}
  \item{initn}{Vector of number of observations from each group or fixed factor
         subset to draw for primary stage of step 2}
  \item{inc}{Logical. TRUE causes relaxation of lmeControl}
  \item{finalm}{List of expanding subset observation numbers}
  \item{fbg}{List of observation numbers by factor subgroup}
  \item{b.d}{Indicator of place in code to begin diagnostic printouts}
}
\details{Support function, usually not called independently}
\value{List of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
