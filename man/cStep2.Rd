\name{cStep2}
\alias{cStep2}
\title{Update Observation Set in Step 2}
\description{
Derives the next set of observation numbers for forsearch in Cox proportional hazard models
}
\usage{
cStep2(df1, rim, formula.elements)
}
\arguments{
  \item{df1}{Complete data frame}
  \item{rim}{Vector of Observation numbers already in dataset}
  \item{formula.elements}{Character vector of names of independent variables in model}
}
\details{Support function, usually not called independently}
\value{Vector of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
