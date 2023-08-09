\name{cStep2}
\alias{cStep2}
\title{Update Observation Numbers in Step 2}
\description{
Derives the next set of observation numbers for forsearch in Cox proportional hazard models
}
\usage{
cStep2(df1, rim, formula.elements, r2)
}
\arguments{
  \item{df1}{Complete data frame}
  \item{rim}{Vector of Observation numbers already in dataset}
  \item{formula.elements}{Character vector of names of independent variables in model}
  \item{r2}{Value of correlation required by Hmisc::redun to declare variables to be
                 redundant}
}
\details{Support function, usually not called independently}
\value{Vector of expanding number sets corresponding to observation numbers}
\author{William R. Fairweather}
\keyword{ manip }
