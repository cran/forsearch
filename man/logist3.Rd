\name{logist3}
\alias{logist3}
\title{Three Parameter Logistic}
\description{
Calculates three parameter logistic function
}
\usage{
logist3(x, a, b, c)
}
\arguments{
  \item{x}{Independent variable}
  \item{a}{First parameter. See details}
  \item{b}{Second parameter}
  \item{c}{Third parameter}
}
\details{Calculates y = a/(1 + xx), where xx = exp(-(x - b)/c)
}
\value{
Numerical output of the function described in details
}
\author{
William R. Fairweather
}
\keyword{ manip }
