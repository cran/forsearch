\name{aStep2}
\alias{aStep2}
\title{Update Observation Set in Step 2
}
\description{
Derives the next set of observation numbers for forsearch in linear models
}
\usage{
aStep2(thislm, data, ycol, thisi)
}
\arguments{
  \item{thislm}{A lm object 
}
  \item{data}{Data frame being analyzed by forward search
}
  \item{ycol}{Column number of response variable
}
  \item{thisi}{Iteration number
}
}
\details{Support function, usually not called independently
}
\value{Vector of integers corresponding to observation numbers}
\author{William R. Fairweather
}
\keyword{ manip }
