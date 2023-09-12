\name{plotdiag.loglik}
\alias{plotdiag.loglik}
\title{Plot Diagnostic Statistics of LOGLIK Output of COXPH Function 
}
\description{
Plot output from forsearch_cph to show change in loglik pairs as the number of 
observations in the forward search procedure increases. Save plot in 
folder containing working directory. 
}
\usage{
plotdiag.loglik(forn,  
     maintitle= "Put main title here", 
     subtitle= "Put subtitle here" , 
     caption="Put caption here",
     wmf = "Put_stored_name_here", 
     Cairo=TRUE,
     printgraph = TRUE,
     verbose=TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_cph}
  \item{maintitle}{Main title of plot}
  \item{subtitle}{Subtitle of plot}
  \item{caption}{Content of caption}
  \item{wmf}{File name of stored plot; omit ".wmf"}
  \item{Cairo}{TRUE causes use of Cairo graphics}
  \item{printgraph}{TRUE causes graph to print to file and
          closes device}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\value{{Process and plot Wald Test statistics from forsearch_cph}
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\keyword{ hplot }
