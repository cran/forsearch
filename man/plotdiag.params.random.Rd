\name{plotdiag.params.random}
\alias{plotdiag.params.random}
\title{Plot Diagnostic Statistics Of Random Coefficients 
}
\description{
Plot output from forsearch_lme to show change in random coefficients as the number of 
observations in the forward search procedure increases. Save plot in folder containing 
working directory. 
}
\usage{
plotdiag.params.random(forn, coeff.codenums, maintitle = "Put maintitle here", 
subtitle = "Put subtitle here", caption = "Put caption here",wmf = "Put_stored_name_here",
Cairo=TRUE,printgraph = TRUE,
legend = "Dummy legend name", diagnose = FALSE, verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_lme
}
  \item{coeff.codenums}{Numeric vector of random coefficients to include together on the 
  plot. Codes are output by identifyCoeffs function
}
  \item{maintitle}{Main title of plot
}
  \item{subtitle}{Subtitle of plot
}
  \item{caption}{Content of caption
}
  \item{wmf}{File name of stored plot; omit ".wmf"  
}
  \item{Cairo}{TRUE causes use of Cairo graphics
}
  \item{printgraph}{TRUE causes graph to print to file and
          closes device
}
  \item{legend}{Name of legend
}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors
}
  \item{verbose}{If TRUE, indicates beginning and end of function
}
}
\value{Process and plot Cook distance statistics from forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
info3 <- system.file("extdata","MOf.R",package="forsearch");
Machines.O.forlme <- source(info3);
Machines.O.forlme <-Machines.O.forlme[[1]];
plotdiag.params.random(Machines.O.forlme,coeff.codenums=c(3,5,6), 
wmf="Machines_Random_Coefficients",Cairo=FALSE,printgraph=FALSE)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }
