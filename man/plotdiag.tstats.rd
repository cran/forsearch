\name{plotdiag.tstats}
\alias{plotdiag.tstats}
\title{Plot Diagnostic T Statistics
}
\description{
Plot output from forsearch_lm or forsearch_lme to show change in t statistics as the 
  number of observations in the forward search procedure increases. Save plot in folder 
  containing working directory. 
}
\usage{
plotdiag.tstats(forn, coeff.codenums=NULL, maintitle = "Put main title here", 
subtitle = "Put subtitle here", caption="Put caption here", wmf = "Put_stored_name_here", 
Cairo=TRUE, printgraph=TRUE,legend = "Dummy legend name", subdiag=FALSE, subverb=FALSE, 
diagnose = FALSE,verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_lm or forsearch_lme
}
  \item{coeff.codenums}{Numeric vector of coefficients to include together on the plot. 
  Codes are output by identifyFixedCoeffs (for lm files) or by identifyCoeffs function 
  (for lme files)
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
  \item{subdiag}{If TRUE, displays code to help diagnose subfunction errors
}
  \item{subverb}{If TRUE, indicates beginning and end of subfunction
}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors
}
  \item{verbose}{If TRUE, indicates beginning and end of function
}
}
\value{Process and plot t statistics of fixed coefficients from forsearch_lm or 
   forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
info3 <- system.file("extdata","MOf.R",package="forsearch");
Machines.O.forlme <- source(info3)[[1]];
plotdiag.tstats(Machines.O.forlme,coeff.codenums=NULL, wmf="Machines_t_statistics",
Cairo=FALSE,printgraph=FALSE
)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }