\name{plotdiag.residuals}
\alias{plotdiag.residuals}
\title{Plot Diagnostic Statistics Of Residuals Or Squared Residuals
}
\description{
Plot output from forsearch_lm or forsearch_lme to show change in residuals or squared 
residuals as the number of observations in the forward search procedure increases. Save 
plot in folder containing working directory. 
}
\usage{
plotdiag.residuals(forn, squared = FALSE, hilos = c(1, 0), maintitle, subtitle, 
caption, wmf, Cairo=TRUE,printgraph=TRUE,
legend = "Dummy legend name", diagnose = FALSE, verbose = TRUE)
}
\arguments{
  \item{forn}{Name of forward search output file
}
  \item{squared}{TRUE causes residuals to be squared before plotting
}
  \item{hilos}{Number of observations having high and number having low values of 
  residuals to identify. No low values are identified for squared residual plot. 
}
  \item{maintitle}{Main title of plot
}
  \item{subtitle}{Subtitle of plot
}
  \item{caption}{Caption of plot}
  
  \item{wmf}{File name of stored plot; omit ".wmf" 
}
  \item{Cairo}{TRUE causes use of Cairo graphics
}
  \item{printgraph}{TRUE causes graph to print to file and
          closes device
}
  \item{legend}{Legend title
}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors
}
  \item{verbose}{If TRUE, indicates beginning and end of function
}
}
\value{
Process and plot changes in residuals or squared residuals from forsearch_lm or 
forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
\testonly{
info3 <- system.file("extdata","MOf.R",package="forsearch");
Machines.O.forlme <- source(info3);
Machines.O.forlme <- Machines.O.forlme[[1]];
plotdiag.residuals(Machines.O.forlme, hilos=c(1,1),wmf="Machines Residuals",Cairo=FALSE,
printgraph=FALSE)
}
}
 \keyword{ hplot }
