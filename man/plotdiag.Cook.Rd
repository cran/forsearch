\name{plotdiag.Cook}
\alias{plotdiag.Cook}
\title{Plot Diagnostic Statistics of Modified Cook's Distance 
}
\description{
Plot output from forsearch_lm or forsearch_lme to show change in Modified Cook's distance 
as the number of observations in the forward search procedure increases. Save plot in 
folder containing working directory. 
}
\usage{
plotdiag.Cook(forn, maintitle = "Put main title here", subtitle = "Put subtitle here", 
wmf = "Put_plot_file_title_here", Cairo=TRUE,printgraph=TRUE, loess = FALSE, 
subdiag = FALSE,subverb = FALSE,diagnose = FALSE, verbose = TRUE)
}
\arguments{
  \item{forn}{
Name of forward search output file
}
  \item{maintitle}{
Main title of plot
}
  \item{subtitle}{
Subtitle of plot
}
  \item{wmf}{
File name of stored plot; omit ".wmf"  
}
\item{Cairo}{TRUE causes use of Cairo graphics
}
\item{printgraph}{TRuE causes graph to print to file and
      closes device
}
\item{loess}{
If TRUE, adds loess curve to plot, otherwise, straight line
}
  \item{subdiag}{
If TRUE, displays code to help diagnose subfunction errors
}
  \item{subverb}{
If TRUE, indicates beginning and end of function
}
  \item{diagnose}{
If TRUE, displays code to help diagnose main function errors
}
  \item{verbose}{
If TRUE, indicates beginning and end of function
}
}
\value{
Process and plot Cook distance statistics from forsearch_lm or forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
%print("Cook");
info3 <- system.file("extdata","Alfalfa.O.forlme.R",package="forsearch");
%print(info3);
Alfalfa.O.forlme <- source(info3);
Alfalfa.O.forlme <- Alfalfa.O.forlme[[1]];
plotdiag.Cook(Alfalfa.O.forlme, wmf="Alfalfa_Cook", Cairo=FALSE,
printgraph=FALSE)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }
