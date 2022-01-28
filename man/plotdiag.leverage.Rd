\name{plotdiag.leverage}
\alias{plotdiag.leverage}
\title{Plot Diagnostic Statistics Of Leverage
}
\description{
Plot output from forsearch_lm or forsearch_lme to show change in leverage of each 
observation as the number of observations in the forward search procedure increases. 
Save plot in folder containing working directory.
}
\usage{
plotdiag.leverage(forn, hilos = c(1, 0), maintitle = "Put main title here", 
subtitle = "Put subtitle here", wmf = "Put_graph_title_here", Cairo=TRUE,
printgraph = TRUE,subdiag = FALSE, subverb = FALSE, diagnose = FALSE, 
verbose = TRUE)
}
\arguments{
  \item{forn}{Name of forward search output file
}
  \item{hilos}{Vector with number of highest observations and number of lowest obserations 
  on graph to identify
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
%print("leverage");
info3 <- system.file("extdata","MOf.R",package="forsearch");
%print(info3);
Machines.O.forlme <- source(info3);
Machines.O.forlme <-Machines.O.forlme[[1]];
plotdiag.leverage(Machines.O.forlme,wmf="Machines_Leverage",Cairo=FALSE,printgraph=FALSE
)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }
