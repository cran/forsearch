\name{plotdiag.params.random}
\alias{plotdiag.params.random}
\title{Plot Diagnostic Statistics Of Random Coefficients 
}
\description{
Plot output from forsearch_lme to show change in root mean squares of random coefficients 
as the number of observations in the forward search procedure increases. Save plot in 
folder containing working directory. 
}
\usage{
plotdiag.params.random(forn, coeff.codenums=NULL, 
asfacets=FALSE, facetdir=c("h","v"),
maintitle = "Put maintitle here", subtitle = "Put subtitle here", 
caption = "Put caption here",wmf = "Put_stored_name_here",
Cairo=TRUE,printgraph = TRUE,
legend = "Dummy legend name", diagnose = FALSE, verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_lme
}
  \item{coeff.codenums}{columns of output file to be included in graph}
  \item{asfacets}{TRUE causes printing in facets}
  \item{facetdir}{"v" lays out the facets vertically, "h" lays them out 
  horizontally}
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
\value{Process and plot RMS of random coefficients from forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
\testonly{
info3 <- system.file("extdata","Machines.O.forlme2.R",package="forsearch");
Machines.O.forlme2 <- source(info3);
Machines.O.forlme2 <- Machines.O.forlme2[[1]];
%print(Machines.O.forlme2$"Random parameter estimates")
plotdiag.params.random(Machines.O.forlme2,asfacets=FALSE, 
wmf="Machines_Random_Coefficients",Cairo=FALSE,printgraph=FALSE)
}
}
 \keyword{ hplot }
