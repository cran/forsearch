\name{plotdiag.phihatx}
\alias{plotdiag.phihatx}
\title{Plot Diagnostic PhiHat Statistics
}
\description{
Plot output from forsearch_glm to show change in phiHat statistics as the 
  number of observations in the forward search procedure increases. Save plot in folder 
  containing working directory. 
}
\usage{
plotdiag.phihatx(forn, maintitle = "Put main title here", 
subtitle = "Put subtitle here", caption="Put caption here", 
wmf = "Put_graph_filename_here", 
Cairo=TRUE, printgraph=TRUE, addline=c("none","loess","straight"),  
verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_glm}
  \item{maintitle}{Main title of plot}
  \item{subtitle}{Subtitle of plot}
  \item{caption}{Content of caption}
  \item{wmf}{File name of stored plot; omit ".wmf"}
  \item{Cairo}{TRUE causes use of Cairo graphics}
  \item{addline}{add a line to the graph; abbreviation allowed}
  \item{printgraph}{TRUE causes graph to print to file and closes device}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\value{Process and plot phiHat statistics from forsearch_glm}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
%\examples{
%\testonly{
%info3 <- system.file("extdata","micem1.for.R",package="forsearch");
%info3 <- source(info3);
%info3 <- info3[[1]];
%plotdiag.phihatx(info3,wmf="phiHat statistics", Cairo=FALSE, printgraph=FALSE,
%      addline="n")
%}
%}
 \keyword{ hplot }
