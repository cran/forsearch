\name{plotdiag.s2}
\alias{plotdiag.s2}
\title{Plot Diagnostic Statistics Of Residual Variation 
}
\description{
Plot output from forsearch_lm to show change in residual variation as the number of 
observations in the forward search procedure increases. Save plot in folder containing 
working directory. 
}
\usage{
plotdiag.s2(forn, maintitle = "Put main title here", subtitle = "Put subtitle here", 
caption = "Put caption here", wmf = "Put_graph_filename_here", 
Cairo=TRUE,printgraph=TRUE, addline = c("none","loess","straight"), 
verbose = TRUE)
}
\arguments{
  \item{forn}{
Name of output file from forsearch_lm
}
  \item{maintitle}{
Main title of plot
}
  \item{subtitle}{
Subtitle of plot
}
  \item{caption}{
Content of caption
}
  \item{wmf}{
File name of stored plot; omit ".wmf"  
}
  \item{Cairo}{TRUE causes use of Cairo graphics
}
  \item{printgraph}{TRUE causes graph to print to file and
          closes device
}
   \item{addline}{add a line to the graph; abbreviation allowed
}
  \item{verbose}{
If TRUE, indicates beginning and end of function
}
}
\value{
Process and plot residual variation statistics from forsearch_lm
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
%\examples{
%\testonly{
%info3 <- system.file("extdata","crossdata.for1.R",package="forsearch");
%crossdata.for1 <- source(info3);
%crossdata.for1 <-crossdata.for1[[1]];
%plotdiag.s2(forn=crossdata.for1, wmf="Crossdata_s2", Cairo=FALSE,
%     printgraph=FALSE, addline="n")
%}
%}
 \keyword{ hplot }
