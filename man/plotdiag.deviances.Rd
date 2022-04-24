\name{plotdiag.deviances}
\alias{plotdiag.deviances}
\title{Plot Diagnostic Deviance Statistics
}
\description{
Plot output from forsearch_glm to show change in deviances as the 
  number of observations in the forward search procedure increases. Save plot in folder 
  containing working directory. 
}
\usage{
plotdiag.deviances(forn, devtype, maintitle = "Put main title here", 
subtitle = "Put subtitle here", caption="Put caption here", 
wmf = "Put_plot_file_title_here", 
Cairo=TRUE, printgraph=TRUE,loess=FALSE, 
diagnose = FALSE,verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_glm}
  \item{devtype}{Type of deviance: "R" or "N" for Residual deviance or Null deviance}
  \item{maintitle}{Main title of plot}
  \item{subtitle}{Subtitle of plot}
  \item{caption}{Content of caption}
  \item{wmf}{File name of stored plot; omit ".wmf"}
  \item{Cairo}{TRUE causes use of Cairo graphics}
  \item{printgraph}{TRUE causes graph to print to file and closes device}
  \item{loess}{If TRUE, loess line is drawn through points, otherwise straight line}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\value{Process and plot deviances from forsearch_glm}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
\testonly{
info3 <- system.file("extdata","micem1.for.R",package="forsearch");
info3 <- source(info3);
info3 <- info3[[1]];
plotdiag.deviances(info3, devtype="R", wmf="Deviance statistics: Residual", 
Cairo=FALSE, printgraph=FALSE)
}
}
 \keyword{ hplot }
