\name{plotdiag.deviance.residuals}
\alias{plotdiag.deviance.residuals}
\title{Plot Diagnostic Statistics Of Deviance Residuals}
\description{
Plot output from forsearch_glm to show change in deviance residuals or augmented 
deviance residuals, either of which can be squared, as the number of observations 
in the forward search procedure increases. Save plot in folder containing working 
directory. 
}
\usage{
plotdiag.deviance.residuals(forn, squared = FALSE, augmented=TRUE, hilos = c(1, 0), 
maintitle="Put main title here", subtitle="Put subtitle here", caption="Put caption here", 
wmf= "Put_graph_title_here", Cairo=TRUE,printgraph=TRUE,
legend = "Dummy legend name", subdiag = FALSE, subverb = FALSE, diagnose = FALSE, 
verbose = TRUE)
}
\arguments{
  \item{forn}{Name of forward search output file}
  \item{squared}{TRUE causes residuals to be squared before plotting}
  \item{augmented}{TRUE causes graphing of augmented deviance residuals, see Details}
  \item{hilos}{Number of observations having high and number having low values of 
  residuals to identify. No low values are identified for squared residual plot}
  \item{maintitle}{Main title of plot}
  \item{subtitle}{Subtitle of plot}
  \item{caption}{Caption of plot}
  \item{wmf}{File name of stored plot; omit ".wmf"}
  \item{Cairo}{TRUE causes use of Cairo graphics}
  \item{printgraph}{TRUE causes graph to print to file and
          closes device}
  \item{legend}{Legend title}
  \item{subdiag}{If TRUE, displays code to help diagnose subfunction errors}
  \item{subverb}{If TRUE, indicates beginning and end of function}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\details{
We reserve the use of the term 'Deviance residuals' to deviance residuals of the
observations that were used to create the model fit, and use the term 'Augmented
deviance residuals' to refer to deviance residuals of all available observations. The 
latter are created by predicting the fit of the model to all observations.
}
\value{
Process and plot changes in deviance residuals or squared deviance residuals from forsearch_glm}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
\testonly{
info3 <- system.file("extdata","micem1.for.R",package="forsearch");
info3 <- source(info3);
info3 <-info3[[1]];
plotdiag.deviance.residuals(info3,hilos=c(1,1),wmf="Micem1_for_DR",Cairo=FALSE,
printgraph=FALSE)
}
}
 \keyword{ hplot }
