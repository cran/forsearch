\name{plotdiag.fit3}
\alias{plotdiag.fit3}
\title{
Plot Diagnostic Statistics of AIC, BIC, and Log Likelihood
}
\description{
Plot output from forsearch_lme to show change in AIC, BIC, and log likelihood as the 
number of observations in the forward search procedure increases. Save plot in folder 
containing working directory. 
}
\usage{
plotdiag.fit3(forn, maintitle = "Put main title here", subtitle = "Put subtitle here", 
    caption = "Put caption here", wmf = "Put_stored_name_here", 
    Cairo=TRUE,printgraph=TRUE, legend="Dummy legend name",
    diagnose = FALSE, verbose = TRUE)
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
  \item{legend}{Legend name
}
  \item{diagnose}{
If TRUE, displays code to help diagnose main function errors
}
  \item{verbose}{
If TRUE, indicates beginning and end of function
}
}
\value{
Process and plot trends of AIC, BIC, and log likelihood statistics from 
forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
\testonly{
info3 <- system.file("extdata","Alfalfa.O.forlme2.R",package="forsearch");
crossdata.for1 <- source(info3);
crossdata.for1 <-crossdata.for1[[1]];
plotdiag.fit3(forn=crossdata.for1, wmf="Alfalfadata_fit3",Cairo=FALSE,
    printgraph=FALSE)
}
}
 \keyword{ hplot }
