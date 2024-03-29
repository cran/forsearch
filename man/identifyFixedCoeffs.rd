\name{identifyFixedCoeffs}
\alias{identifyFixedCoeffs}
\title{Index To Identify Fixed Coefficients To Appear Together on Plot
}
\description{
 Runs the defined linear (lm) model. Displays the resulting coefficients. 
    Attaches codes for identifying them to the plotting functions of this package.
}
\usage{identifyFixedCoeffs(formula, data, verbose = TRUE)
}
\arguments{
  \item{formula}{2-sided formula for fixed effects}
  \item{data}{Name of file (to be) run by forsearch_lm}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\details{
Plotting functions cannot plot more than a few coefficients on one graph. This function 
prepares an index of the coefficients so that the user can more easily identify which ones 
should appear together in a plot.
}
\value{Index of coefficients from forsearch_lm.
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
info3 <- system.file("extdata", "crossdata.R", package="forsearch");
crossdata <- source(info3);
crossdata <- crossdata[[1]];
identifyFixedCoeffs(formula=y~x1*x2, data=crossdata)
}
 \keyword{ manip }
