\name{identifyCoeffs}
\alias{identifyCoeffs}
\title{Index To Identify Fixed and Random Coefficients To Appear Together on 
     Plot}
\description{
   Runs the defined, grouped linear mixed effects (lme) model. Displays the 
   resulting fixed and random coefficients. Attaches codes for identifying them 
   to the plotting functions of this package.}

\usage{identifyCoeffs(fixed, data, random, 
    XmaxIter = 1000, XmsMaxIter = 1000, 
    Xtolerance = 0.01, XniterEM = 1000, XmsMaxEval = 400, XmsTol = 1e-05, 
    Xopt = "optim",
    diagnose = FALSE, verbose = TRUE)
}
\arguments{
  \item{fixed}{2-sided formula for fixed effects}
  \item{data}{Name of file (to be) run by forsearch_lme}
  \item{random}{1-sided formula for random effects}
  \item{XmaxIter}{lme control parameter}
  \item{XmsMaxIter}{lme control parameter} 
  \item{Xtolerance}{lme control parameter}
  \item{XniterEM}{lme control parameter}
  \item{XmsMaxEval}{lme control parameter}
  \item{XmsTol}{lme control parameter}
  \item{Xopt}{lme control parameter}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\details{Plotting functions cannot plot more than a few coefficients on one graph. 
This function prepares an index of the coefficients so that the user can more easily 
identify which ones should appear together in a plot.
}
\value{Index of fixed and random coefficients from forsearch_lme.
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
\testonly{
info3 <- system.file("extdata","Machines.O.R",package="forsearch");
info3 <- source(info3);
info3 <- info3[[1]];
identifyCoeffs(fixed=score~1, data=info3, random= ~1 | Worker)
}
}
 \keyword{ manip }

