\name{forsearch-package}
\alias{forsearch-package}
\alias{forsearch}
\docType{package}
\title{Diagnostic Analysis Using Forward Search Procedure for Various Models
\packageTitle{forsearch}
}
\description{
\packageDescription{forsearch}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{forsearch}
\packageIndices{forsearch}
Ensure that data frame has a leading column of observation numbers. Run forsearch_xxx to 
create a file of diagnostic statistics to be used as input to such 
plotting functions as plotdiag.residuals, plotdiag.params.fixed, plotdiag.params.random, 
plotdiag.s2, plotdiag,leverage, and plotdiag.Cook. The file of diagnostic statistics can 
be voluminous, and the utility function showme displays the output more succinctly. 
Plotting of statistics for fixed and for random coefficients is limited by graphical 
restraints in some cases. The function identifyCoeffs provides a set of indexing codes 
so that plotdiag.params.random can display diagnostics for selected fixed or random model parameters. The function identifyFixedCoeffs does the same for lm models.
}

\author{William R. Fairweather, Flower Valley Consulting, Inc., Silver Spring MD USA
\packageAuthor{forsearch}

Maintainer: \packageMaintainer{forsearch}
William R. Fairweather <wrf343 at flowervalleyconsulting.com>
}
\references{
  Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
  Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.
}
\keyword{ package }
