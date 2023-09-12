\name{forsearch-package}
\alias{forsearch-package}
\alias{forsearch}
\docType{package}
\title{Diagnostic Analysis Using Forward Search Procedure for Various Models}
\description{
The DESCRIPTION file:
%\packageDESCRIPTION{forsearch}
%\packageIndices{forsearch}
%{
Assume that you have a database and that you intend to analyze it with the foo function. The
forsearch_xxx functions of this package do not alter the data, but only reorganize the 
database to identify which observations best fit with the assumed foo model and which ones worst
fit that model. 

Add a leading column of observation numbers to the database. Run foo.out <- forsearch_foo 
to create a file of diagnostic statistics. Run the relevant plotting functions, such as 
plotdiag_residuals, plotdiag_params.fixed, and plotdiag_params.random on foo.out. 
Also run search.history to see the sequence of reorganization datasets. 
Plotting of statistics for fixed and for random coefficients is limited by graphical restraints 
in some cases. The function identifyCoeffs provides a set of indexing codes so that 
plotdiag_params.random can display diagnostics for selected fixed or random model parameters. 
The function identifyFixedCoeffs does the same for lm models.

To perform these evaluations blindly with regard to treatment or another variable, refer to the 
package blindreview.
}
\author{William R. Fairweather}
\references{
  Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
  Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.
  https://CRAN.R-project.org/package=blindreview 
}
\keyword{ package }
