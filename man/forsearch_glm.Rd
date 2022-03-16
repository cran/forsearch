\name{forsearch_glm}
\alias{forsearch_glm}
\title{Create Statistics Of Forward Search in Generalized Linear Model Database}
\description{
Prepares summary statistics at each stage of forward search for subsequent plotting. 
Forward search is conducted in two steps: Step 1 to identify minimal set of observations 
to estimate unknown parameters, and Step 2 to add one observation at each stage such that 
observations in the set are best fitting at that stage.
}
\usage{
forsearch_glm(initial.sample, cobs, response.cols, indep.cols, family, data, 
estimate.phi=TRUE, skip.step1=NULL, diagnose=FALSE, verbose=TRUE)
}
\arguments{
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{cobs}{Number of observations to include in Step 1}
  \item{response.cols}{Column number(s) of response(s)}
  \item{indep.cols}{Column number(s) of independent variables}
  \item{family}{Error distribution and link}
  \item{data}{Name of database}
  \item{estimate.phi}{TRUE causes phi to be estimated; FALSE causes phi to be set = 1}
  \item{skip.step1}{NULL or vector of observation numbers to include at end of Step 1}
  \item{diagnose}{TRUE causes printing of intermediate steps of function}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\value{LIST
  \item{Rows in stage}{Observation numbers of rows included at each stage}
  \item{family}{Family and link}
  \item{Number of model parameters}{Number of fixed effect parameters}
  \item{Fixed parameter estimates}{Matrix of parameter estimates at each stage}  
  \item{Residual deviance}{Vector of deviances}
  \item{Null deviance}{Vector of null deviances}
  \item{PhiHat}{Vector of values of phi parameter}
  \item{Deviance residuals and augments}{Deviance residuals with indication of 
whether each is included in fit}
  \item{AIC}{Vector of AIC values}
  \item{Leverage}{Matrix of leverage of each observation at each stage}
  \item{Call}{Call to this function}
}
\references{ Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather}
\examples{
print("forsearch_glm");
info3 <- system.file("extdata","micem1.R", package="forsearch");
%print(info3);
cross <- source(info3);
cross <- cross[[1]];
%print(cross);
forsearch_glm(initial.sample=5,cobs=5,response.cols=c(4,8),indep.cols=c(9,3), 
   family=binomial("logit"), data=cross, skip.step1=NULL)
}   
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }