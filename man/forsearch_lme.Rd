\name{forsearch_lme}
\alias{forsearch_lme}
\title{Create Statistics Of Forward Search In Database To Be Analyzed by LME Function 
}
\description{
Prepares summary statistics at each stage of forward search for subsequent plotting. Forward seach is conducted in three steps: Step 0 to set up accounting for group structure, Step 1 to identify minimal set of observations to estimate unknown parameters, and Step 2 to add one observation at each stage such that observations in the set are best fitting at that stage.
}
\usage{
forsearch_lme(fixed,  data, random, formula,  response.column,  initial.sample, robs, 
    skip.step1=NULL, XmaxIter = 1000, XmsMaxIter = 1000, 
    Xtolerance = 0.01, XniterEM = 1000, XmsMaxEval = 400, XmsTol = 1e-05, 
    Xopt = "optim", diagnose = FALSE, verbose = TRUE) 
}
\arguments{
  \item{fixed}{2-sided formula for fixed effects}
  \item{data}{data frame, first column of which must be "Observation"}
  \item{random}{1-sided formula for random effects}
  \item{formula}{a formula of the form resp ~ cov | group where resp is the response, cov is the primary covariate, and group is the grouping factor}
  \item{response.column}{Column number of response variable}
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{robs}{Number of observations to include in Step 1 of forward search from each subgroup}
  \item{skip.step1}{NULL or a vector of integers for rows to be included in Step 1}
  \item{XmaxIter}{lme control parameter}
  \item{XmsMaxIter}{lme control parameter}
  \item{Xtolerance}{lme control parameter}
  \item{XniterEM}{lme control parameter}
  \item{XmsMaxEval}{lme control parameter}
  \item{XmsTol}{lme control parameter}
  \item{Xopt}{lme control parameter}
  \item{diagnose}{TRUE causes printing of intermediate steps of function}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\details{
Group structure is ignored in calculating errors of fit. That is, predictions derive from lm fit and not lme fit. Diagnostic statistics are obtained from lme fits
}
\value{LIST
  \item{Rows in stage}{Observation numbers of rows included at each stage}
  \item{Standardized residuals }{Matrix of errors at each stage}
  \item{Number of model parameters}{Rank of model}
  \item{Sigma}{Estimate of random error at final stage; used to standardize all residuals}
  \item{Fixed parameter estimates}{Matrix of parameter estimates at each stage}
  \item{s^2}{Estimate of random error at each stage}
  \item{Leverage}{Matrix of leverage of each observation at each stage}
  \item{Modified Cook distance}{Estimate of sum of squared changes in parameter estimates at each stage}
  \item{Call}{Call to this function}
}
\references{ Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
             Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.
             https://CRAN.R-project.org/package=nlme
}
\author{William R. Fairweather
}
\examples{
%print("forsearch_lme");
info3 <- system.file("extdata","Alfalfa.O.R",package="forsearch");
%print(info3);
Alfalfa.O <- source(info3);
Alfalfa.O <- Alfalfa.O[[1]];
%print(Alfalfa.O);

firstrun <- (1:72)[c(-3,-47)]
Alfalfa.O.forlme<-forsearch_lme(
fixed=Yield~Date,
data=Alfalfa.O,
random= ~1 | Block/Variety,
formula=Yield~1|Block/Variety,
response.column=5,
initial.sample=4,robs=2,skip.step1=firstrun
)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }
