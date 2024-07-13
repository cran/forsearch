\name{forsearch_glm}
\alias{forsearch_glm}
\title{Create Statistics of Forward Search in a Generalized Linear Model Database}
\description{
Prepares summary statistics at each stage of forward search for subsequent plotting. 
Forward search is conducted in three steps: Step 1 to identify minimal set of observations 
to estimate unknown parameters, and Step 2 to add one observation at each stage such that 
observations in the set are best fitting at that stage. A preliminary step (Step 0) 
contains code for pre-processing of the data.
}
\usage{
forsearch_glm(initial.sample=1000, response.cols, indep.cols, family,  
   formula=NULL, binomialrhs=NULL, formula.cont.rhs, data, n.obs.per.level = 1,
   estimate.phi = TRUE, skip.step1=NULL, unblinded=TRUE, begin.diagnose=100, 
   verbose=TRUE)
}
\arguments{
  \item{initial.sample}{Number of random sets of observations in Step 1 of forward 
  search}
  \item{response.cols}{Vector of column numbers (1 or 2) of responses and nonresponses (if binomial)}
  \item{indep.cols}{Column number(s) of independent variables}
  \item{family}{Error distribution and link}
  \item{formula}{Formula relating response to independent variables. Required
       except for family=binomial}
  \item{binomialrhs}{Quoted character.Right-hand side of formula. Required for 
       family=binomial}
  \item{formula.cont.rhs}{Quoted character.Right-hand side of formula, omitting
            factor variables. Required for all families}       
  \item{data}{Name of database}
  \item{n.obs.per.level}{Number of observations per level of (possibly crossed) 
           factor levels}
  \item{estimate.phi}{TRUE causes phi to be estimated; FALSE causes phi to be set = 1}
  \item{skip.step1}{NULL, or vector of observation numbers to include at end of Step 1}
  \item{unblinded}{TRUE allows print of formula of analysis function}
  \item{begin.diagnose}{Numeric. Indicates where in code to begin printing diagnostics.
                0 prints all; 100 prints none}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\value{LIST
  \item{Rows in stage}{Observation numbers of rows included at each stage}
  \item{Family}{Family and link}
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
\details{Step 2 is determined by the results of Step 1, which itself is random. 
     So, it is possible to reproduce the entire run by using the skip.step1 
     argument. Inner subgroups are produced by presence of categorical
     variables. Current version assumes independent variables are all continuous.
}
\references{ Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, 
     Springer, New York, 2000.
}
\author{William R. Fairweather}
\examples{
# Train deaths (Atkinson and Riani, 2000) with Rolling Stock as a factor
Observation<-1:67
Month<-c(9,8,3,1,10,6,7,1,8,4,3,3,12,11,10,9,9,4,3,12,12,10,7,2,12,2,12,12,12,
    11,3,10,4,2,12,12,9,11,1,10,8,6,1,10,6,12,8,4,9,6,12,10,7,2,5,12,5,5,4,3,1,
    9,11,9,7,3,2)
Year<-c(97,96,96,95,94,94,91,91,90,89,89,89,88,88,87,86,86,86,86,84,84,84,84,84,
    83,83,82,81,81,80,80,79,79,79,78,78,77,76,76,75,75,75,75,74,74,73,73,73,72,
    72,71,71,71,71,70,69,69,69,69,69,69,68,67,67,67,67,67)
RollingStock<-c(2,2,3,2,1,1,1,1,2,3,1,1,1,2,1,2,1,3,2,2,1,2,2,3,1,2,1,1,2,3,1,
    1,1,1,1,1,1,3,3,2,3,1,2,3,1,1,1,3,3,1,3,3,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1)
RollingStock <- as.factor(RollingStock)    
Traffic<-c(0.436,0.424,0.424,0.426,0.419,0.419,0.439,0.439,0.431,0.436,0.436,
    0.436,0.443,0.443,0.397,0.414,0.414,0.414,0.414,0.389,0.389,0.389,0.389,
    0.389,0.401,0.401,0.372,0.417,0.417,0.43,0.43,0.426,0.426,0.426,0.43,0.43,
    0.425,0.426,0.426,0.436,0.436,0.436,0.436,0.452,0.452,0.433,0.433,0.433,
    0.431,0.431,0.444,0.444,0.444,0.444,0.452,0.447,0.447,0.447,0.447,0.447,
    0.447,0.449,0.459,0.459,0.459,0.459,0.459)
Deaths<-c(7,1,1,1,5,2,4,2,1,1,2,5,35,1,4,1,2,1,1,3,1,3,13,2,1,1,1,4,1,2,1,5,7,
    1,1,3,2,1,2,1,2,6,1,1,1,10,5,1,1,6,3,1,2,1,2,1,1,6,2,2,4,2,49,1,7,5,9)
train2022 <- data.frame(Observation, Year, RollingStock, Traffic, Deaths)
forsearch_glm(initial.sample = 100, response.cols = 5, 
    indep.cols = 2:4, formula=Deaths~Year + RollingStock + Traffic,
    formula.cont.rhs="Year + Traffic", 
    family = poisson("log"), data = train2022, 
    n.obs.per.level = 1, estimate.phi = TRUE, skip.step1 = NULL, 
    unblinded = TRUE, begin.diagnose=100)
}
\keyword{ datagen }
