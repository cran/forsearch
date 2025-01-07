\name{forsearch_cph}
\alias{forsearch_cph}
\title{Create Statistics of Forward Search in a Cox Proportional Hazard Database}
\description{Prepares summary statistics at each stage of forward search for 
          subsequent plotting.}
\usage{
forsearch_cph(alldata, formula.rhs, nofactform, initial.sample=1000, 
      skip.step1=NULL, ties = "efron", maxdisturb=.01, proportion=TRUE, 
      unblinded=TRUE, begin.diagnose= 100, verbose=TRUE)
}
\arguments{
  \item{alldata}{Data frame whose first 3 columns are Observation, event.time, 
      and status, and whose last columns are independent variables. Factor 
      variables must be defined in advance.}
  \item{formula.rhs}{The right hand side of a formula object (omit tilde '~')}
  \item{nofactform}{Formula, omitting all factor variables}
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
       in Step 1}
  \item{ties}{Method for handling ties in event time, "efron", "breslow" or 
        "exact"; see survival::coxph}
  \item{maxdisturb}{Amount of maximum disturbance to add to each event time to
        prevent ties}
  \item{proportion}{Logical. TRUE causes running of survival::cox.zph on each 
      stage}
  \item{unblinded}{TRUE causes printing of presumed analysis structure}
  \item{begin.diagnose}{Numeric. Indicates where in code to begin printing 
      diagnostics. 0 prints all; 100 prints none}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\details{Step 1 requires one additional observation to be included in database to 
     prevent failure of Wald statistic in Step 2 calculation of coxph function}
\value{LIST
   \item{Rows in stage }{Observation numbers of rows included at each stage}
   \item{Number of model parameters }{Number of fixed coefficients in Cox model}
   \item{Fixed parameter estimates }{Vector of parameter estimates at each stage}
   \item{Wald Test }{Vector of Wald tests at each stage}
   \item{LogLikelihood }{Vector of null and overall coefficients log likelihoods at
        each stage}
   \item{Likelihood ratio test }{Vector of LRTs at each stage}
   \item{Leverage }{Matrix of leverage of each observation at each stage}
   \item{Call }{Call to this function}
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.}
\author{William R. Fairweather}
\examples{
%\dontrun{
Observation <- 1:60
xx1 <- runif(60)
xx2 <- runif(60)
treat <- c(rep("A",30),rep("B",20),rep("C",10))
treat <- as.factor(treat)
event.time <- trunc(100*runif(60))
statusx <- c(rep(1,22),rep(0,8),rep(1,17),rep(0,3),rep(1,6),rep(0,4))
index3 <- sample(1:60,60)
xindep <- data.frame(Observation,event.time,status=statusx,xx1,xx2,treat)
formula.rhs <- "xx1 + xx2 + treat"
coxph.out05 <- forsearch_cph(alldata=xindep, formula.rhs=formula.rhs, 
    nofactform="xx1 + xx2", skip.step1=NULL, ties="efron", initial.sample=500, 
    begin.diagnose=100)
}
%}
\keyword{  datagen  }
