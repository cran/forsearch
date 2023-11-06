\name{forsearch_cph}
\alias{forsearch_cph}
\title{Create Statistics Of Forward Search in a Cox Proportional Hazard Database}
\description{Prepares summary statistics at each stage of forward search for subsequent
        plotting.}
\usage{
forsearch_cph(formula.elements, event.time, status, x, initial.sample=1000, 
n.obs.per.level=1, skip.step1=NULL, ties = c("efron", "breslow", "exact"), unblinded=TRUE, 
begin.diagnose=100, verbose=TRUE)
}
\arguments{
  \item{formula.elements}{Character vector of names of independent variables in model}
  \item{event.time}{Vector of event times, censored or not}
  \item{status}{Vector indicator of event or censoring: 1 = event, 0 = censored, 
        same length as event.time}
  \item{x}{Data frame of independent variables, with number of rows = length of event.time,
       First column is Observation. Factor variables must be defined in advance. Does not
       include event time or status.}
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{n.obs.per.level}{Number of observations per level of (possibly crossed) 
           factor levels to include in Step 1}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
           in Step 1}
  \item{ties}{Method for handling ties in event time; = "efron", "breslow", or 
          "exact"; see survival::coxph}
  \item{unblinded}{TRUE causes printing of presumed analysis structure}
  \item{begin.diagnose}{Numeric. Indicates where in code to begin printing diagnostics.
                0 prints all; 100 prints none}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\details{
Step 1 of a forward search function endeavors to find the smallest, most central 
subset of observations that is capable of estimating the fixed parameters of the 
model. For example, for a linear regression model with r continuous independent 
variables (excluding the intercept), subsets of size r+1 would be considered. 
This allows estimation of the fixed parameters, but not the random one. The 
underlying analytical function is then run repeatedly on each subset of 
observations in order to identify the next larger subset. This general approach 
must be modified if the underlying analytical function presents some problems. 
These problems are not encountered when the sample size is of the usual magnitude. 
In the case of the Cox proportional hazards analysis function, coxph, every call 
to this function involves calculation of the Wald analysis, and this analysis fails 
if it cannot calculate the variation. It is not currently an option not to run this 
analysis. Instead, we increase the number in the Step 1 subsample to ensure that the 
Wald analysis will not fail.

It is also required that a cubic spline is fit to any sample when coxph is run. This 
can also fail (for different reasons) in small samples. To prevent this, the user can 
increase the value of n.obs.per.level from the default (2) to some larger value.
}
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
Observation <- 1:60       
xx1 <- runif(60)
xx2 <- runif(60)
treat <- c(rep("A",30),rep("B",20),rep("C",10))
treat <- as.factor(treat)
event.tim <- trunc(100*runif(60))
statusx <- c(rep(1,22),rep(0,8),rep(1,17),rep(0,3),rep(1,6),rep(0,4))
index3 <- sample(1:60,60)
xindep <- data.frame(xx1,xx2,treat,event.tim, statusx)
xindep <- xindep[index3,]
xindep <- xindep[order(xindep[,5],xindep[,3]),]
xindep <- data.frame(Observation,xindep)
event.tim <- xindep[,5]
statusx <- xindep[,6]
xindep <- xindep[,1:4]
elems <- c("xx1","xx2","treat")
## Not run: 
\dontrun{
coxph.out22 <- forsearch_cph(formula.elements=elems, event.time=event.tim, status=statusx, 
    x=xindep, skip.step1=NULL, ties="efron", initial.sample=500, n.obs.per.level=4,
    begin.diagnose=100)
}
coxph.out22 <- forsearch_cph(formula.elements=elems, event.time=event.tim, status=statusx, 
    x=xindep, skip.step1=c(12,36,42,11,31,45,17,29,44,18,30,41), ties="efron", initial.sample=500, 
    n.obs.per.level=4, begin.diagnose=100)
}
\keyword{  datagen  }
