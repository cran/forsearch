\name{forsearch_lm}
\alias{forsearch_lm}
\title{Create Statistics Of Forward Search in a Linear Model Database
}
\description{Prepares summary statistics at each stage of forward search for subsequent plotting. 
   Forward search is conducted in two steps: Step 1 to identify minimal set of observations 
   to estimate unknown parameters, and Step 2 to add one observation at each stage such that 
   observations in the set are best fitting at that stage.
}
\usage{forsearch_lm(formula, data, initial.sample=1000, n.obs.per.level = 1,
                   skip.step1 = NULL, unblinded = TRUE, diagnose = FALSE,
                   verbose = TRUE)
}
\arguments{
  \item{formula}{Fixed effects formula as described in stats::lm}
  \item{data}{Name of database}
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{n.obs.per.level}{Number of observations per level of (possibly crossed) 
    factor levels}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
   in Step 1}
  \item{unblinded}{TRUE causes printing of presumed analysis structure}
  \item{diagnose}{TRUE causes printing of intermediate steps of function}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\details{Step 2 is determined by the results of Step 1, which itself is random. So,
    it is possible to reproduce the entire run by using the skip.step1 argument.}
\value{LIST
  \item{Rows in stage }{Observation numbers of rows included at each stage}
  \item{Standardized residuals }{Matrix of errors at each stage}
  \item{Number of model parameters }{Rank of model}
  \item{Sigma }{Estimate of random error at final stage; used to standardize all residuals}
  \item{Fixed parameter estimates }{Matrix of parameter estimates at each stage}
  \item{s^2 }{Estimate of random error at each stage}
  \item{Leverage }{Matrix of leverage of each observation at each stage}
  \item{Modified Cook distance }{Estimate of sum of squared changes in parameter estimates at each stage}
  \item{Call }{Call to this function}
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather}
\examples{
# Multiple regression
Observation <- 1:16
y <- runif(16)
x1 <- runif(16)
x2 <- runif(16)
x3 <- runif(16)
lmtest1 <- data.frame(Observation,y,x1,x2,x3)
forsearch_lm(formula=y~x1+x2+x3, data=lmtest1, initial.sample=100,diagnose=FALSE)
\dontrun{
# Analysis of variance 
Observation <- 1:30
y <- runif(30)
AN1 <- as.factor(c(rep("A1",5),rep("A2",5),rep("A3",5)))
AN1 <- c(AN1,AN1)
AN2 <- as.factor(c(rep("B1",15),rep("B2",15)))
lmtest2 <- data.frame(Observation,y,AN1,AN2)
forsearch_lm(formula=y~AN1*AN2, data=lmtest2, initial.sample=100,diagnose=TRUE)

# Analysis of covariance
Observation <- 1:60
y <- runif(60)
AN1 <- as.factor(c(rep("A1",10),rep("A2",10),rep("A3",10)))
AN1 <- c(AN1,AN1)
AN2 <- as.factor(c(rep("B1",30),rep("B2",30)))
COV <- runif(60)
lmtest3 <- data.frame(Observation,y,AN1,AN2,COV)
forsearch_lm(formula=y~AN1*AN2+COV, data=lmtest3, initial.sample=100,diagnose=FALSE)
}
}
\keyword{  datagen  }
