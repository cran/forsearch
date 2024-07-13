\name{forsearch_nls}
\alias{forsearch_nls}
\title{Create Statistics Of Forward Search in a Nonlinear Model Database
}
\description{Prepares summary statistics at each stage of forward search for 
   subsequent plotting. Forward search is conducted in two steps: Step 1 to 
   identify minimal set of observations to estimate unknown parameters, and 
   Step 2 to add one observation at each stage such that observations in the set
   are best fitting at that stage.
}
\usage{forsearch_nls(phaselist, data, poolstart, poolformula, algorithm=
  "default", controlarg=NULL, initial.sample=1000, skip.step1=NULL, 
  begin.diagnose=100, verbose=TRUE)
}      
\arguments{
  \item{phaselist}{LIST of formula, formulacont, start, nopl for each phase}
  \item{data}{Name of database. First 2 variables are Observation and Phases
           (both mandatory)}
  \item{poolstart}{List Start values for Step 2}
  \item{poolformula}{Formula for pooled data from all phases for Step 2}
  \item{algorithm}{algorithm for nls function.}
  \item{controlarg}{nls control. Default is NULL to use preset nls.control}
  \item{initial.sample}{Number of observation sets in Step 1 of forward search}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
      in Step 1}
  \item{begin.diagnose}{Numeric. Indicates where in code to begin printing 
      diagnostics. 0 prints all; 100 prints none}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\details{All datasets are considered to be in phases. See vignette for 
      definition and discussion. There is a phaselist for each phase and an 
      element for each phaselist input variable. In addition, there is a
      (pool)start and a (pool)formula input variable for the pooled dataset.
}
\value{LIST
  \item{Rows in stage }{Observation numbers of rows included at each stage}
  \item{Standardized residuals }{Matrix of errors at each stage}
  \item{Number of model parameters }{Same as number of levels of poolstart input
      variable}
  \item{Sigma }{Estimate of random error at final stage; used to standardize all
      residuals}
  \item{Fixed parameter estimates }{Vector of parameter estimates at each stage}
  \item{s^2 }{Estimate of random error at each stage}
  \item{Call }{Call to this function}
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, 
     New York, 2000.
Pinheiro, JC and DM Bates. Mixed Effects Models in S and S-PLUS, Springer, 
     New York, 2000.
Example from nlstools package     
}
\author{William R. Fairweather}
\examples{
\dontrun{
t<-(0:35)/3
VO2<-c(377.1111,333.3333,352.1429,328.7500,369.8750,394.4000,352.6667,337.3333,
  366.4286,364.0000,293.8889,387.0000,364.8889,342.2222,400.3000,375.1111,
  320.5556,385.1667,527.0714,688.6364,890.8182,1145.1538,1254.9091,1327.5000,
  1463.9000,1487.8333,1586.6667,1619.1000,1494.4167,1640.4545,1643.3750,
  1583.6364,1610.8000,1568.5000,1464.5833,1652.8000)
Observation <- 1:36
Phases <- as.factor(c(rep("REST",18), rep("EXERCISE",18)))
test01 <- data.frame(Observation,Phases,t,VO2)

formula.1 <-as.formula(VO2~VO2rest)
formulacont.1 <- as.formula(VO2~VO2rest)
start.1 <- list(VO2rest = 400)
nopl.1 <- 1

formula.2<-
  as.formula(VO2~(VO2rest+(VO2peak-VO2rest)*(1-exp(-(t-5.883)*I(1/mu)))))
formulacont.2<-
  as.formula(VO2~(VO2rest+(VO2peak-VO2rest)*(1-exp(-(t-5.883)*I(1/mu)))))
start.2 <- list(VO2rest = 400, VO2peak = 1600, mu = 1)
nopl.2 <- 6

phaselist <- list(
             REST=
 list(formula=formula.1,formulacont=formulacont.1,start=start.1,nopp=nopl.1),
             EXERCISE=
 list(formula=formula.2,formulacont=formulacont.2,start=start.2,nopp=nopl.2))

pstart <- list(VO2rest=400, VO2peak = 1600, mu = 1)
pformula <- as.formula(VO2~(t<=5.883)*(VO2rest)+          
            (t>5.883)*(VO2rest+(VO2peak-VO2rest)*
            (1-exp(-(t-5.883)*I(1/mu)))))
forsearch_nls(phaselist=phaselist, data=test01, 
   poolstart=pstart, poolformula=pformula, algorithm="default", 
   controlarg=nls.control(maxiter=50,warnOnly=TRUE), initial.sample = 155)
}
}
 \keyword{ datagen }
