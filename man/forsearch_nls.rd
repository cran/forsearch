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
\usage{forsearch_nls(nlsform, data, start, algorithm="default", 
   nls.control=FALSE, initial.sample=1000, skip.step1=NULL, begin.diagnose=100, 
   verbose=TRUE)
}  
\arguments{
  \item{nlsform}{Formula for nls function}
  \item{data}{Name of database. First 2 variables are Observation and Section
           (both mandatory)}
  \item{start}{LIST of starting values for nls}
  \item{algorithm}{algorithm for nls function}
  \item{nls.control}{Logical. TRUE makes nls controls more liberal}
  \item{initial.sample}{Number of observation sets in Step 1 of forward search}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
      in Step 1}
  \item{begin.diagnose}{Numeric. Indicates where in code to begin printing 
      diagnostics. 0 prints all; 100 prints none}
  \item{verbose}{TRUE causes function identifier to display before and after run}
}
\value{LIST
  \item{Rows in stage }{Observation numbers of rows included at each stage}
  \item{Standardized residuals }{Matrix of errors at each stage}
  \item{Number of model parameters }{Same as number of levels of start
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
https://cran.r-project.org/web/packages/nlstools/vignettes/vignetteJSS.pdf
}
\author{William R. Fairweather}
\examples{
\dontrun{
%Circumference of prune trees (hypothetical)
Observation <- 1:70
Section <- rep(c(1,1,1,1,1,2,2,2,2,2), times=7)
Tree <- rep(1:7, each=10)
age <- rep(c(1,3,5,7,9,11,13,15,17,19), times = 7) * 100
circum <- forsearch::logist3(age, a=170, b=7, c=500) + rnorm(70)*0.1
test02 <- data.frame(Observation, Section, Tree, age, circum) 
startPru <- list(Asym=170, xmid=7, scal=500)
formulaPru  <- circum ~ I(logist3(x=age, a=Asym, b=xmid, c=scal))
forsearch_nls(nlsform = formulaPru,  
           data=test02, start=startPru, nls.control=TRUE,  
           initial.sample = 179, skip.step1=NULL, begin.diagnose=100, verbose=TRUE)
           
% 6-minute walk test
t <- (0:35)/3
VO2<- c(377.1111, 333.3333, 352.1429, 328.7500, 369.8750 ,394.4000, 352.6667, 
       337.3333, 366.4286, 364.0000, 293.8889, 387.0000, 364.8889, 342.2222, 
       400.3000, 375.1111 ,320.5556, 385.1667)
VO2<- c(VO2,527.0714,688.6364,890.8182,1145.1538, 1254.9091, 1327.5000,1463.9000,
       1487.8333 ,1586.6667, 1619.1000, 1494.4167 ,1640.4545, 1643.3750,
       1583.6364, 1610.8000 ,1568.5000, 1464.5833, 1652.8000) 
test01 <- data.frame(t,VO2)       
Observation <- 1:36
Section <- c(rep(1,20),rep(2,8),rep(3,8))
test01 <- cbind(Observation, Section, test01)
pstart <- list(VO2rest=400, VO2peak = 1600, mu = 1)
pformula <- as.formula(VO2~(t<=5.883)*(VO2rest)+          
            (t>5.883)*(VO2rest+(VO2peak-VO2rest)*
            (1-exp(-(t-5.883)*I(1/mu)))))
test2 <- forsearch_nls(nlsform=pformula, data = test01, 
           start=pstart, nls.control=FALSE, initial.sample = 300, skip.step1=NULL, 
           begin.diagnose=100, verbose=TRUE)

}
}
 \keyword{ datagen }
