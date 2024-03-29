\name{forsearch_lme}
\alias{forsearch_lme}
\title{Create Statistics Of Forward Search For a Linear Mixed Effects Database
}
\description{Prepares summary statistics at each stage of forward search for subsequent
   plotting. Forward search is conducted in four steps: Step 0 to set up accounting for 
   group structure, Step 1 to identify minimal set of observations to estimate unknown 
   fixed parameters, Step 2 to identify the order of the remaining observations, and a
   final stage to extract the intermediate statistics based on increasing sample size.
}
\usage{
forsearch_lme(fixedform, alldata, randomform, initial.sample=1000, n.obs.per.level=1, 
   skip.step1=NULL, unblinded=TRUE, begin.diagnose = 100, verbose = TRUE)
}
\arguments{
  \item{fixedform}{2-sided formula for fixed effects}
  \item{alldata}{data frame, first column of which must be "Observation"}
  \item{randomform}{1-sided formula for random effects}
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{n.obs.per.level}{Number of observations per level of (possibly crossed) 
         factor levels}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
         in Step 1}
  \item{unblinded}{TRUE causes printing of presumed analysis structure}
  \item{begin.diagnose}{Numeric indicator of place in coding to begin printing 
         diagnostic information. 0 prints all information, 100 prints none.}
  \item{verbose}{TRUE causes function identifier to display before and after 
         run}
}
\details{data will be grouped within the function, regardless of initial 
    layout. Step 2 is determined by the results of Step 1, which itself is random.
    So, it is possible to reproduce the entire run by using the skip.step1 argument.
    Variables in the randomform formula must be character variables, but *not* factors}
\value{LIST
   \item{Number of observations in Step 1 }{Number of observations included in Step 1}
   \item{Step 1 observation numbers}{Observation numbers useful in skipping step 1}
   \item{Rows by outer subgroup}{List of row numbers, by outer subgroup}
   \item{Rows by outer-inner subgroups}{List of row numbers, by outer-inner subgroup}
   \item{Rows in stage }{Observation numbers of rows included at each stage}
   \item{Sigma }{Estimate of random error at final stage; used to standardize all residuals}
   \item{Standardized residuals}{Matrix of errors at each stage}
   \item{Fixed parameter estimates}{Matrix of parameter estimates at each stage}
   \item{Random parameter estimates}{Matrix of parameter estimates at each stage}
   \item{Leverage }{Matrix of leverage of each observation at each stage}
   \item{Modified Cook distance}{Estimate of sum of squared changes in parameter estimates at each stage}
   \item{Dims }{Dims from fit of lme function}
   \item{t statistics }{t statistics for each fixed parameter}
   \item{Fit statistics }{AIC, BIC, and log likelihood}
   \item{Call }{Call to this function}
}
\references{
{Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.}
{Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.}
{https://CRAN.R-project.org/package=nlme}
}
\author{William R. Fairweather
}
\examples{
# Multiple regression in grouped data
Observation <- 1:16
y <- runif(16)
x1 <- runif(16)
x2 <- runif(16)
x3 <- runif(16)
group <- as.factor(rep(c("G1","G2"),each=8))
lmetest1 <- data.frame(Observation,y,x1,x2,x3,group)
%lmetest1 <- nlme::groupedData(formula = y ~ 1 | group, data = lmetest1)
forsearch_lme(fixedform=y~x1+x2+x3, alldata=lmetest1, randomform= ~1|group, 
   n.obs.per.level=1, initial.sample=200)
\dontrun{

# Analysis of variance in grouped data
Observation <- 1:60
y <- runif(60)
AN1 <- as.factor(c(rep("A1",5),rep("A2",5),rep("A3",5)))
AN1 <- c(AN1,AN1,AN1,AN1)
AN2 <- as.factor(c(rep("B1",15),rep("B2",15)))
AN2 <- c(AN2,AN2)
group <- as.factor(rep(c("G1","G2"),each=30))
lmetest2 <- data.frame(Observation,y,AN1,AN2,group)
%lmetest2 <- nlme::groupedData(formula = y ~ 1|group, alldata = lmetest2)
forsearch_lme(fixedform=y~AN1*AN2, alldata=lmetest2, randomform= ~1|group,
             initial.sample=500)

# Analysis of covariance in grouped data

Observation <- 1:120
y <- runif(120)
AN1 <- as.factor(c(rep("A1",10),rep("A2",10),rep("A3",10)))
AN1 <- c(AN1,AN1,AN1,AN1)
AN2 <- as.factor(c(rep("B1",10),rep("B2",10)))
AN2 <- c(AN2,AN2,AN2,AN2,AN2,AN2)
COV <- runif(120)
group <- as.factor(rep(c("G1","G2"),each=30))
group <- c(group,group)
lmetest3 <- data.frame(Observation,y,AN1,AN2,COV,group)
%lmetest3 <- nlme::groupedData(formula = y ~ 1 | group, data = lmetest3)
forsearch_lme(fixedform=y~AN1*AN2+COV, alldata=lmetest3, randomform= ~ 1 | group,
        initial.sample=500)
}
}
\keyword{ datagen }
