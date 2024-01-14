\name{forsearch_cph}
\alias{forsearch_cph}
\title{Create Statistics Of Forward Search in a Cox Proportional Hazard Database}
\description{Prepares summary statistics at each stage of forward search for subsequent
        plotting.}
\usage{forsearch_cph(alldata, formula.rhs, initial.sample=1000, 
        n.obs.per.level=1, skip.step1=NULL, ties = "efron", proportion=TRUE, 
        unblinded=TRUE, begin.diagnose= 100, verbose=TRUE)
}
\arguments{
  \item{alldata}{Data frame containing Observation number, event time, status, 
        and independent variables}
  \item{formula.rhs}{Character vector of names of independent variables in 
        model}
  \item{initial.sample}{Number of observations in Step 1 of forward search}
  \item{n.obs.per.level}{Number of observations per level of (possibly crossed) 
           factor levels to include in Step 1}
  \item{skip.step1}{NULL or a vector of integers for observations to be included
           in Step 1}
  \item{ties}{Method for handling ties in event time; = "efron", "breslow", or 
          "exact"; see survival::coxph}
  \item{proportion}{TRUE causes evaluation of proportionality of Cox regression
  }
  \item{unblinded}{TRUE causes printing of presumed analysis structure}
  \item{begin.diagnose}{Numeric. Indicates where in code to begin printing 
           diagnostics. 0 prints all; 100 prints none}
  \item{verbose}{TRUE causes function identifier display before and after run}
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
\dontrun{
{# Forsearch for Cox Proportional Survival, including Step 1
veteran <- survival::veteran
veteran <- veteran[order(veteran$celltype),]
veteranx <- veteran[,c(3,4,1,2)]
veteranx$trt <- as.factor(veteranx$trt)
dimv <- dim(veteran)[1]
Observation <- 1:dimv
veteranx <- data.frame(Observation,veteranx)
names(veteranx)[2] <- "event.time"
form.1 <- "trt + celltype"  
forskip <- NULL
# forskip <- c(12,  23,  38,  71,  91, 104, 116, 130,  31,  73,  62,  76)
cphtest1a.out <- forsearch_cph(alldata=veteranx, formula.rhs=form.1, 
       n.obs.per.level=2, skip.step1=forskip, ties="efron", unblinded=TRUE, 
       initial.sample=467, begin.diagnose = 100, verbose = TRUE)
}
{# Same, but skipping Step 1.
forskip <- c(12,  23,  38,  71,  91, 104, 116, 130,  31,  73,  62,  76)
cphtest1b.out <- forsearch_cph(alldata=veteranx, formula.rhs=form.1, 
      n.obs.per.level=2, skip.step1=forskip, ties="efron", unblinded=TRUE, 
      initial.sample=467, begin.diagnose = 100, verbose = TRUE) 
}
}
}
\keyword{  datagen  }
