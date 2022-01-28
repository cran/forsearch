\name{forsearch_lm}
\alias{forsearch_lm}
\title{Create Statistics Of Forward Search in Database To Be Analyzed by LM Function 
}
\description{
Prepares summary statistics at each stage of forward search for subsequent plotting. Forward seach is conducted in two steps: Step 1 to identify minimal set of observations to estimate unknown parameters, and Step 2 to add one observation at each stage such that observations in the set are best fitting at that stage.
}
\usage{
forsearch_lm(formula, data, initial.sample, diagnose = FALSE, verbose = TRUE)
}
\arguments{
  \item{formula}{Fixed effects formula as described in stats::lm
}
  \item{data}{Name of database
}
  \item{initial.sample}{Number of observations in Step 1 of forward search
}
  \item{diagnose}{TRUE causes printing of intermediate steps of function
}
  \item{verbose}{TRuE causes function identifier to display before and after run}
}
\value{LIST
\item{Rows in stage}{Observation numbers of rows included at each stage}
\item{Standardized residuals}{Matrix of errors at each stage}
\item{Number of model parameters}{Rank of model}
\item{Sigma}{Estimate of random error at final stage; used to standardize all residuals}
\item{Fixed parameter estimates}{Matrix of parameter estimates at each stage}
\item{s^2}{Estimate of random error at each stage}
\item{Leverage}{Matrix of leverage of each observation at each stage}
\item{Modified Cook distance}{Estimate of sum of squared changes in parameter estimates at each stage}
\item{Call}{Call to this function}
}
\references{ Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather}

\examples{
%print("forsearch_lm");
info3 <- system.file("extdata","crossdata.R",package="forsearch");
%print(info3);
cross <- source(info3);
cross <- cross[[1]];
%print(cross);
crossdata.for1 <- forsearch_lm(formula=y~x1*x2, data=cross, initial.sample = 100)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }
