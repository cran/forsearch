\name{plotdiag.AICX}
\alias{plotdiag.AICX}
\title{Plot Diagnostic AIC Statistics
}
\description{
Plot output from forsearch_glm to show change in AIC statistics as the 
  number of observations in the forward search procedure increases. Save plot in folder 
  containing working directory. 
}
\usage{
plotdiag.AICX(forn, maintitle = "Put main title here", 
subtitle = "Put subtitle here", caption="Put caption title here", 
wmf = "Put_plot_file_title_here",
Cairo=TRUE, printgraph=TRUE,loess = FALSE, subdiag=FALSE, subverb=FALSE, 
diagnose = FALSE,verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_glm}
  \item{maintitle}{Main title of plot}
  \item{subtitle}{Subtitle of plot}
  \item{caption}{Content of caption}
  \item{wmf}{File name of stored plot; omit '.wmf'}
  \item{Cairo}{TRUE causes use of Cairo graphics}
  \item{printgraph}{TRUE causes graph to print to file and closes device}
  \item{loess}{TRUE causes plot of loess line, otherwise straight line}
  \item{subdiag}{If TRUE, displays code to help diagnose subfunction errors}
  \item{subverb}{If TRUE, indicates beginning and end of subfunction}
  \item{diagnose}{If TRUE, displays code to help diagnose main function errors}
  \item{verbose}{If TRUE, indicates beginning and end of function}
}
\value{Process and plot AIC statistics from forsearch_glm}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\examples{
info3 <- system.file("extdata","micem1.for.R",package="forsearch");
info3 <- source(info3);
info3 <- info3[[1]]
plotdiag.AICX(forn=info3,wmf="micem1_for_AIC", Cairo=FALSE, printgraph=FALSE)
}
 \keyword{ attribute }
 \keyword{ debugging }
 \keyword{ optimize }
 \concept{ outliers }