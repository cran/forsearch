\name{plotdiag.Cook}
\alias{plotdiag.Cook}
\title{Plot Diagnostic Statistics of Modified Cook's Distance 
}
\description{
Plot output from forsearch_lm or forsearch_lme to show change in Modified Cook's distance 
as the number of observations in the forward search procedure increases. Save plot in 
folder containing working directory. 
}
\usage{
plotdiag.Cook(forn, maintitle = "Put main title here", subtitle = "Put subtitle here", 
caption = "Put caption here", wmf = "Put_plot_file_title_here", 
Cairo=TRUE, printgraph=TRUE, addline = "none", verbose = TRUE)
}
\arguments{
  \item{forn}{Name of forward search output file
}
  \item{maintitle}{Main title of plot
}
  \item{subtitle}{Subtitle of plot
}
  \item{caption}{Content of caption
}
  \item{wmf}{File name of stored plot; omit ".wmf"  
}
\item{Cairo}{TRUE causes use of Cairo graphics
}
\item{printgraph}{TRUE causes graph to print to file and
      closes device
}
   \item{addline}{Character variable to add a line to the graph; options: "none", "loess",
         and "straight"; abbreviation allowed
}
  \item{verbose}{If TRUE, indicates beginning and end of function
}
}
\value{
Process and plot Cook distance statistics from forsearch_lm or forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
%\examples{
%\testonly{
%info3 <- system.file("extdata","Alfalfa.O.forlme.R",package="forsearch");
%Alfalfa.O.forlme <- source(info3);
%Alfalfa.O.forlme <- Alfalfa.O.forlme[[1]];
%plotdiag.Cook(Alfalfa.O.forlme, wmf="Alfalfa_Cook", Cairo=FALSE,
%printgraph=FALSE,addline="n")
%}
%}
 \keyword{ hplot }
