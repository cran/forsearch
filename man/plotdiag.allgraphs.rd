\name{plotdiag.allgraphs}
\alias{plotdiag.allgraphs}
\title{Execute All Plotting Functions For a Select Forsearch Object
}
\description{
Executes all the plotting functions for a select analytical function
   such as lm or glm; default omits titles and subtitles and attempts
   to plot all fixed and random coefficients.
}
\usage{
plotdiag.allgraphs(object, mt=" ", st=" ", cpt=" ", cc=NULL, ccrand = NULL,Cairo=TRUE)
}
\arguments{
  \item{object}{Name of forsearch object file}
  \item{mt}{Maintitle of graph}
  \item{st}{Subtitle of graph}
  \item{cpt}{Caption on the graph}
  \item{cc}{Fixed variable code numbers of coefficients to be included in graph}
  \item{ccrand}{Random variable code numbers of parameters to be included in 
      graph}
  \item{Cairo}{TRUE causes use of Cairo graphics}    
}
\value{
Prints search history and creates graphical files in current subdirectory
}
\author{William R. Fairweather
}
%\examples{
%\dontrun{
%info3 <- system.file("extdata", "train.for3.R", package="forsearch");
%info3 <- source(info3);
%info3 <- info3[[1]];
%plotdiag.allgraphs(object=info3, mt=" ", st=" ", cpt=" ", cc=NULL, ccrand = NULL
%   Cairo=FALSE)
%}
%}
 \keyword{ hplot }
