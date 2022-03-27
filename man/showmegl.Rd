\name{showmegl}
\alias{showmegl}
\title{Display Abbreviated Output Of FORSEARCH_GLM Function 
}
\description{Output of forsearch_glm function can be voluminous. This function displays 
the output in an abbreviated format. Primarily for programmer use.
}
\usage{showmegl(x, verbose = TRUE)
}
\arguments{
  \item{x}{Name of forsearch_glm output file
}
  \item{verbose}{
If TRUE, indicates the beginning and end of function run
}
}
\value{Abbreviated printout of output of forsearch_glm function
}
\author{William R. Fairweather}
\examples{
\testonly{
info3 <- system.file("extdata","micem1.for.R",package="forsearch");
invo4 <- source(info3);
invo4 <-invo4[[1]];
showmegl(invo4)
}
}
 \keyword{ manip }
