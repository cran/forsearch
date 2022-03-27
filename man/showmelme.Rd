\name{showmelme}
\alias{showmelme} 
\title{Display Abbreviated Output Of FORSEARCH_LME Function 
}
\description{Output of forsearch_lme function can be voluminous. This function displays 
the output in an abbreviated format. Primarily for programmer use.
}

\usage{showmelme(x, verbose = TRUE)
}
\arguments{
  \item{x}{Name of forsearch_lme output file
}
  \item{verbose}{
If TRUE, indicates the beginning and end of function run
}
}
\value{Abbreviated printout of output of forsearch_lme function
}
\author{William R. Fairweather}
\examples{
\testonly{
info3 <- system.file("extdata","MOf.R",package="forsearch");
Machines.O.forlme <- source(info3);
Machines.O.forlme <-Machines.O.forlme[[1]];
showmelme(Machines.O.forlme)
}
}
 \keyword{ manip }
