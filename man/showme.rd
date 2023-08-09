\name{showme}
\alias{showme} 
\title{Display Abbreviated Output Of FORSEARCH_xxx Function 
}
\description{Output of forsearch_xxx function can be voluminous. This function displays the 
output in an abbreviated format. Primarily for programmer use.
}

\usage{showme(x, verbose = TRUE)
}
\arguments{
  \item{x}{Name of forsearch_xxx output file
}
  \item{verbose}{If TRUE, indicates the beginning and end of function run
}
}
\value{Abbreviated printout of output of forsearch_lm function
}
\author{William R. Fairweather 
}
\examples{
\testonly{
info3 <- system.file("extdata","crossdata.for1.R",package="forsearch");
crossdata.for1 <- source(info3);
crossdata.for1 <- crossdata.for1[[1]];
showme(crossdata.for1)
}
}
 \keyword{ manip }
