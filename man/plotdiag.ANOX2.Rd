\name{plotdiag.ANOX2}
\alias{plotdiag.ANOX2}
\title{Plot Diagnostic Statistics of Analysis of Variance Tables
}
\description{
Plot output from forsearch_xxx to show change in anova p-values as the number of 
observations in the forward search procedure increases. Save plot in folder containing 
working directory. 
}
\usage{
plotdiag.ANOX2(forn, anova.rows=NULL, ylab.extend=c("proportionality","variance"), 
maintitle = "Put main title here", 
subtitle = "Put subtitle here", caption="Put caption here",wmf = "Put_stored_name_here", 
Cairo=TRUE, printgraph=TRUE,legend = "Dummy legend name",  
verbose = TRUE)
}
\arguments{
  \item{forn}{Name of output file from forsearch_xxx
}
  \item{anova.rows}{Row numbers of p values to include together on the plot 
}
  \item{ylab.extend}{Type of anova table. "proportionality" is a test of proportionality
  for a coxph analysis; "variance" is a test of null hypothesis of a lm or lme test
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
  \item{legend}{Name of legend
}
  \item{verbose}{If TRUE, indicates beginning and end of function
}
}
\value{Process and plot anova test p values from forsearch_lm or forsearch_lme
}
\references{
Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
}
\author{William R. Fairweather
}
\keyword{ hplot }
