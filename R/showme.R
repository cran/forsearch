#' @export
showme <-
function(x, verbose=TRUE)
{
     #                          showme
     #
     # VALUE     Abbreviated output from forsearch_lm function.  Largely, a support for programming efforts
     #
     # INPUT    x            lm diagnostics object
     #
     #          verbose      Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running showme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     Hmisc::prn(names(x))
     print(search.history(x)[[1]])
     Hmisc::prn(utils::head(x$"Standardized residuals"))
     Hmisc::prn(utils::tail(x$"Standardized residuals"))
     Hmisc::prn(x$"Number of model parameters")
     Hmisc::prn(x$Sigma)
     Hmisc::prn(x$"Fixed parameter estimates")
     Hmisc::prn(x$"s^2")
     Hmisc::prn(utils::head(x$Leverage))
     Hmisc::prn(utils::tail(x$Leverage))
     Hmisc::prn(x$"Modified Cook distance")

#[1] "Rows in stage"              "Standardized residuals"    
#[3] "Number of model parameters" "Sigma"                     
#[5] "Fixed parameter estimates"  "s^2"                       
#[7] "Leverage"                   "Modified Cook distance"    
#[9] "Call"                      
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running showme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
