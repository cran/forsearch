#' @export
showme <-
function(x, verbose=TRUE)
{
     #                          showme
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
     Hmisc::prn(x$"t statistics")
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running showme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return()
}
