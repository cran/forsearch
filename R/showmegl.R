#' @export
showmegl <-
function(x, verbose=TRUE)
{
     #                          showmegl
     #
     # VALUE     Abbreviated output from forsearch_glm function.  Largely, a support for programming efforts
     #
     # INPUT    x            glm diagnostics object
     #
     #          verbose      Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running showmegl", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     Hmisc::prn(names(x))
     print(search.history(x)[[1]])
     print(x$Family)
     Hmisc::prn(x$"Number of model parameters")

     Hmisc::prn(utils::head(x$"Fixed parameter estimates"))
     Hmisc::prn(utils::tail(x$"Fixed parameter estimates"))

#     Hmisc::prn(utils::head(x$"Studentized deviance residuals"))
#     Hmisc::prn(utils::tail(x$"Studentized deviance residuals"))

     Hmisc::prn(utils::head(x$"Deviance residuals and augments"))
     Hmisc::prn(utils::tail(x$"Deviance residuals and augments"))

     Hmisc::prn(x$"Residual deviance")
     Hmisc::prn(x$"Null deviance")

     Hmisc::prn(x$PhiHat)
     Hmisc::prn(x$AIC)

     Hmisc::prn(utils::head(x$Leverage))
     Hmisc::prn(utils::tail(x$Leverage))

#    Hmisc::prn(utils::head(x$"Modified Cook distance"))
#    Hmisc::prn(utils::tail(x$"Modified Cook distance"))

     Hmisc::prn(utils::head(x$"t statistics"))
     Hmisc::prn(utils::tail(x$"t statistics"))

     Hmisc::prn(x$Call)

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running showmegl", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
