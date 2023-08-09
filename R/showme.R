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
     uu <- as.character(x$Call[[1]])

     if(substr(uu,1,13)=="forsearch_cph"){
          print("This is a forsearch_cph output file", quote=FALSE)
          print("", quote = FALSE)
          Hmisc::prn(names(x))
          print(search.history(x)[[1]])
          Hmisc::prn(x$"Number of model parameters")

          Hmisc::prn(x$"Fixed parameter estimates")

          Hmisc::prn(utils::head(x$Leverage))
          Hmisc::prn(utils::tail(x$Leverage))

          Hmisc::prn(x$"Wald Test")
          Hmisc::prn(x$LogLikelihood)    
          Hmisc::prn(x$"Likelihood ratio test")
     }        # cph
     # 
     if(substr(uu,1,13)=="forsearch_glm"){
          print("This is a forsearch_glm output file", quote=FALSE)
          print("", quote = FALSE)
          Hmisc::prn(names(x))
          print(search.history(x)[[1]])
          print(x$Family)
          Hmisc::prn(x$"Number of model parameters")

          Hmisc::prn(utils::head(x$"Fixed parameter estimates"))
          Hmisc::prn(utils::tail(x$"Fixed parameter estimates"))

          Hmisc::prn(utils::head(x$"Deviance residuals and augments"))
          Hmisc::prn(utils::tail(x$"Deviance residuals and augments"))

          Hmisc::prn(x$"Residual deviance")
          Hmisc::prn(x$"Null deviance")

          Hmisc::prn(x$PhiHat)
          Hmisc::prn(x$AIC)

          Hmisc::prn(utils::head(x$Leverage))
          Hmisc::prn(utils::tail(x$Leverage))

#          Hmisc::prn(utils::head(x$"Modified Cook distance"))
#          Hmisc::prn(utils::tail(x$"Modified Cook distance"))

          Hmisc::prn(utils::head(x$"t statistics"))
          Hmisc::prn(utils::tail(x$"t statistics"))

          Hmisc::prn(x$Call)

     }          # glm
     #
     if(substr(uu,1,12)=="forsearch_lm"){
          if(substr(uu,13,13)=="e"){
               print("This is an forsearch_lme output file", quote=FALSE)
               print("", quote = FALSE)

               Hmisc::prn(names(x))
               Hmisc::prn(x$Dims)
               Hmisc::prn(x$"Number of rows included in Step 1")
               Hmisc::prn(x$Subgroups)
               Hmisc::prn(x$"Rows by subgroup")
               print(search.history(x)[[1]])
               Hmisc::prn(x$Sigma)

               print("", quote=FALSE)
               print("Head and tail of standardized residuals arranged in columns", quote=FALSE)
               Hmisc::prn(utils::head(x$"Standardized residuals"))
               Hmisc::prn(utils::tail(x$"Standardized residuals"))

               print("", quote=FALSE)
               print("Head and tail of fixed parameter estimates arranged in columns", quote=FALSE)
               Hmisc::prn(utils::head(x$"Fixed parameter estimates"))
               Hmisc::prn(utils::tail(x$"Fixed parameter estimates"))

               print("", quote=FALSE)
               print("Head and tail of random parameter estimates arranged in columns", quote=FALSE)
               Hmisc::prn(utils::head(x$"Random parameter estimates"))
               Hmisc::prn(utils::tail(x$"Random parameter estimates"))

               print("", quote=FALSE)
               print("Head and tail of leverage estimates", quote=FALSE)
               Hmisc::prn(utils::head(x$Leverage))
               Hmisc::prn(utils::tail(x$Leverage))

               print("", quote=FALSE)
               Hmisc::prn(x$"Modified Cook distance")
               Hmisc::prn(x$"t statistics")

               print("", quote=FALSE)
               print("Head and tail of fit statistics", quote=FALSE)
               Hmisc::prn(utils::head(x$"Fit statistics"))
               Hmisc::prn(utils::tail(x$"Fit statistics"))
          }        #    lme
          else{
               print("This is a forsearch_lm output file", quote=FALSE)
               print("", quote = FALSE)
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
          }      # else
     }
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
