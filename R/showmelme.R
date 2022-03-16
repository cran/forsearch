#' @export
showmelme <-
function(x, verbose=TRUE)
{
     #                          showmelme
     #
     # VALUE     Shortened output from forsearch_lme function.  Largely, a support for programming efforts
     #
     # INPUT    x            lme diagnostics object
     #
     #          verbose      Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running showmelme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     Hmisc::prn(names(x))
     Hmisc::prn(x$Dims)
     Hmisc::prn(x$"Number of rows included in Step 1")
     Hmisc::prn(x$Subgroups)
     Hmisc::prn(x$"Rows by subgroup")
     print(search.history(x)[[1]])
     Hmisc::prn(x$Sigma)
     Hmisc::prn(utils::head(x$"Standardized residuals"))

     print("", quote=FALSE)
     print("Head and tail of fixed parameter estimates arranged in columns", quote=FALSE)
     Hmisc::prn(utils::head(x$"Fixed parameter estimates"))
     Hmisc::prn(utils::tail(x$"Fixed parameter estimates"))

     print("", quote=FALSE)
     print("Head and tail of random parameter estimates arranged in columns", quote=FALSE)
     vv<- data.frame(x$"Random parameter estimates")
     Hmisc::prn(utils::head(vv))
     Hmisc::prn(utils::tail(vv))

     Hmisc::prn(utils::head(x$Leverage))
     Hmisc::prn(x$"Modified Cook distance")
     Hmisc::prn(x$"t statistics")



#     list("Number of rows included in Step 1"=mstart-1,
#           Subgroups=                         groups.df1, 
#          "Rows by subgroup"=                 df1.by.group,
#          "Rows in stage"=                    rows.in.set,
#           Sigma=                             sigma,
#          "Standardized residuals"=           hold.residuals,            
#          "Fixed parameter estimates"=        param.est,
#          "Random parameter estimates"=       out1,
#           Leverage=                          leverage[-1,],
#          "Modified Cook distance"=           modCook,
#           Dims=                              zholdlme$dims,
#          "t statistics"=                     t.set,
#           Call= MC)

     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running showmelme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
