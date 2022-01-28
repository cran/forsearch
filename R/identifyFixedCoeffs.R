#' @export
identifyFixedCoeffs <-
function(formula, data, diagnose=FALSE, verbose=TRUE)
{
     #                                                       identifyFixedCoeffs
     #
     # VALUE    Runs the linear (lm) model. Displays the resulting coefficients. 
     #          Supplies the codes for identifying them to the plotting functions of this package.
     #
     # INPUT    formula             2-sided formula for fixed effects
     #          data                Name of grouped or ungrouped data file
     #
     #          diagnose            Logical. TRUE causes printing of diagnostic content
     #          verbose             Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running identifyFixedCoeffs", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     ##########################################################################
     # Run the lm function on the input provided and extract the coefficients #
     ##########################################################################
     zzzz <- data
     zzzz <<- zzzz
     on.exit(rm(zzzz, pos=1))

     lmrun <- stats::lm(formula, data) 
     coeffs <- lmrun$coefficients
     fixd <- coeffs
     print("Typical fixed coefficient estimates from a run of lm( ) on data", quote=FALSE)
     print(fixd)
     print("", quote = FALSE)
     print("", quote = FALSE)
     #
     ###############################################
     # Develop the codes for selected coefficients #
     ###############################################
     ncoeffs <- length(fixd)
     code <- 1:ncoeffs
     coefficient <- names(fixd)
     codeset <- data.frame(coefficient,code)

     print("When only a few coefficients can be graphed, specify them using the following codes:", quote=FALSE)
     print(codeset)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running identifyFixedCoeffs", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     list("Fixed coefficients"=fixd, "Use these fixed codes"=codeset, Call=MC)
}
