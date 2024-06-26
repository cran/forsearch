#' @export
identifyCoeffs <-
function(fixed, data, random, 

    XmaxIter=1000,
    XmsMaxIter=1000, 
    Xtolerance=.01,
    XniterEM=1000,
    XmsMaxEval=400,
    XmsTol=.00001, 
    Xopt='optim',
    verbose=TRUE)
{
     #                                                       identifyCoeffs
     #
     # VALUE    Runs the defined, grouped linear mixed effects (lme) model. Displays the resulting fixed and random coefficients. 
     #          Supplies the codes for identifying them to the plotting functions of this package.
     #
     # INPUT    fixed               2-sided formula for fixed effects
     #          data                Name of grouped or ungrouped data file
     #          random              1-sided formula for random effects
     #
     #          Xmaxiter, XmsMaxIter, Xtolerance, XniterEM, XmsMaxEval, XmsTol, Xopt 
     #                             control variates for lme function
     #
     #          verbose             Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running identifyCoeffs", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     ###########################################################################
     # Run the lme function on the input provided and extract the coefficients #
     ###########################################################################
     newcontrol <- nlme::lmeControl(maxIter=XmaxIter,  msMaxIter=XmsMaxIter, tolerance=Xtolerance, 
                      niterEM=XniterEM, msMaxEval=XmsMaxEval, msTol=XmsTol, opt=Xopt)
     newcontrol <<- newcontrol
     zzzz <- data
     zzzz <<- zzzz
     on.exit(rm(zzzz, pos=1))
     on.exit(rm(newcontrol, pos=1))

     lmerun <- nlme::lme(fixed,data,random, control=newcontrol) 
     coeffs <- lmerun$coefficients
     fixd <- coeffs[[1]]
     rndm <- coeffs[[2]]

     print("Typical fixed coefficient estimates from a run of lme( )", quote=FALSE)
     print(fixd)
     print("", quote = FALSE)
     print("", quote = FALSE)
     print("Typical random coefficient estimates from a run of lme( )", quote=FALSE)
     print(rndm)
     print("", quote = FALSE)
     #
     ######################################################
     # Develop the codes for selected random coefficients #
     ######################################################
     namesrndm <- names(rndm)
     dimsN <- lmerun$dims$N
     dimsQ <- lmerun$dims$Q
     IDs <- NULL
     for(iq in 1:dimsQ){
          uu <- rndm[[iq]]
          vv <- apply(uu^2,2,mean)
          namesvv <- names(vv)
          IDs1 <- namesrndm[iq]
          IDs2 <- namesvv
          IDs <- c(IDs, paste(rep(IDs1,times=length(IDs2)), rep(IDs2,each =length(IDs1)),sep="--"))
          IDs <- c(IDs)
     }     #  iq


     vectfixd <- 1:length(fixd)
     names(vectfixd) <- names(fixd)
     vectrandm <- 1:length(IDs)
     names(vectrandm) <- IDs
     vectrandm <- as.matrix(vectrandm, col=1)

     print("When only a few coefficients can be graphed, specify them using the following codes:", quote=FALSE)
     print("", quote = FALSE)
     print("Fixed coefficient estimates", quote=FALSE)
     print(vectfixd, quote = FALSE)
     print("", quote = FALSE)
     print("", quote = FALSE)

     print("Root mean squares of random coefficient estimates", quote=FALSE)
     print(vectrandm)
     print("", quote = FALSE)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running identifyCoeffs", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     list("Fixed coefficients codes"=vectfixd, "Random RMS codes"=vectrandm)
}
