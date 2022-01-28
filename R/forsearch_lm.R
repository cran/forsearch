#' @export
forsearch_lm <-
function(formula, data, initial.sample, diagnose=FALSE, verbose=TRUE)
{
     #                          forsearch_lm
     #
     # VALUE         Prepares input for diagnostic plotting of database to be analyzed by lm function.
     #
     # INPUT 
     #             formula            Include -1 to remove assumed constant independent variable (results in cell means analysis)
     #             data               Name of data frame
     #             initial.sample     Number of reorderings of 1:n
     #
     #             diagnose           Logical. TRUE causes printing of diagnostic content
     #             verbose            Logical. TRUE causes printing of program ID before and after running.
     #
     # EXAMPLE:  
     #
     # REF:  Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running forsearch_lm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     ##################################################################
     # Ensure that first independent variable is Observation          #
     # Make x1 the matrix of independent variables without obs number #
     # x1 will represent the formula for independent obs              #
     ##################################################################
     uu <- names(data)
     if(uu[1] != "Observation")stop("First column of data must be 'Observation'")

     lmAlldata <- stats::lm(formula, data, singular.ok=TRUE, x=TRUE, y=TRUE)                                         # lm
     coeffnames <- names(lmAlldata$coefficients)
     z1 <- lmAlldata$x
     y1 <- lmAlldata$y

     #########################################################
     # Add a little random difference to avoid singularities #
     #########################################################
     dimz1x <- dim(z1)
     rantimes <- stats::runif(prod(dimz1x),0,min(abs(z1))/100)
     rantimes <- matrix(rantimes,dimz1x[1],dimz1x[2])
     x1 <- z1 + rantimes
                                      if(diagnose) {Hmisc::prn(x1); Hmisc::prn(y1)}
     #
     ############################################################################
     # Resampling algorithm for Step 1 of forward search     p31                #
     # p is the dimension of the independent variables, omitting 'Observation'  #
     # rows.in.model is a list whose elements are vectors of increasing length  #
     # residuals contains model residual for each term at each stage of the fit #
     # The residuals for i-th observation are in i-th row                       #
     ############################################################################
 
     p <- dim(lmAlldata$x)[2]
     randset <- 1:initial.sample
     dimx <- dim(data)
     dimx1 <- dimx[1]
     OBS <- data[,1]
     Z <- cbind(OBS,x1,y1)
     zlist <- vector("list",initial.sample)
     result <- matrix(0, nrow=dimx1, ncol=2)
     rows.in.model <- vector("list",dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=p, ncol=dimx1)




     xtemp.list <- vector("list",dimx1)
     modCook <- rep(0,dimx1)
     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)
     #
     ################################################################################
     # Step 1 of the procedure follows.                                             #
     # By creating an index, randomly reorder the rows of matrix Z = (X,y).         #
     # Within each element of zlist, calculate the several estimates of b from the  #
     # first p of these rows and calculate the median of the residuals in each set. #
     # Don't run lm( ) in Step 1                                                    # 
     # Don't enter errors into residuals in Step 1                                  #
     ################################################################################

     print("ENTERING STEP 1", quote=FALSE)

     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(1:dimx1,dimx1)                        #    sample permutation
     }      #   i
                                            if(diagnose)Hmisc::prn(zlist)    
     medaugx <- matrix(1, nrow=initial.sample, ncol=2)
     for(i in 1:initial.sample){
               rowz <- zlist[[i]]
               xtemp <- as.matrix(x1[rowz,])
               ytemp <- as.matrix(c(y1[rowz]))
                                            if(diagnose){
                                               Hmisc::prn(i)
                                               Hmisc::prn(rowz)
                                               Hmisc::prn(xtemp)
                                               Hmisc::prn(ytemp)
                                            }
               transtemp <- t(xtemp)
               cross <- transtemp %*% xtemp

               crossinv <- solve(cross)
               betahat <- (crossinv %*% transtemp) %*% ytemp
               errors <- ytemp - xtemp %*% betahat
                                            if(diagnose)Hmisc::prn(errors)
               errors2 <- errors * errors
                                # rank order the squared errors and extract the median
               augx <- matrix(rowz[,1],nrow=dimx1,ncol=1)    
               augx <- cbind(augx,errors2)
               augx <- augx[order(augx[,2]),]
               medaugx[i,] <- augx[floor((dimx1+p+1)/2),]
                                            if(diagnose)Hmisc::prn(augx)
     }   #  i in 1:initial.sample
     randset <- 1:initial.sample
     medaugx <- cbind(medaugx,randset)
     minmed <- min(medaugx[,2])
     locatemin <- medaugx[medaugx[,2]==minmed,]
     if(is.matrix(locatemin)) locatemin <- locatemin[1,]
     zliststar <- zlist[[locatemin[3]]]
                                           if(diagnose){
                                              Hmisc::prn(medaugx)
                                              Hmisc::prn(minmed)
                                              Hmisc::prn(locatemin)
                                              Hmisc::prn(locatemin)
                                              Hmisc::prn(zliststar)
                                              Hmisc::prn(zliststar[1:p,1])
                                           }
     rows.in.model[[p]] <- c(zliststar[1:p,1])
     mstart <- p + 1
     betaMMinus1 <- betahat
                                           if(diagnose){
                                              Hmisc::prn(rows.in.model[[p]])
                                              print("",quote=F)
                                              print("Step 2",quote=F)
                                              print("",quote=F)
                                           }
     #
     ##################################################################
     # Step 2 of the procedure follows.                               #
     # Adding observations to the initial set                         #
     # Outer i loop adds 1 each time to list of observations in model #
     # May not be the same set of observations                        #
     ##################################################################

     print("ENTERING STEP 2", quote=FALSE)
     betahatset <- matrix(0,nrow=dimx1,ncol=p)
     for(i in mstart:(dimx1+1)){                  # mstart is the step after the original p obs entered
          rim <- rows.in.model[[i-1]]         # picks up rows for previous step
          Zlatest <- Z[rim,]
                                            if(diagnose){Hmisc::prn(i); Hmisc::prn(Zlatest)}
          xtemp <- x1[rim,]
          xtemp.list[[i-1]] <- xtemp                                              # used in modified Cook
          transtemp <- t(xtemp)
          cross <- transtemp %*% xtemp
          crossinv <- solve(cross)                               # inverts a matrix
          getmatch <- data[,1] %in% rim
          subdata <- data[getmatch,]

          getthislm <- stats::lm(formula, subdata, singular.ok=TRUE)                                  #   lm
                                               if(diagnose) Hmisc::prn(getthislm)
          betahat <- getthislm$coefficients
          # Allowinig singular results means we have to replace NAs with 0 #
          if(any(is.na(betahat))){
               nabetahat <- is.na(betahat)
               betahat[nabetahat] <- 0
          }

          model.resids <- getthislm$residuals
          betahatset[i-1,] <- betahat
          medaugx <- matrix(1:dimx1, nrow=dimx1, ncol=2, byrow=FALSE)      # initialize medaugx
          medaugx[,2] <- 0 
          param.est[,i-1] <- betahat
          if(i > p) s.2[i-1] <- sum(model.resids * model.resids)/(i-p)
                                            if(diagnose) {Hmisc::prn(model.resids); Hmisc::prn(s.2[i-1])}
          #
          ############
          # Leverage #
          ############
          thisleverage <- 1
          if(is.matrix(x1)){
               for(j in 1:(i-1)){
                   Zlatest2 <- data.frame(Zlatest)
                   if(dim(x1)[2]==1){
                        thisleverage <- c(c(matrix(xtemp[j],nrow=1) %*% crossinv %*% matrix(xtemp[j],ncol=1)))
                        thisleverage <- c(i-1,Zlatest2[j,1],thisleverage)
                   }else{
                        thisleverage <- c(c(matrix(xtemp[j,],nrow=1) %*% crossinv %*% matrix(xtemp[j,],ncol=1)))
                        thisleverage <- c(i-1,Zlatest2[j,1],thisleverage)
                   }
                   leverage <- rbind(leverage,thisleverage)
               }   # j 1:i
          }        #   x1 is matrix
                                             if(diagnose) Hmisc::prn(leverage)
          #
          #####################################################################
          # Calculate the error for each observation using this model betahat #
          # Keep the original order of the observations                       #
          #####################################################################
          errors <- rep(-999,dimx1)
          for(j in 1:dimx1){
               errors[j] <- y1[j] - sum((betahat) * x1[j,]) 
          }         # j
                                            if(diagnose)Hmisc::prn(errors)
          residuals[,i-1] <- errors
          errors2 <- errors * errors
          medaugx[,2] <- errors2             
                                            if(diagnose){Hmisc::prn(medaugx); Hmisc::prn(residuals)}
          medaugx <- medaugx[order(medaugx[,2]),]
          if(i <= dimx1){
               rows.in.model[[i]] <- medaugx[1:i,1]
          }
     }            # i in mstart ...
     #
     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m", "Observation", "leverage"))
     param.est <- as.data.frame(t(param.est))
     names(param.est) <- coeffnames
     m <- 1:dimx1
     param.est <- cbind(m,param.est)
     # 
     ############################################
     # Estimate sigma and standardize residuals #
     ############################################
     sigma <- sqrt(sum(medaugx[,2])/(dimx1-p))
     residuals <- residuals/sigma
     #
     ##########################
     # Modified Cook distance #
     ##########################
     nms <- dim(param.est)[1]
     param.est.current <- param.est[-1,]
     param.est.prev <- param.est[-nms,]
     param.diff <- param.est.prev - param.est.current
                         if(diagnose){
                             Hmisc::prn(param.est)
                             Hmisc::prn(param.est.prev)
                             Hmisc::prn(param.est.current)
                             Hmisc::prn(param.diff)
                             Hmisc::prn(xtemp.list)     
                         }
     for(i in mstart:dimx1){
         aa <- param.diff[i-1,]
         aa <- as.numeric(aa[,-1])
         aa <- matrix(aa,nrow=1)
         bb <- xtemp.list[[i]]
                       if(diagnose){
                             Hmisc::prn(aa)
                             Hmisc::prn(bb)
                       }
         www <- aa %*% t(bb)
                       if(diagnose){Hmisc::prn(www)}
         modCook[i] <- (www %*% t(www))/(p*s.2[i])
     }       #    i
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_lm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     list(
          "Rows in stage"=                     rows.in.model,
          "Standardized residuals"=            residuals, 
          "Number of model parameters"=        p, 
           Sigma=                              sigma, 
          "Fixed parameter estimates"=         param.est, 
          "s^2"=                               s.2, 
           Leverage=                           leverage, 
          "Modified Cook distance"=            modCook, 
           Call=                               MC)
}
