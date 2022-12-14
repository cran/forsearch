#' @export
forsearch_lm <-
function(formula, data, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL, unblinded=TRUE, diagnose=FALSE, verbose=TRUE)
{
     #                          forsearch_lm
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
     #
     #
     ####################################################
     # Ensure that there is replication in the database #
     ####################################################
     varlist <- variablelist(data, verbose=FALSE)
     if(length(varlist)==dim(data)[1]){
          print("",quote=FALSE)
          print("There is no replication in this dataset. All observations are defined as a combination of factors.", quote=FALSE)
          print("This function does not support such datasets.", quote=FALSE)
          print("",quote=FALSE)
          stop("Try eliminating one or more of the factors in the database and in the formulas of the call.")
     }
     #
     #######################################################################
     # Print the structure of the analysis that will be done on these data #
     # unless the treatment group is blinded                               #
     #######################################################################
     options(warn = -1)
     on.exit(options(warn = 0))
     nulldata <- data
     lhsformula <- formula.tools::lhs(formula)
     ncols <- dim(nulldata)[2]
     ycol <- (1:ncols)[uu==lhsformula]
     nulldata[,ycol] <- 0
     print("", quote = FALSE)
     lmAlldata <- stats::lm(formula, nulldata, singular.ok=TRUE, x=TRUE, y=TRUE)                          # lm
     coeffnames <- names(lmAlldata$coefficients)
     lmAlldata$df.residual <- 99999
     if(unblinded){
          print("**************************************", quote = FALSE)
          print("", quote = FALSE)
          print("Check the following anova paradigm for source of variation and associated degrees of freedom (Df)", quote=FALSE)
          print("", quote = FALSE)
          print(stats::anova(lmAlldata))
          print("", quote = FALSE)
          print("All categorical variables must be defined to be factors in the database, not in the formula.", quote=FALSE)
          print("", quote = FALSE)
          print("**************************************", quote = FALSE)
     }
     y1 <- data[,ycol]
     #
     ######################################
     # Check for factor status of dataset #
     ######################################
     dimdata <- dim(data)[2]
     ufactor <- rep(TRUE, dimdata)
     for(m in 1:dimdata) ufactor[m] <- is.factor(data[,m])
     yesfactor <- any(ufactor)     
     lmAlldata <- stats::lm(formula, data, singular.ok=TRUE, x=TRUE, y=TRUE)                            # lm
     p <- rnk <- lmAlldata$rank
     nopl <- n.obs.per.level 
     if(yesfactor){
          # there are factors in the dataset
          ss77 <- variablelist(datadf = data, verbose=FALSE)
          pickm <- picksome(subsetlist=ss77, nobs=dim(data)[1], initial.sample = initial.sample, n.obs.per.level=nopl, rank=rnk, verbose = FALSE)
          dimpickm <- dim(pickm)[2]  
     }    # yesfactor
     #
     #################################################
     # Set up all output and intermediate structures #
     #################################################
     randset <- 1:initial.sample
     dimx <- dim(data)
     dimx1 <- dimx[1]

     x1 <- data[,-c(1,ycol)]
     if(!is.data.frame(x1))x1 <- data.frame(x1)
     OBS <- data[,1]
     Z <- data                         #  ??????
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result
     result <- matrix(0, nrow=dimx1, ncol=2)
     rows.in.model <- vector("list", dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=p, ncol=dimx1)
     t.set <- param.est

     xtemp.list <- vector("list",dimx1)
     modCook <- rep(0,dimx1)
     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)
     #
     #################################################################################
     # Step 1 of the procedure follows.                                              #
     # By creating an index, randomly reorder the rows of matrix Z = (X,y).          #
     # Within each element of zlist, calculate the several estimates of b from the   #
     # first p of these rows and calculate the median of the residuals in each set.  #
     #################################################################################
     medaugx <- matrix(1, nrow=initial.sample, ncol=2)                       #  median of augx
     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(1:dimx1,dimx1)                        #    sample permutation
          if(yesfactor)zlist[[i]][1:dimpickm,1] <- pickm[i,]
     }      #   i
                                            if(diagnose)Hmisc::prn(zlist)    
     if(is.null(skip.step1)){
          print("ENTERING STEP 1", quote=FALSE)
          inner.rank <- lmAlldata$rank
          mstart <- inner.rank + 1
          firstrim <- aStep1(yesfactor, data, inner.rank=lmAlldata$rank, initial.sample, formula, ycol, nopl=n.obs.per.level)
          rows.in.model[[p]] <- firstrim
     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          p <- length(skip.step1)
          rows.in.model[[p]] <- skip.step1
          datastep1 <- data[as.vector(rows.in.model[[p]]),]
          lmchosen <- stats::lm(formula, datastep1, singular.ok=TRUE)                                  #  lm 
          betahat <- lmchosen$coefficients
          mstart <- p + 1
          betaMMinus1 <- betahat
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
          mstart <- p + 1
          betaMMinus1 <- betahat
     }       # skipping step 1
     #
     #############################################
     # Step 2 of the procedure follows.          #
     # Adding observations to the initial set    #
     # First, get all stats for Step 1 subset.   #
     # Then calculate next subset at end of loop #
     #############################################
     print("ENTERING STEP 2", quote=FALSE)
     betahatset <- matrix(0,nrow=dimx1,ncol=p)
     for(i in mstart:(dimx1+1)){                  # mstart is the step after the original p obs entered
          rim <- rows.in.model[[i-1]]             # picks up rows for previous step
          Zlatest <- Z[rim,]
                                       if(diagnose){Hmisc::prn(i); Hmisc::prn(Zlatest)}
          if(is.data.frame(x1))xtemp <- x1[rim,]
          else xtemp <- data.frame(x1[rim])
          subdata <- data[rim,]
          getthislm <- stats::lm(formula, subdata, singular.ok=TRUE, x=TRUE, y=TRUE)                                  #   lm
          getthislmx <- getthislm$x
          xtemp.list[[i-1]] <- getthislmx                          # used in modified Cook
          betahat <- getthislm$coefficients                                                             #    correct????
          # Allowing singular results means we have to replace NAs with some value, here 0 
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
          #
          #################################
          # Extract t values from summary #
          #################################
          if(i > mstart){
               t.setbase <- c(summary(getthislm)$coefficients[,3], rep(0,p))
               t.set[,i-1] <- t.setbase[1:p]
          }
          #
          ############
          # Leverage #
          ############
          thisleverage <- 1
          x1 <- data[,-c(1,ycol)]
          xtemp <- getthislmx
          crossinv <- solve(t(xtemp) %*% xtemp)
          for(j in 1:(i-1)){
              Zlatest2 <- data.frame(Zlatest)
              if(dim(xtemp)[2]==1){
                   thisleverage <- c(c(matrix(xtemp[j],nrow=1) %*% crossinv %*% matrix(xtemp[j],ncol=1)))
                   thisleverage <- c(i-1,Zlatest2[j,1],thisleverage)
              }else{
                   thisleverage <- c(c(matrix(xtemp[j,],nrow=1) %*% crossinv %*% matrix(xtemp[j,],ncol=1)))
                   thisleverage <- c(i-1,Zlatest2[j,1],thisleverage)
              }
              leverage <- rbind(leverage,thisleverage)
               }   # j 1:i
                                             if(diagnose) Hmisc::prn(leverage)
          #
          #####################################################################
          # Calculate the error for each observation using this model betahat #
          # Select next subset.                                               #
          #####################################################################
          preds <- stats::predict(getthislm, data, pred.var=1)                                              #   predict
          errors <- y1 - preds
          residuals[,i-1] <- errors
          errors2 <- errors * errors
          medaugx[,2] <- errors2             
          medaugx <- medaugx[order(medaugx[,2]),]
          if(i < dimx1+1){
               thisrim <- aStep2(getthislm, data, ycol, thisi=i)                                      #   aStep2 returns a list
               rows.in.model[[i]] <- thisrim[[1]]
          }        #  if i < dimx1
     }            # i in mstart ...
     #
     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m", "Observation", "leverage"))
     param.est <- as.data.frame(t(param.est))
     names(param.est) <- coeffnames
     m <- 1:dimx1
     param.est <- cbind(m,param.est)                                 #  here we add the m column
     # 
     t.set <- as.data.frame(t(t.set))
     names(t.set) <- coeffnames
     t.set <- cbind(m,t.set)
     #
     ############################################
     # Estimate sigma and standardize residuals #
     ############################################
     sigma <- sqrt(sum(medaugx[,2])/(dimx1-p))
     residuals <- residuals/sigma
     ##########################
     # Modified Cook distance #
     ##########################
     nms <- dim(param.est)[1]
     param.est.prev <- param.est[-nms,]
     param.est.current <- param.est[-1,]          # remove first and last rows
     param.diff <- param.est.prev - param.est.current
                         if(diagnose){
                             Hmisc::prn(param.est)
                             Hmisc::prn(param.est.prev)
                             Hmisc::prn(param.est.current)
                             Hmisc::prn(param.diff)
                             Hmisc::prn(xtemp.list)     
                         }
     for(i in mstart:dimx1){
         aa <- param.diff[i-1,]         # select row i-1
         aa <- as.numeric(aa[,-1])      # remove col 1
         aa <- matrix(aa,nrow=1)        # turn this into a matrix with 1 row
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
     listout <- list(
          "Rows in stage"=                     rows.in.model,
          "Standardized residuals"=            residuals, 
          "Number of model parameters"=        p, 
          "Sigma"=                             sigma, 
          "Fixed parameter estimates"=         param.est, 
          "s^2"=                               s.2, 
          "Leverage"=                          leverage, 
          "Modified Cook distance"=            modCook, 
          "t statistics"=                      t.set,
          "Call"=                              MC)
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_lm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(listout)
}
