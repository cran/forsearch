#' @export
forsearch_glm <-
function(
          initial.sample=1000, 
          response.cols,
          indep.cols,
          family,
          formula=NULL,
          binomialrhs=NULL,
          data,  
          n.obs.per.level=1,
          estimate.phi= TRUE, 
          skip.step1=   NULL,   
          unblinded = TRUE,

         begin.diagnose=  100,          verbose=    TRUE)
{
#                                                                      forsearch_glm
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running forsearch_glm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     # Doesn't follow general plan for diagnostics because doesn't call xStepy #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                   "     # used for begin.diagnose prints

     ############################
     # Identity matrix function #
     ############################
     my.identity <- function(n){
          matrix(rep(c(1, rep(0, n)), times = n)[1:(n * n)], nrow = n, ncol = n, byrow = T)
     }
     ######################
     # 2 Helper functions #
     ######################

     ########################################################################################
     # Helper function for calculaton of deviance code for each observation based on family #
     ########################################################################################
     devianceCode <- function (obs, pred, ni=NULL, fam) 
     {
          if(obs/pred < 0) {
               out <- 0
          }
          else{
               if(fam=="binomial"){
                    if(obs==0){
                         out <- ni*log(ni/(ni-pred))                   # expects whole numbers, not proportions
                    }
                    else if(obs==ni){
                         out <- obs*log(obs/pred) 
                    }
                    else{
                         out <- obs*log(obs/pred) + (ni-obs)*log((ni-obs)/(ni-pred))  
                    }
               }
               else if(fam=="Gamma"){
                   out <- -log(obs/pred) + (obs-pred)/pred
               }
               else if(fam=="poisson"){
                    if(obs==0){
                        out <- pred
                    }
                    else{
                        out <- obs*log(obs/pred) - obs + pred
                    }
               }
               else if(fam=="exponential"){
                    out <- -log(obs/pred) + obs*pred - 1
               }
               else{
                    stop("family name not recognized")
               }
               if(is.na(out)){
                    out <- 0
               }
               else{
                    if(abs(out)<= 10^(-12)){
                         out <- 0
                    }
                    else{
                         out <- 2 * out
                         vv <- obs - pred
                         vv <- vv/abs(vv)            #  1 or -1
                         out <- sqrt(out)*vv
                    }
               }
          }
          return(out)
     }                                   # end of devianceCode function
     #

                                                                         #############################
                                                                         # MAIN FUNCTION STARTS HERE #
                                                                         #############################
     options(warn = -1)
     on.exit(options(warn = 0))

     nobs <- dim(data)[1]
 
                                 if(begin.diagnose <= 1){print(paste(spacer,"Section 1",sep=" "),quote=FALSE);Hmisc::prn(utils::head(data));
                                        Hmisc::prn(formula);Hmisc::prn(nobs);Hmisc::prn(response.cols);Hmisc::prn(indep.cols)       }

     ##################################################################
     # Get parts of formula to enable subsetting                      #
     # Ensure that first independent variable is Observation          #
     ##################################################################
     uu <- names(data)
     if(uu[1] != "Observation")stop("First column of data must be 'Observation'")
     #
     if(family[[1]]=="binomial"){
          name.response <- uu[response.cols]
          inresponse <- data[,response.cols]
          bin.wts <- apply(data[,response.cols],1,sum)
          proportion1 <- data[,response.cols[1]]
          proportion1 <- proportion1/bin.wts
          indata <- data.frame(data,proportion1,bin.wts)                             # indata has proportion1 and bin.wts
          bin.wts <- NULL                                                            # insist that these come from indata
          genformula <- formula(paste("proportion1", binomialrhs, sep=" ~ ")) 
    }      # binomial  
     else{
          name.response <- uu[response.cols]
          inresponse <- data[,response.cols]
          indata <- data
          genformula <- formula
     }     # not binomial

     print("The formula applied in this analysis was:", quote=FALSE)
     print("", quote=FALSE)
     print(genformula, quote=FALSE)
     print("", quote=FALSE)
     print(paste("response.cols = ", name.response, sep=""), quote=FALSE)
     print("", quote=FALSE)

     ############################################################
     # Print structure of analysis by blinding of real response #
     ############################################################
     indataXX <- indata
     indataXX[,response.cols] <- 1
     if(family[[1]]=="binomial"){
          lmAlldata <- stats::glm(formula=genformula, family=family, data=indataXX, weights=bin.wts, singular.ok=TRUE, x=TRUE, y=TRUE)         # glm
     }
     else if(family[[1]]=="Gamma"){
          lmAlldata <- stats::glm(formula=genformula, family=family, data=indataXX, singular.ok=TRUE, x=TRUE, y=TRUE)                          # glm
     }
     else if(family[[1]]=="poisson"){
          lmAlldata <- stats::glm(formula=genformula, family=family, data=indataXX, singular.ok=TRUE, x=TRUE, y=TRUE)                          # glm
     }
     else{
          stop("family name not recognized")
     }
     lmAlldata$df.residual <- 99999
     lmAlldata$df.null     <- 99999

     print("", quote = FALSE)
     if(unblinded){
          print("Check the following paradigm structure to be sure that the assumed analysis for these data is correct:", quote=FALSE)
          print("", quote = FALSE)
          print(lmAlldata)
          print("",quote=FALSE)
     }
     ############################################################################
     # Check for factor status of dataset                                       #
     # Define number of factors and continuous variables in the independent set #
     ############################################################################
     ndatacols <- length(indep.cols)
     ufactor <- rep(TRUE, ndatacols)
     for(m in 1:ndatacols){
          index <- indep.cols[m]
          ufactor[m] <- is.factor(data[,index])
     }
     yesfactor <- any(ufactor)
     lmAlldata <- stats::lm(formula=genformula, data, singular.ok=TRUE, x=TRUE, y=TRUE)                                    # lm
     ncoeffs <- length(lmAlldata$coefficients)
     rnk <- lmAlldata$rank

                                 if(begin.diagnose <= 2){print(paste(spacer,"Section 2",sep=" "),quote=FALSE);Hmisc::prn(ncoeffs);
                                       Hmisc::prn(rnk)       }
 
     nopl <- n.obs.per.level
     ######################################################
     # Define p = number of continuous variables in X and #
     #   nfacts = number of factor variables in X         #
     ######################################################

     p <- length(indep.cols)       # temporary
     nfacts <- 0
     if(yesfactor){
          nfacts <- sum(ufactor)
          p <- p - nfacts
          ss77all <- variablelist(datadf = data, prank = p)
          nss77all <- length(ss77all)
     }
 
                                 if(begin.diagnose <= 3){print(paste(spacer,"Section 3",sep=" "),quote=FALSE);Hmisc::prn(ndatacols);
                                       Hmisc::prn(nfacts);Hmisc::prn(p)       }
 
     if(yesfactor){
          print("There are factors in the design", quote=FALSE)
          ss77 <- variablelist(datadf = data, prank=p)

                                if(begin.diagnose <= 3.5){print(paste(spacer,"Section 3.5",sep=" "),quote=FALSE);Hmisc::prn(ss77);
                                       Hmisc::prn(nfacts);Hmisc::prn(p)       }

          pickm <- picksome(subsetlist=ss77, nobs=dim(data)[1], initial.sample = initial.sample, n.obs.per.level=nopl, rank=p+1)
          dimpickm <- dim(pickm)[2]  
     }
     else{ print("There are no factors in the design", quote=FALSE) }

                                if(begin.diagnose <= 4){print(paste(spacer,"Section 4",sep=" "),quote=FALSE);Hmisc::prn(nfacts);Hmisc::prn(p)       }

     #
     ############################################
     # Restore to previous version of lmAlldata #
     ############################################
     if(family[[1]]=="binomial"){
          lmAlldata <- stats::glm(formula=genformula, family=family, data=indata, weights=bin.wts, singular.ok=TRUE, x=TRUE, y=TRUE)         # glm
     }
     else if(family[[1]]=="Gamma"){
          lmAlldata <- stats::glm(formula=genformula, family=family, data=indata, singular.ok=TRUE, x=TRUE, y=TRUE)                          # glm
     }
     else if(family[[1]]=="poisson"){
          lmAlldata <- stats::glm(formula=genformula, family=family, data=data, singular.ok=TRUE, x=TRUE, y=TRUE)                            # glm
     }
     else{
          stop("family name not recognized")
     }
     coeffnames <- names(lmAlldata$coefficients)
     ncoeffs <- length(coeffnames)

                                 if(begin.diagnose <= 5){print(paste(spacer,"Section 5",sep=" "),quote=FALSE);Hmisc::prn(coeffnames);
                                       Hmisc::prn(ncoeffs)       }

     x1 <- z1 <- lmAlldata$x
     y1 <- lmAlldata$y
     dimdata <- dim(data)
      #
     ############################################################################
     # Resampling algorithm for Step 1 of forward search     p31                #
     # p is the dimension of the continuous independent variables               #
     # rows.in.model is a list whose elements are vectors of increasing length  #
     # residuals contains model residual for each term at each stage of the fit #
     # The residuals for i-th observation are in i-th row                       #
     ############################################################################
     randset <- 1:initial.sample
     dimx <- dim(data)
     dimx1 <- dimx[1]
     OBS <- data[,1]
     Z <- cbind(OBS,x1,y1)
     rownums <- rep(FALSE,dimx1)
     rownums <- data.frame(rownums)
     zlist <- vector("list",initial.sample)                                           # Step 1
     result <- rep(0,p)                                                               # step 1
     deviance.matrix <- matrix(1:dimx1, nrow=dimx1, ncol=4, byrow=FALSE)              # step 1

     zlist.inner <- list(result=result,devmat=deviance.matrix,weights=NULL)           # step 1
     rows.in.model <- vector("list",dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=ncoeffs, ncol=dimx1)
     t.set <- param.est                                        # t statistics
     xtemp.list <- vector("list",dimx1)          
     modCook <- rep(0,dimx1)                                                                                              #CD
     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)

                                 if(begin.diagnose <= 6){print(paste(spacer,"Section 6",sep=" "),quote=FALSE);Hmisc::prn(utils::head(param.est));
                                       Hmisc::prn(utils::tail(param.est));Hmisc::prn(s.2)       }
     #
####################################################################################################################################################
# Step 1
     ################################################################################
     # Step 1 of the procedure follows.                                             #
     # By creating an index, randomly reorder the rows of matrix Z = (X,y).         #
     # Within each element of zlist, calculate the several estimates of b from the  #
     # first p of these rows and calculate the median of the residuals in each set. #
     # Don't run lm( ) in Step 1                                                    # 
     ################################################################################
     if(is.null(skip.step1)){
          print("ENTERING STEP 1", quote=FALSE)                                                                          # Step 1
          for(i in 1:initial.sample){
               if(!yesfactor){
                   this.sample <- sample(1:dimx1,nobs)                                                                   # sample
                   zlist.inner$result <- this.sample
                   this.data <- indata[this.sample,]
               }     # if not yesfactor
               if(yesfactor){
                   this.sample <- pickm[i,]
                   zlist.inner$result <- pickm[i,]
                   this.data <- indata[this.sample,]
               }
               if(family[[1]]=="binomial"){
                    glm.s2 <- stats::glm(formula=genformula, family=family, data=this.data, weights=bin.wts, singular.ok=TRUE)          #   glm
               }
               else if(family[[1]]=="Gamma"){
                    glm.s2 <- stats::glm(formula=genformula, family=family, data=this.data, singular.ok=TRUE)                           #   glm
               }
               else if(family[[1]]=="poisson"){
                    glm.s2 <- stats::glm(formula=genformula, family=family, data=this.data, singular.ok=TRUE)                           #   glm
               }
               else{
                    stop("family name not recognized")
               }
               pred.s2 <- stats::predict.glm(glm.s2, newdata=data, type="response")                                                     #   predict
               zlist.inner$devmat[,2] <- pred.s2      #                                                            # pred.s2 for each candidate ?

                                    if(begin.diagnose <= 20){print(paste(spacer,"Section 20",sep=" "),quote=FALSE);Hmisc::prn(pred.s2)       }

               if(family[[1]]=="binomial"){
                    for(j in 1:dimdata[1]){
                         this.obsy1 <- y1[j]*indata$bin.wts[j]
                         this.pred <- zlist.inner$devmat[j,2]*indata$bin.wts[j]
                         this.weight <- indata$bin.wts[j]
                         if(this.pred==this.obsy1){
                              zlist.inner$devmat[j,3] <- 0
                         }
                         else{
                              zlist.inner$devmat[j,3] <- (devianceCode(this.obsy1, this.pred, this.weight, fam=family[[1]]))^2          # devianceCode
                         }
                    }        # j in 1:dimdata 2
               }
               else if(family[[1]]=="Gamma"){
                    for(j in 1:dimdata[1]){
                         this.obsy1 <- y1[j]
                         this.pred <- zlist.inner$devmat[j,2]
                         zlist.inner$devmat[j,3] <- (devianceCode(this.obsy1, this.pred, fam=family[[1]]))^2                           # devianceCode
                    }        # j in 1:dimdata 2
               }
               else if(family[[1]]=="poisson"){
                    for(j in 1:dimdata[1]){
                         this.obsy1 <- y1[j]
                         this.pred <- zlist.inner$devmat[j,2]
                         zlist.inner$devmat[j,3] <- (devianceCode(this.obsy1, this.pred, fam=family[[1]]))^2                          # devianceCode
                    }        # j in 1:dimdata 2
               }
               else{
                    stop("family name not recognized")
               }
               #
 
                                 if(begin.diagnose <= 22){print(paste(spacer,"Section 22",sep=" "),quote=FALSE);Hmisc::prn(y1);
                                       Hmisc::prn(zlist.inner$devmat)       }
 
 
               #######################
               # get median deviance #
               #######################
               devs <- sort(zlist.inner$devmat[,3])
               median.dev <- devs[floor((dimx1+p+1)/2)]
               zlist.inner$devmat[,4] <- median.dev 
               zlist[[i]] <- zlist.inner
          }       # for i in 1 initial.sample

                                        if(begin.diagnose <= 24){print(paste(spacer,"Section 24",sep=" "),quote=FALSE);Hmisc::prn(zlist[[i]][[1]]);Hmisc::prn(utils::head(zlist[[i]][[2]]));
                                                          Hmisc::prn(utils::tail(zlist[[i]][[2]]))       }

         # select set with min(median deviance)
          candidate.sets<- rep(0,initial.sample)
          for(ib in 1:initial.sample){
               if(is.na(zlist[[ib]]$devmat[1,4])){
                    zlist[[i]]$devmat[1,4] <- zlist[[i]]$devmat[1,4] + 1
               }
               candidate.sets[ib] <- zlist[[ib]]$devmat[1,4]
          }     #  ib
 
          min.candidate <- min(candidate.sets)
          for(id in 1:initial.sample){
               if(zlist[[id]]$devmat[1,4]==min.candidate){
                    got.set <- id
                    break
               }    # if zlist
          }   # id
          chosen.set <- zlist[[got.set]]$result                             # chosen set



          ninStep1 <- nopl*rnk + 1
          chosen.set <- chosen.set[1:ninStep1]        
          lengotset <- length(chosen.set)
          rows.in.model[[lengotset]] <- chosen.set                          # obs chosen for Step 1
          mstart <- lengotset + 1
          SOON <- chosen.set
     }       #   skip.step1 is NULL ?
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          lenskip <- length(skip.step1)
          rows.in.model[[lenskip]] <- skip.step1
          SOON <- skip.step1
          mstart <- lenskip + 1
     }   #   if not skip step 1
                                                  if(begin.diagnose <= 49){print(paste(spacer,"Section 49",sep=" "),quote=FALSE);Hmisc::prn(SOON);Hmisc::prn(mstart)       }
     #
###################################################################################################################################################
# Step2
     ##################################################################
     # Step 2 of the procedure follows.                               #
     # Adding observations to the initial set                         #
     # Outer i loop adds 1 each time to list of observations in model #
     # May not be the same set of observations                        #
     ##################################################################

     print("ENTERING STEP 2", quote=FALSE)

     #####################################################
     # Initialize storage matrices and vectors           #
     # Variables starting "glm" created in this function #
     #####################################################
     betahatset <- matrix(0,nrow=dimx1,ncol=ncoeffs)
     glmcoeff <- matrix(0,nrow=dimx1, ncol=p)
     glmdeviance <- rep(0,dimx1)
     glmphi <- rep(0,dimx1)
     glmnulldeviance <- glmdeviance
     glmaic <- rep(0,dimx1)
     ADdevres <- vector("list",dimx1)
     storeDevianceResiduals <- vector("list",dimx1)         # initially saved without studentization as vector
     glmstudresids <- NULL                                  # after studentization 
     glmHat <- vector("list", dimx1)                        # will store Hat matrices for leverage, Cook, and (last one) studentization of errors
     #

                                             if(begin.diagnose <= 51){print(paste(spacer,"Section 51",sep=" "),quote=FALSE);Hmisc::prn(glmdeviance)       }

     ###################################################################################################
     # The next section increments the number of observations while collecting entermediate statistics #
     ###################################################################################################
     for(i in mstart:(dimx1-1)){                                                  # mstart is the step after the original p obs entered
          rim <- rows.in.model[[i-1]]                                             # picks up rows for previous step
          Zlatest <- Z[rim,]
                                          if(begin.diagnose <= 53){print(paste(spacer,"Section 53",sep=" "),quote=FALSE);Hmisc::prn(Zlatest)       }

          xtemp <- x1[rim,]

          xtemp.list[[i-1]] <- xtemp                                                                                                 #CD
          subdata <- indata[rim,]

          if(family[[1]]=="binomial"){
               getthisglm <- stats::glm(formula=genformula, family=family, data=subdata, weights=bin.wts, singular.ok=TRUE)          #   glm
          }
          else if(family[[1]]=="Gamma"){
               getthisglm <- stats::glm(formula=genformula, family=family, data=subdata, singular.ok=TRUE)                           #   glm
          }
          else if(family[[1]]=="poisson"){
               getthisglm <- stats::glm(formula=genformula, family=family, data=subdata, singular.ok=TRUE)                           #   glm
          }
          else{
               stop("family name not recognized")
          }

                                                if(begin.diagnose <= 55){print(paste(spacer,"Section 55",sep=" "),quote=FALSE);Hmisc::prn(getthisglm)              }

          ##################################################
          # Form weighted hat matrix for leverage and Cook #
          ##################################################
          W <- getthisglm$weights
          lengthW <- length(W)
          matrixW <- c(W[1],rep(0,lengthW))
          for(bb in 2:lengthW){
               matrixW <-c(matrixW,W[bb],rep(0,lengthW))
          } 
          matrixW <- matrixW[1:(lengthW*lengthW)]
          matrixW <- matrix(matrixW,lengthW,lengthW)            # weights in a diagonal matrix  
          transtemp <- t(xtemp)
          cross <- transtemp %*% matrixW %*% xtemp
          uueigen <- unlist(eigen(cross)[1])
#          if(any(uueigen < 10^(-12))){
          if(any(uueigen < 10^(-8))){
               print(uueigen - 10^(-12))
               print("", quote=FALSE)
               print("The X'X matrix of the model is singular. Increase the value", quote=FALSE)
               print("of n.obs.per.level until all eigenvalues (below) are > 10^-12 .", quote=FALSE)
               print("", quote=FALSE)
               print(uueigen)
               print("", quote=FALSE)
               stop(" ")
          }
          crossinv <- solve(cross)                              # inverts a matrix
          hat <- xtemp %*% crossinv %*% transtemp      
          sqrtW <- sqrt(matrixW)
          glmHat[[i-1]] <- sqrtW %*% hat %*% sqrtW        

                                           if(begin.diagnose <= 60){print(paste(spacer,"Section 60",sep=" "),quote=FALSE);Hmisc::prn(cross);Hmisc::prn(glmHat[[i-1]])       }

          #
          ###################################################################
          # Extract values from glm object with current set of observations #
          ###################################################################

          # Coefficients #
          betahat <- getthisglm$coefficients
          # Allowinig singular results means we have to replace NAs with 0 #
          if(any(is.na(betahat))){
               nabetahat <- is.na(betahat)
               betahat[nabetahat] <- 0
          }
          param.est[,i-1] <- betahat
          model.resids <- getthisglm$residuals          # NEEDED????
          betahatset[i-1,] <- betahat
          medaugx <- matrix(1:dimx1, nrow=dimx1, ncol=2, byrow=FALSE)                 # initialize medaugx
          medaugx[,2] <- 0 
          if(i > p) s.2[i-1] <- sum(model.resids * model.resids)/(i-p)

                                            if(begin.diagnose <= 64){print(paste(spacer,"Section 64",sep=" "),quote=FALSE);Hmisc::prn(model.resids);
                                                      Hmisc::prn(s.2[i-1])       }

          # Deviances #
          glmdeviance[i-1] <- getthisglm$deviance
          if(estimate.phi){
               glmphi[i-1] <- glmdeviance[i-1]/(length(rim)-p)
          }
          else{
               glmphi[i-1] <- 1
          }
         # Null Deviance #
          glmnulldeviance[[i-1]] <- getthisglm$null.deviance

          # AIC #
          glmaic[[i-1]] <- getthisglm$aic
          #
          if(i > mstart){
               t.setbase <- c(summary(getthisglm)$coefficients[,3], rep(0,ncoeffs))
               t.set[,i-1] <- t.setbase[1:ncoeffs]
          }
                                   if(begin.diagnose <= 65){print(paste(spacer,"Section 65",sep=" "),quote=FALSE);Hmisc::prn(glmphi[i-1]);
                                               Hmisc::prn(glmnulldeviance[[i-1]]);Hmisc::prn(t.set[,i-1])       }

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
                   
                                                   if(begin.diagnose <= 70){print(paste(spacer,"Section 70",sep=" "),quote=FALSE);Hmisc::prn(leverage)       }
          #
          #####################################################################
          # Determine next set of observations to include in set              #
          # Calculate the deviance for each observation using the current     #
          # model by means of the predict function                            #
          #####################################################################
          if(i <= dimx1){
               predicted <- stats::predict(object=getthisglm, newdata=indata, type="response")                                  # predict proportion
               disparity <- y1 - predicted
               storeDevianceResiduals[[i-1]] <- disparity
               disparity <- disparity^2
               candidates <- data.frame(data,disparity)                                # data has Observation column
               candidates <- candidates[order(candidates$disparity),]
               rows.in.model[[i]] <- sort(candidates[1:i,1])

                                             if(begin.diagnose <= 75){print(paste(spacer,"Section 75",sep=" "),quote=FALSE);Hmisc::prn(candidates)       }

               #
               ######################################
               # Deviance Residuals (unstudentized) #
               ######################################
               getdeviance <- rep(0,dimx1)
               for(tt in 1:dimx1){
                    if(family[[1]]=="binomial"){
                         obs <- y1[tt]
                         pred <- predicted[tt]
                         this.weight <- indata$bin.wts[tt]
                         if(obs==pred){
                             getdeviance[tt] <- 0
                         }
                         else{
                              getdeviance[tt] <- devianceCode(y1[tt]*indata$bin.wts[tt], predicted[tt]*indata$bin.wts[tt], 
                                    indata$bin.wts[tt],fam=family[[1]])                                    # devianceCode
                         }
                    }
                    else if(family[[1]]=="Gamma"){
                         getdeviance[tt] <- devianceCode(y1[tt], predicted[tt], fam=family[[1]])           # devianceCode
                    }
                    else if(family[[1]]=="poisson"){
                         getdeviance[tt] <- devianceCode(y1[tt], predicted[tt], fam=family[[1]])           # devianceCode
                    }
                    else{
                         stop("family name not recognized")
                    }
               }
               obs <- 1:dimx1
               mm1 <- i-1
               uudev <- data.frame(mm1,"A",obs,getdeviance)          # initially label all residuals as "A"ugmented
               names(uudev) <- c("m", "AorD", "obs", "deviance")
               uudev$AorD[rim] <- "D"                                # label the deviance residuals "D"
               ADdevres[[i-1]] <- uudev 
          }
          else{
               getdeviance <- rep(0,dimx1)
               if(family[[1]]=="binomial"){
                    getthisglm <- stats::glm(formula=genformula, family=family, data=indata, weights=bin.wts, singular.ok=TRUE)
               }
               else{
                    getthisglm <- stats::glm(formula=genformula, family=family, data=indata, singular.ok=TRUE)
               }
               predicted <- stats::predict(object=getthisglm, newdata=indata, type="response")                  # predict proportion
               for(tt in 1:dimx1){
                    if(family[[1]]=="binomial"){
                         obs <- y1[tt]
                         pred <- predicted[tt]
                         this.weight <- indata$bin.wts[tt]
                         if(obs==pred){
                              getdeviance[tt] <- 0
                         }
                         else{
                              getdeviance[tt] <- devianceCode(y1[tt]*indata$bin.wts[tt], predicted[tt]*indata$bin.wts[tt], 
                                   indata$bin.wts[tt],fam=family[[1]])                             # devianceCode
                         }
                    }
                    else if(family[[1]]=="Gamma"){
                         getdeviance[tt] <- devianceCode(y1[tt], predicted[tt], fam=family[[1]])   # devianceCode
                    }
                    else if(family[[1]]=="poisson"){
                         getdeviance[tt] <- devianceCode(y1[tt], predicted[tt], fam=family[[1]])   # devianceCode
                    }
                    else{
                         stop("family name not recognized")
                    }
               }                                                                # should use prior.weights? need to collect them in file
               obs <- 1:dimx1
               uudev <- data.frame(dimx1,"A",obs,getdeviance)
               names(uudev) <- c("m", "AorD", "obs", "deviance")
               uudev$AorD[rim] <- "D"
               ADdevres[[i-1]] <- uudev 
          }
     }            # i in mstart ...                                                    END OF STEP 2
     #
     ###############################################
     # Reformat deviance residuals as a data frame #
     ###############################################
     ADdevres.df <- NULL
     for(i in mstart:(dimx1+1)){
          ADdevres.df <- rbind(ADdevres.df, ADdevres[[i-1]])
     }
     #
     ############################################################
     # Set up deviance residuals for plotting                   #
     ############################################################
#     denom <- sqrt(glmphi * (1 - diag(glmHat[[dimx1]])))                            # see p 199 of A&R  
#     obs <- data$Observation
#     for(ii in mstart:(dimx1)){
#          yyy <- data.frame(obs, storeDevianceResiduals[[ii-1]], denom, storeDevianceResiduals[[ii-1]]/denom)
#          yyy <- data.frame(i, yyy)
#          glmstudresids <- rbind(glmstudresids,yyy)
#     }   #  ii
#     names(glmstudresids) <- c("m","obs","devresid")

     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m", "Observation", "leverage"))
     param.est <- as.data.frame(t(param.est))
     names(param.est) <- coeffnames
     m <- 1:dimx1
     param.est <- cbind(m,param.est)
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
     #                                                                                                                   #Cook
if(F){
     ##########################
     # Modified Cook distance #
     ##########################
     nms <- dim(param.est)[1]
     param.est.current <- param.est[-1,]
     param.est.prev <- param.est[-nms,]
     param.diff <- param.est.prev - param.est.current

                                                  if(begin.diagnose <= 85){print(paste(spacer,"Section 85",sep=" "),quote=FALSE);Hmisc::prn(param.est);
                                                            Hmisc::prn(param.est.prev);
                                                            Hmisc::prn(param.est.current);Hmisc::prn(param.diff);Hmisc::prn(xtemp.list)      }    
     getw <- getthisglm$weights
     for(ij in mstart:dimx1){
         aa <- param.diff[ij-1,]
         aa <- as.numeric(aa[,-1])
         aa <- matrix(aa,nrow=1)
         bb <- xtemp.list[[ij]]
                                     if(begin.diagnose <= 88){Hmisc::prn(paste(spacer,"Section 88",sep=" "));Hmisc::prn(aa);Hmisc::prn(bb)    }

         lastHat <- glmHat[[dimx1]]
         lastHat <- lastHat * my.identity(dimx1)    # select diagonal elements
         www <- aa %*% t(bb)

                                 if(begin.diagnose <= 89){print(paste(spacer,"Section 89",sep=" "),quote=FALSE);Hmisc::prn(www)       }

         modCook[i] <- (www %*% t(www))/(p*s.2[ij])
     }       #    ij
}     # if(F)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_glm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     if(!estimate.phi)glmphi <- "Phi not estimated, instead set = 1"

     listout <- list(
          "Step 1 observation numbers"=        SOON,
          "Rows in stage"=                     rows.in.model,
          "Family"=                            family,   
          "Number of model parameters"=        p, 
          "Fixed parameter estimates"=         param.est,          
#         "Studentized deviance residuals"=    glmstudresids,
          "Residual deviance"=                 glmdeviance,  
          "Null deviance"=                     glmnulldeviance,
          "PhiHat"=                            glmphi,
          "Deviance residuals and augments"=   ADdevres.df, 
          "AIC"=                               glmaic,       
          "Leverage"=                          leverage[-1,],
#         "Modified Cook distance"=            modCook, 
          "t statistics"=                      t.set,
          "Call"=                              MC)

     return(listout)
}
