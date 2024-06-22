#' @export
forsearch_glm <-
function(
          initial.sample=1000, 
          response.cols,
          indep.cols,
          family,

          formula=NULL,
          binomialrhs=NULL,
          formula.cont.rhs,      

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

     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXX           forsearch_glm  "

     ############################
     # Identity matrix function #  Used only in Cook calculation
     ############################
     my.identity <- function(n){
          matrix(rep(c(1, rep(0, n)), times = n)[1:(n * n)], nrow = n, ncol = n, byrow = T)
     }
     ########################################################################################
     # Helper function for calculaton of deviance code for each observation based on family #
     ########################################################################################
     devianceCode <- function (obs, pred, ni=NULL, fam) 
     {
          if(obs/pred < 0) { out <- 0  }
          else{
               if(fam=="binomial"){
                    if(obs==0){ out <- ni*log(ni/(ni-pred))     }      # expects whole numbers, not proportions
                    else if(obs==ni){ out <- obs*log(obs/pred) }
                    else{ out <- obs*log(obs/pred) + (ni-obs)*log((ni-obs)/(ni-pred)) }
               }
               else if(fam=="Gamma"){ out <- -log(obs/pred) + (obs-pred)/pred  }
               else if(fam=="poisson"){ if(obs==0){ out <- pred }
                    else{ out <- obs*log(obs/pred) - obs + pred  }
               }
               else if(fam=="exponential"){ out <- -log(obs/pred) + obs*pred - 1 }
               else{
                    stop("family name not recognized")
               }
               if(is.na(out)){ out <- 0 }
               else{
                    if(abs(out)<= 10^(-12)){ out <- 0 }
                    else{ out <- 2 * out; vv <- obs - pred; vv <- vv/abs(vv)            #  1 or -1
                         out <- sqrt(out)*vv }
               }
          }
          return(out)
     }                                   # end of devianceCode function
     #
     ##########################################################################################################
     # MAIN FUNCTION STARTS HERE #
     options(warn = -1)
     on.exit(options(warn = 0))

     alldata <- data
     nobs <- dim(data)[1]
     nalldata2 <- dim(data)[2] 
                                 if(begin.diagnose <= 1){print(paste(spacer,"Section 1",sep=" "),quote=FALSE);
                                     Hmisc::prn(utils::head(data));Hmisc::prn(formula.cont.rhs);Hmisc::prn(nobs);
                                     Hmisc::prn(response.cols);Hmisc::prn(indep.cols)       }

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
          indata <- data.frame(data,proportion1,bin.wts)        # indata has proportion1 and bin.wts
          bin.wts <- NULL                                       # insist that these come from indata, not argument
          genformula <- formula(paste("proportion1", binomialrhs, sep=" ~ ")) 
     }      # binomial  
     else{
          name.response <- uu[response.cols]
          inresponse <- data[,response.cols]
          indata <- data
          genformula <- formula
    }       # not binomial

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
     #
if(F){
     ################################################################
     # Check for constructed variables in formula, ie, use of I()   #
     # First convert formula to a vector of character pairs. Then   #
     # recode the I( letters as  I(A) and the test accordingly.     #      
     # Determine whether any of these is 'I(A)'.  If so, count them #
     ################################################################
     nAsIs <- 0
     charform <- as.character(formula.cont.rhs)
     nstrs <- nchar(charform)
     output <- rep("S", nstrs)
     for(i in 1:nstrs){
          output[i] <- substr(charform,start=i, stop=i+1)
     }
     formpairs <- output
     formpairs <- paste(formpairs, "A)", sep="")
     jj <- "I(A)" %in% formpairs
     if(jj){
          kk <- grep(as.character("I(A)"), as.character(formpairs))
          nAsIs <- length(kk)
     }

                                 if(begin.diagnose <= 3){print("", quote = FALSE);print(paste(spacer,"Section 3",sep=" "),quote=FALSE);
                                      Hmisc::prn(charform);Hmisc::prn(formpairs);Hmisc::prn(jj);Hmisc::prn(nAsIs)       }
#prn(nAsIs)
 nAsIs <- 0
}


     ########################################
     # Check for factor status of dataset   #
     # Get rank of analysis without factors #
     ########################################
     datacontrank <- NULL
     ndatacols <- dim(data)[2]
     ufactor <- rep(TRUE, ndatacols)
     for(m in 1:ndatacols){
          ufactor[m] <- is.factor(data[,m])
     }
     yesfactor <- any(ufactor)
     lmAlldata <- stats::lm(formula=genformula, data, singular.ok=TRUE, x=TRUE, y=TRUE)                                    # lm
     ncoeffs <- length(lmAlldata$coefficients)
     rnk <- lmAlldata$rank
                                 if(begin.diagnose <= 4){print(paste(spacer,"Section 4",sep=" "),quote=FALSE);
                                      Hmisc::prn(ncoeffs);Hmisc::prn(rnk)       }
 
     nopl <- n.obs.per.level
     ##########################################################
     # Ensure that if yesfactor, formula.cont.rhs is not null #
     ##########################################################
     if(yesfactor & is.null(formula.cont.rhs))stop("There is a factor variable in the data, but formula.cont.rhs is null.")

     x1 <- z1 <- lmAlldata$x
     y1 <- lmAlldata$y
     dimdata <- dim(data)
     #
     ############################################################################
     # Resampling algorithm for Step 1 of forward search     p31                #
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
     deviance.matrix <- matrix(1:dimx1, nrow=dimx1, ncol=4, byrow=FALSE)              # step 1

     rows.in.model <- vector("list",dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=ncoeffs, ncol=dimx1)
     t.set <- param.est                                        # t statistics
     xtemp.list <- vector("list",dimx1)          
     modCook <- rep(0,dimx1)                                                                                      
     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)

                                 if(begin.diagnose <= 6){print(paste(spacer,"Section 6",sep=" "),quote=FALSE);
                                      Hmisc::prn(utils::head(param.est));Hmisc::prn(utils::tail(param.est));Hmisc::prn(s.2)       }
     #

     #####################################################################################################
     # Create three files needed below:                                                                  #
     #    fixdat.df,   a data frame with factor level indicator and containing all independent variables #
     #    fixdat.list, a list with the same variables as alldata but no factor variables by factor level #
     #    datacont, a data frame with no factor variables for determining the rank of regression         #
     #####################################################################################################

     ############################################################
     # Check for factor status of alldata and get factor names  #
     ############################################################
     datacontrank <- NULL
     ufactor <- rep(TRUE, nalldata2)
     for(m in 1:nalldata2) ufactor[m] <- is.factor(alldata[,m])
     alldataNames <- names(alldata)
     yesfactor <- any(ufactor)     

     # Add factor subset indicator #
     factorNames <- alldataNames[ufactor]

     ############################################
     # Add the factor grouping code to the data #
     # If binomial, add proportion1 to the data #
     # Define ycol in any case                  #
     ############################################
     fixdat.df <- alldata
     isfactor <- rep(TRUE,nalldata2)
     for(nn in 1:nalldata2){
          isfactor[nn] <- is.factor(fixdat.df[,nn])
     }
     justfactors <- fixdat.df[isfactor]
     holdISG <- apply(justfactors, 1, paste,collapse="/")
     holdISG <- paste("_",holdISG, sep="")
     fixdat.df <- data.frame(fixdat.df, holdISG)                                                 # fixdat.df    building
     if(family[[1]]=="binomial"){
          fixdat.df <- data.frame(fixdat.df, proportion1); ycol <- dim(fixdat.df)[2]
          responseName <- "proportion1"
     }
     else{responseName <- names(data)[response.cols]     }

                                               if(begin.diagnose <=17){ print(paste(spacer,"Section 17",sep=" "),quote=FALSE);
                                                      Hmisc::prn(utils::head(fixdat.df));Hmisc::prn(utils::tail(fixdat.df));
                                                      Hmisc::prn(dim(fixdat.df));Hmisc::prn(responseName)   }

     ##########################################################################################
     # Create a list by factor subset levels, each of which does not contain factor variables #
     ##########################################################################################
     ufixdatISG <- unique(fixdat.df$holdISG)
     nlevfixdat <- length(ufixdatISG)
     fixdat.list <- vector("list", nlevfixdat)

     innerfact <- isfactor
     innerfact[1:3] <- FALSE
     for(i in 1:nlevfixdat){
          fixdat.sub <- fixdat.df[fixdat.df$holdISG==ufixdatISG[i],]
          fixdat.sub <- fixdat.sub[,!innerfact]
          fixdat.list[[i]] <- fixdat.sub
     }                                                                           # fixdat.list   Step 2
     names(fixdat.list) <- ufixdatISG
     nfacts <- nlevfixdat                                        # used in extracting statistics
     #

     ########################################
     # Remove factor variables and get rank #
     # Call this  inner.rnk                 #
     ###################################### #
     this.form2 <- NULL
     lmnofactor <- NULL
     datacont <- data[,!ufactor]
     ufactornames <- names(datacont)          # includes Observation and response
     ufactornames <- ufactornames[-1]         # includes response

     if(length(ufactornames)==1){             # there are no continuous independent variables
          if(family[[1]] == "binomial"){
               this.form2 <-paste("proportion1", "1", sep=" ~ ")
               this.form2 <- stats::formula(this.form2)
          }
          else{
               nameresponse <- names(data)[response.cols]
               formulacont <- paste(nameresponse, "1", sep=" ~ ")
               formulacont <- stats::formula(formulacont)
          }
          inner.rnk <- 1                      # just the intercept
     }
     else{                                    # there are some continuous variables; maybe all of them
          if(family[[1]] == "binomial"){
               this.form2 <-paste("proportion1", formula.cont.rhs, sep=" ~ ")
               this.form2 <- stats::formula(this.form2)
          }
          else{
               nameresponse <- names(data)[response.cols]
               this.form2 <- paste(nameresponse, formula.cont.rhs, sep=" ~ ")
               this.form2 <- stats::formula(this.form2)
         }
          lmnofactor <- stats::lm(formula = this.form2, data=datacont)
          inner.rnk <- lmnofactor$rank
     }
                                 if(begin.diagnose <= 18){print(paste(spacer,"Section 18",sep=" "),quote=FALSE);
                                      Hmisc::prn(ufactornames);Hmisc::prn(this.form2);Hmisc::prn(lmnofactor);
                                      Hmisc::prn(inner.rnk)       }
     #
# stop("xnorm")
##############################################################################################################################
# Step 1
     if(is.null(skip.step1)){
          print("BEGINNING STEP 1", quote=FALSE)                   
          print("",quote=FALSE)

          ##################################################
          # Create formula string, and send this to dStep1 #
          ##################################################
          this.form <- paste(formula.tools::lhs(genformula),formula.cont.rhs, sep = " ~ ") 
          if(yesfactor){
               yk <- names(fixdat.list[[1]])
               yk <- match(responseName,yk)
               firstrim <- dStep1(yesfactor,df1=fixdat.list, inner.rank=inner.rnk, initial.sample=initial.sample, 
                         formuladStep=this.form, fam=family, ycol=yk, nopl=nopl, b.d=begin.diagnose)                      # dStep1
         }
          else{
               yk <- names(fixdat.df)
               yk <- match(responseName,yk)
               firstrim <- dStep1(yesfactor, df1=fixdat.df, inner.rank=inner.rnk, initial.sample=initial.sample,
                         formuladStep=this.form, fam=family, ycol=yk, nopl=nopl, b.d=begin.diagnose)                      # dStep1)
          }
          alength <- length(firstrim)
          rows.in.model[[alength]] <- firstrim
          mstart <- alength + 1
          SOON <- firstrim
     }       #   skip.step1 is NULL ?
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          print("",quote=FALSE)
          lenskip <- length(skip.step1)
          rows.in.model[[lenskip]] <- skip.step1
          SOON <- skip.step1
          mstart <- lenskip + 1
     }   #   if not skip step 1
                                                  if(begin.diagnose <= 49){print(paste(spacer,"Section 49",sep=" "),quote=FALSE);
                                                        Hmisc::prn(SOON);Hmisc::prn(mstart)       }
     #
# stop("end of Step 1")
###################################################################################################################################################
# Step2

     print("BEGINNING STEP 2", quote=FALSE)
     print("",quote=FALSE)

     yk <- names(fixdat.df)
     yk <- match(responseName, yk)
     zzzz <- dStep2(f2=genformula, dfa2=fixdat.df, ms=mstart, finalm=rows.in.model, fbg=fixdat.list, b.d=begin.diagnose, 
                rnk2=inner.rnk, ycol=yk, fam=family)                                                                       # dStep2

     rows.in.model <- zzzz[[1]]
     rows.in.model[[nobs]] <- 1:nobs
                                                 if(begin.diagnose <= 52){print(paste(spacer,"Section 52",sep=" "),quote=FALSE);
                                                        Hmisc::prn(rows.in.model)       }
     #
#prn(rows.in.model)
# stop("before extracting intermediate stats")
###################################################################################################################################################
     print("BEGINNING EXTRACTION OF INTERMEDIATE STATISTICS", quote=FALSE)
     print("",quote=FALSE)

     #####################################################
     # Initialize storage matrices and vectors           #
     #####################################################
     betahatset <- matrix(0,nrow=dimx1,ncol=ncoeffs)
     glmdeviance <- rep(0,dimx1)
     glmphi <- rep(0,dimx1)
     glmnulldeviance <- glmdeviance
     glmaic <- rep(0,dimx1)
     ADdevres <- vector("list",dimx1)
     storeDevianceResiduals <- vector("list",dimx1)         # initially saved without studentization as vector
     glmstudresids <- NULL                                  # after studentization 
     glmHat <- vector("list", dimx1)                        # will store Hat matrices for leverage, Cook, and (last one) studentization of errors
     #

                                             if(begin.diagnose <= 53){print(paste(spacer,"Section 53",sep=" "),quote=FALSE);
                                                  Hmisc::prn(glmdeviance)       }
# Insert result for last level of fooResult
     fooResult <- zzzz[[2]]
          if(family[[1]]=="binomial"){
               fooResult[[nobs]] <- stats::glm(formula=genformula, family=family, data=alldata, y=TRUE, x=TRUE)          #   glm
          }
          else if(family[[1]]=="Gamma"){
               fooResult[[nobs]] <- stats::glm(formula=genformula, family=family, data=alldata,  y=TRUE, x=TRUE)          #   glm
          }
          else if(family[[1]]=="poisson"){
               fooResult[[nobs]] <- stats::glm(formula=genformula, family=family, data=alldata,  y=TRUE, x=TRUE)          #   glm
          }

     for(i in mstart:nobs){                                
          getthisglm <- fooResult[[i]]



          rim <- rows.in.model[[i]] 
          Zlatest <- Z[rim,]
                                          if(begin.diagnose <= 53){print(paste(spacer,"Section 53",sep=" "),quote=FALSE);
                                                Hmisc::prn(Zlatest)       }

          xtemp <- getthisglm$x
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

          an.error.occurred <- FALSE
          tryCatch( crossinv <- solve(cross), error=function(e) {an.error.occurred <<- TRUE})
          if(an.error.occurred){
               print("Error occurred in matrx inversion.", quote=FALSE)
               glmHat[[i-1]] <- NA
          }
          else{
               hat <- xtemp %*% crossinv %*% transtemp      
               sqrtW <- sqrt(matrixW)
               glmHat[[i-1]] <- sqrtW %*% hat %*% sqrtW        
          }
                                           if(begin.diagnose <= 60){print(paste(spacer,"Section 60",sep=" "),quote=FALSE);
                                                 Hmisc::prn(cross);Hmisc::prn(glmHat[[i-1]])       }

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
#          if(i > p) s.2[i-1] <- sum(model.resids * model.resids)/(i-p)

                                            if(begin.diagnose <= 64){print(paste(spacer,"Section 64",sep=" "),quote=FALSE);
                                                            Hmisc::prn(model.resids)       }

          # Deviances #
          glmdeviance[i-1] <- getthisglm$deviance
          if(estimate.phi){
               glmphi[i-1] <- glmdeviance[i-1]/(length(rim)-rnk)
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
                                   if(begin.diagnose <= 65){print(paste(spacer,"Section 65",sep=" "),quote=FALSE);
                                         Hmisc::prn(glmphi[i-1]);Hmisc::prn(glmnulldeviance[[i-1]]);
                                         Hmisc::prn(t.set[,i-1])       }

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
                   
                                        if(begin.diagnose <= 70){print(paste(spacer,"Section 70",sep=" "),quote=FALSE);
                                              Hmisc::prn(leverage)       }
          #
          #####################################################################
          # Determine next set of observations to include in set              #
          # Calculate the deviance for each observation using the current     #
          # model by means of the predict function                            #
          #####################################################################
          if(i < dimx1){
               predicted <- stats::predict(object=getthisglm, newdata=indata, type="response")                                  # predict proportion
               disparity <- y1 - predicted
               storeDevianceResiduals[[i-1]] <- disparity
               disparity <- disparity^2
               candidates <- data.frame(data,disparity)                                # data has Observation column
               candidates <- candidates[order(candidates$disparity),]
               rows.in.model[[i]] <- sort(candidates[1:i,1])

                                             if(begin.diagnose <= 75){print(paste(spacer,"Section 75",sep=" "),quote=FALSE);
                                                   Hmisc::prn(candidates)       }

               #
               ######################################
               # Deviance Residuals (unstudentized) #
               ######################################
               getdeviance <- rep(0,dimx1)
               for(tt in 1:(dimx1-1)){
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
               }    # tt
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
               }                                                                # should use prior.weights? need to collect them in file
               obs <- 1:dimx1
               uudev <- data.frame(dimx1,"A",obs,getdeviance)
               names(uudev) <- c("m", "AorD", "obs", "deviance")
               uudev$AorD[rim] <- "D"
               ADdevres[[i-1]] <- uudev 
          }
     }            # i in mstart ...  
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
     coeffnames <- names(lmAlldata$coefficients)
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
     sigma <- sqrt(sum(medaugx[,2])/(dimx1-inner.rnk))
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

                                         if(begin.diagnose <= 85){print(paste(spacer,"Section 85",sep=" "),quote=FALSE);
                                               Hmisc::prn(param.est);Hmisc::prn(param.est.prev);
                                               Hmisc::prn(param.est.current);Hmisc::prn(param.diff);Hmisc::prn(xtemp.list)      }    
     getw <- getthisglm$weights
     for(ij in mstart:dimx1){
         aa <- param.diff[ij-1,]
         aa <- as.numeric(aa[,-1])
         aa <- matrix(aa,nrow=1)
         bb <- xtemp.list[[ij]]
                                     if(begin.diagnose <= 88){Hmisc::prn(paste(spacer,"Section 88",sep=" "),quote=FALSE);
                                            Hmisc::prn(aa);Hmisc::prn(bb)    }

         lastHat <- glmHat[[dimx1]]
         lastHat <- lastHat * my.identity(dimx1)    # select diagonal elements
         www <- aa %*% t(bb)

                                 if(begin.diagnose <= 89){print(paste(spacer,"Section 89",sep=" "),quote=FALSE);
                                          Hmisc::prn(www)       }

         modCook[i] <- (www %*% t(www))/(inner.rnk*s.2[ij])
     }       #    ij
}     # if(F)


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
          "Number of model parameters"=        rnk, 
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
