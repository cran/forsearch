#' @export
forsearch_cph <-
function(alldata, formula.rhs, nofactform, initial.sample=1000, skip.step1=NULL,
    ties = "efron", maxdisturb=.01, proportion=TRUE, wiggle=1, unblinded=TRUE, begin.diagnose= 100, verbose=TRUE)
{
     #                                           forsearch_cph    
     #
     # VALUE    List of datasets and statistics for plotting in forward search procedure to diagnose coxph observations. Currently, fits only right-censored data.
     #
     # INPUT    alldata            Data frame whose first 3 columns are Observation, event.time, and status, and whose last columns are independent variables
     #          formula.rhs        Right hand side of formula (omitting ~)
     #          nofactform         Right hand side of formula (omitting ~ and factor variables)
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          skip.step1         NULL or a list, each element of which is a vector of integers for observations from 1 subgroup to be included in Step 1
     #          ties               Character string specifying the method for handling ties in event time.  See survival::coxph() for definitions.
     #          proportion         TRUE causes running of survival::cox.zph on each stage
     #          wiggle
     #          unblinded          TRUE permits printing of ultimate coxph analysis, as specified above
     #          begin.diagnose     Numeric. Indicates where in code to begin printing diagnostics. 0 prints all; 100 prints none.
     #          verbose            Logical. TRUE causes printing of function ID before and after running.
     #
     #
     # DETAILS  No need to identify variable for which proportionality would be assessed: proportionality is assessed for EACH independent variable if for any
     #
     MC <- match.call()
     if(verbose) {
          print("", quote=FALSE)
          print("Running forsearch_cph", quote=FALSE)
          print("", quote=FALSE)
          print(date(), quote=FALSE)
          print("", quote=FALSE)
          print("Call:", quote=FALSE)
          print(MC, quote=FALSE)
          print("", quote=FALSE)
     }
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      forsearch_cph   "     # used for begin.diagnose prints
     options(warn=-1)
     on.exit(options(warn=0))
     #
     nalldata1 <- dim(alldata)[1]
     nalldata2 <- dim(alldata)[2]

     statusOK <- sum(alldata[,3])
     nindepvars <- nalldata2 - 3

     OBS <- alldata[,1]
     wiggle <- wiggle*stats::runif(nalldata1)/100
     alldata <- cbind(alldata, wiggle)
     nalldata2 <- dim(alldata)[2]
     alldataNames <- names(alldata)

     #########################################################
     # Ensure that first independent variable is Observation #
     # Ensure that Observation conforms to standard 1:N      #
     # Ensure that there are some uncensored observations    #
     # Ensure that status is only 1's and 0's                #
     #########################################################
     uu <- names(alldata)
     if(uu[1] != "Observation"){Hmisc::prn(uu);stop("First column of alldata must be 'Observation'")}

     xObs <- alldata[,1]
     dimobs <- length(xObs)
     for(i in 1:dimobs){if(xObs[i] != i)stop("Observation must be 1 through N (the number of rows) in row order")     }

     if(sum(alldata[,3])==0)stop("There must be some uncensored observations (ie, status=1)")
     if(statusOK + sum(alldata[,3]==0) - nalldata1 != 0)stop("status must consist only of 0's and 1's")
     #
     ############################################################################
     # Add small disturbance to each event.time in order to ensure against ties #
     ############################################################################
     disturb <-  -maxdisturb + stats::runif(nalldata1, min=0, max=2*maxdisturb)
     alldata[,2] <- alldata[,2] + disturb
     #
     #######################################################################
     # Print the structure of the analysis that will be done on these data #
     # unless the treatment group is blinded                               #
     #######################################################################
     options(warn = -1)
     on.exit(options(warn = 0))

     randevent <- round(100*stats::runif(nalldata1,0))
     nulldata <- alldata
     nulldata[,2] <- randevent
     xform <- paste("survival::Surv(time=randevent, event=nulldata[,3], type='right')", formula.rhs, sep=" ~ ")                                      # Surv
     formulatest <- stats::as.formula(xform)

                                    if(begin.diagnose <=3){ print(paste(spacer,"Section 3",sep=" "),quote=FALSE);
                                                Hmisc::prn(randevent);Hmisc::prn(utils::head(nulldata));Hmisc::prn(formulatest)   }

     coxrandAlldata <- do.call(survival::coxph, list(formula=formulatest, data=nulldata, model=TRUE, x=TRUE, y=TRUE))                 # coxph using do.call
     ncols <- dim(nulldata)[2]
     print("", quote = FALSE)

     coeffnames <- names(coxrandAlldata$coefficients)
     ncoeffs <- length(coeffnames)                      # used to calculate inner.rank; omits intercept
     if(unblinded){
          print("**************************************", quote = FALSE)
          print("", quote = FALSE)
          print("Check the following ANOVA paradigm for source of variation and associated degrees of freedom (Df)", quote=FALSE)
          print("Actual event times have been replaced by random numbers here.", quote=FALSE)
          print("", quote = FALSE)
          print(stats::anova(coxrandAlldata))
          #       
          print("", quote = FALSE)
          print("", quote = FALSE)
          print("All categorical variables must be defined to be factors in the database, not in the formula.", quote=FALSE)
          print("", quote = FALSE)
          print("**************************************", quote = FALSE)
     }
     #####################################################################################################
     # Create all files needed below:                                                                    #
     #    fixdat.df,   a data frame with factor level indicator and containing all independent variables #
     #    fixdat.list, a list with the same variables as alldata but no factor variables by factor level #
     #    uncensored.list,  a list with the same variables as fixdat.list but no censored observations   #
     #    uncenrank.df, a data frame the same as uncensored.list for determining rank                    #
     #####################################################################################################

     ############################################################
     # Check for factor status of alldata and get factor names  #
     ############################################################
     datacontrank <- NULL
     ufactor <- rep(TRUE, nalldata2)
     for(m in 1:nalldata2) ufactor[m] <- is.factor(alldata[,m])
     yesfactor <- any(ufactor)     
     if(!yesfactor)stop("There must be at least 1 independent variable that is a factor (eg, Treatment)")
     #
     # Add factor subset indicator #

     factorNames <- alldataNames[ufactor]

     ############################################
     # Add the factor grouping code to the data #
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

                                               if(begin.diagnose <=11){ print(paste(spacer,"Section 11",sep=" "),quote=FALSE);
                                                      Hmisc::prn(utils::head(fixdat.df));Hmisc::prn(utils::tail(fixdat.df));Hmisc::prn(dim(fixdat.df))   }

     #########################################
     # Create a list by factor subset levels #
     #########################################
     ufixdatISG <- unique(fixdat.df$holdISG)
     nlevfixdat <- length(ufixdatISG)
     fixdat.list <- vector("list", nlevfixdat)

     innerfact <- isfactor
     innerfact[1:3] <- FALSE
     for(i in 1:nlevfixdat){
          fixdat.sub <- fixdat.df[fixdat.df$holdISG==ufixdatISG[i],]
          fixdat.sub <- fixdat.sub[,!innerfact]
          fixdat.list[[i]] <- fixdat.sub
     }                                                                                          # fixdat.list   Step 2
     names(fixdat.list) <- ufixdatISG
     nfacts <- nlevfixdat                                        # used in extracting statistics
     #
     ####################################################
     # Subset this list to remove censored observations #
     ####################################################
     uncensored.list <- fixdat.list
     for(i in 1:nlevfixdat){
          uuu <- uncensored.list[[i]]
          uncensored.list[[i]] <- uuu[uuu[,3]==1,]
     }                                                                                           # uncensored.list  Step 1
     #
     ######################################
     # Collapse this file to a data frame #
     ######################################
     uncenrank.df <- NULL
     for(i in 1:nlevfixdat){
          uncenrank.df <- rbind(uncenrank.df, uncensored.list[[i]])                             # uncenrank.df      rank within factor levels
     }
     #
     ################################################################################### 
     # Calculate datacontrank = vector of observations to pull from each factor subset #
     ################################################################################### 
     ncoeffsP1 <- ncoeffs + 1    # added because ncoeffs ignores intercept
#prn(ncoeffsP1)
#prn(nlevfixdat)
     model.df.source <- floor(ncoeffsP1/nlevfixdat + .00001)
#prn(model.df.source)
     datacontrank <- rep(model.df.source, nlevfixdat)
#prn(datacontrank)
     additional <- ncoeffsP1 - model.df.source * nlevfixdat
#prn(additional)
     datacontrank[1:additional] <- datacontrank[1:additional] + 1
#prn(datacontrank)
     datacontrank[datacontrank==0] <- 1
#prn(datacontrank)
                                               if(begin.diagnose <=17) {print(paste(spacer,"Section 17",sep=" "),quote=FALSE)   
                                                   Hmisc::prn(ncoeffsP1);Hmisc::prn(model.df.source);Hmisc::prn(additional);
                                                   Hmisc::prn(coeffnames);Hmisc::prn(datacontrank)     }

# stop("prior to entering step 1")
######################################################################################################################################################
# Step 1
     #####################################################
     # Use only the uncensored observatons for this step #   
     # ycol is event time, column 2                      #
     #####################################################
     rows.in.model <- vector("list", nalldata1)
     LLL <- vector("list", nalldata1)
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result
     if(is.null(skip.step1)){
          print("BEGINNING STEP 1", quote=FALSE)
          inner.rank <- NULL
          #####################################################################
          # There is at least one factor variable among independent variables #
          #####################################################################
          firstrim <- cStep1(df1=uncenrank.df, df1.ls=uncensored.list, inner.rank=datacontrank, 
                 initial.sample, cphties=ties, f.e = nofactform, ycol=2, b.d=begin.diagnose)                           # cStep1 

          nfirstrim <- length(firstrim)
          rows.in.model[[nfirstrim]] <- firstrim
          SOON <- firstrim
          mstart <- length(firstrim) + 1
     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)       # no guarantees for user-defined skip.step1
          nfirstrim <- length(skip.step1)
          rows.in.model[[nfirstrim]] <- skip.step1
 
                                              if(begin.diagnose <=55) {print(paste(spacer,"Section 55",sep=" "),quote=FALSE);
                                                   Hmisc::prn(rows.in.model);Hmisc::prn(alldata[skip.step1,])   }

          SOON <- skip.step1
          mstart <- nfirstrim + 1
     }       # skipping step 1
     # 
                                               if(begin.diagnose <=56) {print(paste(spacer,"Section 56",sep=" "),quote=FALSE);
                                                   Hmisc::prn(SOON)   }
# stop("prior to entering step 2")
#####################################################################################################################################
# Step 2
     print("", quote=FALSE)
     print("BEGINNING STEP 2", quote = FALSE)
     print("", quote=FALSE)
     pushaside <- nofactform == "1"
     heresStep2 <- cStep2(fe=formula.rhs, finalm=rows.in.model, rimbs=fixdat.list, dfa2=fixdat.df, onlyfactor=pushaside,  
                      ycol=NULL, mstart=mstart, cphties=ties, rnk=NULL, b.d=begin.diagnose)                           # cStep2

     rows.in.model <- heresStep2[[1]]
     rows.in.model[[nalldata1]] <- 1:nalldata1

     LLL <- heresStep2[[2]]
     rows.in.model[[nalldata1]] <- 1:nalldata1
     adet <- alldata$event.time
     adst <- alldata$status
     xform <- paste("survival::Surv(time=adet, event=adst, type='right')", formula.rhs, sep=" ~ ")                    # Surv
     formulalast <- stats::as.formula(xform)

     coxAlldata <- survival::coxph(formula=formulalast, data=alldata, model=TRUE, x=TRUE, y=TRUE)                  # coxph
     LLL[[nalldata1]] <- coxAlldata
     #
#stop("before intermediate data extraction")
#######################################################################################################################################
     print("BEGINNING INTERMEDIATE RESULTS EXTRACTION", quote=FALSE)
     print("", quote=FALSE)

     #################################################    
     # Set up all output and intermediate structures #
     # except that proprtnlty is definced above      #
     #################################################
     dimx1 <- nalldata1
     param.est <- matrix(0,nrow=ncoeffs, ncol=nalldata1) 
     WaldTest <- rep(0,dimx1)
     LL <- matrix(0,nrow=dimx1,ncol=3)

     if(nfacts > 0){
          cox.zph.h1 <- survival::cox.zph(coxAlldata, global=FALSE)
          dn.cox.zph.h1 <- dimnames(cox.zph.h1[[1]])                               # get both sets of dimnames for first element of cox.zph.h1
          proprtnlty <- matrix(0, nrow=length(dn.cox.zph.h1[[1]]), ncol=nalldata1)      # Counts GLOBAL, which we may drop later
     }
     leverage <- matrix(0,nrow=1,ncol=3)

     for(dd in mstart:nalldata1){                                            #    dd
          coxphout <- LLL[[dd]]
          ####################################
          # Extract the following statistcs: #
          # Estimated coefficients           #
          # Wald test                        #
          # Log likelihood                   #
          # Likelihood ratio test            #
          ####################################
          param.est[,dd] <- coxphout$coefficients
          WaldTest[dd] <- coxphout$wald.test
          LLout <- coxphout$loglik
          LL[dd,] <- c(dd,LLout) 
          if(proportion){
               ##############################
               # Proportionality test       #
               # Skip first stage in Step 2 #
               ##############################
               if(dd > mstart){              
                    an.error.occurred <- FALSE
                    thiszph <- data.frame(0,0,"Skipping proportionality test")
                    if(nfacts > 0 & proportion){
                         tryCatch( {thiszph <- survival::cox.zph(coxphout, global=FALSE)[[1]]   }
                             , error = function(e) {an.error.occurred <<- TRUE})
                    }
                    if(an.error.occurred){
                         messprop <- paste("Skipping proportionality test in stage",dd, sep=" ")
                         print(messprop,quote=FALSE)
                         proprtnlty[,dd] <- NA
                    }
                    else{
                         proprtnlty[,dd] <- thiszph[,3]
                    }
               }        # if dd > mstart
          }     # if proportion
          else{
               proprtnlty <- "Not computed"
          }     
          ############
          # Leverage #
          ############
          an.inversion.failed <- FALSE
          xtemp <- coxphout$x
          tryCatch( expr=crossinvtest <- solve(t(xtemp) %*% xtemp)     
              , error = function(e) {an.inversion.failed <<- TRUE})
          if(an.inversion.failed){
               messprop2 <- paste("Skipping leverage calculation in stage",dd, sep="  ")
               print(messprop2, quote=FALSE)
               leverage <- rbind(leverage,NA)
          }
          else{
               crossinv <- crossinvtest
               Zlatest <- rows.in.model[[dd]]
               Zlatest <- alldata[Zlatest,]
               for(j in 1:(dd-1)){
                   Zlatest2 <- data.frame(Zlatest)
                   if(dim(xtemp)[2]==1){
                        thisleverage <- c(c(matrix(xtemp[j],nrow=1) %*% crossinv %*% matrix(xtemp[j],ncol=1)))
                        thisleverage <- c(dd,Zlatest2[j,1],thisleverage)
                   }else{
                        thisleverage <- c(c(matrix(xtemp[j,],nrow=1) %*% crossinv %*% matrix(xtemp[j,],ncol=1)))
                        thisleverage <- c(dd,Zlatest2[j,1],thisleverage)
                   }
                   leverage <- rbind(leverage,thisleverage)
               }   # j 1:dd-1
          }

     }      # dd
     #
     ###################
     # Clean up output #
     ###################
     m <- 1:nalldata1

     param.est <- as.data.frame(t(param.est))
     names(param.est) <- coeffnames
     param.est <- cbind(m,param.est)                                 #  here we add the m column
     param.est <- param.est[-1,]                  # additional row removals below

     if(proportion){
          proprtnlty <- as.data.frame(t(proprtnlty))
          names(proprtnlty) <- dn.cox.zph.h1[[1]]           # already dimnames
          proprtnlty <- cbind(m, proprtnlty)
     }

     Wald <- WaldTest 
     WaldTest <- data.frame(m,Wald) 
     WaldTest <- WaldTest[-1,]                # additional row removals below

     LL <- as.data.frame(LL)
     names(LL) <- c("m", "Null", "Coefficient")

     ########################################
     # Create likelihood ratio test from LL #
     ########################################
     m <- LL[,1]
     Difference <- LL[,2] - LL[,3]
     LRT <- exp(Difference)
     LRT <- data.frame(m,LRT)
     LRT <- LRT[-1,]                # additional row removals below

     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m","Observation","leverage"))
     p <- datacontrank - 1
     #
     listout <- list(
          "Step 1 observation numbers"=                    SOON,
          "Rows in stage"=                                 rows.in.model,
#          "Standardized residuals"=                       residuals,
          "Number of continuous model parameters"=         p, 
#          "Sigma"=                                        sigma,
          "Fixed parameter estimates"=                     param.est[-c(1:p),], 
          "Proportionality Test"=                          proprtnlty,
          "Wald Test"=                                     WaldTest[-c(1:p),],
          "LogLikelihood"=                                 LL,
          "Likelihood ratio test"=                         LRT[-c(1:p),],
#          "s^2"=                                          s.2,  
          "Leverage"=                                      leverage[-1,], 
#          "Modified Cook distance"=                       modCook,     
#          "t statistics"=                                 t.set,     
          "Call"=                                          MC)

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_cph", quote = FALSE)
          print("", quote = FALSE)
          print("REMINDER: Transfer from Step 1 to Step 2 is artificial and likely", quote=FALSE)
          print("                     shows many changes in search criteria. Ignore these.", quote=FALSE) 
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(listout)
}
