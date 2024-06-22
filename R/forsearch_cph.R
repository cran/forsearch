#' @export
forsearch_cph <-
function(alldata, formula.rhs, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL,
    ties = "efron", proportion=TRUE, unblinded=TRUE, begin.diagnose= 100, verbose=TRUE)
{
     #                                           forsearch_cph    
     #
     # VALUE    List of datasets and statistics for plotting in forward search procedure to diagnose coxph observations. Currently, fits only right-censored data.
     #
     # INPUT    alldata            Data frame whose first 3 columns are Observation, event.time, and status, and whose last columns are independent variables
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          n.obs.per.level    Number of observations per level of (possibly crossed) factor levels
     #          skip.step1         NULL or a list, each element of which is a vector of integers for observations from 1 subgroup to be included in Step 1
     #          ties               Character string specifying the method for handling ties in event time.  See survival::coxph() for definitions.
     #          proportion         TRUE causes running of survival::cox.zph on each stage
     #          unblinded          TRUE permits printing of ultimate coxph analysis, as specified above
     #          begin.diagnose     Numeric. Indicates where in code to begin printing diagnostics. 0 prints all; 100 prints none.
     #          verbose            Logical. TRUE causes printing of function ID before and after running.
     #

     #   begin.diagnose      Step 0: 1 - 19        Step 1: 20     -     49     Step 2:    50 - 59           Extraction:     81 - 
     #                                                aStep1:  31 - 39                    cStep2:  60 - 80

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
# STEP 0
     print("", quote=FALSE)
     print("BEGINNING STEP 0", quote=FALSE)
     print("", quote=FALSE)

     nalldata1 <- dim(alldata)[1]
     nalldata2 <- dim(alldata)[2]

     nopl <- n.obs.per.level 
     statusOK <- sum(alldata[,3])
     nindepvars <- nalldata2 - 3

     OBS <- alldata[,1]

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

                                               if(begin.diagnose <=3){ print(paste(spacer,"Section 3",sep=" "),quote=FALSE);Hmisc::prn(randevent);
                                                 Hmisc::prn(utils::head(nulldata));Hmisc::prn(formulatest)   }

     coxrandAlldata <- do.call(survival::coxph, list(formula=formulatest, data=nulldata, model=TRUE, x=TRUE, y=TRUE))                 # coxph using do.call
     ncols <- dim(nulldata)[2]
     print("", quote = FALSE)

     coeffnames <- names(coxrandAlldata$coefficients)
     ncoeffs <- length(coeffnames)
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

     ################################################################
     # Check for constructed variables in formula, ie, use of I()   #
     # First convert formula to a vector of character pairs. Then   #
     # recode the I( letters as  I(A) and the test accordingly.     #      
     # Determine whether any of these is 'I(A)'.  If so, count them #
     ################################################################
     nAsIs <- 0
     charform <- as.character(xform)
     nstrs <- nchar(charform) - 1
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
                                 if(begin.diagnose <= 6){print("", quote = FALSE);print(paste(spacer,"Section 6",sep=" "),quote=FALSE);
                                      Hmisc::prn(charform);Hmisc::prn(formpairs);Hmisc::prn(jj);Hmisc::prn(nAsIs)       }

# stop("xform")
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

                                               if(begin.diagnose <=17){ print(paste(spacer,"Section 17",sep=" "),quote=FALSE);
                                                      Hmisc::prn(utils::head(fixdat.df));Hmisc::prn(utils::tail(fixdat.df));Hmisc::prn(dim(fixdat.df))   }

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
     #######################################################################
     # Remove factor variables and get rank within uncensored observations #
     # and list of uncensored observation numbers ss77                     #
     #######################################################################
     unnames <- names(uncenrank.df)
     unnames <- unnames[-length(unnames)]    # remove factor level indicator
     if(length(unnames)==3){
          datacontrank <- 1
          mini.df <- uncenrank.df[,-(1:3)]
     }
     else{
          unnames <- unnames[-c(1:3)]
          unnames <- paste(unnames, collapse=" + ")
          rank.form <- paste("event.time", unnames, sep=" ~ ")      
          formulacont <- rank.form
          rank.form <- stats::as.formula(rank.form)
          lm4rank <- stats::lm(formula=rank.form, data=uncenrank.df)                    #    lm for rank
          datacontrank <- lm4rank$rank
          mini.df <- uncenrank.df[,-(1:3)]
     }

######################################################################################################################################################
# Step 1
     #####################################################
     # Use only the uncensored observatons for this step #
     # ycol is event time, column 2                      #
     #####################################################
     rows.in.model <- vector("list", nalldata1)
     LLL <- vector("list", nalldata1)
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result

                                               if(begin.diagnose <=22) {print(paste(spacer,"Section 22",sep=" "),quote=FALSE)   }

     if(is.null(skip.step1)){
          print("BEGINNING STEP 1", quote=FALSE)
          inner.rank <- datacontrank

                                               if(begin.diagnose <=23){ print(paste(spacer,"Section 23",sep=" "),quote=FALSE);Hmisc::prn(inner.rank);
                                                 Hmisc::prn(nopl);
                                                 Hmisc::prn(utils::head(uncenrank.df));Hmisc::prn(utils::tail(uncenrank.df))   }


          if(datacontrank==1){
               xform2 <- "event.time ~ 1"
          }else{
               formula.rhs2 <- unnames
               if(length(formula.rhs2)>1)formula.rhs2 <- paste("",formula.rhs2, collapse=" + ")
               xform2 <- paste("event.time", formula.rhs2, sep=" ~ ")
          }
          formulacont <- stats::as.formula(xform2)
          firstrim2 <- aStep1(yesfactor=TRUE, df1=uncensored.list, inner.rank=datacontrank + nAsIs, initial.sample, formula=formulacont, ycol=2, 
                    nopl=nopl, b.d=begin.diagnose)                                                                                               # aStep1 
          firstrim <- firstrim2[[1]]
          nfirstrim <- length(firstrim)
          rows.in.model[[nfirstrim]] <- firstrim2                        # save list with pool and individual subsets
          SOON <- firstrim
          mstart <- nfirstrim + 1
     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)       # no guarantees for user-defined skip.step1
          nfirstrim <- length(skip.step1)
          rows.in.model[[nfirstrim]] <- skip.step1
 
                                              if(begin.diagnose <=41) {print(paste(spacer,"Section 41",sep=" "),quote=FALSE);
                                                   Hmisc::prn(rows.in.model);Hmisc::prn(alldata[skip.step1,])   }

          SOON <- skip.step1
          mstart <- nfirstrim + 1
     }       # skipping step 1
     # 
#                                               if(begin.diagnose <=42) {print(paste(spacer,"Section 42",sep=" "),quote=FALSE);
 #                                                  Hmisc::prn(thiscph);Hmisc::prn(SOON)   }

#
#####################################################################################################################################
# Step 2
     print("", quote=FALSE)
     print("BEGINNING STEP 2", quote = FALSE)
     print("", quote=FALSE)
     heresStep2 <- cStep2(f.e=formula.rhs, finalm=rows.in.model, dfa2=fixdat.df, ms=mstart,  
                         rnk2=datacontrank + nAsIs, ss=skip.step1, b.d=begin.diagnose)                                             # cStep2

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
# stop("before interim stats extraction")
#UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
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

     proprtnlty <- "Not computed"
     
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
          ##############################
          # Proportionality test       #
          # Skip first stage in Step 2 #
          ##############################
          if(dd > mstart){              
               an.error.occurred <- FALSE
               thiszph <- data.frame(0,0,"Skipping proportionality test")
#prn(thiszph)
#stop("thiszph")
               if(nfacts > 0 & proportion){
                    tryCatch( {thiszph <- survival::cox.zph(coxphout, global=FALSE)[[1]]}
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
          } # dd > mstart
          ############
          # Leverage #
          ############
          xtemp <- coxphout$x
          crossinv <- solve(t(xtemp) %*% xtemp)

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

     proprtnlty <- as.data.frame(t(proprtnlty))
     names(proprtnlty) <- dn.cox.zph.h1[[1]]           # already dimnames
     proprtnlty <- cbind(m, proprtnlty)

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
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(listout)
}
