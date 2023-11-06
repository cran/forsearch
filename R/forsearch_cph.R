#' @export
forsearch_cph <-
function(formula.elements, event.time, status, x, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL,
    ties = c("efron", "breslow", "exact"), unblinded=TRUE, begin.diagnose= 100, verbose=TRUE)
{
     #                                           forsearch_cph    
     #
     # VALUE    List of datasets and statistics for plotting in forward search procedure to diagnose coxph observations. Currently, fits only right-censored data.
     #
     # INPUT    formula.elements   The individual names of independent variables in a formula object (omit tilde '~'). 
     #          event.time         Vector of event times, censored or not 
     #          status             Vector indicator of event or censoring: 1 = event, 0 = censored, same length as event.time 
     #          x                  Data frame of independent variables, with number of rows = length of event.time, first column is Observation. Factor variables must
     #                                    be defined in advance. Does not include response time or status.
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          n.obs.per.level    Number of observations per level of (possibly crossed) factor levels
     #          skip.step1         NULL or a list, each element of which is a vector of integers for observations from 1 subgroup to be included in Step 1
     #          ties               Character string specifying the method for handling ties in event time.  See survival::coxph() for definitions.
     #          unblinded          TRUE permits printing of ultimate coxph analysis, as specified above
     #          begin.diagnose     Numeric. Indicates where in code to begin printing diagnostics. 0 prints all; 100 prints none.
     #          verbose            Logical. TRUE causes printing of function ID before and after running.
     #
     # diagnose Step0: 0 - 10    Step1: 13 - 17    Step2: 20-55   >2:  60 - 
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
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      forsearch_cph             "     # used for begin.diagnose prints
     options(warn=-1)
     on.exit(options(warn=0))
     #
# STEP 0
     print("", quote=FALSE)
     print("BEGINNING STEP 0", quote=FALSE)
     print("", quote=FALSE)

     #########################################################
     # Ensure that first independent variable is Observation #
     #########################################################
     uu <- dimnames(x)
     if(uu[[2]][1] != "Observation"){Hmisc::prn(uu);stop("First column of x must be 'Observation'")}
     #
     ####################################################
     # Ensure that Observation conforms to standard 1:N #
     ####################################################
     xObs <- x[,1]
     dimobs <- length(xObs)
     for(i in 1:dimobs){
         if(xObs[i] != i)stop("Observation must be 1 through N (the number of rows) in row order")     }
     #
     #################################
     # Ensure that x is a data frame #
     #################################
     if(!is.data.frame(x))stop("x must be a data frame")
     #
     ######################################################
     # Ensure that there are some uncensored observations #
     ######################################################
     if(sum(status)==0)stop("There must be some uncensored observations (ie, status=1)")
     #
     ###############################################
     # Construct formula.rhs from formula.elements #
     ###############################################
     formula.rhs <- paste(formula.elements, collapse=" + ")
     #
     ##########################################################
     # Reduce x to those variables listed in formula.elements #
     ##########################################################
     namesx <- names(x)
     nnames <- length(namesx)
     namesformula.x <- rep(NA,nnames)   
     for(j in 1:nnames){
          namesformula.x[j] <- charmatch(formula.elements[j], namesx, nomatch = -99)
     }  # for j
     savenames <- namesformula.x != -99
     namesformula.x <- c(1, namesformula.x[savenames])
     x <- x[,namesformula.x]
     #
     ####################################################
     # Ensure that there is replication in the database #
     ####################################################
     varlist <- variablelist(x, verbose=TRUE)
     if(length(varlist)==dim(x)[1]){
          print("",quote=FALSE)
          print("There is no replication in this dataset. All observations are defined as a combination of factors.", quote=FALSE)
          stop("This function does not support such datasets.")
     }
     #
     #######################################################################
     # Print the structure of the analysis that will be done on these data #
     # unless the treatment group is blinded                               #
     #######################################################################
     options(warn = -1)
     on.exit(options(warn = 0))
     randevent <- round(100*stats::runif(length(event.time)),0)
     nulldata <- data.frame(x, event=randevent, status)
     xform <- paste("survival::Surv(time=randevent, event=status, type='right')", formula.rhs, sep=" ~ ")                                      # Surv
     formulatest <- stats::as.formula(xform)

                                               if(begin.diagnose <=3){ print(paste(spacer,"Section 3",sep=" "),quote=FALSE);Hmisc::prn(randevent);Hmisc::prn(status);
                                                 Hmisc::prn(formulatest);Hmisc::prn(utils::head(nulldata));Hmisc::prn(utils::tail(nulldata))   }

     # coxph is run without a preliminary test of variable redundancy
     coxrandAlldata <- survival::coxph(formula=formulatest, data=nulldata, singular.ok=TRUE, x=TRUE, y=TRUE)                          # coxph 
     ncols <- dim(nulldata)[2]
     print("", quote = FALSE)

     coeffnames <- names(coxrandAlldata$coefficients)
     if(unblinded){
          print("**************************************", quote = FALSE)
          print("", quote = FALSE)
          print("Check the following ANOVA paradigm for source of variation and associated degrees of freedom (Df)", quote=FALSE)
          print("Actual event times have been replaced by random numbers here.", quote=FALSE)
          print("", quote = FALSE)
          print(stats::anova(coxrandAlldata))
          print("", quote = FALSE)
          print("All categorical variables must be defined to be factors in the database, not in the formula.", quote=FALSE)
          print("", quote = FALSE)
          print("**************************************", quote = FALSE)
     }
     y1 <- event.time
     #
     ######################################################
     # Run a paradigm test for parallelism on this output #
     # to get structure and names                         #
     ######################################################
     cox.zph.h1 <- survival::cox.zph(coxrandAlldata)
     dn.cox.zph.h1 <- dimnames(cox.zph.h1[[1]])               # get both sets of dimnames for first element of cox.zph.h1
     proprtnlty <- matrix(0,nrow=length(dn.cox.zph.h1[[1]]), ncol=dim(x)[1])
     #
     ######################################
     # Check for factor status of dataset #
     ######################################
     dimdata <- dim(x)[2]
     ufactor <- rep(TRUE, dimdata)
     for(m in 1:dimdata) ufactor[m] <- is.factor(x[,m])
     yesfactor <- any(ufactor)     
     alldata <- data.frame(x, event=event.time, status)                                                                            # alldata = [x,event,status]
     xform <- paste("survival::Surv(time=event.time, event=status, type='right')", formula.rhs, sep=" ~ ")                         # Surv
     formulaStep2 <- stats::as.formula(xform)
                                               if(begin.diagnose <=5){ print(paste(spacer,"Section 5",sep=" "),quote=FALSE);Hmisc::prn(alldata);Hmisc::prn(yesfactor)  }
     ################################################
     # run coxph just to get number of coefficients #
     # without preliminary test of redundancy       #
     ################################################
     coxAlldata <- survival::coxph(formula=formulaStep2, data=alldata, singular.ok=TRUE, x=TRUE, y=TRUE)                          # coxph 

     p <- rnk <- length(coeffnames)   # These values of p and rnk may result in failed Wald test; don't change value, but increment by 2 in their use below
     nopl <- n.obs.per.level 
     if(yesfactor){
          # there are factors in the dataset
          #############################################################################################
          # Ensure that if there are factors, that there are uncensored observations at each sublevel #
          #############################################################################################
          warncensor <- NULL
          ss77 <- variablelist(datadf = alldata[alldata$status==1,], verbose=FALSE)                                                                 # variablelist
          for(nss77 in 1:length(ss77)){
               uu <- ss77[[nss77]][,1]
               uu <- alldata[uu,]
               uu <- uu$status
               if(sum(uu)==0)warncensor <- "Not NULL; issue warning"
          }
          if(!is.null(warncensor))print("Some subsets among the (crossed) factors have no uncensored observations and may be inestimable") 
# Note use of rnk+2 below
          pickm <- picksome(subsetlist=ss77, nobs=dim(x)[1], initial.sample = initial.sample, n.obs.per.level=nopl, rank=rnk+2, 
                      verbose = FALSE)    # picksome
          dimpickm <- dim(pickm)[2]  
     }    # yesfactor
     #
     #################################################
     # Set up all output and intermediate structures #
     # except that proprtnlty is definced above      #
     #################################################
     randset <- 1:initial.sample
     dimx <- dim(x)
     dimx1 <- dimx[1]

     x1 <- x
     OBS <- x[,1]

     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result
     result <- matrix(0, nrow=dimx1, ncol=2)
     rows.in.model <- vector("list", dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=p, ncol=dimx1)    # p is OK here
     WaldTest <- rep(0,dimx1)
     LL <- matrix(0,nrow=dimx1,ncol=3)
#     t.set <- param.est                             # same dimensions

     xtemp.list <- vector("list",dimx1)
#     modCook <- rep(0,dimx1)
#     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)
                                               if(begin.diagnose <=10) {print(paste(spacer,"Section 10",sep=" "),quote=FALSE);Hmisc::prn(rows.in.model);Hmisc::prn(param.est);
                                                          Hmisc::prn(WaldTest);Hmisc::prn(utils::head(LL));Hmisc::prn(utils::tail(LL));Hmisc::prn(utils::head(leverage));
                                                          Hmisc::prn(utils::tail(leverage))    }
     #
######################################################################################################################################################
# Step 1
     #################################################################################
     # Step 1 of the procedure follows.                                              #
     # Use only the uncensored observatons for this step. Create new observation nums#
     # By creating an index, randomly reorder the rows of matrix Z = (X,y).          #    Don't need this
     # Within each element of zlist, calculate the several estimates of b from the   #
     # first p of these rows and calculate the median of the residuals in each set.  #
     #################################################################################
     uncensored <- alldata[alldata$status==1,]
     names(uncensored)[1] <- "OriginalObs"

     dimunobs <- dim(uncensored)[1]
     Observation <- 1:dimunobs
     uncensored <- data.frame(Observation, uncensored)
     medaugx <- matrix(1, nrow=initial.sample, ncol=2)                       #  median of augx
     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(x=1:dimx1, size=dimx1)                        #    sample permutation
          if(yesfactor)zlist[[i]][1:dimpickm,1] <- pickm[i,]
     }      #   i
                                               if(begin.diagnose <=13) {print(paste(spacer,"Section 13",sep=" "),quote=FALSE);Hmisc::prn(zlist[[initial.sample]])    }
     if(is.null(skip.step1)){
          print("ENTERING STEP 1", quote=FALSE)
# increment p here
          inner.rank <- p + 2           #  ????

          mstart <- inner.rank + 1                       # this is just a default
          formlm <-paste("event",formula.rhs,sep=" ~ ") 
          formlm <- stats::as.formula(formlm)

          namunCen <- names(uncensored)
          nnames <- length(namunCen)
          ycol <- c(1:nnames)[namunCen=="event"]
#   increment p here
          firstrim <- aStep1(yesfactor, uncensored, inner.rank=p+2, 
                      initial.sample, formula=formlm, ycol, nopl=n.obs.per.level)          #aStep1
     ############################################################################################################
     # firstrim is the first set of observations according to the "observation" number of the uncensored subset #
     # Delete uncensored$Observation and rename uncensored$OriginalObs back to Observation                      #
     ############################################################################################################
#   increment p here
         pskip <- length(firstrim) 
         rows.in.model[[pskip]] <- firstrim
         redun.rim <- firstrim
         SOON <- firstrim
     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          pskip <- length(skip.step1)
          rows.in.model[[pskip]] <- skip.step1
                                               if(begin.diagnose <=14) {print(paste(spacer,"Section 14",sep=" "),quote=FALSE);Hmisc::prn(rows.in.model);
                                              Hmisc::prn(alldata[skip.step1,])    }
          redun.rim <- skip.step1
          SOON <- skip.step1
     }       # skipping step 1
     mstart <- pskip + 1
     # 
     #  Place error diagnosis here to avoid looping the output #
                                               if(begin.diagnose <=15) {print(paste(spacer,"Section 15",sep=" "),quote=FALSE);Hmisc::prn(yesfactor);
                                                   Hmisc::prn(utils::head(uncensored));Hmisc::prn(utils::tail(uncensored));Hmisc::prn(p);Hmisc::prn(initial.sample);
                                                   Hmisc::prn(SOON)    }
     #
     #####################################################
     # test for redundant variables before running coxph #
     #####################################################
     thisdata <- alldata[redun.rim,]
     #
######################################################################################################################################################
# Step 2
     #############################################
     # Step 2 of the procedure follows.          #
     # Adding observations to the initial set    #
     # First, get all stats for Step 1 subset.   #
     # Then calculate next subset at end of loop #
     #############################################
     print("ENTERING STEP 2", quote=FALSE)
     betahatset <- matrix(0,nrow=dimx1,ncol=p)
     for(i in mstart:(dimx1)){                  # mstart is the step after the original p obs entered
          ############################################################################################################
          # Set up coxph arguments in order to extract statistics along the way. For formula to work,need to subset. #
          ############################################################################################################
          rim <- rows.in.model[[i-1]]             # picks up rows for previous step
          tempdata <- alldata[rim,]
          innertempdata <- length(formula.elements)
          xtempdata <- tempdata[,1+(1:innertempdata)]
          xform <- paste("survival::Surv(time=event, event=status, type='right')", formula.rhs, sep=" ~ ")                                      # Surv
          formulanewrim <- stats::as.formula(xform)
                                               if(begin.diagnose <=20) {print(paste(spacer,"Section 20",sep=" "),quote=FALSE);Hmisc::prn(i);Hmisc::prn(tempdata[,1])    }
          coxph.out07 <- survival::coxph(formula=formulanewrim, data=alldata, subset=tempdata[,1], 
                ties="efron", singular.ok=TRUE, model=TRUE, x=TRUE, y=TRUE)                                                                     # coxph

                                               if(begin.diagnose <=25) {print(paste(spacer,"Section 25",sep=" "),quote=FALSE);Hmisc::prn(i);Hmisc::prn(coxph.out07)    }

          param.est[,i] <- coxph.out07$coefficients
          thiszph <- survival::cox.zph(coxph.out07)[[1]]
          proprtnlty[,i] <- thiszph[,3]
          WaldTest[i] <- coxph.out07$wald.test
          LLout <- coxph.out07$loglik
          LL[i,] <- c(i,LLout) 
          Zlatest <- alldata[rim,]
                                               if(begin.diagnose <=30) {print(paste(spacer,"Section 30",sep=" "),quote=FALSE);Hmisc::prn(param.est[,i]);
                                                                       Hmisc::prn(WaldTest[i]);Hmisc::prn(LL[i,]);Hmisc::prn(Zlatest)     }
          #
          ############
          # Leverage #
          ############
          xtemp <- coxph.out07$x
          ###########################################
          # Set limits on determinant of hat matrix #
          ###########################################
          myM1 <- my.Machine.double.eps  <- 2.22e-14       #-16
#          myM2 <- my.Machine.double.xmax <- 2.23e-306      # -308
          myM3 <- my.Machine.double.xmin <- 1.79e306       # 308

          uuu <- prod(eigen(t(xtemp) %*% xtemp)[[1]])
                                               if(begin.diagnose <=35) {print(paste(spacer,"Section 35",sep=" "),quote=FALSE);Hmisc::prn(xtemp);
                                                                       Hmisc::prn(uuu)     }
          if(uuu > myM1 & uuu < myM3){
               crossinv <- solve(t(xtemp) %*% xtemp)
                                               if(begin.diagnose <=37) {print(paste(spacer,"Section 37",sep=" "),quote=FALSE);Hmisc::prn(xtemp);Hmisc::prn(crossinv)    }
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
          }    # if uuu
                                               if(begin.diagnose <=40) {print(paste(spacer,"Section 40",sep=" "),quote=FALSE);Hmisc::prn(leverage)       }
 
         #
          #####################################################################
          # Select next subset.                                               #
          #####################################################################
          df1 <- data.frame(x, event.time, status)
          if(i <= dimx1){
                                     if(begin.diagnose <=50) {print(paste(spacer,"Section 50",sep=" "),quote=FALSE);Hmisc::prn(i);Hmisc::prn(utils::head(df1));
                                             Hmisc::prn(utils::tail(df1));Hmisc::prn(rim);Hmisc::prn(event.time);Hmisc::prn(status);Hmisc::prn(formula.rhs);
                                          print("End of Section 50")   }

               rows.in.model[[i]] <- cStep2(df1, rim, formula.rhs)                                                           # cStep2  


#                                               if(begin.diagnose <=55) {print(paste(spacer,"Section 55",sep=" "),quote=FALSE);Hmisc::prn(i);Hmisc::prn(rows.in.model[[i]]);
#                      Hmisc::prn(df1[rim,])    }
          }
     }            # i in mstart ...
                                               if(begin.diagnose <=57) {print(paste(spacer,"Section 57",sep=" "),quote=FALSE);Hmisc::prn(utils::head(leverage));
                                                   Hmisc::prn(utils::tail(leverage))    }



     #
     ###################
     # Clean up output #
     ###################
     m <- 1:dimx1

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
     #
     ########################################
     # Create likelihood ratio test from LL #
     ########################################
     m <- LL[,1]
     Difference <- LL[,2] - LL[,3]
     LRT <- exp(Difference)
                                               if(begin.diagnose <=60) {print(paste(spacer,"Section 60",sep=" "),quote=FALSE);Hmisc::prn(Difference);Hmisc::prn(LRT)    }
     LRT <- data.frame(m,LRT)
     LRT <- LRT[-1,]                # additional row removals below

     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m","Observation","leverage"))

                                               if(begin.diagnose <=70) {print(paste(spacer,"Section 70",sep=" "),quote=FALSE);Hmisc::prn(dimnames(leverage[[2]]))    }
     # 

     listout <- list(
          "Step 1 observation numbers"=        SOON,
          "Rows in stage"=                     rows.in.model,
#          "Standardized residuals"=            residuals, 
          "Number of model parameters"=        p, 
#          "Sigma"=                             sigma, 
          "Fixed parameter estimates"=         param.est[-c(1:p),], 
          "Proportionality Test"=              proprtnlty,
          "Wald Test"=                         WaldTest[-c(1:p),],
          "LogLikelihood"=                     LL,
          "Likelihood ratio test"=             LRT[-c(1:p),],
#          "s^2"=                               s.2, 
          "Leverage"=                          leverage[-1,], 
#          "Modified Cook distance"=            modCook, 
#          "t statistics"=                      t.set,
          "Call"=                              MC)

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_cph", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(listout)
}
