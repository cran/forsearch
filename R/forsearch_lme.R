#' @export
forsearch_lme <-
function(fixedform, alldata, randomform, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL, 
    unblinded=TRUE, begin.diagnose= 100, verbose=TRUE)
{
    MC <- match.call()
     if(verbose) {
          print("", quote=FALSE)
          print("Running forsearch_lme", quote=FALSE)
          print("", quote=FALSE)
          print(date(), quote=FALSE)
          print("", quote=FALSE)
          print("Call:", quote=FALSE)
          print(MC, quote=FALSE)
          print("", quote=FALSE)
     }
     #                                           forsearch_lme 
     #
     # VALUE    List of datasets and statistics for plotting in forward search procedure to diagnose lme observations.  These include: scaled residuals, s^2, leverage,
     #                 estimated coefficients, variance of random effects. 
     #
     # INPUT    fixedform          2-sided formula for fixed effects
     #          alldata            Data frame, first column of which must be "Observation".    
     #          randomform         1-sided formula for random effects
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          n.obs.per.level    Number of observations per level of (possibly crossed) factor levels
     #          skip.step1         NULL or a list, each element of which is a vector of integers for observations from 1 subgroup to be included in Step 1
     #
     #          unblinded          TRUE permits printing of ultimate lme analysis, as specified above
     #          begin.diagnose     Numeric. Indicates where in code to begin printing diagnostics. 0 prints all; 100 prints none.
     #          verbose            Logical. TRUE causes printing of function ID before and after running.
     # 

     #   begin.diagnose      Step 0: 1 - 19        Step 1: 20     -     49     Step 2:    50 - 59           Extraction:     81 - 
     #                                                aStep1:  31 - 39                    aStep2:  60 - 80

# Step 0
     spacer <- "                                            forsearch_lme        "
     options(warn=-1)
     on.exit(options(warn=0))
     #
     #####################################################################################################
     # Create all data files that will be used in this function                                          #
     # Includes fixdat.df data frame with all variables, including inner factor and grouping variables   #
     # Includes fixdat.list, all continuous variables within each level of factor and grouping variables # 
     # Includes datacont.df, all continuous variables only                                               #          
     #####################################################################################################
     nalldata1 <- dim(alldata)[1]
     nalldata2 <- dim(alldata)[2]
     alldataNames <- names(alldata)
     nopl <- n.obs.per.level
     #####################################################################################################
     # Create all files needed below:                                                                    #
     #    fixdat.df,   a data frame with factor level indicator and containing all independent variables #
     #    fixdat.list, a list with the same variables as alldata but no factor variables by factor level #
     #    fixdatrank.df, a data frame for determining rank                                               #
     #####################################################################################################

     ############################################################
     # Check for factor status of alldata and get factor names  #
     ############################################################
     datacontrank <- NULL
     ufactor <- rep(TRUE, nalldata2)
     for(m in 1:nalldata2) ufactor[m] <- is.factor(alldata[,m])
     yesfactor <- any(ufactor)     
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
     innerfact[1:2] <- FALSE
     for(i in 1:nlevfixdat){
          fixdat.sub <- fixdat.df[fixdat.df$holdISG==ufixdatISG[i],]
          fixdat.sub <- fixdat.sub[,!innerfact]
          fixdat.list[[i]] <- fixdat.sub
     }                                                                                          # fixdat.list   Step 2
     names(fixdat.list) <- ufixdatISG

                                               if(begin.diagnose <=19){ print(paste(spacer,"Section 19",sep=" "),quote=FALSE);
                                                   Hmisc::prn(fixdat.list);Hmisc::prn(names(fixdat.list))  }

     nfacts <- nlevfixdat                                        # used in extracting statistics
     #
     ######################################
     # collapse this list to a data frame #
     ######################################
     uncenrank.df <- NULL
     for(i in 1:nlevfixdat){
          uncenrank.df <- rbind(uncenrank.df, fixdat.list[[i]])
     }
     #
     #######################################################################
     # Remove factor variables and get rank within continuous observations #
     #######################################################################
     unnames <- names(uncenrank.df)
     unnames <- unnames[-length(unnames)]    # remove factor level indicator
     if(length(unnames)==2){          # contains only Observation and response
           datacontrank <- 1
           mini.df <- fixdat.df[,-(1:2)]
     }
     else{
          unnames <- unnames[-c(1:2)]
          unnames <- paste(unnames, collapse=" + ")
          fixed.lhs <- formula.tools::lhs(fixedform)
          rank.form <- paste(fixed.lhs, unnames, sep=" ~ ")
          rank.form <- stats::as.formula(rank.form)
          lm4rank <- stats::lm(formula=rank.form, data=alldata)                    #    lm for rank
          datacontrank <- lm4rank$rank
          mini.df <- fixdat.df[,-(1:2)]
     }
###########################################################################################################################################
# Step 1
     rows.in.model <- vector("list", nalldata1)
     LLL <- vector("list", nalldata1)
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result

                                               if(begin.diagnose <=22) {print(paste(spacer,"Section 22",sep=" "),quote=FALSE)   }

     if(is.null(skip.step1)){
          print("ENTERING STEP 1", quote=FALSE)
          inner.rank <- datacontrank

                                               if(begin.diagnose <=23){ print(paste(spacer,"Section 23",sep=" "),quote=FALSE);Hmisc::prn(inner.rank);
                                                 Hmisc::prn(nopl);
                                                 Hmisc::prn(utils::head(uncenrank.df));Hmisc::prn(utils::tail(uncenrank.df))   }

          if(datacontrank==1){
               fixed.lhs <- formula.tools::lhs(fixedform)
               xform2 <- paste(fixed.lhs, "1", sep=" ~ ")
          }else{
               formula.rhs2 <- unnames
               if(length(formula.rhs2)>1)formula.rhs2 <- paste("",formula.rhs2, collapse=" + ")
               xform2 <- paste(fixed.lhs, formula.rhs2, sep=" ~ ")
          }
          formulacont <- stats::as.formula(xform2)
          firstrim2 <- aStep1(yesfactor=TRUE, df1=fixdat.list, inner.rank=datacontrank, initial.sample, formula=formulacont, ycol=2, 
                    nopl=nopl, b.d=begin.diagnose)                                                                                     # aStep1 
          firstrim <- firstrim2[[1]]
          nfirstrim <- length(firstrim)
          rows.in.model[[nfirstrim]] <- firstrim                        # save list with pool and individual subsets
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
########################################################################################################################################################################
# Step 2
     print("BEGINNING STEP 2", quote=FALSE)
     print("", quote=FALSE)

     fixedformform <- stats::formula(fixedform)

     zzzz <- bStep2(f2=fixedform, dfa2=fixdat.df, randm2=randomform, ms=mstart,  
                    finalm=rows.in.model, fbg=fixdat.list, b.d=begin.diagnose, rnk2=datacontrank, ss=skip.step1, LLL=LME)                 # bStep2

     rows.in.set <- zzzz[[1]]
     LME <- zzzz[[2]]
     rows.in.set[[nalldata1]] <- 1:nalldata1
     LME[[nalldata1]] <- nlme::lme(fixedform, alldata, randomform)

########################################################################################################################################################################
# Extracting summary statistics
     print("", quote=FALSE)
     print("BEGINNING INTERMEDIATE RESULTS EXTRACTION", quote=FALSE)
     ###############################################################
     # Show analysis of full dataset
     ###############################################################
     print("", quote=FALSE)
     zholdlm <- nlme::lme(fixed=fixedform, data=fixdat.df, random=randomform)                                  #    lme
     zholdcoeffs <- zholdlm$coefficients
     rnkfixed <- length(zholdcoeffs[[1]])                                                # these go into param.est

     if(unblinded){
          print("The assumed analysis for these data will be as follows:", quote=FALSE)
          print("", quote=FALSE)
          print(zholdlm)
          print("", quote=FALSE)
     }
     #################################################################################################
     # Retrieve metadata for plotting from each row of rows.in.set                                   #
     # Set up files to hold residuals, dimensions, sigma, fixed parameter estimates,                 #
     #    random parameter estimates, leverage, Cook distance and list for storage of xtemp matrices #                                                                                            #
     #################################################################################################
     nrowsdf1 <- nalldata1
     dddd <- zholdlm$dims
     ss77 <- variablelist(datadf=alldata, prank=datacontrank)
     hold.residuals <- matrix(0,nrowsdf1,nrowsdf1)          # for standardized residuals across all observations
     hold.subset.residuals <- rep(0,nrowsdf1)               # for use in Cook distance   
     hold.dims <- vector("list", nrowsdf1)
     hold.sigma <- rep(-999, nrowsdf1)   
     param.est <- matrix(0,nrow=rnkfixed, ncol=nrowsdf1)
     t.set <- param.est
     hold.coeffs.fixed <- vector("list",nrowsdf1)

     leverage <- matrix(-999,nrow=1,ncol=3)
     xtemp.list <- vector("list",nrowsdf1)
     modCook <- rep(0,nrowsdf1-1)
     hold.summary.stats <- matrix(0,nrow=nrowsdf1,3)        # unnamed columns: AIC, BIC, log likelihood
     #
     ###################################################################################
     # Set up for holding anova p values and standard deviations for random components #
     ###################################################################################
     templme <- zholdlm
     tempanova <- stats::anova(templme)
     dnAVlme <- dimnames(tempanova)
     dnAVlme <- dnAVlme[[1]]
     anova.pvalues <- matrix(0, nrow=length(dnAVlme), ncol = dim(alldata)[1])
     #
     VC <- nlme::VarCorr(templme)
     VCnames <- dimnames(VC)[[1]]
     VC2 <- VC[,2]
     hold.coeffs.random <- matrix(0,nrow=nrowsdf1, ncol=length(VC2))       # won't need to transpose    
     #########################################################
     # Set up for extraction of leverage and Cook's distance #
     #########################################################
     levAlldata <- stats::lm(formula=fixedform, data=alldata, x=TRUE, y=TRUE)                                  # lm
     x1 <- levAlldata$x

     for(dd in mstart:nalldata1){                                                                   #    for loop starts here
          rim <- rows.in.set[[dd]]             # picks up row numbers for a set of observations
          Zlatest <- fixdat.df[rim,]  
          ############################################################################
          # Extract indep vars of the subset for use in leverage and Cook's distance #
          # lm function allows collection of x matrix and y vector                   #
          ############################################################################
          xtemp <- x1[rim,]
          xtemp.list[[dd]] <- xtemp 
          transtemp <- t(xtemp)
          cross <- transtemp %*% xtemp

          solvecross <- FALSE
          if(abs(det(cross)) > 0){
               crossinv <- solve(cross)
               solvecross <- TRUE          }
          #
          ################################################################
          # Capture results of lme analysis for this set of observations #
          ################################################################
          zholdlme <- LME[[dd]]
          VC <- nlme::VarCorr(zholdlme)
          VC2 <- as.numeric(VC[,2])
          hold.coeffs.random[dd,] <- VC2              # insert extraction into row dd
          temprand <- zholdlme$coefficients[[2]]
          listnames <- names(temprand)
          nlistnames <- length(listnames)
          outnames <- NULL

          for(jj in 1:nlistnames){
               colnames2 <- dimnames(temprand[[jj]])[2][[1]]
               stu <- NULL
               for(kk in 1:length(colnames2)){
                    stu <- c(stu, paste(listnames[jj], colnames2[kk], sep=" - ")  )
               }          #   kk
               outnames <- c(outnames, stu)
          }
          #
          resids <- zholdlme$residuals
          hold.coeffs.fixed[[dd]] <- zholdlme$coefficients[[1]]
          temprand<-zholdlme$coefficients[[2]]

                                               if(begin.diagnose <= 88){print(paste(spacer,"Section 88         dd=",dd));Hmisc::prn(zholdlme)}

          xbar <- mean(c(temprand[[1]]))
          devs <- c(temprand[[1]])-xbar
          sumsq <- mean(devs^2)
          RMS <- sqrt(sumsq)
          logLik <- zholdlme$logLik
          AIC <- summary(zholdlme)$AIC
          BIC <- summary(zholdlme)$BIC
          hold.summary.stats[dd,] <- c(AIC, BIC, logLik)
          param.est[,dd] <- c(zholdlme$coefficients[[1]])                        # same as holdcoeffs.fixed[[i-1]]
          t.set[,dd] <- summary(zholdlme)$tTable[,4]
          hold.dims[[dd]] <- zholdlme$dims
          hold.sigma[dd] <- zholdlme$sigma
          newrows <- NULL
          potential <- NULL
          hold.subset.residuals[dd] <- sum(resids^2)/(dd-rnkfixed)
          #
          #########################################################
          #          "Standardized residuals"=     hold.residuals #
          #########################################################
          td1 <- dim(fixdat.df)[1]
          errors <- rep(-999, td1)
          y1 <- fixdat.df[,2]
          for(j in 1:td1){
               errors[j] <- y1[j] - sum(hold.coeffs.fixed[[dd]] * x1[j,])
          }             #   j
          hold.residuals[,dd] <- errors
          #
          ############################################################
          #          "Fixed parameter estimates"=          param.est #
          ############################################################
          param.est[,dd] <- c(zholdlme$coefficients[[1]])     
                                                              if(begin.diagnose <= 90){print(paste(spacer,"Section 90         dd=",dd));Hmisc::prn(param.est[,dd])}
          #
          ###########################################################
          # ANOVA test of fixed effects               anova.pvalues #
          ###########################################################
          AVlme <- stats::anova(zholdlme)
          AVlmeps <- AVlme[,4]                     # no need to remove p value (NA) for residuals; not in AVlme
          anova.pvalues[,dd] <- AVlmeps
                                                              if(begin.diagnose <= 92){print(paste(spacer,"Section 92         dd=",dd));
                                                                     Hmisc::prn(anova.pvalues[,dd])}
          #
          ############################################################
          #           Leverage=                        leverage[-1,] #
          ############################################################
          thisleverage <- 1
          if(is.matrix(x1)){
               for(j in 1:dd){
                                                              if(begin.diagnose <= 93){print(paste(spacer,"Section 93         dd=",dd));Hmisc::prn(j);
                                                                    #Hmisc::prn(Zlatest);
                                                                    Hmisc::prn(dim(x1)[2]);Hmisc::prn(matrix(xtemp[j,],nrow=1));Hmisc::prn(det(crossinv))       }

                    if(solvecross){
                         Zlatest2 <- data.frame(Zlatest)
                         if(dim(x1)[2]==1){
                               thisleverage <- c(c(matrix(xtemp[j],nrow=1) %*% crossinv %*% matrix(xtemp[j],ncol=1)))
                               thisleverage <- c(dd,Zlatest2[j,1],thisleverage)
                          }
                          else{
                               thisleverage <- c(c(matrix(xtemp[j,],nrow=1) %*% crossinv %*% matrix(xtemp[j,],ncol=1)))
                               thisleverage <- c(dd,Zlatest2[j,1],thisleverage)
                          }
                          leverage <- rbind(leverage,thisleverage)
                    }        # if solvecross
                    else{
                        leverage <- rbind(leverage,thisleverage)     # this will append the last set value of thisleverage
                    }
               }   # j 1:dd
          }        #   x1 is matrix
          #
          ########################################################
          #          "t statistics"=                       t.set #
          ########################################################
          t.set[,dd] <- summary(zholdlme)$tTable[,4]
          ###########################################################
          #          "Fit statistics"=           hold.summary.stats #
          ###########################################################
          AIC <- summary(zholdlme)$AIC
          BIC <- summary(zholdlme)$BIC
          logLik <- zholdlme$logLik
          hold.summary.stats[dd,] <- c(AIC, BIC, logLik)
     }                                                                  #   dd     for loop ends here                                                    

     param.est <- as.data.frame(t(param.est))
     names(param.est) <- paste("b",1:rnkfixed, sep="")
     m <- 1:nalldata1
     xxparam <- t(param.est)
     param.est <- as.data.frame(xxparam)
     names(param.est) <- paste("b",1:rnkfixed, sep="")
     param.est <- rbind(m,param.est)
     #
     hold.summary.stats <- data.frame(m, hold.summary.stats)                          # hold.summary.stats defined
     names(hold.summary.stats) <- c("m", "AIC", "BIC", "logLik")
     # 
     t.set <- as.data.frame(t(t.set))
     nz <- names(zholdcoeffs[[1]])
     dimnames(t.set)[[2]] <- nz
     t.set <- cbind(m,t.set)
     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m","Observation","leverage"))
     #
     ############################################
     # Sigma used to standardize residuals      #  
     # Average squares of last col of residuals #
     # and take square root                     #
     ############################################
     s.2 <- sum((hold.residuals[,nrowsdf1])^2)
     sigma.squared <- s.2/(nrowsdf1-rnkfixed)
     sigma <- sqrt(sigma.squared)
     hold.residuals <- hold.residuals/sigma
     #
     ############################
     # Modified Cook distance #
     ##########################
     param.est <- t(as.matrix(param.est))
     nms <- dim(param.est)[1]
     param.est.current <- param.est[-1,]
     param.est.prev <- param.est[-nms,]
     param.diff <- param.est.prev - param.est.current
     for(i in mstart:nalldata1){
          aa <- param.diff[i-1,]
          aa <- as.numeric(aa[-1])
          aa <- matrix(aa,nrow=1)
          bb <- as.matrix(xtemp.list[[i]])
                                                              if(begin.diagnose <= 95){print(paste(spacer,"Section 95"));
                                                                   Hmisc::prn(aa);Hmisc::prn(bb)}
          www <- aa %*% t(bb)
          modCook[i-1] <- (www %*% t(www))/(rnkfixed * hold.subset.residuals[i-1])

                                                              if(begin.diagnose <= 96){print(paste(spacer,"Section 96"));
                                                                     Hmisc::prn(hold.subset.residuals[i-1])}
     }     #  for i
     #
     #######################################
     # Clean up files written to workspace #
     #######################################
     if(is.null(skip.step1)){
          rm(list=c("zzzz", "Zlatest", "newcontrol"), pos=1)
     }
     else{
          rm(list=c("Zlatest","newcontrol"), pos=1)
     }
     hold.coeffs.random <- as.data.frame(hold.coeffs.random)
     names(hold.coeffs.random) <- VCnames
     dimhcr <- dim(hold.coeffs.random)[1]
     testrow <- hold.coeffs.random[dimhcr,]
     if(any(is.na(testrow))){
          # eliminate columns with any NA
          index <- 1:length(testrow)
          index <- index[!is.na(testrow)]
          hold.coeffs.random <- hold.coeffs.random[,index]
     }
     #
     anova.pvalues <- as.data.frame(t(anova.pvalues))
     names(anova.pvalues) <- dnAVlme
     if(dnAVlme[1]=="(Intercept)"){
          anova.pvalues <- anova.pvalues[,-1]           # remove (Intercept)
     }
     m <- 1:(dim(anova.pvalues)[1])
     anova.pvalues <- cbind(m, anova.pvalues)
     listout <- list(
          "Number of observations in Step 1"=   mstart-1,
          "Step 1 observation numbers"=         SOON,
          "Rows by subgroup"=                   fixdat.list,
          "Rows in stage"=                      rows.in.set,
           Sigma=                               sigma,
          "Standardized residuals"=             hold.residuals,            
          "Fixed parameter estimates"=          param.est,
          "Random parameter estimates"=         hold.coeffs.random,
           Leverage=                            leverage[-1,],
           ANOVA=                               anova.pvalues,
          "Modified Cook distance"=             modCook,
           Dims=                                zholdlme$dims,
          "t statistics"=                       t.set,
          "Fit statistics"=                     hold.summary.stats,
           Call=                                MC )
    #
     if(verbose) {
          print("", quote=FALSE)
          print("Finished running forsearch_lme", quote=FALSE)
          print("", quote=FALSE)
          print(date(), quote=FALSE)
          print("", quote=FALSE)
     }

     return(listout)
}
