#' @export
forsearch_lme <-
function(fixedform, nofactform, alldata, randomform, groupname, randfactnames=NULL, initial.sample=1000,
    skip.step1=NULL, unblinded=TRUE, begin.diagnose= 100, incCont=FALSE, verbose=TRUE)
{
     #                                           forsearch_lme 
     #
     # VALUE    List of datasets and statistics for plotting in forward search procedure to 
     #                 diagnose lme observations.  These include: scaled residuals, s^2, leverage,
     #                 estimated coefficients, variance of random effects. 
     #
     # INPUT    fixedform          2-sided formula for fixed effects
     #          nofactform         2-sided formula for fixed effects, omitting factors but retaining
     #                                           random terms                                           NEEDED?
     #          alldata            Data frame, first column of which must be "Observation".    
     #          randomform         1-sided formula for random effects
     #          groupname          Quoted name of single grouping variable within randomform
     #          randfactnames      Quoted names of random factor variables in randomform, or NULL
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          skip.step1         NULL or a list, each element of which is a vector of integers for 
     #                                 observations from 1 subgroup to be included in Step 1
     #          unblinded          TRUE permits printing of ultimate lme analysis, as specified above
     #          begin.diagnose     Numeric. Indicates where in code to begin printing diagnostics.  
     #                                  0 prints all; 100 prints none.
     #          incCont            Logical TRUE causes increase in allowable iterations and 
     #                                 allowable tolerances for duration of this function
     #          verbose            Logical. TRUE causes printing of function ID before and after running.
     # 

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
     spacer <- "$$$$$$$$$$$$$$$$$$$$$$$$$$                forsearch_lme        "
     options(warn=-1)
     on.exit(options(warn=0))

     ##################
     # Set lmeControl #
     ##################
     if(incCont){
          utils::str(lCtr <- nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, tolerance = 1e-2, 
                 niterEM = 1000, msMaxEval = 1000, msTol = 1e-2, optimMethod = "L-BFGS-B",
                 msVerbose = FALSE, returnObject = FALSE) )
          do.call(nlme::lmeControl, lCtr)
     }
     nalldata1 <- dim(alldata)[1]
     nalldata2 <- dim(alldata)[2]
     alldataNames <- names(alldata)
     #
     #########################################################
     # Ensure that first independent variable is Observation #
     #########################################################
     if(alldataNames[1] != "Observation") stop("First column of data must be 'Observation'")
     #
     ##########################
     # Locate response column #
     ##########################
     respName <- formula.tools::lhs(fixedform)
     ycolfixed <- (1:nalldata2)[respName==alldataNames]
     #
     ###############################################################
     # Locate grouping variable and ensure that it is not a factor #
     ###############################################################
     if(length(groupname) > 1)stop("Only one group name is allowed")
     groupcol <- (1:nalldata2)[groupname==alldataNames]
     if(is.factor(alldata[,groupcol])) stop("groupname must not be a factor in the database")
     #
     ######################################################################################
     # Create all files needed below:                                                     #
     #    fixdat.df,   alldata with factor level indicator and group variable indicator   #
     #                                                                                    #
     #    fixdatouter.list, a 2-layer list with 1st layer group subset levels and second  #
     #        layer factor subset within group layers. If there are no factors, the       #
     #        layer will be a list of 1 level. The name of that level will be 'None'      #                                                              #
     #                                                                                    #
     #    fixdatcombo.list, a list whose elements are data frames by combiation of        #
     #        factor subset and group level                                               #
     ######################################################################################
     fixdat.df <- alldata
     #
     ##############################################################
     # Check for factor status of fixedform and get factor names  #
     # Add factor subset indicator if there are any factors       #
     ##############################################################
     ufactor <- rep(TRUE, nalldata2)
     for(m in 1:nalldata2) ufactor[m] <- is.factor(alldata[,m])
     yesfactor <- any(ufactor)     
     #
     fixedISG <- "_None"                                       # default if no factors
     if(yesfactor){
          factorNames <- alldataNames[ufactor]  
          #############################################
          # Append the factor subset code to the data #
          #############################################
          isfactor <- rep(TRUE,nalldata2)
          for(nn in 1:nalldata2){
               isfactor[nn] <- is.factor(fixdat.df[,nn])
          }
          factorcount <- sum(isfactor)
          justfactors <- fixdat.df[isfactor]                        
          fixedISG <- apply(justfactors, 1, paste,collapse="/")
          fixedISG <- paste("_", fixedISG, sep="")
     }     #   yesfactor                                         # does this work with F1 and F2?
     fixdat.df <- data.frame(fixdat.df, fixedISG)   
     namesSubsets <- unique(fixedISG)
     nsubsets <- length(namesSubsets)                                  # number of model factor levels
     #
     ###############################################################
     # Identify factor-group cells                                 #
     # Append a combined ISG and determine number of combo subsets #
     ###############################################################
     comboISG <- paste(fixdat.df$fixedISG, fixdat.df[,groupcol], sep="_F,G ")
     fixdat.df <- data.frame(fixdat.df, comboISG)
     ucombo <- unique(comboISG)
     nucombo <- length(ucombo)
     #
     ugroups <- unique(fixdat.df[,groupcol])
     ngrouplevels <- length(ugroups)
     ###################################################
     # Structure fixdatouter.list and then populate it #
     ###################################################
     fixdatouter.list <- vector("list", ngrouplevels)
     fixdatinner.list <- vector("list", nsubsets)
     for(k in 1:ngrouplevels){
          fixdatouter.list[[k]] <- fixdatinner.list
     }
     grouplevels <- unique(alldata[,groupcol])
     ngrouplevels <- length(grouplevels)
     for(k in 1:ngrouplevels){
          for(j in 1:nsubsets){
              uu <- fixdat.df[    fixdat.df$fixedISG==namesSubsets[j] ,]
              vv <- uu[    uu[,groupcol] == grouplevels[k]    ,     ]
              fixdatouter.list[[k]][[j]] <- vv
          }     #   j
     }          #   k         
     #
     ###################################################
     # Structure fixdatcombo.list and then populate it #
     ###################################################
     fixdatcombo.list <- vector("list", nucombo)
     for(i in 1:nucombo){
          uu <- fixdat.df[ fixdat.df$comboISG==ucombo[i], ]
          fixdatcombo.list[[i]] <- uu
     }
#
###########################################################################################################################################
# Step 1
     holdrim <- NULL 
     rows.in.model <- vector("list", nalldata1)
     LLL <- vector("list", nalldata1)
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result

     if(is.null(skip.step1)){
          print("BEGINNING STEP 1", quote=FALSE)
          print(" ", quote=FALSE)

          ###############################################################################
          # Determine the number of coefficients in the entire database analyzed by lme #
          ###############################################################################
          lmeAll <- nlme::lme(fixed=fixedform, data=alldata, random=randomform)                        # lme no control
          coeffAll <- stats::coef(lmeAll)
          coeffAll <- coeffAll[1,]
          ncoeffs <- length(coeffAll)
          #
          if(yesfactor){
               #####################################################
               # Determine the structure of the database           #
               # number of group, factor, and continuous variables #
               #####################################################
               nG <- length(groupname)
               nF <- sum(isfactor)
               nC <- nalldata2 - nF - nG - 2    #excludes, Observation, response
               ##################################################
               # Calculate inner.r and source matrix (yf TRUE   #
               # nG determines which part of step1.dist is run  #
               ##################################################
               ufixed <- unique(fixedISG)

               ##############################################
               # Determine now many observations are needed #
               ##############################################
               N <- max(nucombo,ncoeffs) + nC
               ########################################
               # Now determine how to distribute them #
               ########################################
               basenumber <- max(1,floor(N/nucombo + .00001)   )
               additional <- rep(0, nucombo)
               nadd <- N - basenumber * nucombo
               if(nadd > 0)additional[1:nadd] <- 1
               inner.r <- rep(basenumber, times=nucombo) + additional
               inner.r <- inner.r[1:nucombo]
               #
               ################################################################################
               # Define number of source slots and define the number of observations per slot #
               ################################################################################
               firstrim <- bStep1(yesfactor, df1=fixdat.df, df1.ls=fixdatouter.list, groups=grouplevels, 
                            inner.rank=inner.r, source=comboISG, 
                            initial.sample=initial.sample, nofactform=nofactform,formulaA=fixedform, 
                            randform=randomform, inc=incCont, ycol=ycolfixed, b.d=begin.diagnose)          #  bStep1 
          SOON <- sort(firstrim)
          nfirstrim <- length(firstrim) 
          mstart <- nfirstrim + 1
          rows.in.model[[nfirstrim]] <- sort(firstrim)
          }     # End of yesfactor
          else{
               nG <- 1
               nC <- nalldata2 - nG - 2    #excludes, Observation, response
               #################################################
               # Calculate inner.r and source matrix (yf FALSE #
               #################################################
               ufixed <- unique(fixedISG)
               ##############################################
               # Determine now many observations are needed #
               ##############################################
               N <- max(nucombo,ncoeffs) + nC
               ########################################
               # Now determine how to distribute them #
               ########################################
               basenumber <- max(1,floor(N/nucombo + .00001)   )
               additional <- rep(0, nucombo)
               nadd <- N - basenumber * nucombo
               if(nadd > 0)additional[1:nadd] <- 1
               inner.r <- rep(basenumber, times=nucombo) + additional
               inner.r <- inner.r[1:nucombo]
               #
               firstrim <- bStep1(yesfactor, df1=fixdat.df, df1.ls=fixdatouter.list, groups=grouplevels,
                            inner.rank=inner.r, source=comboISG, initial.sample=initial.sample, nofactform=nofactform, 
                            formulaA=fixedform, randform=randomform, inc=incCont, ycol=ycolfixed, b.d=begin.diagnose)     #  bStep1 

          }     # no factor present
          SOON <- sort(firstrim)
          nfirstrim <- length(firstrim) 
          mstart <- nfirstrim + 1
          rows.in.model[[nfirstrim]] <- sort(firstrim)
     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)       # no guarantees for user-defined skip.step1
          nfirstrim <- length(skip.step1)
          rows.in.model[[nfirstrim]] <- sort(skip.step1)
 
                               if(begin.diagnose <=42) {print(paste(spacer,"Section 42",sep=" "),quote=FALSE);
                                                   Hmisc::prn(rows.in.model);Hmisc::prn(alldata[skip.step1,])   }

          SOON <- sort(skip.step1)
          mstart <- nfirstrim
     }       # skipping step 1
     # 
                                               if(begin.diagnose <=43) {print(paste(spacer,"Section 43",sep=" "),quote=FALSE);
                                                  Hmisc::prn(SOON)   }
# stop("before entering step 2")
########################################################################################################################################################################
# Step 2
     print("BEGINNING STEP 2", quote=FALSE)
     fixedform <- stats::formula(fixedform)

     #############################################################################
     # Bind a combination index to fixdat.df and add a column for squared errors #
     #############################################################################
     diffs2 <- 0
     pushaside <- formula.tools::rhs(randomform)=="1"
     if(yesfactor){inner.r <- rep(inner.r, length(SOON))}
     zzzz <- bStep2(yf=yesfactor, f2=fixedform, dfa2=fixdat.df, onlyfactor=pushaside, randm2=randomform, ms=mstart, 
                    ycol=ycolfixed, initn=inner.r, inc=incCont, finalm=rows.in.model, fbg=fixdatcombo.list, b.d=begin.diagnose)     # bStep2

     rows.in.set <- zzzz[[1]]
     LME <- zzzz[[2]]

     rows.in.set[[nalldata1]] <- 1:nalldata1
     LME[[nalldata1]] <- nlme::lme(fixedform, alldata, randomform)                                # lme

# stop("before extraction")
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
          print("The assumed analysis for these data is as follows:", quote=FALSE)
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
     levAlldata <- stats::lm(formula=fixedform, data=alldata, x=TRUE, y=TRUE)         # lm ????
     x1 <- levAlldata$x

     for(dd in mstart:nalldata1){                               #    for loop starts here
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

                                               if(begin.diagnose <= 88){print(paste(spacer,"Section 88         dd=",dd));
                                                              Hmisc::prn(zholdlme)    }

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
          y1 <- fixdat.df[,ycolfixed]
          for(j in 1:td1){
               errors[j] <- y1[j] - sum(hold.coeffs.fixed[[dd]] * x1[j,])
          }             #   j
          hold.residuals[,dd] <- errors
          #
          ############################################################
          #          "Fixed parameter estimates"=          param.est #
          ############################################################
          param.est[,dd] <- c(zholdlme$coefficients[[1]])     
                                     if(begin.diagnose <= 90){print(paste(spacer,"Section 90         dd=",dd));
                                              Hmisc::prn(param.est[,dd])}
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
 #                                       if(begin.diagnose <= 93){print(paste(spacer,"Section 93         dd=",dd));Hmisc::prn(j);
 #                                                  Hmisc::prn(dim(x1)[2]);Hmisc::prn(matrix(xtemp[j,],nrow=1));Hmisc::prn(det(crossinv))       }

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
     hold.summary.stats <- data.frame(m, hold.summary.stats)         # hold.summary.stats defined
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
     m <- 1:(dim(anova.pvalues)[1])
     anova.pvalues <- cbind(m, anova.pvalues)
     if(dnAVlme[1]=="(Intercept)"){
          anova.pvalues <- anova.pvalues[,-1]           # remove (Intercept)
     }

                                       if(begin.diagnose <= 99){print(paste(spacer,"Section 99"));
                                                   Hmisc::prn(inner.r)     }

    listout <- list(
          "Number of observations in Step 1"=   mstart-1,
          "Step 1 observation numbers"=         SOON,
          "Rows by subgroup"=                   fixdatouter.list,
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
