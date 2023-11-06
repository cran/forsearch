#' @export
forsearch_lme <-
function(fixedform, data, randomform, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL, 

    XmaxIter=1000,
    XmsMaxIter=1000, 
    Xtolerance=.01,
    XniterEM=10,
    XmsMaxEval=400,
    XmsTol=.00001, 
    Xopt='optim',
    unblinded=TRUE,

    begin.diagnose= 100, verbose=TRUE)
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
     #          data               Grouped or ungrouped data frame, first column of which must be "Observation".    
     #          randomform         1-sided formula for random effects
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          n.obs.per.level    Number of observations per level of (possibly crossed) factor levels
     #          skip.step1         NULL or a list, each element of which is a vector of integers for observations from 1 subgroup to be included in Step 1
     #
     #          Xmaxiter, XmsMaxIter, Xtolerance, XniterEM, XmsMaxEval, XmsTol, Xopt 
     #                             control variates for lme function
     #          unblinded          TRUE permits printing of ultimate lme analysis, as specified above
     #          begin.diagnose     Numeric. Indicates where in code to begin printing diagnostics. 0 prints all; 100 prints none.
     #          verbose            Logical. TRUE causes printing of function ID before and after running.
     # 
     # NOTE: The lme function will be run to extract the intermediate results. However, implementing the increase in subset size will modify the Atkinson
     #       and Riani approach because of the possibility of crashing the program due to a bad configuration of observations. Simple increase in subset size 
     #       is never a problem in this regard, but substitution of a superior set can be a problem if a subset were completely removed. In that case, the 
     #       prediction procedure would find a new parameter in the full data set that wasn't present in the predictor set. Accordingly, we determine candidate
     #       sets of additional observations within each subgroup and pick the best one, thus ensuring that any removals will be replaced only by members of 
     #       that subgroup. Also, in this function all of the sets of observation subsets is defined before any of the extractions are made.
     #
     # REF:  Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
     #       Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.
     #       https://CRAN.R-project.org/package=nlme
     #
     #      diagnostic sections--  Step 0:   1 - 25      Step 1:  29 - 29     Step 2: 30 - 50 (in bStep2)     Extract: 
     # 
# Step 0
     print("BEGINNING STEP 0", quote=FALSE)
#     print("", quote=FALSE)
     spacer <- "                                            forsearch_lme        "
     options(warn=-1)
     on.exit(options(warn=0))
     newcontrol <- nlme::lmeControl(maxIter=XmaxIter,  msMaxIter=XmsMaxIter, tolerance=Xtolerance, 
                      niterEM=XniterEM, msMaxEval=XmsMaxEval, msTol=XmsTol, opt=Xopt)
     newcontrol <<- newcontrol
     #
     #########################################################
     # Ensure that first independent variable is Observation #
     #########################################################
     uu <- dimnames(data)
     if(uu[[2]][1] != "Observation"){Hmisc::prn(uu);stop("First column of dataset must be 'Observation'")}
     #
     ####################################################
     # Ensure that there is replication in the database #
     ####################################################
     varlist <- variablelist(data, verbose=FALSE)
     if(length(varlist)==dim(data)[1]){
          print("", quote=FALSE)
          print("There is no replication in this dataset. All observations are defined as a combination of factors.", quote=FALSE)
          print("This function does not support such datasets.", quote=FALSE)
          print("", quote=FALSE)
          stop("Try eliminating one or more of the factors in the database and in the formulas of the call.")
     }
     #
     ################################################
     # Determine column number of response variable #
     ################################################
     response.variable <- as.character(formula.tools::lhs(fixedform))
     response.colnum <- names(data)==response.variable
     response.colnum <- (1:length(response.colnum))[response.colnum]
     #
     ###################################################################
     # Define list to hold lme objects and insert object with all data #
     ###################################################################
     nobs <- dim(data)[1]
     LME <- vector("list",nobs)
     preliminary.lme <- nlme::lme(fixed=fixedform, data=data, random=randomform)                              # lme
     p <- preliminary.lme$dims$ncol[2]
     LME[[nobs]] <- preliminary.lme

     gG <- nlme::getGroups(preliminary.lme, sep=".")
     rank.prelim <- preliminary.lme$coefficients
     rank.prelim <- length(rank.prelim[[1]])
     #################################################
     # Add variable to input data for outer subgroup # 
     #################################################
     namesgG <- names(gG)
     nnamesgG <- length(namesgG)
     namesgG.1 <- gG
     if(nnamesgG > 1){
          namesgG.1 <- gG[,1]
     }
     else {
          namesgG.1 <- gG
     }                              # single vector of groups 
     OuterSubgroup <- paste("SG-", namesgG.1,sep="")
     fixdat <- data.frame(data,OuterSubgroup)
     #
     print("", quote=FALSE)
     print("BEGINNING STEP 0", quote=FALSE)

     namesfixdat <- names(fixdat)
     nnamesfixdat <- length(namesfixdat)
     response.variable <- as.character(formula.tools::lhs(fixedform))
     indep.vars <- formula.tools::rhs(fixedform)
     indep.vars <- formula.tools::get.vars(indep.vars)

     if(length(indep.vars) > 0){
          incols <- c(outer(namesfixdat,indep.vars, "=="))
          colnums <- 1:nnamesfixdat
          nindepvars <- length(indep.vars)
          colnums <- rep(colnums,times=nindepvars)
          indep.columns <- colnums[incols]
     }
     else{
          indep.vars <-"NoIndepVars"
     }
     obslist <- c("Observation", response.variable, indep.vars, "OuterSubgroup")               # this is just a character vector IGNORES GROUPING VARIABLES

                                                              if(begin.diagnose <= 1){print(paste(spacer,"Section 1"), quote=FALSE);Hmisc::prn(obslist)}

     ###########################################################
     # Check for factor status of each variable within dataset #
     # Assumes same structure for every group and subgroup     #
     # If there are factors, add a variable for inner groups   #
     # Assumes that only independent variables can be factors  #
     ###########################################################
     yesfactor <- FALSE
     if(indep.vars[1] != "NoIndepVars"){
          dimdata <- length(indep.columns)
          ufactor <- rep(TRUE, dimdata)
          for(m in 1:dimdata) ufactor[m] <- is.factor(data[,indep.columns[m]])
          yesfactor <- any(ufactor)
     }
     if(yesfactor){
          print("There are factors in the design", quote=FALSE)
     }
     else{
          print("There are no factors in the design", quote=FALSE)
     }

                                                              if(begin.diagnose <= 1){print(paste(spacer,"Section 1"), quote=FALSE);Hmisc::prn(obslist)}

     dimdata2 <- dim(fixdat)[2]
     isfactor <- rep(TRUE,dimdata2)
     if(yesfactor){
          for(nn in 1:dimdata2){
               isfactor[nn] <- is.factor(fixdat[,nn])
          }
          justfactors <- fixdat[isfactor]
          holdISG <- apply(justfactors, 1, paste, collapse="/")
          OuterInner <- data.frame(OuterSubgroup,holdISG)
          OuterInner <- apply(OuterInner,1,paste,collapse="_")
          fixdat <- data.frame(fixdat, OuterInner)         #   replace outer variable with OuterInner
     }    # yesfactor
     #
     crostab <- outer(obslist, namesfixdat, FUN="==") 
     crostab <- apply(crostab,2,any)
     varnums <- (1:length(crostab))[crostab]
     fixdat2 <- fixdat[,varnums]
     #
     #######################################################################################################################
     # Set up lists to hold observations by group (fixdat.by.Outer and fixdat.by.OuterInner).  These are permanent sets    #
     # Renumber all Observations and save original observation numbers in saved.obsnums.Outer and saved.obsnums.OuterInner # 
     #######################################################################################################################
     fixdat.by.OuterInner <- NULL
     nufixdatOuterInner <- NULL
     if(yesfactor){
          ufixdatOuterInner <- unique(fixdat$OuterInner)
          nufixdatOuterInner <- length(ufixdatOuterInner)
          fixdat.by.OuterInner <- vector("list", nufixdatOuterInner)                                        #  list
          saved.obsnums.OuterInner <- vector("list",nufixdatOuterInner  )                                         #  list

                                                              if(begin.diagnose <= 5){print(paste(spacer,"Section 5"));Hmisc::prn(saved.obsnums.OuterInner)}
                                                              if(begin.diagnose <= 10){print(paste(spacer,"Section 10"));Hmisc::prn(obslist)}
          n.fixdat <- rep(0, nufixdatOuterInner)
          for(jj in 1:nufixdatOuterInner){
               fixdat.by.OuterInner[[jj]] <- fixdat[fixdat$OuterInner==ufixdatOuterInner[jj],]
               extractObs <- fixdat.by.OuterInner[[jj]]
               saved.obsnums.OuterInner[[jj]] <- extractObs[,1]
               dimextract <- dim(extractObs)[1]
               extractObs[,1] <- 1:dimextract
               fixdat.by.OuterInner[[jj]] <- extractObs 
               n.fixdat[jj] <- dim(fixdat.by.OuterInner[[jj]])[1]
          }        #  jj
                                                              if(begin.diagnose <= 15){print(paste(spacer,"Section 15"), quote=FALSE);Hmisc::prn(saved.obsnums.OuterInner)}
          maxnfixdat <- max(n.fixdat)
     }    # there are inner factors
     #
     # Set up fixdat.by.Outer even if there are also inner factors
     ufixdatOuter <- unique(fixdat$OuterSubgroup)
     nufixdatOuter <- length(ufixdatOuter)
     fixdat.by.Outer <- vector("list", nufixdatOuter)                                        #  list
     saved.obsnums.Outer <- vector("list",nufixdatOuter  )                                         #  list

                                                              if(begin.diagnose <= 15){print(paste(spacer,"Section 15"));Hmisc::prn(saved.obsnums.Outer)}
                                                              if(begin.diagnose <= 20){print(paste(spacer,"Section 20"));Hmisc::prn(obslist)}
     n.fixdat <- rep(0, nufixdatOuter)
     for(jj in 1:nufixdatOuter){
          fixdat.by.Outer[[jj]] <- fixdat[fixdat$OuterSubgroup==ufixdatOuter[jj],]
          extractObs <- fixdat.by.Outer[[jj]]
          saved.obsnums.Outer[[jj]] <- extractObs[,1]
          dimextract <- dim(extractObs)[1]
          extractObs[,1] <- 1:dimextract
          fixdat.by.Outer[[jj]] <- extractObs 
          n.fixdat[jj] <- dim(fixdat.by.Outer[[jj]])[1]
     }        #  jj
                                                              if(begin.diagnose <= 25){print(paste(spacer,"Section 25"), quote=FALSE);Hmisc::prn(saved.obsnums.Outer)}
          maxnfixdat <- max(n.fixdat)
     #
     ####################################################################
     # Set up variables to hold partial results. This is augmenting set #
     # rim.all.subgroups is indexed by [[kk]] then [[dd]]               #
     ####################################################################
     rim.all.subgroups <- vector("list", nufixdatOuter)                          # holds vectors of final determination of each rim in each subgroup kk runs over columns within a row
     current.obs.by.group <- vector("list", nobs)                                # contains rim for each stage within each subgroup          dd runs down column, over rows
     for(kk in 1:nufixdatOuter){
          rim.all.subgroups[[kk]] <- current.obs.by.group                        # this only sets up a list within a list structure, not data  NOT USED BEYOND HERE
     }
     newrims <- rim.all.subgroups                                                # for translation back to original observation numbers
     #
     ########################################################################
     # Renumber any set of observation numbers being used instead of step 1 #     Let's see if we really need this
     ########################################################################
     if(!is.null(skip.step1)){
          if(!is.matrix(skip.step1))stop("skip.step1 must be a matrix")
          if(yesfactor){
               if(!dim(skip.step1)[1]==nufixdatOuterInner)stop("The number of rows in skip.step1 must equal the number of subgroups")
          }
          else{
               if(!dim(skip.step1)[1]==nufixdatOuter)stop("The number of rows in skip.step1 must equal the number of subgroups")
          }
     }    # if !is null skip.step1
     #
#######################################################################################################################################################################
# Step 1
     print("", quote=FALSE)
     if(yesfactor){
          print("Observations to enter Step 1 must be selected from each OuterInner subgroup:", quote=FALSE)
          print("", quote=FALSE)
          Hmisc::prn(ufixdatOuterInner)
     }
     else{
          print("Observations to enter Step 1 must be selected from each OuterSubgroup:", quote=FALSE)
          print("", quote=FALSE)
          uOuter <- unique(OuterSubgroup)
          Hmisc::prn(uOuter)
     }
     print("", quote=FALSE)
     ###############################################################
     # Begin Step 1 or insert initial observation numbers manually #
     ###############################################################
     rows.in.set <- vector("list",nobs)
     ##########################################################
     # Retrieve independent and dependent columns from fixdat #
     ##########################################################
     lmAlldata <- stats::lm(formula=fixedform, data=fixdat, singular.ok=TRUE, x=TRUE, y=TRUE)                               #   lm
     y1 <- lmAlldata$y
     x1 <- lmAlldata$x                                              # takes care of more complex fixed effect formulas
     coeffnames <- dimnames(x1)[2]
     p <- dim(x1)[2]
     #
     inner.rank.all <- lmAlldata$rank     -1
     if(is.null(skip.step1)){
          print("", quote=FALSE)
          print("BEGINNING STEP 1", quote=FALSE)
          print("", quote=FALSE)

          inner.rim <- vector("list", nufixdatOuter)
          for(i in 1:nufixdatOuter){
               zz <- aStep1(yesfactor=yesfactor, data=fixdat.by.Outer[[i]], inner.rank=inner.rank.all, initial.sample=initial.sample, 
                       formula=fixedform, ycol=response.colnum, nopl=n.obs.per.level)                                                #  aStep1
               inner.rim[[i]] <- saved.obsnums.Outer[[i]][zz]
          }    #   i
         inner.rim <- unlist(inner.rim)
         rim.all.subgroups <- sort(inner.rim)

         mstart.step1 <- length(rim.all.subgroups)
     }   # skip.step1 is NULL
     else{
          print("", quote=FALSE)
          print("SKIPPING STEP 1", quote=FALSE)
          print("", quote=FALSE)
          skip.step1b <- c(skip.step1)
          mstart.step1 <- length(skip.step1b)
          rim.all.subgroups <- sort(skip.step1b)
     }   # skip.step1 is not null
     rows.in.set[[nobs]] <- 1:nobs
     rows.in.set[[mstart.step1]] <- rim.all.subgroups     # drop the first set of observation numbers into rows.in.set
     mstart <- mstart.step1 + 1 
                                                              if(begin.diagnose <= 30){print(paste(spacer,"Section 30"), quote=FALSE);Hmisc::prn(rim.all.subgroups) }
     SOON <- rim.all.subgroups
     #
########################################################################################################################################################################
# Step 2
     print("BEGINNING STEP 2", quote=FALSE)
     print("", quote=FALSE)
     if(yesfactor){
          zzzz <- bStep2(fixed2=fixedform, fulldata=fixdat, random2=randomform, yf=yesfactor, nOuter=nufixdatOuterInner, mstart=mstart, nobs=nobs, yobs=response.colnum,
                    fbg=fixdat.by.OuterInner, n.f=n.fixdat, s.o=saved.obsnums.OuterInner, ufixdat=ufixdatOuterInner,
                    ras=rows.in.set, b.d=begin.diagnose, LLL=LME, verbose=TRUE)                                   # bStep2
     }
     else{
          zzzz <- bStep2(fixed2=fixedform, fulldata=fixdat, random2=randomform, yf=yesfactor, nOuter=nufixdatOuter, mstart=mstart, nobs=nobs, yobs=response.colnum,
                    fbg=fixdat.by.Outer, n.f=n.fixdat, s.o=saved.obsnums.Outer, ufixdat=ufixdatOuter,
                    ras=rows.in.set, b.d=begin.diagnose, LLL=LME, verbose=FALSE)                                   # bStep2
     }
     #

     rows.in.set <- zzzz[[1]]
     LME <- zzzz[[2]]
     ###############################################################
     # Show assumed analysis if call is for unblinded presentation #
     ###############################################################
     print("", quote=FALSE)
     zholdlm <- nlme::lme(fixed=fixedform, data=fixdat, random=randomform)                                  #    lme
     if(unblinded){
          print("The assumed analysis for these data will be as follows:", quote=FALSE)
          print("", quote=FALSE)
          print(zholdlm)
          print("", quote=FALSE)
     }
########################################################################################################################################################################
# Extracting summary statistics
     print("", quote=FALSE)
     print("BEGINNING INTERMEDIATE RESULTS EXTRACTION", quote=FALSE)
     #################################################################################################
     # Retrieve metadata for plotting from each row of rows.in.set                                   #
     # Set up files to hold residuals, dimensions, sigma, fixed parameter estimates,                 #
     #    random parameter estimates, leverage, Cook distance and list for storage of xtemp matrices #                                                                                            #
     #################################################################################################
     nrowsdf1 <- nobs
     dddd <- zholdlm$dims
     ss77 <- variablelist(datadf=data, verbose=FALSE)
     hold.residuals <- matrix(0,nrowsdf1,nrowsdf1)          # for standardized residuals across all observations
     hold.subset.residuals <- rep(0,nrowsdf1)               # for use in Cook distance   
     hold.dims <- vector("list", nrowsdf1)
     hold.sigma <- rep(-999, nrowsdf1)   
     param.est <- matrix(0,nrow=p, ncol=nrowsdf1)
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
     templme <- preliminary.lme
     tempanova <- stats::anova(templme)
     dnAVlme <- dimnames(tempanova)
     dnAVlme <- dnAVlme[[1]]
     anova.pvalues <- matrix(0, nrow=length(dnAVlme), ncol = dim(data)[1])
     #
     VC <- nlme::VarCorr(templme)
     VCnames <- dimnames(VC)[[1]]
     VC2 <- VC[,2]
     hold.coeffs.random <- matrix(0,nrow=nrowsdf1, ncol=length(VC2))       # won't need to transpose    
     for(dd in mstart:nobs){                                                        #    for loop starts here
          rim <- rows.in.set[[dd]]             # picks up row numbers for a set of observations
          Zlatest <- fixdat[rim,]  
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

                                                              if(begin.diagnose <= 35){print(paste(spacer,"Section 35         dd=",dd));Hmisc::prn(zholdlme)}

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
          hold.subset.residuals[dd] <- sum(resids^2)/(dd-p)
          #
          #########################################################
          #          "Standardized residuals"=     hold.residuals #
          #########################################################
          td1 <- dim(fixdat)[1]
          errors <- rep(-999, td1)
          for(j in 1:td1){
               errors[j] <- y1[j] - sum(hold.coeffs.fixed[[dd]] * x1[j,])
          }             #   j
          hold.residuals[,dd] <- errors
          #
          ############################################################
          #          "Fixed parameter estimates"=          param.est #
          ############################################################
          param.est[,dd] <- c(zholdlme$coefficients[[1]])     
                                                              if(begin.diagnose <= 40){print(paste(spacer,"Section 40         dd=",dd));Hmisc::prn(param.est[,dd])}
          #
          ###########################################################
          # ANOVA test of fixed effects               anova.pvalues #
          ###########################################################
          AVlme <- stats::anova(zholdlme)
          AVlmeps <- AVlme[,4]                     # no need to remove p value (NA) for residuals; not in AVlme
          anova.pvalues[,dd] <- AVlmeps
                                                              if(begin.diagnose <= 45){print(paste(spacer,"Section 45         dd=",dd));Hmisc::prn(anova.pvalues[,dd])}
          #
          ############################################################
          #           Leverage=                        leverage[-1,] #
          ############################################################
          thisleverage <- 1
          if(is.matrix(x1)){
               for(j in 1:dd){
                                                              if(begin.diagnose <= 47){print(paste(spacer,"Section 47         dd=",dd));Hmisc::prn(j);
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
     names(param.est) <- paste("b",1:p, sep="")
     m <- 1:nobs
     xxparam <- t(param.est)
     param.est <- as.data.frame(xxparam)
     names(param.est) <- paste("b",1:p, sep="")
     param.est <- rbind(m,param.est)
     #
     hold.summary.stats <- data.frame(m, hold.summary.stats)                          # hold.summary.stats defined
     names(hold.summary.stats) <- c("m", "AIC", "BIC", "logLik")
     # 
     t.set <- as.data.frame(t(t.set))
     dimnames(t.set)[2] <- coeffnames
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
     sigma.squared <- s.2/(nrowsdf1-p)
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
     for(i in mstart:nobs){
          aa <- param.diff[i-1,]
          aa <- as.numeric(aa[-1])
          aa <- matrix(aa,nrow=1)
          bb <- as.matrix(xtemp.list[[i]])
                                                              if(begin.diagnose <= 50){print(paste(spacer,"Section 50"));Hmisc::prn(aa);Hmisc::prn(bb)}
          www <- aa %*% t(bb)
          modCook[i-1] <- (www %*% t(www))/(p * hold.subset.residuals[i-1])

                                                              if(begin.diagnose <= 53){print(paste(spacer,"Section 53"));Hmisc::prn(hold.subset.residuals[i-1])}
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
          "Rows by outer subgroup"=             fixdat.by.Outer,
          "Rows by outer-inner subgroups"=      fixdat.by.OuterInner,
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
