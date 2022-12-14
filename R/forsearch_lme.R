#' @export
forsearch_lme <-
function(fixed, data, random, formula, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL, 

    XmaxIter=1000,
    XmsMaxIter=1000, 
    Xtolerance=.01,
    XniterEM=10,
    XmsMaxEval=400,
    XmsTol=.00001, 
    Xopt='optim',
    unblinded=TRUE,

    begin.diagnose= 1000, verbose=TRUE)
{
     #                                           forsearch_lme 
     #
     # VALUE    List of datasets and statistics for plotting in forward search procedure to diagnose lme observations.  These include: scaled residuals, s^2, leverage,
     #                 estimated coefficients, variance of random effects. 
     #
     # INPUT    fixed              2-sided formula for fixed effects
     #          data               Grouped or ungrouped data frame, first column of which must be "Observation".    
     #          random             1-sided formula for random effects
     #          formula            a formula of the form resp ~ cov | group where resp is the response, cov is the primary covariate, and group is the grouping factor, 
     #                                   as in nlme::groupedData.
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
     # NOTE: Data will be grouped, a la Pihneiro and Bates.  The approach here is to form outer groups and renumber the observations, keeping track of the
     #       original values.  All groups will be represented in the rim of Step 1.  Then the rim will be augmented by 1 observation only. We keep track
     #       of the rim of each subgroup at each stage and combine them, reporting out the results in terms of the original observation numbers.
     #  
     # REF:  Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
     #       Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.
     #       https://CRAN.R-project.org/package=nlme
     #
     MC <- match.call()
     if(verbose) {
          print("")
          print("Running forsearch_lme")
          print("")
          print(date())
          print("")
          print("Call:")
          print(MC)
          print("")
     }
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
          print("")
          print("There is no replication in this dataset. All observations are defined as a combination of factors.")
          print("This function does not support such datasets.")
          print("")
          stop("Try eliminating one or more of the factors in the database and in the formulas of the call.")
     }
     #
     nobs <- dim(data)[1]
     gG <- nlme::getGroups(data)
     namesgG <- names(gG)
     nnamesgG <- length(namesgG)
     namesgG.1 <- gG
     if(nnamesgG > 1){namesgG.1 <- gG[,1] }
     else {namesgG.1 <- gG}                              # single vector of groups 
     OuterSubgroup <- paste("SG-", namesgG.1,sep="")
     fixdat <- data.frame(data,OuterSubgroup)
     #
     ############################################################################
     # Initialization of data frame and of indices of fixed and random effects  #
     # data might be a groupedData object; regroup as fixdat with known formula #
     ############################################################################
     print("")
     print("BEGINNING STEP 0")

     ############################################################################
     # Print a summary of the assumed analysis so that the user can be sure     #
     #           that the forward search is done on the relevant analysis.      #
     # Remove variables not included in Observation, response, indep variables, #
     #          OuterSubgroup from fixdat                                       #
     ############################################################################
     namesfixdat <- names(fixdat)
     nnamesfixdat <- length(namesfixdat)
     response.variable <- as.character(formula.tools::lhs(fixed))
     indep.vars <- formula.tools::rhs(fixed)
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
     ###########################################
     # Ensure that fixdat is a grouped dataset #
     ###########################################
     fixdat <- nlme::groupedData(formula, fixdat)

     obslist <- c("Observation", response.variable, indep.vars, "OuterSubgroup")
                                                              spacer <- rep(" ", 20)
                                                              if(begin.diagnose <= 0){print(c(spacer,"Section 0"), quote=FALSE);Hmisc::prn(obslist)}

     crostab <- outer(obslist, namesfixdat, FUN="==") 
     crostab <- apply(crostab,2,any)
     varnums <- (1:length(crostab))[crostab]
     fixdat2 <- fixdat[,varnums]
     print("")
     print("Observations will be reordered based on the following data (partial listings):")
     print("")
     print(utils::head(fixdat2))
     print("")
     print(utils::tail(fixdat2))
     print("")
     #
     ###########################################################
     # Check for factor status of each variable within dataset #
     # Assumes same structure for every group and subgroup     #
     ###########################################################
     yesfactor <- FALSE
     if(indep.vars[1] != "NoIndepVars"){
          dimdata <- length(indep.columns)
          ufactor <- rep(TRUE, dimdata)
          for(m in 1:dimdata) ufactor[m] <- is.factor(data[,indep.columns[m]])
          yesfactor <- any(ufactor)
     }
     if(yesfactor){print("There are factors in the design")}else{print("There are no factors in the design")}
     #
     #######################################################################################
     # Set up list to hold observations by group (fixdat.by.group).  This is permanent set #
     # Renumber all Observations and save original observation numbers in saved.obsnums    #
     #######################################################################################
     ufixdatOuter <- unique(fixdat$OuterSubgroup)
     nufixdatOuter <- length(ufixdatOuter)
     fixdat.by.group <- vector("list", nufixdatOuter)                                        #  list
     saved.obsnums <- vector("list",nufixdatOuter  )                                         #  list
                                                              if(begin.diagnose <= 5){print(c(spacer,"Section 5"));Hmisc::prn(saved.obsnums)}
                                                              if(begin.diagnose <= 10){print(c(spacer,"Section 10"));Hmisc::prn(obslist)}
     n.fixdat <- rep(0, nufixdatOuter)
     for(jj in 1:nufixdatOuter){
          fixdat.by.group[[jj]] <- fixdat[fixdat$OuterSubgroup==ufixdatOuter[jj],]
          extractObs <- fixdat.by.group[[jj]]
          saved.obsnums[[jj]] <- extractObs[,1]
          dimextract <- dim(extractObs)[1]
          extractObs[,1] <- 1:dimextract
          fixdat.by.group[[jj]] <- extractObs 
          n.fixdat[jj] <- dim(fixdat.by.group[[jj]])[1]
     }        #  jj
                                                              if(begin.diagnose <= 15){print(c(spacer,"Section 15"), quote=FALSE);Hmisc::prn(saved.obsnums)}
     maxnfixdat <- max(n.fixdat)
     #
     ####################################################################
     # Set up variables to hold partial results. This is augmenting set #
     # rim.all.subgroups is indexed by [[kk]] then [[dd]]               #
     ####################################################################
     rim.all.subgroups <- vector("list", nufixdatOuter)                          # holds vectors of final determination of each rim in each subgroup kk runs over columns within a row
     current.obs.by.group <- vector("list", nobs)                                # contains rim for each stage within each subgroup          dd runs down column, over rows
     for(kk in 1:nufixdatOuter){
          rim.all.subgroups[[kk]] <- current.obs.by.group                        # this only sets up a list within a list structure, not data
     }
     newrims <- rim.all.subgroups                                                # for translation back to original observation numbers
     #
     ########################################################################
     # Renumber any set of observation numbers being used instead of step 1 #
     ########################################################################
     if(!is.null(skip.step1)){
          if(!is.matrix(skip.step1))stop("skip.step1 must be a matrix with the number of rows equal to the number of subgroups")
          vv <- skip.step1
          for(i in 1:nufixdatOuter){
               rr <- outer(vv[i,], saved.obsnums[[i]], "==")
               rr <- apply(rr,2,any)
               vv[i,] <- (1:n.fixdat[i])[rr]
          }    # for i
          skip.step1 <- vv
     }    # if !is null skip.step1
     #
     ##################################################################
     # Perform Step 1 or insert initial observation numbers manually? #
     ##################################################################
     if(is.null(skip.step1)){
          print("")
          print("BEGINNING STEP 1")
          print("")
          zz <- bStep1(fixed=fixed, yf=yesfactor, mnf=maxnfixdat, nOuter=nufixdatOuter, 
                     s.o=saved.obsnums, nopl=n.obs.per.level, nobs=nobs, i.s=initial.sample, 
                     fbg=fixdat.by.group, b.d=begin.diagnose,verbose=FALSE)                                   # bStep1
          mstart <- zz[[1]]
          rim.all.subgroups <- zz[[2]]
     }   # skip.step1 is NULL
     else{
          print("")
          print("SKIPPING STEP 1")
          print("")
          mstart <- dim(skip.step1)[2]
          for(i in 1:nufixdatOuter){
               rim.all.subgroups[[i]][[mstart]] <- skip.step1[i,]
          }
     }   # skip.step1 is not null
          print("")

     print("BEGINNING STEP 2")
     print("")
     rows.in.set <- bStep2(fixed=fixed, nOuter=nufixdatOuter, mnf=maxnfixdat, mstart=mstart, nobs=nobs, 
                    fbg=fixdat.by.group, n.f=n.fixdat, s.o=saved.obsnums, 
                    ras=rim.all.subgroups, b.d=begin.diagnose, verbose=FALSE)                                   # bStep2
     #
     ######################################################################################################
     # Multiple subgroups were combined to form rows.in.set. mstart must be increased to account for this #
     ######################################################################################################
     mstart <- mstart * nufixdatOuter
                                                              if(begin.diagnose <= 30){print(c(spacer,"Section 30"), quote=FALSE);Hmisc::prn(mstart)}
     #
     ###############################################################
     # Show assumed analysis if call is for unblinded presentation #
     ###############################################################
     print("")
     zholdlm <- nlme::lme(fixed=fixed, data=fixdat, random=random, control=newcontrol)                                  #    lme
     if(unblinded){
          print("The assumed analysis for these data will be as follows:")
          print("")
          print(zholdlm)
          print("")
     }
     ##########################################################
     # Retrieve independent and dependent columns from fixdat #
     ##########################################################
     lmAlldata <- stats::lm(formula=fixed, data=fixdat, singular.ok=TRUE, x=TRUE, y=TRUE)                               #   lm
     y1 <- lmAlldata$y
     x1 <- lmAlldata$x                                              # takes care of more complex fixed effect formulas
     coeffnames <- dimnames(x1)[2]
     p <- dim(x1)[2]
     #
     #################################################################################################
     # Retrieve metadata for plotting from each row of rows.in.set                                   #
     # Set up files to hold residuals, dimensions, sigma, fixed parameter estimates,                 #
     #    random parameter estimates, leverage, Cook distance and list for storage of xtemp matrices #                                                                                            #
     #################################################################################################
     nrowsdf1 <- nobs

     dddd <- zholdlm$dims
     totalvars <- sum(unlist(dddd$ncol[1:dddd$Q]))
     ss77 <- variablelist(datadf=data, verbose=FALSE)
     hold.residuals <- matrix(0,nrowsdf1,nrowsdf1)          # for standardized residuals across all observations

     hold.subset.residuals <- rep(0,nrowsdf1)               # for use in Cook distance
     hold.dims <- vector("list", nrowsdf1)
     hold.sigma <- rep(-999, nrowsdf1)   
     param.est <- matrix(0,nrow=p, ncol=nrowsdf1)
     t.set <- param.est
     hold.coeffs.fixed <- vector("list",nrowsdf1)

     hold.coeffs.random <- matrix(0,nrow=nrowsdf1, ncol=totalvars) 
     leverage <- matrix(-999,nrow=1,ncol=3)
     xtemp.list <- vector("list",nrowsdf1)
     modCook <- rep(0,nrowsdf1-1)
     hold.summary.stats <- matrix(0,nrow=nrowsdf1,3)        # unnamed columns: AIC, BIC, log likelihood

     print("BEGINNING INTERMEDIATE RESULTS EXTRACTION")

     for(dd in mstart:nobs){
          rim <- rows.in.set[[dd]]             # picks up row numbers for a set of observations
          Zlatest <- fixdat[rim,]  
          Zlatest <<- Zlatest                  # allows lme to find file.    Needed?  directly tested and the answer is YES
          ############################################################################
          # Extract indep vars of the subset for use in leverage and Cook's distance #
          # lm function allows collection of x matrix and y vector                   #
          ############################################################################
          xtemp <- x1[rim,]
          xtemp.list[[dd]] <- xtemp 
          transtemp <- t(xtemp)
          cross <- transtemp %*% xtemp
          crossinv <- solve(cross)
          ################################################################
          # Capture results of lme analysis for this set of observations #
          ################################################################
          zholdlme <- nlme::lme(fixed=fixed, data=Zlatest, random=random, control=newcontrol)               #   lme
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
          ###########################################################
          # There is no predict function for lme, so use predict.lm #
          ###########################################################
          LMzholdlme <- stats::lm(formula=fixed, data=data.frame(Zlatest))                                  #   lm
          resids <- zholdlme$residuals
          hold.coeffs.fixed[[dd]] <- zholdlme$coefficients[[1]]
          temprand<-zholdlme$coefficients[[2]]
          gotvars <- rep(0,totalvars)
                                                              if(begin.diagnose <= 35){print(c(spacer, "Section 35         dd = ",dd), quote=FALSE);Hmisc::prn(LMzholdlme)}
          for(level in 1:totalvars){
               uvars <- c(temprand[[level]])
               gotvars[level] <- sqrt(sum(uvars*uvars)/length(uvars))
          }
          hold.coeffs.random[dd,] <- gotvars
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
               predlmeall <- stats::predict(object=LMzholdlme, newdata=fixdat)
               errors[j] <- y1[j] - sum(hold.coeffs.fixed[[dd]] * x1[j,])
          }             #   j
          hold.residuals[,dd] <- errors
          #
          ############################################################
          #          "Fixed parameter estimates"=          param.est #
          ############################################################
          param.est[,dd] <- c(zholdlme$coefficients[[1]])     
                                                              if(begin.diagnose <= 40){print(paste(c(spacer,"Section 40         dd="),dd));Hmisc::prn(param.est[,dd])}
          #
          #####################################################################
          #          "Random parameter estimates"=         hold.coeffs.random #
          #####################################################################
          temprand<-zholdlme$coefficients[[2]]
          gotvars <- rep(0,totalvars)
          for(level in 1:totalvars){
               uvars <- c(temprand[[level]])
               gotvars[level] <- sqrt(sum(uvars*uvars)/length(uvars))
          }
          hold.coeffs.random[dd,] <- gotvars
          ############################################################
          #           Leverage=                        leverage[-1,] #
          ############################################################
          thisleverage <- 1
          if(is.matrix(x1)){
               for(j in 1:dd){
                   Zlatest2 <- data.frame(Zlatest)
                   if(dim(x1)[2]==1){
                        thisleverage <- c(c(matrix(xtemp[j],nrow=1) %*% crossinv %*% matrix(xtemp[j],ncol=1)))
                        thisleverage <- c(dd,Zlatest2[j,1],thisleverage)
                   }else{
                        thisleverage <- c(c(matrix(xtemp[j,],nrow=1) %*% crossinv %*% matrix(xtemp[j,],ncol=1)))
                        thisleverage <- c(dd,Zlatest2[j,1],thisleverage)
                   }
                   leverage <- rbind(leverage,thisleverage)
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
     }                                                                  #   dd                                                     

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
                                                              if(begin.diagnose <= 50){print(c(spacer,"Section 50"));Hmisc::prn(aa);Hmisc::prn(bb)}
          www <- aa %*% t(bb)
          modCook[i-1] <- (www %*% t(www))/(p * hold.subset.residuals[i-1])
                                                              if(begin.diagnose <= 53){print(c(spacer,"Section 53"));Hmisc::prn(hold.subset.residuals[i-1])}
     }
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
     names(hold.coeffs.random) <- outnames
    #
     listout <- list(
          "Number of rows included in Step 1"=  mstart-1,
          "Rows by subgroup"=                   fixdat.by.group,
          "Rows in stage"=                      rows.in.set,
           Sigma=                               sigma,
          "Standardized residuals"=             hold.residuals,            
          "Fixed parameter estimates"=          param.est,
          "Random parameter estimates"=         hold.coeffs.random,
           Leverage=                            leverage[-1,],
          "Modified Cook distance"=             modCook,
           Dims=                                zholdlme$dims,
          "t statistics"=                       t.set,
          "Fit statistics"=                     hold.summary.stats,
           Call=                                MC )
    #
     if(verbose) {
          print("")
          print("Finished running forsearch_lme")
          print("")
          print(date())
          print("")
     }

     return(listout)
}
