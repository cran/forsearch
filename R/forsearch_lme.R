#' @export
forsearch_lme <-
function(
    fixed,  
    data,
    random,
    formula,
    response.column,
    initial.sample, 
    robs, 
    skip.step1=NULL, 

    XmaxIter=1000,
    XmsMaxIter=1000, 
    Xtolerance=.01,
    XniterEM=1000,
    XmsMaxEval=400,
    XmsTol=.00001, 
    Xopt='optim',

    diagnose=FALSE, verbose=TRUE)
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
     #          response.column    Column number of response variable
     #          indep.cols         Column numbers of independent variables
     #          initial.sample     Number of reorderings of observations (= m in Atkinson and Riani)
     #          robs               Number of observations per subgroup in initial sample 
     #          skip.step1         NULL or a vector of integers for rows to be included in Step 1
     #
     #          Xmaxiter, XmsMaxIter, Xtolerance, XniterEM, XmsMaxEval, XmsTol, Xopt 
     #                             control variates for lme function
     #
     #          diagnose           Logical. TRUE causes printing of diagnostic content
     #          verbose            Logical. TRUE causes printing of program ID before and after running.
     #  
     # REF:  Atkinson, A and M Riani. Robust Diagnostic Regression Analysis, Springer, New York, 2000.
     #       Pinheiro, JC and DM Bates. Mixed-Effects Models in S and S-Plus, Springer, New York, 2000.
     #       https://CRAN.R-project.org/package=nlme
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running forsearch_lme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     #########################################################
     # Ensure that first independent variable is Observation #
     #########################################################
     uu <- dimnames(data)
     if(uu[[2]][1] != "Observation"){Hmisc::prn(uu);stop("First column of dataset must be 'Observation'")}
     #
     ############################################################################
     # Initialization of data frame and of indices of fixed and random effects  #
     # data might be a groupedData object; regroup as fixdat with known formula #
     ############################################################################

     print("BEGINNING STEP 0", quote=FALSE)
     ########################################################################
     # Print a summary of the assumed analysis so that the user can be sure #
     # that the forwards search is done on the relevant analysis            #
     ########################################################################
     fixed2 <- fixed
     random2 <- random
     print("", quote = FALSE)
     print("The assumed analysis for these data will be as follows:", quote=FALSE)
     print("", quote = FALSE)
     zholdlm <- nlme::lme(fixed=fixed2, data=data, random=random2)                      #    lme
     print(summary(zholdlm))
     print("****************", quote = FALSE)
     fixdat <- nlme::groupedData(formula, data)
     #     
     ##################################
     # Reassign internal name of file #
     # Add small random amount        #
     ##################################
     df1 <- fixdat                #   fixdat not used subsequently; df1 is a groupedData object
     dimdf1 <- dim(df1)
     for(j in 2:dimdf1[2]){
          if(is.numeric(df1[,j])){
               thiscol <- df1[,j]
               absthiscol <- abs(thiscol)
               absthiscol <- absthiscol[absthiscol>0]
               mindf1 <- min(absthiscol)
               newcol <- thiscol +  stats::runif(dimdf1[1], min=-.001, max=.001)/mindf1
               df1[,j] <- newcol 
          }
     }
     data2 <- data.frame(df1)                        #  data2 is ungrouped data frame with random additions                      
     #######################################################
     # Retrieve independent and dependent columns from df1 #
     #######################################################
                          if(diagnose) print("Step 0 before lm")
     lmAlldata <- stats::lm(formula=fixed2, data=df1, singular.ok=TRUE, x=TRUE, y=TRUE)                               #   lm
     y1 <- lmAlldata$y
     x1 <- lmAlldata$x                                              # takes care of more complex fixed effect formulas
     coeffnames <- dimnames(x1)[2]
#prn(coeffnames)
     p <- dim(x1)[2]

     nrowsdf1 <- dim(df1)[1]
     rows.in.set <- vector("list", nrowsdf1)
     rows.in.set[[nrowsdf1]] <- -999

     ###################################################################################################
     # Gather data on subgroups as defined in (simplified) formula as applied to df1                   #
     # The subgroups don't need to be those of the analysis; it is only a diversified start for Step 1 #
     ###################################################################################################
     groups.df1 <- unique(nlme::getGroups(df1))
                                   if(diagnose){Hmisc::prn(groups.df1)}
     ngroups.df1 <- length(groups.df1)                                            # a vector of length 27
     df1.by.group <- vector("list", ngroups.df1)                                 # receptical list
     step1 <- vector("list", ngroups.df1)
     for(i in 1:ngroups.df1){
          uu <- df1[nlme::getGroups(df1)==groups.df1[i],]
          df1.by.group[[i]] <- uu
     }    # i
     #
     ############################################################################
     # STEP 1                                                                   #
     # p is the number of the independent variables, omitting 'Observation'     #
     # rows.in.set is a list whose elements are vectors of increasing length    #
     # residuals contains model residual for each term at each stage of the fit #
     # The residuals for i-th observation are in i-th row                       #
     # This step is conducted within groups with results left in group step1    #
     ############################################################################
     newcontrol <- nlme::lmeControl(maxIter=XmaxIter,  msMaxIter=XmsMaxIter, tolerance=Xtolerance, 
                      niterEM=XniterEM, msMaxEval=XmsMaxEval, msTol=XmsTol, opt=Xopt)
     newcontrol <<- newcontrol

     if(is.null(skip.step1)){
          print("BEGINNING STEP 1", quote=FALSE)
          randset <- 1:initial.sample
          step2 <- df1[1:2,]                 # data frame paradygm; first 2 rows will later be deleted
          #
          ############################################################
          # For each subgroup generate initial.sample randomizations #
          ############################################################
          for(k in 1:ngroups.df1){
               uu <- df1.by.group[[k]]
               dimuu1 <- dim(uu)[1]
                                        if(diagnose){Hmisc::prn(k);Hmisc::prn(dimuu1)} 
               zholdlist <- vector("list",initial.sample)
               for(i in 1:initial.sample){
                                        if(diagnose) Hmisc::prn(c(k,i))
                    zholdlist[[i]] <- matrix(0, nrow=dimuu1, ncol=3)
                    zholdlist[[i]][,1] <- sample(uu[,1],dimuu1)                               # sample permutation
                                        if(diagnose) Hmisc::prn(df1[zholdlist[[i]][1:robs,],])
               }    #   i
               step1[[k]] <- zholdlist
          }      # k
          #
          ########################################################
          # Combine these across groups into initial.sample sets #
          ########################################################
          topSamples <- vector("list",initial.sample)
          botSamples <- vector("list",initial.sample)
          for(i in 1:initial.sample){
               combtopSamples <- NULL
               combbotSamples <- NULL
               for(k in 1:ngroups.df1){
                    vv <- step1[[k]][[i]]
                    sizevv <- dim(vv)[1]
                    combtopSamples <- rbind(combtopSamples, vv[1:robs,])
                    combbotSamples <- rbind(combbotSamples, vv[(robs+1):sizevv,])
               }    #  k
               topSamples[[i]] <- combtopSamples
               botSamples[[i]] <- combbotSamples
          }         #  i
          #
          #######################################################
          # Determine which one to use as the result for Step 1 #
          #######################################################
          catchsamp <- vector("list", initial.sample)
          for(i in 1:initial.sample){
                                   if(diagnose) Hmisc::prn(i)
               uu <- topSamples[[i]]
               dimuu1 <- dim(uu)[1]
               if(2*round(dimuu1/2)==dimuu1){
                    amedian <- round(dimuu1/2) + 1
               }
               else{
                    amedian <- round(dimuu1/2)
               }
               zzzz <- df1[uu[,1],]
               zzzz <<- zzzz  
                                            # needed? directly tested and the answer is YES
                                               if(diagnose) {print("Step 1 before lme"); print(i); Hmisc::prn(zzzz)}
#               attach(zzzz)
               zholdlm <- nlme::lme(fixed=fixed2, data=zzzz, random=random2, control=newcontrol)                      #    lme
#               detach()
                                        if(diagnose) Hmisc::prn(zholdlm)
               dddd <- zholdlm$dims
               totalvars <- sum(unlist(dddd$ncol[1:dddd$Q]))
               temprand <- zholdlm$coefficients[[2]]
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

               allrows <- rbind(topSamples[[i]],botSamples[[i]])
               df2 <- df1[allrows[,1],]
               LMzholdlm <- stats::lm(formula=fixed2, data=zzzz)
               preds <- stats::predict(object=LMzholdlm, newdata=data2)   
               yinput <- df2[,response.column]
               err2 <- (preds - yinput)^2
               catchsamp[[i]] <- data.frame(allrows[,1],preds,yinput,round(err2,4),-999)     
               #
               ##################################
               # Determine median squared error #
               ##################################
               catchsamp[[i]] <- catchsamp[[i]][order(catchsamp[[i]][,4]),]
               catchsamp[[i]][,5] <- catchsamp[[i]][amedian,4]
          }      #    i
                             if(diagnose) Hmisc::prn(catchsamp)
          #
          catch <- catchsamp[[1]] 
          for(ii in 2:initial.sample){
               if(catchsamp[[ii]][1,5] < catch[1,5]){
                    catch <- catchsamp[[ii]]
               }
          }    #   ii
                           if(diagnose) Hmisc::prn(catch)
          ngo <- dim(topSamples[[1]])[1]
                             if(diagnose) Hmisc::prn(ngo)
          step2 <- df1[catch[1:ngo,1],]
          step2 <- step2[order(step2[,1]),]
                                  if(diagnose) Hmisc::prn(step2)  
          mstart <- dim(step2)[1] + 1
          rows.in.set[[mstart-1]] <- step2[,1]
                                  if(diagnose) Hmisc::prn(rows.in.set)
     }          # skip.step1 is null
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          stage <- length(skip.step1)
          rows.in.set[[stage]] <- skip.step1
          mstart <- stage + 1






               dddd <- zholdlm$dims
               totalvars <- sum(unlist(dddd$ncol[1:dddd$Q]))
               temprand <- zholdlm$coefficients[[2]]
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

























     }
     #
     ##################################################################
     # STEP 2                                                         #
     # Adding observations to the initial set                         #
     # Outer i loop adds 1 each time to list of observations in model #
     # May not be the same set of observations                        #
     ##################################################################

     print("BEGINNING STEP 2", quote=FALSE)

     ##################################################################################################################################### 
     # Set up files to hold residuals, dimensions, sigma, fixed parameter estimates, random parameter estimates, leverage, Cook distance #
     # and list for storage of xtemp matrices                                                                                            #
     ##################################################################################################################################### 
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

     for(i in mstart:(nrowsdf1+1)){                  # mstart is the step after the original obs entered
                                             if(diagnose) Hmisc::prn(i)
          rim <- rows.in.set[[i-1]]            # picks up rows for previous set of observations
          Zlatest <- df1[rim,]  
          Zlatest <<- Zlatest                  # allows lme to find file.    Needed?  directly tested and the answer is YES

          ############################################################################
          # Extract indep vars of the subset for use in leverage and Cook's distance #
          # lm function allows collection of x matrix and y vector                   #
          ############################################################################
          xtemp <- x1[rim,]
                                                 if(diagnose) {Hmisc::prn(xtemp)}
          xtemp.list[[i-1]] <- xtemp     
          transtemp <- t(xtemp)
          cross <- transtemp %*% xtemp
          crossinv <- solve(cross)
          getmatch <- df1[,1] %in% rim
          subdata <- data.frame(df1[getmatch,])     # removes grouping structure
          getthislm <- stats::lm(formula=fixed2, data=subdata, singular.ok=TRUE)                                  #   lm
          betahat <- getthislm$coefficients
          # Allowinig singular results means we have to replace NAs with 0 #
          if(any(is.na(betahat))){
               nabetahat <- is.na(betahat)
               betahat[nabetahat] <- 0
          }
          ################################################################
          # Capture results of lme analysis for this set of observations #
          ################################################################

                                    if(diagnose){print("before lme"); print(i); Hmisc::prn(Zlatest)}
          zholdlme <- nlme::lme(fixed=fixed2, data=Zlatest, random=random2, control=newcontrol)                   #   lme
                                    if(diagnose){ Hmisc::prn(Zlatest); Hmisc::prn(zholdlme)}
          LMzholdlme <- stats::lm(formula=fixed2, data=data.frame(Zlatest))                                       # lm

          resids <- zholdlme$residuals
          hold.coeffs.fixed[[i-1]] <- zholdlme$coefficients[[1]]

          temprand<-zholdlme$coefficients[[2]]
          gotvars <- rep(0,totalvars)
          for(level in 1:totalvars){
               uvars <- c(temprand[[level]])
               gotvars[level] <- sqrt(sum(uvars*uvars)/length(uvars))
          }
          hold.coeffs.random[i-1,] <- gotvars
          xbar <- mean(c(temprand[[1]]))
          devs <- c(temprand[[1]])-xbar
          sumsq <- mean(devs^2)
          RMS <- sqrt(sumsq)

          logLik <- zholdlme$logLik
          AIC <- summary(zholdlme)$AIC
          BIC <- summary(zholdlme)$BIC
          hold.summary.stats[i-1,] <- c(AIC, BIC, logLik)


          param.est[,i-1] <- c(zholdlme$coefficients[[1]])                        # same as holdcoeffs.fixed[[i-1]]
          t.set[,i-1] <- summary(zholdlme)$tTable[,4]
          hold.dims[[i-1]] <- zholdlme$dims
          hold.sigma[i-1] <- zholdlme$sigma
          newrows <- NULL
          potential <- NULL
          hold.subset.residuals[i-1] <- sum(resids^2)/(i-p)
          #
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
                                             if(diagnose) Hmisc::prn(leverage)
          #
          ###############################################################################
          # Manipulate remaining observations to determine the ones to go into next set #
          ###############################################################################
          for(j in 1:ngroups.df1){
                                                 if(diagnose) {Hmisc::prn(c(i,j)); Hmisc::prn(df1.by.group[[j]])}
               predlme <- stats::predict(object=LMzholdlme, newdata=df1.by.group[[j]])
               errs2 <- (predlme - df1.by.group[[j]][,response.column])^2 
               vv <- data.frame(df1.by.group[[j]]$Observation, predlme, df1.by.group[[j]][,response.column], errs2)
               vv <- vv[order(vv[,4]),]
               newrowsj <- vv[c(1,2),1]            # take first 2 rows of vv (smallest errs2) to put into newrows
               newrows <- c(newrows, newrowsj) 
               vv <- vv[c(-1,-2),]                 # drop first 2 rows of vv
               potential <- rbind(potential, vv)   #  concatenate remainder into a pool of potential observations for further selection across groups
          }                  #    j
          potential <- potential[order(potential[,4]),]
          addnew <- i - length(newrows)
          newrows <- c(newrows, potential[1:addnew,1])
          if(i < nrowsdf1 + 1) {rows.in.set[[i]] <- newrows}
          #
          ##########################################################
          # Gather errors from predictions to all data in database #
          ##########################################################
          td1 <- dim(df1)[1]
          errors <- rep(-999, td1)
          for(j in 1:td1){
               predlmeall <- stats::predict(object=LMzholdlme, newdata=df1)
               errors[j] <- y1[j] - sum(hold.coeffs.fixed[[i-1]] * x1[j,])
          }             #   j
                                             if(diagnose){Hmisc::prn(errors)}
          hold.residuals[,i-1] <- errors
     }                       # i in mstart ...
     #
     param.est <- as.data.frame(t(param.est))
     names(param.est) <- paste("b",1:p, sep="")
     m <- 1:nrowsdf1
     param.est <- cbind(m,param.est)
     #
     hold.summary.stats <- data.frame(m, hold.summary.stats) 
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
     nms <- dim(param.est)[1]
     param.est.current <- param.est[-1,]
     param.est.prev <- param.est[-nms,]
     param.diff <- param.est.prev - param.est.current
                         if(diagnose){
                             Hmisc::prn(param.est)
                             Hmisc::prn(param.est.prev)
                             Hmisc::prn(param.est.current)
                             Hmisc::prn(param.diff)
                             Hmisc::prn(xtemp.list)     
                         }
     for(i in mstart:nrowsdf1){
          aa <- param.diff[i-1,]
          aa <- as.numeric(aa[-1])
          aa <- matrix(aa,nrow=1)
          bb <- matrix(xtemp.list[[i-1]], nrow=p)
          www <- aa %*% bb
          modCook[i-1] <- (www %*% t(www))/(p * hold.subset.residuals[i-1])
     }
     #
     #######################################
     # Clean up files written to workspace #
     #######################################
     if(is.null(skip.step1))rm(list=c("zzzz", "Zlatest", "newcontrol"), pos=1)
     else rm(list=c("Zlatest","newcontrol"), pos=1)

     hold.coeffs.random <- as.data.frame(hold.coeffs.random)
     names(hold.coeffs.random) <- outnames
    #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_lme", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     list("Number of rows included in Step 1"=mstart-1,
           Subgroups=                         groups.df1, 
          "Rows by subgroup"=                 df1.by.group,
          "Rows in stage"=                    rows.in.set,
           Sigma=                             sigma,
          "Standardized residuals"=           hold.residuals,            
          "Fixed parameter estimates"=        param.est,
          "Random parameter estimates"=       hold.coeffs.random,
           Leverage=                          leverage[-1,],
          "Modified Cook distance"=           modCook,
           Dims=                              zholdlme$dims,
          "t statistics"=                     t.set,
          "Fit statistics"=                   hold.summary.stats,
           Call=                              MC )
}
