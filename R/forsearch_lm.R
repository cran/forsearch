#' @export
forsearch_lm <-
function(formula, data, initial.sample=1000, n.obs.per.level=1, skip.step1=NULL, unblinded=TRUE, begin.diagnose=100, verbose=TRUE)
{
     #                          forsearch_lm
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running forsearch_lm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }

     #   begin.diagnose      Step 0: 1 - 19        Step 1: 20     -     49     Step 2:    50 - 59           Extraction:     81 - 
     #                                                aStep1:  31 - 39                    aStep2:  60 - 80

     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        forsearch_lm               "     # used for begin.diagnose prints

     ##################################################################
     # Ensure that first independent variable is Observation          #
     # Make x1 the matrix of independent variables without obs number #
     # x1 will represent the formula for independent obs              #
     ##################################################################
     uu <- names(data)
     if(uu[1] != "Observation")stop("First column of data must be 'Observation'")
     #
     #
     #######################################################################
     # Print the structure of the analysis that will be done on these data #
     # unless the treatment group is blinded                               #
     #######################################################################
     options(warn = -1)
     on.exit(options(warn = 0))

     nulldata <- data
     lhsformula <- formula.tools::lhs(formula)
     rhsformula <- formula.tools::rhs(formula)
     ncols <- dim(nulldata)[2]
     ycol <- (1:ncols)[uu==lhsformula]
     nulldata[,ycol] <- 0
     print("", quote = FALSE)
     lmAlldata <- stats::lm(formula, nulldata, singular.ok=TRUE, x=TRUE, y=TRUE)                          # lm
     coeffnames <- names(lmAlldata$coefficients)
     lmAlldata$df.residual <- 99999
     ########################################
     # Get the number of terms in the ANOVA #
     # and set up the holding matrix        #
     ########################################
     AVlm <- stats::anova(lmAlldata)
     dnAVlm <- dimnames(AVlm)                               # level 1 is names of sources of variation, incl residuals
     dnAVlmM1 <- dnAVlm[[1]]
     dnAVlmM1 <- dnAVlmM1[-length(dnAVlmM1)]           # chop off last one 
     anova.pvals <- matrix(0, nrow=length(dnAVlmM1), ncol=dim(data)[1])
     #
     if(unblinded){
          print("**************************************", quote = FALSE)
          print("", quote = FALSE)
          print("Check the following ANOVA paradigm for source of variation and associated degrees of freedom (Df)", quote=FALSE)
          print("", quote = FALSE)
          ux2 <- stats::anova(lmAlldata)
          nux2 <- dim(ux2)[1]
          dimnamesux2 <- dimnames(ux2)
          dn1 <- dimnamesux2[[1]]
          dn1 <- dn1[-nux2]
          ux2 <- ux2[-nux2,]
          Df <- as.character(ux2[,1])
          ux3 <- data.frame(Df, row.names=dn1)
          print(ux3, , quote=FALSE)
          print("", quote = FALSE)
          print("All categorical variables must be defined to be factors in the database, not in the formula.", quote=FALSE)
          print("", quote = FALSE)
          print("**************************************", quote = FALSE)
     }
     y1 <- data[,ycol]
     #
     ################################################################
     # Check for constructed variables in formula, ie, use of I()   #
     # First convert formula to a vector of character pairs. Then   #
     # recode the I( letters as  I(A) and the test accordingly.     #      
     # Determine whether any of these is 'I(A)'.  If so, count them #
     ################################################################
     nAsIs <- 0
     charform <- as.character(formula)
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
     #
     #############################################################
     # Check for factor status of dataset, redefining lmAlldata  #
     # Remove factor variables and get rank of linear regression #
     # Locate response column in reduced data frame              #
     #############################################################
     datacontrank <- NULL
     ndatacols <- dim(data)[2]
     ufactor <- rep(TRUE, ndatacols)
     for(m in 1:ndatacols) ufactor[m] <- is.factor(data[,m])
     yesfactor <- any(ufactor)     
     ncoeffs <- length(lmAlldata$coefficients)
     rnk <- lmAlldata$rank
     nopl <- n.obs.per.level 

     ######################################
     # remove factor variables and get rank
     ######################################
     datacont <- data[,!ufactor]
     ufactornames <- names(datacont)          # includes Observation and response
     ufactornames <- ufactornames[-1]         # includes response
     contform.rhs <- ufactornames
     contform.rhs <- contform.rhs[contform.rhs != lhsformula]
     if(length(ufactornames)==1){             # there are no continuous independent variables
          formulacont <- paste(lhsformula, "1", sep=" ~ ")
          formulacont <- stats::formula(formulacont)
          datacontrank <- 1
     }
     else{
          if(length(contform.rhs) > 1){
               contform.rhs <- paste(contform.rhs, collapse=" + ")
          }
          formulacont <- paste(lhsformula, contform.rhs, sep= " ~ ")
          formulacont <- stats::formula(formulacont)
          lmdatacont <- stats::lm(formulacont,datacont)                                                # lm on data without factor variables
          datacontrank <- lmdatacont$rank
     }
                                 if(begin.diagnose <= 8){print("", quote = FALSE);print(paste(spacer,"Section 8",sep=" "),quote=FALSE);
                                      Hmisc::prn(utils::head(datacont));Hmisc::prn(ufactornames);Hmisc::prn(formulacont);
                                      Hmisc::prn(datacontrank);Hmisc::prn(nAsIs)       }

     if(yesfactor){                               # why pick some here and not in aStep1 ?  Seems to need this
          # there are factors in the dataset
          ss77 <- variablelist(datadf = data, prank=ndatacols)
          pickm <- picksome(subsetlist=ss77, nobs=dim(data)[1], initial.sample = initial.sample, n.obs.per.level=nopl, rank=datacontrank)
          dimpickm <- dim(pickm)[2]  
     }    # yesfactor
     #
     #################################################
     # Set up all output and intermediate structures #
     #################################################
     randset <- 1:initial.sample
     dimx <- dim(data)
     dimx1 <- dimx[1]

     x1 <- data[,-c(1,ycol)]
     if(!is.data.frame(x1))x1 <- data.frame(x1)
     OBS <- data[,1]
     Z <- data                         #  ??????
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result
     result <- matrix(0, nrow=dimx1, ncol=2)
     rows.in.model <- vector("list", dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=ncoeffs, ncol=dimx1)
     t.set <- param.est

     xtemp.list <- vector("list",dimx1)
     modCook <- rep(0,dimx1)
     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)
######################################################################################################################################
# Step 1
     #################################################################################
     # Step 1 of the procedure follows.                                              #
     # By creating an index, randomly reorder the rows of matrix Z = (X,y).          #
     # Within each element of zlist, calculate the several estimates of b from the   #
     # first p of these rows and calculate the median of the residuals in each set.  #
     #################################################################################
     medaugx <- matrix(1, nrow=initial.sample, ncol=2)                       #  median of augx
     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(x=1:dimx1, size=dimx1)                        #    sample permutation
     }      #   i
                                 if(begin.diagnose <= 20){print("", quote = FALSE);print(paste(spacer,"Section 20",sep=" "),quote=FALSE);
                                     Hmisc::prn(zlist);       }

     if(is.null(skip.step1)){
          print("ENTERING STEP 1", quote=FALSE)
          inner.rank <- datacontrank
          ycolcont <- names(datacont)==lhsformula
          ycolcont <- (1:length(ycolcont))[ycolcont]
 
          if(yesfactor){
               ss77C <- vector("list", length(ss77))
               for(nn in 1:length(ss77)){
                    uv <- ss77[[nn]][,1]
                    ss77C[[nn]] <- datacont[uv,]
               }
                                 if(begin.diagnose <= 23){print("", quote = FALSE);print(paste(spacer,"Section 23a",sep=" "),quote=FALSE);
                                     Hmisc::prn(utils::head(datacont));Hmisc::prn(datacontrank);Hmisc::prn(ycolcont)       }

              firstrim <- aStep1(yesfactor, df1=ss77C, inner.rank=datacontrank + nAsIs, initial.sample, formula=formulacont, 
                       ycol=ycolcont, nopl=n.obs.per.level, b.d = begin.diagnose)                                                         #aStep1
              firstrim <- firstrim[[1]]
              lenfr <- length(firstrim)
              rows.in.model[[lenfr]] <- firstrim  
              SOON <- firstrim
              mstart <- length(SOON) + 1
          }          # factors present
          else{
                                 if(begin.diagnose <= 23){print("", quote = FALSE);print(paste(spacer,"Section 23b",sep=" "),quote=FALSE);
                                     Hmisc::prn(utils::head(data));Hmisc::prn(formula);Hmisc::prn(datacontrank);Hmisc::prn(ycol)       }

              firstrim <- aStep1(yesfactor, df1=data, inner.rank=datacontrank + nAsIs, initial.sample, formula=formula, 
                       ycol=ycol, nopl=n.obs.per.level, b.d = begin.diagnose)                                                             #aStep1
               lenfr <- length(firstrim)
               rows.in.model[[lenfr]] <- firstrim  
               SOON <- firstrim
               mstart <- lenfr + 1
          }        # no factor present

     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          nfirst <- length(skip.step1)
          rows.in.model[[nfirst]] <- skip.step1
          datastep1 <- data[as.vector(rows.in.model[[nfirst]]),]
          lmchosen <- stats::lm(formula, datastep1, singular.ok=TRUE)                                  #  lm 
          betahat <- lmchosen$coefficients
          betaMMinus1 <- betahat
          randset <- 1:initial.sample
          medaugx <- cbind(medaugx,randset)
          minmed <- min(medaugx[,2])
          locatemin <- medaugx[medaugx[,2]==minmed,]
          if(is.matrix(locatemin)) locatemin <- locatemin[1,]
          zliststar <- zlist[[locatemin[3]]]

                                 if(begin.diagnose <= 46){print("", quote = FALSE);print(paste(spacer,"Section 46",sep=" "), quote=FALSE);
                                      Hmisc::prn(medaugx);Hmisc::prn(minmed);
                                      Hmisc::prn(locatemin);Hmisc::prn(zliststar);Hmisc::prn(zliststar[1:nfirst,1])       }

          mstart <- nfirst + 1
          betaMMinus1 <- betahat
          SOON <- skip.step1
     }       # skipping step 1
                                 if(begin.diagnose <= 48){print("", quote = FALSE);print(paste(spacer,"Section 48",sep=" "),quote=FALSE);
                                            Hmisc::prn(rows.in.model);Hmisc::prn(SOON)      }
    #
# stop("before step 2")
#############################################################################################################################
# Step 2
     print("ENTERING STEP 2", quote=FALSE)
     thisy <- ycol
     if(yesfactor){
          aStep2out <- aStep2(yesfactor, form.A2=formula, finalm=rows.in.model, rimbs=ss77, dfa2=data, ycol=thisy, mstart=mstart, 
                        rnk=datacontrank+nAsIs, b.d=begin.diagnose)                                                    #   aStep2
     }       #  factors are present
     else{
          nofactrank <- rnk
          aStep2out <- aStep2(yesfactor, form.A2=formula, finalm=rows.in.model, rimbs=NULL, dfa2=data, ycol=thisy, mstart=mstart, 
                        rnk=nofactrank+nAsIs, b.d=begin.diagnose)                                                    #   aStep2
     }      # no factors present

     rows.in.model <- aStep2out[[1]]
     fooResult <- aStep2out[[2]]
     resids2 <- aStep2out[[3]]
     sigma <- aStep2out[[4]]
                                 if(begin.diagnose <= 80){print("", quote = FALSE);print(paste(spacer,"Section 80",sep=" "),quote=FALSE);
                                            Hmisc::prn(rows.in.model);Hmisc::prn(fooResult)      }
#
# stop("before extract stats")
##############################################################################################################################
# EXTRACT STATS
     print("", quote=FALSE)
     print("EXTRACTING INTERMEDIATE STATISTICS", quote=FALSE)
     print("", quote=FALSE)
     
     betahatset <- matrix(0,nrow=dimx1, ncol=ncoeffs)

     ##########################################
     # Extract statistics from each lm object #
     ##########################################
     for(i in mstart:dimx1){
          getthislm <- fooResult[[i]]
          if(i > mstart & i <= dimx1){
               AVlm <- stats::anova(getthislm)
               AVlmps <- AVlm[,5]
               AVlmps <- AVlmps[-length(AVlmps)]
               anova.pvals[,i] <- AVlmps
          }  # i > mstart::nobs 

          getthislmx <- getthislm$x
          xtemp.list[[i]] <- getthislmx                          # used in modified Cook
          betahat <- c(getthislm$coefficients, rep(0,ncoeffs))
          betahat <- betahat[1:ncoeffs]
          model.resids <- getthislm$residuals

          betahatset[i-1,] <- betahat
          medaugx <- matrix(1:dimx1, nrow=dimx1, ncol=2, byrow=FALSE)      # initialize medaugx
          medaugx[,2] <- 0 
          param.est[,i-1] <- betahat
          if(i > rnk) s.2[i-1] <- sum(model.resids * model.resids)/(i-rnk)
          #
          #################################
          # Extract t values from summary #
          #################################
          if(i > mstart){
               t.setbase <- c(summary(getthislm)$coefficients[,3], rep(0,rnk))
               t.set[,i-1] <- t.setbase[1:rnk]
          }
          #
          ############
          # Leverage #
          ############
          thisleverage <- 1
          x1 <- data[,-c(1,ycol)]
          xtemp <- getthislmx
                                 if(begin.diagnose <= 81){print("", quote = FALSE);print(paste(spacer,"Section 81",sep=" "),quote=FALSE);
                                       Hmisc::prn(dimx1);Hmisc::prn(i);Hmisc::prn(xtemp)       }

          vvv <- eigen(t(xtemp) %*% xtemp)[[1]]
          uuu <- prod(vvv)
                                 if(begin.diagnose <= 82){print("", quote = FALSE);print(paste(spacer,"Section 82",sep=" "),quote=FALSE);
                                           Hmisc::prn(i);Hmisc::prn(vvv);Hmisc::prn(uuu)       }

          if(uuu > .Machine$double.eps*10000){
               crossinv <- solve(t(xtemp) %*% xtemp)
               for(j in 1:(i-1)){
                   Zlatest2 <- data[rows.in.model[[i-1]],]
                   if(dim(xtemp)[2]==1){
                        thisleverage <- c(c(matrix(xtemp[j],nrow=1) %*% crossinv %*% matrix(xtemp[j],ncol=1)))
                        thisleverage <- c(i-1,Zlatest2[j,1],thisleverage)
                   }else{
                        thisleverage <- c(c(matrix(xtemp[j,],nrow=1) %*% crossinv %*% matrix(xtemp[j,],ncol=1)))
                        thisleverage <- c(i-1,Zlatest2[j,1],thisleverage)
                   }
                   leverage <- rbind(leverage,thisleverage)
               }   # j 1:i

                                 if(begin.diagnose <= 83){print("", quote = FALSE);print(paste(spacer,"Section 83",sep=" "),quote=FALSE);
                                        Hmisc::prn(i);Hmisc::prn(leverage)       }

          }    # if uuu
          else{print("Leverage not calculated because x matrix singular", quote=FALSE)}
          #
     }        # i extracting statistics now completed
     #
     ###########################################
     # Cleanup after all extractions completed #
     ###########################################
     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m", "Observation", "leverage"))
     param.est <- as.data.frame(t(param.est))
     names(param.est) <- coeffnames
     m <- 1:dimx1
     param.est <- cbind(m,param.est)                                 #  here we add the m column
     # 
     anova.pvals <- as.data.frame(t(anova.pvals))
     names(anova.pvals) <- dnAVlmM1
     anova.pvals <- cbind(m,anova.pvals)
     #
     t.set <- as.data.frame(t(t.set))
     names(t.set) <- coeffnames
     t.set <- cbind(m,t.set)
     #
     ############################################
     # Estimate sigma and standardize residuals #
     ############################################
     residuals <- resids2/sigma
     #
     ##########################
     # Modified Cook distance #
     ##########################
     nms <- dim(param.est)[1]
     param.est.prev <- param.est[-nms,]
     param.est.current <- param.est[-1,]          # remove first and last rows
     param.diff <- param.est.prev - param.est.current

                                 if(begin.diagnose <= 85){print("", quote = FALSE);print(paste(spacer,"Section 85",sep=" "),quote=FALSE);Hmisc::prn(param.est.prev);
                                       Hmisc::prn(param.est.current);Hmisc::prn(param.diff);Hmisc::prn(xtemp.list)   }

     for(ir in mstart:(dimx1-1)){
         aa <- param.diff[ir-1,]         # select row ir-1
         aa <- as.numeric(aa[,-1])      # remove col 1
         aa <- matrix(aa,nrow=1)        # turn this into a matrix with 1 row
         bb <- xtemp.list[[ir]]

                                 if(begin.diagnose <= 86){print("", quote = FALSE);print(paste(spacer,"Section 86",sep=" "),quote=FALSE);
                                       Hmisc::prn(ir);Hmisc::prn(aa);Hmisc::prn(bb)     }
         www <- aa %*% t(bb)
         modCook[ir] <- (www %*% t(www))/(rnk*s.2[ir])

                                 if(begin.diagnose <= 89){print("", quote = FALSE);print(paste(spacer,"Section 89",sep=" "),quote=FALSE);
                                                   Hmisc::prn(ir);Hmisc::prn(www);Hmisc::prn(modCook[ir])       }

     }       #    ir
      modCook <- modCook[-length(modCook)]
     #
     listout <- list(

          "Step 1 observation numbers"=        SOON,
          "Rows in stage"=                     rows.in.model,
          "Standardized residuals"=            residuals, 
          "Number of model parameters"=        rnk,
          "Sigma"=                             sigma, 
          "Fixed parameter estimates"=         param.est, 
          "s^2"=                               s.2, 
          "Leverage"=                          leverage, 
          "ANOVA"=                             anova.pvals,
          "Modified Cook distance"=            modCook, 
          "t statistics"=                      t.set,
          "Call"=                              MC)

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_lm", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(listout)
}
