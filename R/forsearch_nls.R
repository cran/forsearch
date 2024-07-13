#' @export
forsearch_nls <-
function(phaselist, data, poolstart, poolformula, algorithm="default", controlarg=NULL, initial.sample=1000,  
              skip.step1=NULL, begin.diagnose=100, verbose=TRUE)
{
     #                          forsearch_nls
     #
     # VALUE      forsearch_nls does not use covariates besides those included in formula as nonlinear coefficients. See nlsList for a version
     #                with additional covariates. Considers all nls datasets to be in phases, possibly a single phase
     #
     # INPUT
     #     phaselist       Named list, each element of which contains 4 elements that describe 1 phase.
     #          formula           Formula for phase i
     #          formulacont       Formula omitting all factor variables for phase i  
     #          start             Start for each phase i
     #          nopp              nopp for the phase i
     #
     #     data              Dataset. First 2 variables are Observation and Phases (both mandatory)
     #     poolstart         Start for Step 2
     #     algorithm         See stats::nls
     #     control           See stats::nls
     #     initial.sample
     #     skip.step1
     #     unblinded
     #     begin.diagnose
     #     verbose

     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running forsearch_nls", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }

     #   begin.diagnose      Step 0: 1 - 19        Step 1: 20     -     49     Step 2:    50 - 59           Extraction:     81 - 
     #                                                eStep1:  31 - 39                    eStep2:  60 - 80

     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        forsearch_nls               "     # used for begin.diagnose prints
     nlsc <- nls.control()
     on.exit(expr=nls.control <- nlsc)
     if(!is.null(controlarg))nls.control <- controlarg
     carrycontrol <- nls.control
     #
     #####################################################################
     # Ensure that first variable of data is Observation and that second #
     # variable is Phases. Ensure that phaselist has correct elements    #
     # Make x1 the matrix of independent variables without obs number    #
     # x1 will represent the formula for independent obs                 #
     #####################################################################
     uu <- names(data)
     if(uu[1] != "Observation")stop("First column of data must be 'Observation'")
     #
     if(uu[2] != "Phases")stop("Second column of data must be 'Phases'")
     #
     if(!is.list(phaselist))stop("phaselist must be a list")
     #
     argnames <- names(phaselist[[1]])
     if(!{argnames[1]=="formula" &
          argnames[2]=="formulacont" &
          argnames[3]=="start" &
          argnames[4]=="nopp"} )stop("phaselist is incorrectly named or structured")
     #
     ###################################################################
     # Determine the name and variable number of the response variable #
     ###################################################################
     vv <- phaselist[[1]][[1]]
     lhsformula <- formula.tools::lhs(vv)
     nphases <- length(phaselist)
     ncoeffs <- vector("list", nphases)

                                               if(begin.diagnose <=1){ print(paste(spacer,"Section 1",sep=" "),quote=FALSE);
                                                      Hmisc::prn(phaselist);Hmisc::prn(dim(data));Hmisc::prn(lhsformula)   }

     dimdata <- dim(data)
     dataNames <- names(data)
     respnum <- lhsformula==dataNames
     ycolw <- (1:dimdata[2])[respnum]
     #
     #############################################################################
     # Create pooled variables for Step 2 from phaselist for single phase models #
     #############################################################################
     if(nphases==1){
          poolstart <- phaselist[[1]][[3]]
          poolformula <- phaselist[[1]][[1]]
     }
     ncoeffs <- length(poolstart)
     #
     ###################################################################################################
     # Create all files needed below:                                                                  #
     #    fixdat.df,   a data frame (data) with factor level indicator and containing all variables           #
     #    fixdat.list, a list with  #
     ###################################################################################################

     #########################################################
     # Check for factor status of data and get factor names  #
     #########################################################
     datacontrank <- NULL
     ufactor <- rep(TRUE, dimdata[2])
     for(m in 1:dimdata[2]) ufactor[m] <- is.factor(data[,m])
     yesfactor <- any(ufactor)     
     if(!yesfactor)stop("There must be at least 1 variable that is a factor (eg, Phases)")
     #
     # Add factor subset indicator #

     factorNames <- dataNames[ufactor]

     ############################################
     # Add the factor grouping code to the data #
     ############################################
     fixdat.df <- data
     isfactor <- rep(TRUE,dimdata[2])
     for(nn in 1:dimdata[2]){
          isfactor[nn] <- is.factor(fixdat.df[,nn])
     }
     justfactors <- fixdat.df[isfactor]
     holdISG <- apply(justfactors, 1, paste,collapse="/")
     holdISG <- paste("_",holdISG, sep="")
     fixdat.df <- data.frame(fixdat.df, holdISG)                                                 # fixdat.df    building

                                               if(begin.diagnose <=5){ print(paste(spacer,"Section 5",sep=" "),quote=FALSE);
                                                      Hmisc::prn(utils::head(fixdat.df));Hmisc::prn(utils::tail(fixdat.df));Hmisc::prn(dim(fixdat.df))   }

     #########################################################################
     # Create a list of the dataset by factor subset levels and a data frame #
     # of number of observations needed in Step 1 for each subset            #
     #########################################################################
     ufixdatISG <- unique(fixdat.df$holdISG)
     nlevfixdat <- length(ufixdatISG)
     fixdat.list <- vector("list", nlevfixdat)
     nobs <- rep(-9,nlevfixdat)
     nobs.df <- data.frame(ufixdatISG, nobs)
     for(jj in 1:nlevfixdat){
          fixdat.list[[jj]] <- fixdat.df[fixdat.df$holdISG==ufixdatISG[jj],]
     }
                                               if(begin.diagnose <=7){ print(paste(spacer,"Section 7",sep=" "),quote=FALSE);
                                                      Hmisc::prn(fixdat.list)   }
     #
     #############################################################
     # Determine how many obs are needed to get convergence with #
     # and add these to fixdat.df                                #
     #############################################################
     np <- names(phaselist)
     phaseind <- 1:nphases
     np.df <- data.frame(np,phaseind)
     for(i in 1:nlevfixdat){
          thissubset <- fixdat.list[[i]]
          dimtss <- dim(thissubset)[2]
          thisphase <- thissubset[1,2]
          mat <- match(thisphase,np)
          uu <- phaselist[[mat]]
          thisformulacont <- uu[[2]]
          thisstart <- uu[[3]]  
          thisnopp <- uu[[4]]
          for(tryn in max(3,thisnopp):dim(thissubset)[1]){

                                               if(begin.diagnose <= 9){ print(paste(spacer,"Section 9",sep=" "),quote=FALSE);
                                                      Hmisc::prn(thisphase);Hmisc::prn(thisformulacont);Hmisc::prn(tryn)   }
               an.error.occurred <- FALSE
               OKnls <- list(convInfo <- FALSE)

               tryCatch(
                   expr=(OKnls <- stats::nls(formula=thisformulacont, data=thissubset[1:tryn,], start=thisstart, 
                       algorithm = algorithm))
                   , error=function(e) {an.error.occurred <- TRUE}
                        )
               if(an.error.occurred){
                    messprop <- "nls failed to converge; increasing nobs in estimation of Step 1"
                    print(messprop, quote=FALSE)
               }
               else{
                    if(OKnls$convInfo[[1]])break
               }

          }    # tryn
          nobs.df[i,1] <- thissubset[1,dimtss]
          nobs.df[i,2] <- tryn
     }         # i
     maxfirst <- 0
     fixdat.df <- data.frame(fixdat.df,maxfirst)

     for(i in 1:dim(fixdat.df)[1]){
          index <- nobs.df[,1]==fixdat.df$holdISG[i]
          fixdat.df$maxfirst[i] <- nobs.df$nobs[index]      
     }         # i
                                               if(begin.diagnose <=11){ print(paste(spacer,"Section 11",sep=" "),quote=FALSE);
                                                      Hmisc::prn(nobs.df);Hmisc::prn(OKnls)   }

     dimx <- dim(data)
     dimx1 <- dimx[1]
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
     #
######################################################################################################################################
# Step 1
     #################################################################################
     # Repeat the general structure for Step 1 for each factor subset of the data.   #
     # Then combine the output into a single set of observation numbers.             #
     #################################################################################

     medaugx <- matrix(1, nrow=initial.sample, ncol=2)                       #  median of augx
     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(x=1:dimx1, size=dimx1)                        #    sample permutation
     }      #   i
     if(is.null(skip.step1)){
          print(" ",quote=FALSE)
          print("ENTERING STEP 1", quote=FALSE)
          rim <- NULL
          for(ip in 1:nlevfixdat){ 
               thisdata <- fixdat.list[[ip]]
               thisphase <- thisdata[1,2]
               thissubset <- thisdata$holdISG[1]     # name of subset with leading _
               index <- nobs.df[,1]==thissubset
               nobs <- nobs.df[index,2]
               index2 <- np.df$phaseind[np.df$np==thisphase]
               formulacont <- phaselist[[index2]][[2]] 
               thisstart <- phaselist[[index2]][[3]]  
               namesthisdata <- names(thisdata)
               ww <- lhsformula==namesthisdata
               nnames <- dim(thisdata)[2]
               ycolw <- (1:nnames)[ww]    

                                 if(begin.diagnose <= 21){print("", quote = FALSE);print(paste(spacer,"Section 21",sep=" "),quote=FALSE);
                                     Hmisc::prn(thisdata);Hmisc::prn(thisphase);Hmisc::prn(thissubset);Hmisc::prn(nobs);
                                     Hmisc::prn(formulacont);Hmisc::prn(thisstart);Hmisc::prn(ycolw)     }

              firstrim <- eStep1(df1=thisdata, inner.rank=nobs, initial.sample=initial.sample, formula=formulacont, start = thisstart,
                   contR=carrycontrol, algo=algorithm, ycol=ycolw, b.d = begin.diagnose)                                                  #  eStep1

              rim <- c(rim, firstrim)
          }     # ip

          lenrim <- length(rim)
          rows.in.model[[lenrim]] <- sort(rim)
          SOON <- rim
          mstart <- length(SOON) + 1

     }         # is.null skip.step1
     else{
          print("SKIPPING STEP 1", quote=FALSE)
          nfirst <- length(skip.step1)
          rows.in.model[[nfirst]] <- sort(skip.step1)
          datastep1 <- data[as.vector(rows.in.model[[nfirst]]),]
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
          SOON <- sort(skip.step1)
     }       # skipping step 1
                                 if(begin.diagnose <= 48){print("", quote = FALSE);print(paste(spacer,"Section 48",sep=" "),quote=FALSE);
                                            Hmisc::prn(rows.in.model);Hmisc::prn(SOON)      }
    #
#############################################################################################################################
# Step 2
     # Step 2 is conducted entirely within eStep2, but with setup here #
     #     phaselist       Named list, each element of which is a list describing 1 phase. Each element contains a list with elements 
     #          formula           Formula for phase
     #          formulacont       Formula omitting all factor variables for phase  
     #          start             Start for each phase
     #          nopp              nopp for the phase

     print(" ",quote=FALSE)
     print("ENTERING STEP 2", quote=FALSE)

                                 if(begin.diagnose <= 50){print("", quote = FALSE);print(paste(spacer,"Section 50",sep=" "),quote=FALSE);
                                            Hmisc::prn(poolformula)      }

     eStep2out <- eStep2(mstart=mstart, finalm=rows.in.model, start=poolstart, data2=data, pformula=poolformula, algo2=algorithm,    
              fixdb=fixdat.df, fixlist=fixdat.list, contR=carrycontrol, ycol=ycolw, b.d= begin.diagnose)                                                      #  eStep 2

     rows.in.model <- eStep2out[[1]]
     fooResult <- eStep2out[[2]]
     resids2 <- eStep2out[[3]]
     sigma <- eStep2out[[4]]       # Note that this is the next-to-last estimate in eStep2
                                 if(begin.diagnose <= 80){print("", quote = FALSE);print(paste(spacer,"Section 80",sep=" "),quote=FALSE);
                                            Hmisc::prn(rows.in.model);Hmisc::prn(fooResult)      }
     #
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
          getthisnls <- fooResult[[i]]
          sumnls <- summary(getthisnls)
          resgetthis <- residuals(getthisnls)
          s.2[i] <- sum(resgetthis^2)
          ###################################################
          # Get coefficients from summary of each fooResult #
          ###################################################
          sumnls <- summary(getthisnls)
          nlscoeffs <- sumnls[[10]][,1]
          param.est[,i] <- nlscoeffs
          t.set[,i] <- sumnls[[10]][,3]
     }        # i extracting statistics now completed
     #
     ###########################################
     # Cleanup after all extractions completed #
     ###########################################
     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m", "Observation", "leverage"))
     param.est <- as.data.frame(t(param.est))
     names(param.est) <- names(poolstart)
     m <- 1:dimx1
     param.est <- cbind(m,param.est)                                 #  here we add the m column
     # 
     t.set <- as.data.frame(t(t.set))
     names(t.set) <- names(poolstart)
     t.set <- cbind(m,t.set)
     #
     ############################################
     # Estimate sigma and standardize residuals #
     ############################################
     residuals <- resids2/sigma
     residuals <- residuals[,-dimx1]
     #
     listout <- list(

          "Step 1 observation numbers"=        SOON,
          "Rows in stage"=                     rows.in.model,
          "Standardized residuals"=            residuals, 
          "Number of model parameters"=        dim(t.set)[2]-1,
          "Sigma"=                             sigma, 
          "Fixed parameter estimates"=         param.est, 
          "s^2"=                               s.2, 
#          "Leverage"=                          leverage, 
          "t statistics"=                      t.set,
          "Call"=                              MC)

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running forsearch_nls", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(listout)
}
