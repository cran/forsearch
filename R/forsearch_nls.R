#' @export
forsearch_nls <-
function(nlsform, data, start, algorithm="default", 
           nls.control=FALSE, initial.sample=1000, skip.step1=NULL, begin.diagnose=100, verbose=TRUE)
{
     #                                       forsearch_nls
     #
     # VALUE      forsearch_nls does not use covariates besides those included in formula as nonlinear coefficients.
     #                 See nlsList for a version with additional covariates. 
     #
     # INPUT
     #     nlsform           Formula
     #     nofactform        Formula omitting all factor variables, NULL if no factors
     #     data              Dataset. First  variable is Observation; second variable is Section
     #     start             Start
     #     algorithm         See stats::nls
     #     nls.control       Logical See stats::nls
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
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX        forsearch_nls       "

     carrycontrol <- nls.control
     #
     # Following 2 lines will be used in later version
     nofactform <- NULL
     nofactstart <- NULL
     #####################################################
     # Ensure that first variable of data is Observation #
     # Ensure that second variable of data is Section    #
     #####################################################
     uu <- names(data)
     if(uu[1] != "Observation")stop("First column of data must be 'Observation'")
     if(uu[2] != "Section") stop("Second column of data must be 'Section'")

                                               if(begin.diagnose <=1){ print(paste(spacer,"Section 1",sep=" "),quote=FALSE);
                                                      Hmisc::prn(nlsform);Hmisc::prn(dim(data));Hmisc::prn(start)   }

     #######################################################################################################
     # Create all files needed below:                                                                      #
     #    fixdat.df,        a data frame (data) with factor level indicator and section indicator          #
     #    fixdatfact.list   a 1-layer list with the same data as fixdat.df but by factor subset            #
     #    fixdatouter.list, a 2-layer list with the same data as fixdat.df but with [[factor]] [[setion]]  #
     #                      If no factors, the one 'factor' is NONE                                        #
     #    fixdatcombo.list, a 1-layer list with elemnts that combine section and factor subset             #
     #######################################################################################################
     #
     ###################################################################
     # Determine the name and variable number of the response variable #
     ###################################################################
     fixdat.df <- data
     dimdata <- dim(fixdat.df)
     namesfixdat <- names(fixdat.df)
     thisenv <- rlang::env()
     xform <- stats::formula(nlsform, env=thisenv)
                    # nlstest05 <- nls(formula=xform, data=fixdat.df, start=start)
                    # prn(nlstest05)
                    # This all worked. Leave the code here
     lhsformula <- formula.tools::lhs(xform)
     ycolw <- namesfixdat==lhsformula
     ycol <- (1:dimdata[2])[ycolw]      
     #
     ####################################
     # Add section code to the data     #
     # There is only 1 section variable #
     ####################################
     sectISG <- paste("_",fixdat.df$Section, sep="")
     usect <- unique(sectISG)
     nusect <- length(usect)
     fixdat.df <- data.frame(fixdat.df,sectISG)
     #
     #########################################################
     # Check for factor status of data and get factor names  #
     #########################################################
     datacontrank <- NULL
     ufactor <- rep(TRUE, dimdata[2])
     for(m in 1:dimdata[2]) ufactor[m] <- is.factor(data[,m])
     yesfactor <- any(ufactor)     
     nufactor <- sum(ufactor)        # This is the number of factor variables, not the number of levels
     #
     ######################################
     # Add factor subset indicator        #
     # regardless whether factors present #
     ######################################
     if(yesfactor){

          stop("nls does not recognize factor variables. Program aborted.")

          factorNames <- namesfixdat[ufactor]
          ###################################
          # Add the factor code to the data #
          ###################################
          isfactor <- rep(TRUE,dimdata[2])
          for(nn in 1:dimdata[2]){
               isfactor[nn] <- is.factor(fixdat.df[,nn])
          }
          justfactors <- fixdat.df[isfactor]  
          holdISG <- apply(justfactors, 1, paste,collapse="/")
          holdISG <- paste("_",holdISG, sep="")
          uholdfactors <- unique(holdISG)
          nuholdfactors <- length(uholdfactors)
     }
     else{
          factorNames <- "NONE"
          holdISG <- "_NONE"
          uholdfactors <- holdISG
          nuholdfactors <- 1
     }
     fixdat.df <- data.frame(fixdat.df, holdISG)                                                 # fixdat.df    building

                          if(begin.diagnose <= 5){ print(paste(spacer,"Section 5",sep=" "),quote=FALSE);
                                    Hmisc::prn(utils::head(fixdat.df));Hmisc::prn(utils::tail(fixdat.df));
                                    Hmisc::prn(dim(fixdat.df));print("end of  _nls    bd 5")   }
     #
     ##########################
     # Create fixdatfact.list #
     ##########################    
     fixdatfact.list <- vector("list",nuholdfactors)
     for(i in 1:nuholdfactors){
          fixdatfact.list[[i]] <- fixdat.df[fixdat.df$holdISG==uholdfactors[i],]
     }
     #
     #########################################################
     # Add a code for the combination of section and factor) #
     #########################################################
     comboISG <- paste(fixdat.df$sectISG, fixdat.df$holdISG, sep="_S/F")
     fixdat.df <- data.frame(fixdat.df, comboISG)
     ucombolevels <- unique(comboISG)
     nucombolevels <- length(ucombolevels)
     #
     ###########################################
     # Create fixdatouter.list and populate it #
     ###########################################
     fixdatouter.list <- vector("list", nuholdfactors)
     dummylist <- vector("list", nusect)
     for(i in 1:nuholdfactors){
          fixdatouter.list[[i]] <- dummylist
     }    # i
     for(k in 1:nuholdfactors){
          for(j in 1:nusect){
               uu <- fixdat.df[fixdat.df$holdISG==uholdfactors[k], ]
               thesesectlevels <- unique(uu$sectISG)
               uu <- uu[uu$sectISG==usect[j],]
               fixdatouter.list[[k]][[j]] <- uu
          }   #    j
     }        #    k
     #

     #########################################################
     # Determine number of observations to take in Step 1 by #
     # factor, if present and overall in every case          #
     #########################################################
     bignls <- stats::nls(formula=nlsform, data=fixdat.df, start=start)      #   nls
     ncoeffs.outer <- length(stats::coef(bignls))
     if(yesfactor){
          basenumber <- max(1,floor(ncoeffs.outer/nuholdfactors + .00001)   )
          additional <- rep(0, nuholdfactors)
          nadd <- ncoeffs.outer - basenumber * nuholdfactors
          if(nadd > 0)additional[1:nadd] <- 1
          inner.r <- rep(basenumber, times=nuholdfactors) + additional
          inner.r <- inner.r[1:nuholdfactors]
     }    #  if yessfactor
     dimx <- dim(data)
     dimx1 <- dimx[1]
     zlist <- vector("list",initial.sample)         # zlist elements start with matrix result

     result <- matrix(0, nrow=dimx1, ncol=2)
     rows.in.model <- vector("list", dimx1)
     residuals <- matrix(0,nrow=dimx1,ncol=dimx1)
     param.est <- matrix(0,nrow=dimx1, ncol=ncoeffs.outer)
     t.set <- param.est

     xtemp.list <- vector("list",dimx1)
     modCook <- rep(0,dimx1)
     s.2 <- rep(0,dimx1)
     leverage <- matrix(0,nrow=1,ncol=3)
     #

#stop("before step 1")
     #
######################################################################################################################################
# Step 1
     medaugx <- matrix(1, nrow=initial.sample, ncol=2)                       #  median of augx
     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(x=1:dimx1, size=dimx1)                        #    sample permutation
     }      #   i
     if(is.null(skip.step1)){
          print(" ",quote=FALSE)
          print("ENTERING STEP 1", quote=FALSE)

          if(yesfactor){
               firstrim <- eStep1(yf=yesfactor, df1=fixdat.df, 
                             df1.ls=fixdatouter.list, 
                             df1fact.list <- fixdatfact.list,
                             inner.rank=ncoeffs.outer, initial.sample=initial.sample, formulaE=nlsform, nofactform=nofactform, 
                             start=start, nofactstart=nofactstart, inner.r.vector=inner.r, contR=nls.control, 
                             algo=algorithm, ycol=ycol, b.d=begin.diagnose)                                            #   eStep1
          }    # yesfactor
          else{
               firstrim <- eStep1(yf=yesfactor, df1=fixdat.df, 
                             df1.ls=fixdatouter.list, 
                             df1fact.list <- NULL,
                             inner.rank=ncoeffs.outer, initial.sample=initial.sample, formulaE=nlsform, nofactform=NULL, 
                             start=start, contR=nls.control, 
                             algo=algorithm, ycol=ycol, b.d=begin.diagnose)                                               # eStep1
          }     # no factors
          rim <- firstrim
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
#prn(rows.in.model)
#stop("_nls    end of step 1")
#############################################################################################################################
# Step 2
     print(" ",quote=FALSE)
     print("ENTERING STEP 2", quote=FALSE)

     eStep2out <- eStep2(yf=yesfactor, f2=nlsform, dfa2=fixdat.df, start=start, algo=algorithm, ms=mstart, ycol=ycol, 
             initn=length(SOON), inc=nls.control, finalm=rows.in.model, b.d=begin.diagnose)                             # eStep2

     rows.in.model <- eStep2out[[1]]
     fooResult <- eStep2out[[2]]
     resids2 <- eStep2out[[3]]
     sigma <- eStep2out[[4]]     

     rows.in.model[[dimx1]] <- 1:dimx1
     fooResult[[dimx1]] <- bignls

                                 if(begin.diagnose <= 80){print("", quote = FALSE);print(paste(spacer,"Section 80",sep=" "),quote=FALSE);
                                            Hmisc::prn(rows.in.model);Hmisc::prn(fooResult)      }

#stop("before extracting statistics")
     #
##############################################################################################################################
# EXTRACT STATS
     print("", quote=FALSE)
     print("EXTRACTING INTERMEDIATE STATISTICS", quote=FALSE)
     print("", quote=FALSE)
     betahatset <- matrix(0,nrow=dimx1, ncol=ncoeffs.outer)

     ###########################################
     # Extract statistics from each nls object #
     ###########################################
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
          param.est[i,] <- nlscoeffs
          t.set[i,] <- sumnls[[10]][,3]
     }        # i extracting statistics now completed
     #
     ###########################################
     # Cleanup after all extractions completed #
     ###########################################
     dimleverage <- dim(leverage)
     dimnames(leverage) <- list(rep("",dimleverage[1]),c("m", "Observation", "leverage"))
     param.est <- as.data.frame(param.est)
     names(param.est) <- names(start)
     m <- 1:dimx1
     param.est <- cbind(m,param.est)                                 #  here we add the m column
     # 
     t.set <- as.data.frame(t.set)
     names(t.set) <- names(start)
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
          "Number of model parameters"=        dim(t.set)[2] - 1,
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
