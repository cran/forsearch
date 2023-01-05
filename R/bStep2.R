#' @export
bStep2 <-
function (fixed, nOuter, mnf, mstart, nobs, yobs, fbg, n.f, s.o, ras, b.d, verbose=FALSE) 
{
#                                                             bStep2
#
# VALUE		Sets of observation numbers after the first for use in forsearch_lme
#
# INPUT
#
#         begin.diagnose	Numeric indicator of first diagnostic to print for this function. 0 causes all diagnostics to print.
#         verbose             Logical. TRUE causes print of function identifiers at start and end of function
#
# NOTE: Calls aStep2
#
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running bStep2", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
#prn(yobs)
#stop("bStep2     yobs")
     nufixdatOuter <- nOuter
     maxnfixdat<- mnf
     fixdat.by.group <- fbg
     n.fixdat <- n.f
     saved.obsnums <- s.o
     rim.all.subgroups <- ras
     begin.diagnose <- b.d
     spacer <- rep(" ", 20)
     ###########################################################
     # Calculate Step 2 for each subgroup and do this by stage #
     ###########################################################
     DD <- rep(1:maxnfixdat, each=nufixdatOuter)
     KK <- rep(1:nufixdatOuter, times=maxnfixdat)
     largescore <- 10^10
     Score <- rep(largescore,nufixdatOuter * maxnfixdat)    # initialized as 10^10. Below, no value indicated by 10^8 
     RScore <-round(Score,4)
     HoldScore <- data.frame(DD, KK, Score,RScore)
#print("bStep 2      OK to 10")
     for(kk in 1:nufixdatOuter){                                    # run across the rows
          for(dd in (mstart+1):maxnfixdat){                                # run down the columns
               if(dd <= n.fixdat[kk]){
                    fixdatkk <- fixdat.by.group[[kk]]
                    ras <- rim.all.subgroups[[kk]][[dd-1]]
                    fixdatsub <- fixdatkk[ras,]
                    lm.inner2 <- stats::lm(fixed, data=fixdatsub, singular.ok=TRUE)                         #   lm
                    justarim <- aStep2(thislm=lm.inner2, data=fixdatkk, ycol=yobs, dd)
                    rim.all.subgroups[[kk]][[dd]] <- justarim[[1]]
                    index <- (HoldScore[1]==dd) & (HoldScore[2]==kk)
                    HoldScore[index,3] <- justarim[[2]]
               }      #   if dd <=#
          }           #   dd
     }                #   kk
#prn(HoldScore)
#stop("bStep2       HS")
     # Add dummy rows at the bottom of HoldScore #
     HSadd <- HoldScore[HoldScore[,1]==1,]
     HSadd[,1] <- maxnfixdat + 1
     HSadd[,c(3,4)] <- largescore
     HSadd[,4] <- round(HSadd[,3]) 
     HoldScore <- rbind(HoldScore, HSadd)
     HoldScore[,4] <- round(HoldScore[,3],4)
                                                             if(begin.diagnose <= 10){print(c(spacer,"Section 10"), quote=FALSE);Hmisc::prn(HoldScore)}
     #
     ################################################### 
     # Translate into original observation numbers and #
     # combine results over all subgroups              #
     ################################################### 
     rim.all.translated <- rim.all.subgroups
     for(dd in mstart:nobs){
          for(kk in 1:nufixdatOuter){
               index <- rim.all.subgroups[[kk]][[dd]]
               rim.all.translated[[kk]][[dd]] <- saved.obsnums[[kk]][index]
          }    #  kk
     }     #   dd
     #################################################
     # Combine subgroups according to overall scores #
     # Define BB as saturated list of lists          #
     #################################################
     vectint <- integer(nobs)
     innerlist <- vector("list", nufixdatOuter)
     for(kk in 1:nufixdatOuter){
          innerlist[[kk]] <- vectint
     }
     BB <- vector("list", nobs)                                                  #   define BB
     for(dd in 1:nobs){
          BB[[dd]] <- innerlist
     }

     for(kk in 1:nufixdatOuter){
          BB[[mstart]][[kk]] <- rim.all.translated[[kk]][[mstart]]    
     }    #   kk

     Subgroup.counter <- rep(0,nufixdatOuter)
                                                        if(begin.diagnose <= 15){print(c(spacer,"Section 15"),quote=FALSE);for(mm in 1:nufixdatOuter){Hmisc::prn(BB[[mstart]][[mm]])}}
     #
     ###########################################################################
     # For each line of BB, determine which group will add another observation #
     # dd is the index of the final listing BB of observations. dd indexes BB  #
     # Start by placing the last set into the next set of dd                   #
     # The result for BB is all lines from mstart to nobs are filled in but    #
     # still in nufixdatOuter lists.                                           #
     ###########################################################################
     for(dd in (mstart+1):nobs){
          SubgroupID <- as.integer(1:nufixdatOuter)
          Subscore <- rep(0,nufixdatOuter)
          candidates <- data.frame(SubgroupID, Subscore)
          for(kk in 1:nufixdatOuter){
#prn(c(dd,kk))
               BBdummy <- BB[[dd-1]][[kk]]
#prn(BBdummy)
               BB[[dd]][[kk]] <- BBdummy
               Subgroup.counter[kk] <- length(BB[[dd]][[kk]])                   # temporary count
               ###################################################
               # fill in the scores for candidates and sort them #
               ###################################################
               indHS <- (HoldScore[,1]==Subgroup.counter[kk] + 1) & (HoldScore[,2]==kk)
#prn(indHS)
              if(!any(indHS)){
                    print(paste("No viable candidates in (dd, kk) = ",   paste(dd,kk,sep=", "), sep=""), quote=FALSE)
               }
               else{
                    candidates[kk,2] <- HoldScore[indHS,3]
               }   # any indHS
          }     #   kk
 
          candidates <- candidates[order(candidates[,2]),]
#prn(candidates)
          selec <- candidates[1,1] 
          # How many observations in the subgroup counter? #
          picknext <- Subgroup.counter[selec] + 1                                 # equivalent to dd
          BB[[dd]][[selec]] <- rim.all.translated[[selec]][[picknext]]            # replace group selec in stage dd

          BBcount <- length(unlist(BB[[dd]]))
          if(BBcount==nobs)break                     #  jump out if all nobs observations have been assigned
     }       #   dd
     BB2 <- BB                                   # clean up list by reassignment
     BB <- NULL                                # save space
                                                        if(begin.diagnose <= 20){print(c(spacer,"Section 20"),quote=FALSE);Hmisc::prn(BB2[[nobs-1]]);Hmisc::prn(BB2[[nobs]])}
     #
     ############################################################################################
     # Combine subgroups of BB2 and slide the levels of BB2 down to account for multiple groups #
     ############################################################################################
     realstart <- mstart * nufixdatOuter
     CC <- vector("list", nobs)
     for(dd in mstart:nobs){
          CC[[(nufixdatOuter-1)*mstart + dd]] <- unlist(BB2[[dd]])
     }
     if(length(CC) > nobs){
     for(dd in length(CC):(nobs+1)) CC[[dd]] <- NULL
     }

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running bStep2", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }

     return(CC)
}
