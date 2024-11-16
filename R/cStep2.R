cStep2 <-
function (fe, finalm, rimbs, dfa2, onlyfactor=FALSE, ycol, cphties, mstart, rnk, b.d) 
{

     #                                            cStep2 
     #
     # VALUE        A list of 4 elements. First element contains an updated list of the rim during Step 2 
     #                 Second element of the primary list is saved lm output for subsequent extraction of statistics. 
     #
     # INPUT fe                Formula for analysis of entire dataset
     #       finalm            See VALUE above. finalm argument is the same but only for Step 1 values  NO
     #       rimbs             List, each element is a complete matrix of obs numbers and corresponding subset codes
     #       dfa2              Data frame being analyzed by forward search. 
     #       ycol              Response column number, including 1 for Observation
     #       cphties           ties
     #       mstart            First subset to be defined
     #       rnk               Rank of X matrix. For factors, this is rank with factors removed.
     #       b.d               Number at which to begin diagnostic listings
     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        cStep2               "
# print("in cStep2")
     yf <- TRUE        # Continue with yesfactor even though it is now always TRUE
     nobs <- dim(dfa2)[1]
     nfacts <- length(rimbs)
     namesrimbs <- names(rimbs)
     countSubs <- rep(-99, nfacts)

                                 if(b.d <= 60){print("", quote = FALSE);print(paste(spacer,"Section 60",sep=" "),quote=FALSE);
                                      Hmisc::prn(yf);Hmisc::prn(fe);Hmisc::prn(rimbs);Hmisc::prn(finalm);Hmisc::prn(mstart)       }

     fooResult <- vector("list", nobs)
     #
     #############################################
     # Save last fooResult by direct computation #
     #############################################
     thisdf1 <- dfa2
     td.et <- thisdf1$event.time
     td.st <- thisdf1$status

     xform <- paste("survival::Surv(time=td.et, event=td.st)", fe, sep=" ~ ")
     xform <- stats::as.formula(xform)

     thiscph <- NULL

     thiscph <- do.call(what=survival::coxph, args=list(formula=xform, data=thisdf1, 
                 ties=cphties, model=TRUE, singular.ok=TRUE, x=TRUE, y=TRUE))                    # coxph                                                                 
     fooResult[[nobs]] <- thiscph

     finalm[[nobs]] <- 1:nobs
     predictions.base <- data.frame(Observation <- 1:nobs, Diffs <- rep(-999, nobs))
     names(predictions.base) <- c("Observation", "Diffs2")
     xtemp.list <- vector("list", nobs)
     modCook <- rep(0,nobs)
     residuals2 <- matrix(0,nobs,nobs)
     #
     if(yf){
          # Get number of observations in primary to gather for each factor subset #
          temprim <- finalm[[mstart-1]]
          uu <- dfa2[temprim,]
          for(i in 1:nfacts){
               countSubs[i] <- sum(uu$holdISG==namesrimbs[i])
          }    #  i
          ##############################
          # Begin loop creating step 2 #
          ##############################
          for(i in mstart:(nobs-1)){
               dfa2aug <- cbind(-99,dfa2)
               rim <- finalm[[i-1]]
               thisdf1 <- dfa2[rim,]
               td.et <- thisdf1$event.time
               td.st <- thisdf1$status

               xform <- paste("survival::Surv(time=td.et, event=td.st)", fe, sep=" ~ ")
               xform <- stats::as.formula(xform)

               thiscph <- NULL
               thiscph <- do.call(what=survival::coxph, args=list(formula=xform, data=thisdf1, 
                       ties="efron", model=TRUE, singular.ok=TRUE, x=TRUE, y=TRUE))                                         # coxph
               fooResult[[i]] <- thiscph

                           if(b.d <=66 ){print("",quote=FALSE); print(paste(spacer,"Section 66",sep=" "),quote=FALSE);
                                Hmisc::prn(c(mstart,i));Hmisc::prn(thiscph$coefficients);Hmisc::prn(thiscph$x);Hmisc::prn(thiscph$y)    }

                thispredict <- stats::predict(thiscph, dfa2, type="lp")                                         #  predict

                residuals2[,i] <- thispredict - dfa2[,2]
                dfa2aug[,1] <- (thispredict - dfa2[,2])^2 + dfa2$wiggle
                dfa2aug <- dfa2aug[order(dfa2aug[,1]),]
               #
               ##########################################################################
               # Get obs numbers for initial set of countSubs obs in each factor subset #
               ##########################################################################
               collect.final <- NULL
               finalStage <- NULL
               for(j in 1:nfacts){
                    uu <- dfa2aug[dfa2aug$holdISG==namesrimbs[j],]
# prn(uu)
                    uu <- uu[1:countSubs[j],2]
#prn(uu) 
                    finalStage <- c(finalStage,uu)
#prn(finalStage)
               }    #   j
               collect.final <- c(collect.final, finalStage)
               #
 
#                           if(b.d <= 69 ){ print("",quote=FALSE);print(paste(spacer,"Section 69",sep=" "),quote=FALSE);
#                                Hmisc::prn(galaxy);Hmisc::prn(rnk2);print(" ");print("rim entering dStep2:");
#                                Hmisc::prn(thisrim); print(" ");print("Step2:");Hmisc::prn(cumrim);     stop("in dStep2 b.d 69")   }

                          if(b.d <= 69 ){ print("",quote=FALSE);print(paste(spacer,"Section 69",sep=" "),quote=FALSE);
                                 Hmisc::prn(thispredict);Hmisc::prn(dfa2[rim,]);Hmisc::prn(rim);
                                 print(" ");print("rim entering dStep2:");Hmisc::prn(dfa2[1:(mstart-1),]);Hmisc::prn(sort(dfa2[1:(mstart-1),1]));          stop("in cStep2 b.d = 69")   }

               ###################################################
               # Add observation numbers to bring number up to i #
               ###################################################
               dfa2aug <- dfa2aug[order(dfa2aug$Observation),]
               dfa2aug <- dfa2aug[-collect.final,]
               dfa2aug <- dfa2aug[order(dfa2aug[,1]),]
               nfinal <- length(collect.final)
               needed <- i - nfinal
               finalm[[i]] <-sort(c(collect.final, dfa2aug[1:needed,2]))

          }      #  i
     }         # factors present   
     else{
          stop("In cStep2 with yf=FALSE")
     }               # no factors present  

#stop("end of cStep2")
     outlist <- list(finalm, fooResult, residuals2)

     return(outlist)
}
