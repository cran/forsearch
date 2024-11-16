bStep2 <-
function (yf, f2, dfa2, randm2, onlyfactor=FALSE, ms, ycol, initn, finalm, fbg, b.d) 
{
#                                                         bStep2

     # yf                yesfactor
     # f2                2-sided formula for fixed effects
     # dfa2              data frame being analyzed
     # randm2            1-sided formula for random variables and groups
     # ms                mstart
     # ycol              column number of response variable
     # initn             inner.rnk vector of number of observations by source
     # finalm            list of observations at each stage 
     # fbg               df2a as list by subfactors and groups
     # b.d               begin diagnosis          60 - 80
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      bStep1           "    
print("in bStep2")
     ###################################################################
     # Determine how many to pull (reinit) from each comboISG to match #
     #     initial rim from Step 1                                     #
     ###################################################################
     if(yf){
          rim <- finalm[[ms-1]]
          dfa2sub <- dfa2[rim,]
          udfa2sub <- unique(dfa2sub$comboISG)  
          ndfa2sub <- length(udfa2sub)            # number of unique comboISG in rim 
          subISG <- sort(udfa2sub)
          nISG <- length(rim)
          obscount <- rep(-99, ndfa2sub)
          source <- data.frame(subISG, obscount)

          for(i in 1:ndfa2sub){
               rock <- subISG[i]
               rimsub <- dfa2sub[dfa2sub$comboISG==rock,]
               xtotal <- dim(rimsub)[1]
               source[i,2] <- xtotal
          }    #  i

     }   # if yf

     nobs <- dim(dfa2)[1] 

     namesGroups <- unique(dfa2$groupISG)
     namesFixed  <- unique(dfa2$fixedISG)
     namesCombo <- sort(unique(dfa2$comboISG))
     rim0 <- finalm[[ms-1]]
     ncombo <- length(namesCombo)
     reinit <- rep(0,ncombo)
     dfa2sub <- dfa2[rim0,]
     for(i in 1:ncombo){
          bas <- dfa2sub[dfa2sub$comboISG==namesCombo[i],]
          reinit[i] <- dim(bas)[1]
     }     # i
     #

     fooResult <- vector("list", length(finalm))

     ftlhs <- formula.tools::lhs(f2)
     finalm[[ms-1]] <- sort(finalm[[ms-1]])

     for(i in ms:(nobs-1)){                         # loop with mstart
          rim <- finalm[[i-1]]
          thisdata <- dfa2[rim,]
          #
          ######################################################################################
          # Use thisdata to predict to all dfa2 and add squared error column to dfa2 as diffs2 #
          ######################################################################################

          thislme <- do.call(what=nlme::lme,args=list(fixed=f2, data=thisdata, random=randm2)    )         # lme
          fooResult[[i]] <- thislme
          thispredict <- stats::predict(thislme, dfa2)                                                 # predict
 
          thesediffs2 <- (dfa2[,ycol] - thispredict)^2
          dfa2$diffs2 <- thesediffs2 + dfa2$wiggle[i]                        # df2error is in Observation order
          dfa2z <- dfa2[(order(dfa2$diffs2)),]
          #
          ###################################
          # Collect primary observation set #
          ###################################
          newprimary <- NULL
          dfa2z <- dfa2z[order(dfa2z$Observation),]
          if(yf){
               for(ii in 1:ndfa2sub){
                    rock <- source[ii,1]
                    rockcount <- source[ii,2]
                    dfa2zcand <- dfa2z[dfa2z$comboISG==rock,]
                    xx <- dfa2zcand[,1]
                    newprimary <- c(newprimary,xx[1:rockcount] )
               }   # ii

          }   #       yf=TRUE
          else{
               for(j in 1:ncombo){
                    dfa2zcand <- dfa2z[dfa2z$comboISG==namesCombo[j],] 
                    xx <- dfa2zcand[,1]
                    newprimary <- c(newprimary,xx[      1:reinit[j]]   )
               }    #   j
          }   # no factors

                           if(b.d <= 69 ){ print("",quote=FALSE);print(paste(spacehere,"Section 69",sep=" "),quote=FALSE);
                                 Hmisc::prn(thispredict);Hmisc::prn(dfa2[rim0,]);Hmisc::prn(rim0);
                                 Hmisc::prn(dfa2z[1:(ms-1),]);Hmisc::prn(sort(dfa2z[1:(ms-1),1]));          stop("in bStep2 b.d 69")   }

          ############################################
          # Bring the number of observations up to i #
          ############################################
          needed <- i - length(newprimary)
          sourced <- dfa2$Observation[-newprimary] 
          needed <- sourced[1:needed]
          finalm[[i]] <- sort(c(newprimary,needed))
     }   #  i

     return(list(finalm, fooResult)) 
}
