bStep2 <-
function (yf, f2, dfa2, randm2, onlyfactor=FALSE, ms, ycol, initn, inc, finalm, fbg, b.d) 
{
#                                            bStep2

     # yf                yesfactor
     # f2                2-sided formula for fixed effects
     # dfa2              data frame being analyzed
     # randm2            1-sided formula for random variables and groups
     # ms                mstart
     # ycol              column number of response variable
     # initn             inner.rnk vector of number of observations by source
     # inc               Logical, if TRuE causes relaxation of lmeControl
     # finalm            list of observations at each stage 
     # fbg               df2a as list by subfactors and groups
     # b.d               begin diagnosis          60 - 80
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      bStep2           "    
# print("in bStep2")

                      if(b.d <= 60 ){ print("",quote=FALSE);print(paste(spacehere,"Section 60",sep=" "),quote=FALSE);
                            Hmisc::prn(dfa2);Hmisc::prn(dim(dfa2));Hmisc::prn(randm2);Hmisc::prn(ms);Hmisc::prn(initn);
                            Hmisc::prn(inc);print("end bStep2 b.d = 60")   }

     #####################
     # Set up lmeControl #
     #####################
     if(inc){
          utils::str(lCtr <- nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, tolerance = 1e-2, 
                 niterEM = 1000, msMaxEval = 1000, msTol = 1e-2, optimMethod = "L-BFGS-B",
                 msVerbose = FALSE, returnObject = FALSE) )
          do.call(nlme::lmeControl, lCtr)
     }

     ##############################
     # Define the primary subsets #
     ##############################
     nobs <- dim(dfa2)[1]
     rim0 <- finalm[[ms-1]]         
     primary.df <- dfa2[rim0,]
     uprimaryISG <- unique(primary.df$comboISG)
     if(yf){
          nprimary <- length(rim0)
     }else{
          nprimary <- length(initn)
     }
     fooResult <- vector("list", length(finalm))

     for( ii in ms:(nobs - 1)){
          ##################################################################################
          # Calculate squared errors of prediction by the current rim on the complete dfa2 #
          ##################################################################################
          rim <- finalm[[ii-1]]
          this.data <- dfa2[rim,]
          if(inc){
               thislme <- do.call(what=nlme::lme,args=list(fixed=f2, data=this.data, random=randm2)    )
          }                                                                                                          # lme
          else{
               thislme <- do.call(what=nlme::lme,args=list(fixed=f2, data=this.data, random=randm2)    )   
          }
          thispredict <- stats::predict(thislme, dfa2)                                                              # predict
 
          fooResult[[ii]] <- thislme
          thesediffs2 <- (dfa2[,ycol] - thispredict)^2
          dfa2$diffs2 <- thesediffs2                        # df2error is in Observation order
          dfa2z <- dfa2[(order(dfa2$diffs2)),]

                      if(b.d <= 66 ){ print("",quote=FALSE);print(paste(spacehere,"Section 66",sep=" "),quote=FALSE);
                            Hmisc::prn(primary.df);Hmisc::prn(rim);Hmisc::prn(thislme);Hmisc::prn(thispredict)   }

          ##############################################################
          # Select the single observations from each primary subset by #
          # ordering dfa2 based on squared errors and Observation      #
          ##############################################################
          newrim <- NULL
          for(j in 1:nprimary){
               index <- uprimaryISG[j]
               uu <- dfa2z[dfa2z$comboISG==index,]
               uu <- uu[order(uu$diffs2,uu$Observation),]
               newrim<- sort(c(newrim, uu$Observation[1:initn[j]] ))
          }   # j
                      if(b.d <= 68 ){ print("",quote=FALSE);print(paste(spacehere,"Section 68",sep=" "),quote=FALSE);
                            Hmisc::prn(this.data);Hmisc::prn(nprimary);Hmisc::prn(thispredict);Hmisc::prn(newrim)   }

          match.found <- rim0 %in% newrim
          compare.primary <- data.frame(rim0, newrim, match.found)
          number.match <- sum(match.found)

                      if(b.d <= 69 ){ print("",quote=FALSE);print(paste(spacehere,"Section 69",sep=" "),quote=FALSE);
                            Hmisc::prn(compare.primary);Hmisc::prn(number.match);  stop("in bStep2 b.d 69")   }

          ###########################################################################################
          # Remove the selected observations from dfa2 and pick the next ii - nprimary observations #
          ###########################################################################################
          vv <- dfa2[order(dfa2$Observation),]
          vv <- dfa2[-newrim,]
          additional <- vv$Observation[1:(ii - nprimary)]
          newrim <- c(newrim, additional)
          finalm[[ii]] <- sort(newrim)
     }    # ii

     return(list(finalm, fooResult)) 
}
