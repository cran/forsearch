bStep2 <-
function (fixed2, fulldata, random2, yf, nOuter, mstart, nobs, yobs, fbg, n.f, s.o, ras, b.d, ufixdat, LLL, verbose=FALSE) 
{
#                                                             bStep2
#
# VALUE		Sets of observation numbers after the first for use in forsearch_lme. Also an archive of lme runs on the current subset
#
#     yf <- yesfactor 

#     nOuter   <- nufixdatOuter       or nufixdatOuterInner
#     fbg      <- fixdat.by.Outer     or fixdat.by.OuterInner
#     s.o      <- saved.obsnums.Outer or saved.obsnums.OuterInner
#     fixdat   <- ufixdatOuter        or ufixdatOuterInner
#     

#     n.fixdat <-          n.f    ?????
#     ras <- rim.all.subgroups 
#     b.d <- begin.diagnose 
#     yobs <-  response.colnum 

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
     begin.diagnose <- b.d
     spacer <- "                                    bStep2   "
     #################################################################################
     # Result of this function is list fulldata with all levels filled in for Step 2 #
     #################################################################################
     tempSS <- 0
     fulldata <- data.frame(fulldata,tempSS)
                                                               if(begin.diagnose <= 35){print(paste(spacer,"Section 35"), quote=FALSE);Hmisc::prn(s.o);;Hmisc::prn(mstart);Hmisc::prn(nobs-1)}

     for(pp in mstart:(nobs-1)){                                     # runs over levels of fulldata
  
                                                             if(begin.diagnose <= 36){print(paste(spacer, "Section 36     pp = ", pp), quote=FALSE)     }

          ###########################################################################################
          # Calculate yhat using lme and predict to entire data base starting with pp observations. # 
          ###########################################################################################
          rim.all.subgroups <- ras[[pp-1]]
          thisdata <- fulldata[rim.all.subgroups,]
          thislme <- do.call(nlme::lme,list(fixed=fixed2, data=thisdata, random=random2))                       # lme
          LLL[[pp]] <- thislme              # store the lme object for later extractions
          fulldata$tempSS <- 0
          yhat <- stats::predict(thislme, newdata=fulldata)  
          diff <- fulldata[,yobs] - yhat
          fulldata$tempSS <- diff^2
          #
                                                             if(begin.diagnose <= 38){print(paste(spacer, "Section 38     pp = ", pp), quote=FALSE)     }

          ############################################################################################## 
          # Group observations by OuterInner and select the one from each group with the lowest tempSS #
          # Add each one to build.rim.initial                                                          # 
          # Below, bigData will always be the same length because it is defined by removal of the top  #
          #      of the grouped data.                                                                  #
          ############################################################################################## 
          fixdat.by.group2 <- vector("list", nOuter)                                        #  list
          n.fixdat2 <- rep(0, nOuter)
          build.rim.input <- NULL
          if(yf){
               for(ww in 1:nOuter){
                    fixdat.by.group2[[ww]] <- fulldata[fulldata$OuterInner==ufixdat[ww],]
                    n.fixdat2[ww] <- dim(fixdat.by.group2[[ww]])[1]
               }    
                                                             if(begin.diagnose <= 39.1){print(paste(spacer, "Section 39.1     pp = ", pp), quote=FALSE)     }

          }    #  there are factors
         else{
               for(ww in 1:nOuter){
                    fixdat.by.group2[[ww]] <- fulldata[fulldata$OuterSubgroup==ufixdat[ww],]
                    n.fixdat2[ww] <- dim(fixdat.by.group2[[ww]])[1]
               }        #  ww
                                                             if(begin.diagnose <= 39.2){print(paste(spacer, "Section 39.2     pp = ", pp), quote=FALSE)     }

          }    # there are no factors 
          bigData <- NULL
          for(ww in 1:nOuter){
               theseObs <- fixdat.by.group2[[ww]]
               theseObs <- theseObs[order(theseObs$tempSS),]                     # put them in order
               fixdat.by.group2[[ww]] <- theseObs
               build.rim.input <- c(build.rim.input, theseObs$Observation[1])        #   concatenates over the first observation of each OuterInner groups
               fixdat.by.group2[[ww]] <- fixdat.by.group2[[ww]][-1,]              # remove the one added 
               bigData <- rbind(bigData, fixdat.by.group2[[ww]])
          }  # ww
          #
                                                             if(begin.diagnose <= 40){print(paste(spacer, "Section 40     pp = ", pp), quote=FALSE)     }

          ################################################################## 
          # Pool the remainder of the observations and sort them by tempSS #
          # Then add the next observation                                  #
          ################################################################## 
          bigData <- bigData[order(bigData$tempSS),]
          indexNext <- pp - nOuter
          nextObs <-  bigData$Observation[1:indexNext]
          tempras <- c(build.rim.input, nextObs) 
          tempras <- sort(tempras)
          ras[[pp]] <- tempras
                                                              if(begin.diagnose <= 41){print(paste(spacer,"Section 41"), quote=FALSE);Hmisc::prn(ras[[pp]])}

     }    #   pp
     rasLLL <- list(ras, LLL)
                                                               if(begin.diagnose <= 42){print(paste(spacer,"Section 42"), quote=FALSE);Hmisc::prn(rasLLL)}

     if(verbose) {
          print("", quote = FALSE)
          print("Finished running bStep2", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(rasLLL)
}
