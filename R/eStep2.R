eStep2 <-
function (yf, f2, dfa2, start, algo, ms, ycol, initn, inc, finalm, b.d) 
{
     #                                            eStep2
     #
     # VALUE        A list of 4 elements. First element contains an updated list of the rim during Step 2 
     #                 Second element of the primary list is saved nls output for subsequent extraction of statistics. 
     #
     # INPUT 
     #        yf         yesfactor
     #        f2         formula
     #        dfa2       Complete database of the analysis
     #        start
     #        algo       algorithm
     #        ms         mstart
     #        ycol       Response column number
     #        initn      Number of observations in rim0 if no factors
     #        inc        nls.Control
     #        finalm
     #        b.d 



     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX        eStep2               "
     nrim <- length(finalm)

                                 if(b.d <= 61){print("", quote = FALSE);print(paste(spacer,"Section 61",sep=" "),quote=FALSE);
                                      Hmisc::prn(yf);Hmisc::prn(f2);Hmisc::prn(dfa2);Hmisc::prn(start);Hmisc::prn(ycol);
                                      Hmisc::prn(initn);Hmisc::prn(inc);print("end of eStep2 b.d = 61")       }

     fooResult <- vector("list", nrim)

#stop("initial eStep2")

     ##############################
     # Define the primary subsets #
     ##############################
     nobs <- dim(dfa2)[1]
     rim0 <- finalm[[ms-1]]         
     primary.df <- dfa2[rim0,]
     primaryISG <- primary.df$comboISG
     uprimaryISG <- unique(primary.df$comboISG)
     nuprim <- length(uprimaryISG)    
     uprim <- factor(primaryISG)
     primtab <- tabulate(uprim, nuprim)
     primtab.df <- data.frame(uprimaryISG, primtab)
     fooResult <- vector("list", length(finalm))
     residuals2 <- matrix(0,nrim,nrim)

    for( ii in ms:(nobs - 1)){
          ##################################################################################
          # Calculate squared errors of prediction by the current rim on the complete dfa2 #
          ##################################################################################
          rim <- finalm[[ii-1]]
          this.data <- dfa2[rim,]

          if(inc){
#               print("loosening controls")
               thisnls <- do.call(what=stats::nls,args=list
                         (formula=f2, data=this.data, start=start, algorithm=algo,
                         control=list(maxiter=1000, tol=.001, warnOnly=TRUE))
)                                          # nls
          }                                                                                                          
          else{
               thisnls <- do.call(what=stats::nls,args=list(formula=f2, data=this.data, start=start, algorithm=algo))   # nls
          }
          thispredict <- stats::predict(thisnls, dfa2)                                                              # predict
 
          fooResult[[ii]] <- thisnls
          newpreds <- thispredict - dfa2[,ycol]
          thesediffs2 <- newpreds^2

          dfa2$diffs2 <- thesediffs2                        # df2error is in Observation order
          dfa2z <- dfa2[(order(dfa2$diffs2, dfa2$Observation)),]
          residuals2[,ii] <- thispredict - dfa2z[,ycol]

                      if(b.d <= 66 ){ print("",quote=FALSE);print(paste(spacer,"Section 66",sep=" "),quote=FALSE);
                            Hmisc::prn(primary.df);Hmisc::prn(rim);Hmisc::prn(thisnls);Hmisc::prn(thispredict)   }

          ##############################################################
          # Select the observations from each primary subset by #
          # ordering dfa2 based on squared errors and Observation      #
          ##############################################################
          newrim <- NULL
          for(j in 1:nuprim){
               index <- primtab.df[j,1]
               icount <- primtab.df[j,2]
               uu <- dfa2z[dfa2z$comboISG==index,]
               uu <- uu[order(uu$diffs2,uu$Observation),]
               newrim<- c(newrim, uu$Observation[1:icount] )
          }   # j
          nprimary <- length(newrim)
                      if(b.d <= 68 ){ print("",quote=FALSE);print(paste(spacer,"Section 68",sep=" "),quote=FALSE);
                            Hmisc::prn(this.data);Hmisc::prn(nprimary);Hmisc::prn(thispredict);Hmisc::prn(newrim)   }
          newrim <- sort(newrim)
          rim0.match.found <- rim0 %in% newrim
          compare.primary <- data.frame(rim0, newrim, rim0.match.found)
          number.match <- sum(rim0.match.found)

                      if(b.d <= 69 ){ print("",quote=FALSE);print(paste(spacer,"Section 69",sep=" "),quote=FALSE);
                            Hmisc::prn(compare.primary);Hmisc::prn(number.match);Hmisc::prn(dfa2[rim0,]);Hmisc::prn(dfa2[newrim,]);  
                            stop("in bStep2 b.d 69")   }

          ###########################################################################################
          # Remove the selected observations from dfa2 and pick the next ii - nprimary observations #
          ###########################################################################################
          vv <- dfa2[order(dfa2$Observation),]
          vv <- dfa2[-newrim,]
          additional <- vv$Observation[1:(ii - nprimary)]
          newrim <- c(newrim, additional)
          finalm[[ii]] <- sort(newrim)
     }    # ii
     sigma <- sqrt(sum(thesediffs2)/(nrim-length(start)))      # here, we're using the last set of defferences
     outlist <- list(finalm, fooResult, residuals2, sigma)

     return(outlist)
}
