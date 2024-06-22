eStep2 <-
function (mstart, finalm, start, data2, pformula, algo2, cont, fixdb, fixlist, ycol, b.d)
{
     #                                            eStep2
     #
     # VALUE        A list of 4 elements. First element contains an updated list of the rim during Step 2 
     #                 Second element of the primary list is saved lm output for subsequent extraction of statistics. 
     #
     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        eStep2               "
     nrim <- length(finalm)

                                 if(b.d <= 61){print("", quote = FALSE);print(paste(spacer,"Section 61",sep=" "),quote=FALSE);
                                      Hmisc::prn(utils::head(data2));Hmisc::prn(ycol)       }

     fooResult <- vector("list", nrim)
     ##########################################
     # Create fooResult for complete database #
     ##########################################
     fooResult[[nrim]] <- stats::nls(formula=pformula, start=start, data=data2, control=cont, algorithm=algo2)       # nls
     finalm[[nrim]] <- 1:nrim

     Observation <- 1:nrim
     Diffs <- rep(-999,nrim)
     maxfirst <- fixdb$maxfirst
     holdISG <- fixdb$holdISG
     predictions.base <- data.frame(Observation, maxfirst, Diffs, holdISG)
     uholdISG <- unique(holdISG)  
     nsubsets <- length(uholdISG)
     xtemp.list <- vector("list", nrim)
     modCook <- rep(0,nrim)
     residuals2 <- matrix(0,nrim,nrim)
     #
#prn(finalm)
#prn(fooResult)
# stop("step2")
     ##############################
     # Begin loop creating step 2 #
     ##############################
     for(i in mstart:(nrim-1)){
          predictions <- predictions.base
          rim <- finalm[[i-1]]
          thisdf1 <- data2[rim,]

          nlsthis <- stats::nls(formula=pformula, start=start, data=thisdf1, control=cont, algorithm=algo2)              # nls

          fooResult[[i]] <- nlsthis
          preds <- stats::predict(nlsthis, data2)
          residuals2[,i] <- preds - data2[,ycol]
          predictions[,3] <- (preds - data2[,ycol])^2
          predictions <- predictions[order(predictions[,3]),]
          firstrim <- NULL
          for(j in 1:nsubsets){
               zzz <- predictions[predictions$holdISG==uholdISG[j],]
               ntake <- zzz$maxfirst[1]
               zzz <- zzz[1:ntake,]
               firstrim <- rbind(firstrim, zzz)
          }    #  j
                         if(b.d <= 64){print("", quote = FALSE);print(paste(spacer,"Section 64",sep=" "),quote=FALSE);
                                Hmisc::prn(i);Hmisc::prn(uholdISG);Hmisc::prn(firstrim)       }
          #
          ###################################################
          # Add observation numbers to bring number up to i #
          ###################################################
          predictions <- predictions[order(predictions$Observation),]
          predictions <- predictions[-c(firstrim$Observation),]
          predictions <- predictions[order(predictions$Diffs),]
          getnew <- i-dim(firstrim)[1]
          firstrim <- rbind(firstrim, predictions[1:getnew,])
          finalm[[i]] <- sort(firstrim$Observation)
     }     # i
                   if(b.d <= 70){print("", quote = FALSE);print(paste(spacer,"Section 70",sep=" "),quote=FALSE);
                                      Hmisc::prn(finalm);Hmisc::prn(predictions)       }

     sigma <- sqrt(sum(predictions[,3])/(nrim-length(start)))      # here, we're using an overall estimate, not by factor subsets
#     param.est <- as.data.frame(t(param.est))
     #
#     ###################################
#     # Finish off finalm and fooResult #
#     ###################################
#     finalm[[nrim]] <- 1:nrim
#     fooResult[[nrim]] <- stats::nls(formula=pformula, start=start, data=data2, control=cont, algorithm=algo2)              # nls

     outlist <- list(finalm, fooResult, residuals2, sigma)
     return(outlist)
}
