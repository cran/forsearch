aStep2 <-
function (yesfactor, form.A2, finalm, rimbs, dfa2, ycol, mstart, rnk, b.d) 
{
     #                                            aStep2
     #
     # VALUE        A list of 4 elements. First element contains an updated list of the rim during Step 2 
     #                 Second element of the primary list is saved lm output for subsequent extraction of statistics. 
     #                 If there are no factors, the first list contains only the first element above. 
     #
     # INPUT yesfactor         True or False for presence of factors
     #       form.A2           Formula for analysis of entire dataset
     #       finalm            See VALUE above. finalm argument is the same but only for Step 1 values
     #       rimbs             List, each element is a matrix of obs numbers and corresponding subset codes
     #       dfa2              Data frame being analyzed by forward search. Presence of Observation column has 
     #                             no effect on output
     #       ycol              Response column number, including 1 for Observation
     #       mstart            First subset to be defined
     #       rnk               Rank of X matrix. For factors, this is rank with factors removed.
     #       b.d               Number at which to begin diagnostic listings
     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        aStep2               "
     nobs <- dim(dfa2)[1]

     fooResult <- vector("list", nobs)
     fooResult[[nobs]] <- stats::lm(formula=form.A2, data=dfa2, singular.OK=TRUE, x=TRUE, y=TRUE)                                   # lm
     finalm[[nobs]] <- 1:nobs

     predictions.base <- data.frame(Observation <- 1:nobs, Diffs <- rep(-999, nobs))
     names(predictions.base) <- c("Observation", "Diffs2")
     xtemp.list <- vector("list", nobs)
     modCook <- rep(0,nobs)
     param.est <- matrix(0,nrow=rnk, ncol=nobs)
     residuals2 <- matrix(0,nobs,nobs)
     #
# prn(finalm[[mstart-1]])
     if(yesfactor){
          nlevels <- length(rimbs)
          predictions.base <- data.frame(Observation <- 1:nobs, Subset <- rep("A",nobs), Diffs <- rep(-999, nobs))
          names(predictions.base) <- c("Observation", "Subset", "Diffs2")
          matrimbs <- NULL
          for(ni in 1:length(rimbs)){
               matrimbs <- rbind(matrimbs, rimbs[[ni]])
          }
          matrimbs <- matrimbs[order(matrimbs[,1]),]
          predictions.base$Subset <- matrimbs[,2]

                                 if(b.d <= 62){print("", quote = FALSE);print(paste(spacer,"Section 62",sep=" "),quote=FALSE);
                                      Hmisc::prn(matrimbs);Hmisc::prn(predictions.base)       }

          ##############################
          # Begin loop creating step 2 #
          ##############################
          for(i in mstart:(nobs-1)){
               predictions <- predictions.base
               rim <- finalm[[i-1]]
               thisdf1 <- dfa2[rim,]

               lmthis <- stats::lm(formula=form.A2, data=thisdf1, singular.OK=TRUE, x=TRUE, y=TRUE)              # lm
               fooResult[[i]] <- lmthis

               preds <- stats::predict(lmthis, dfa2)
               residuals2[,i] <- preds - dfa2[,ycol]
               predictions[,3] <- (preds - dfa2[,ycol])^2
               predictions <- predictions[order(predictions[,3]),]
# prn(i)
# temppred <- predictions[1:i,]
# temppred <- temppred[order(temppred[,1]),]
# prn(temppred)
               ####################################################################
               # Get obs numbers for initial set of rnk obs in each factor subset #
               ####################################################################
               collect.final <- NULL
               for(j in 1:nlevels){
                    finalStage <- NULL
                    uu <- predictions[predictions[,2]==names(rimbs)[j],] 
                    firstrnk <- uu[,1]
                    finalStage <- rbind(finalStage,firstrnk[1:rnk])
                    collect.final <- c(collect.final, finalStage)
               }      # j
 
                                 if(b.d <= 66){print("", quote = FALSE);print(paste(spacer,"Section 66",sep=" "),quote=FALSE);
                                      Hmisc::prn(i);Hmisc::prn(rim);print("Observations in Initial group");Hmisc::prn(collect.final);
                                      Hmisc::prn(predictions[,1])       }
               #
               ###################################################
               # Add observation numbers to bring number up to i #
               ###################################################
               remainder <- match(collect.final,predictions[,1])
               remainder.mat <- predictions[-remainder,]

                                 if(b.d <= 68){print("", quote = FALSE);print(paste(spacer,"Section 68",sep=" "),quote=FALSE);
                                      Hmisc::prn(remainder.mat[order(remainder.mat$Subset,remainder.mat$Diffs2),])       }

               nfinal <- length(collect.final)
               needed <- remainder.mat[1:(i - nfinal),]
               needed.1 <- needed[,1]
               finalm[[i]] <- sort(c(collect.final, needed.1))
# prn(finalm[i])
# if(i > mstart + 2)stop("first 3 ")
          }     # i   

                                 if(b.d <= 70){print("", quote = FALSE);print(paste(spacer,"Section 70",sep=" "),quote=FALSE);
                                      Hmisc::prn(finalm)       }

          sigma <- sqrt(sum(predictions[,3])/(nobs-rnk))      # here, we're using an overall estimate, not by factor subsets

     }               # factors present
     else{
          for(i in mstart:(nobs-1)){
               predictions <- predictions.base
               rim <- finalm[[i-1]]
               thisdf1 <- dfa2[rim,]

               lmthis <- stats::lm(formula=form.A2, data=thisdf1, singular.OK=TRUE, x=TRUE, y=TRUE)                         # lm
               fooResult[[i]] <- lmthis

               preds <- stats::predict(lmthis, dfa2)
               residuals2[,i] <- preds - dfa2[,ycol]
               predictions[,2] <- (preds - dfa2[,ycol])^2                   # this used to be medaugx[,2]
               predictions <- predictions[order(predictions[,2]),]
# prn(i)
# temppred <- predictions[1:i,]
# prn(temppred)
               finalm[[i]] <- predictions[1:i,1]
# prn(finalm[[i]])
# if(i > mstart + 2)stop("first 3 no factors")

                                 if(b.d <= 71){print("", quote = FALSE);print(paste(spacer,"Section 71",sep=" "),quote=FALSE);
                                      Hmisc::prn(i);Hmisc::prn(thisdf1);Hmisc::prn(predictions);Hmisc::prn(finalm[[i]])       }
          }    #   i
          sigma <- sqrt(sum(predictions[,2])/(nobs-rnk))
     }               # no factors present  

     outlist <- list(finalm, fooResult, residuals2, sigma)
     return(outlist)
}
