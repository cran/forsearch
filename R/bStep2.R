bStep2 <-
function (f2, dfa2, randm2, ms, finalm, fbg, b.d, rnk2, ycol) 
{
     nobs <- dim(dfa2)[1] 
     namesSubsets <- unique(dfa2$holdISG)
     nlevels <- length(namesSubsets)
     fooResult <- vector("list", length(finalm))

     names.dfa2 <- names(dfa2)
     ftlhs <- formula.tools::lhs(f2)
     ycol <- pmatch(ftlhs, names.dfa2)
     for(i in ms:(nobs-1)){
          rim <- finalm[[i-1]]
          thisdata <- dfa2[rim,]

          # Use thisdata to predict to all dfa2 and add squared error column to dfa2 as df2error

          thislme <- do.call(what=nlme::lme,args=list(fixed=f2, data=thisdata, random=randm2)    )         # lme
          fooResult[[i]] <- thislme
          thispredict <- stats::predict(thislme, dfa2, type="response", pred.var=1)                        # predict

          diffs2 <- (dfa2[,ycol] - thispredict)^2
          df2error <- cbind(dfa2,diffs2)         # df2error is in Observation order
          #################################################################################
          # Extract primary subset of rnk from each element of df2error using indices and #
          # order each one in turn. Pool the primary set of observations                  #
          #################################################################################
          primaryHold <- NULL
          for(j in 1:nlevels){
               index <- df2error$holdISG==namesSubsets[j]
               thisSub <- df2error[index,]
               thisSub <- thisSub[order(thisSub$diffs),]
               thisObs <- thisSub$Observation[1:rnk2]
               primaryHold <- c(primaryHold, thisObs)
          }   # j
          ###################################################################
          # Construct the reduced set of df2error and order this data frame #
          ###################################################################
          remainder.df <- df2error[-primaryHold,]
          remainder.df <- remainder.df[order(remainder.df$diffs),]
          #
          ############################################
          # Bring the number of observations up to i #
          ############################################
          needed <- i - length(primaryHold)
          needed <- remainder.df$Observation[1:needed]
          finalm[[i]] <- sort(c(primaryHold,needed))
     }   #  i

     return(list(finalm, fooResult)) 
}
