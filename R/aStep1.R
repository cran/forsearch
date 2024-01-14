aStep1 <-
function (yesfactor, df1, inner.rank, initial.sample, formula, ycol, nopl, b.d) 
{
     #                                    aStep1   
     # 
     # VALUE      Produces rim for Step 1. If there are no factors, random choice 
     #            of inner.rank rows from entire dataset. If there are factors, use 
     #            variablelist and picksome to get a constrained choice of innersample.rows.
# following is no longer true 
    #            output is a list of 2 levels, the first is all obs in Step 1, the second is the
     #            same but is a list of length = number of subsets each level of with is obs by level.
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            df1            Data frame being analyzed by forward search. Presence of
     #                                 Observation column has no effect on output, but numbers 
     #                                 in Observation may not be 1 through n. Data may be full dataset
     #                                 (no factors) or list of subsets (factors present)
     #            inner.rank     Rank of lm analysis on dataset with or without factor variables, depending on yesfactor     
     #            initial.sample Number of random samples from which to take rim
     #            formula        Formula for fixed effects in lm always without factors
     #            ycol           Response column number
     #            nopl           n.obs.per.level multiplier
     #            b.d            begin.diagnose Ranges from 25 to 45
     #
     # NOTES
     #     This function does not involve any inner or outer factors such as grouping variables.
     #     Form inner groups available in data using variablelist.  rim will have nopl 
     #     observation from each of these inner groups plus enough randomly selected other 
     #     observations to make up a set of inner.rank observations, if any more are needed.
     #     If there are no inner groups, need to do the augmentation on the full dataset after
     #     taking inner.rank * nobs.per.level observations.
     #
     
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      aStep1           "    
     
     #######################################
     # The following is a support function #
     #######################################
     supfun <- function(df2, initial.sample, inner.rank, ycol, nopl, formula ){

                          if(b.d <=32 ){ print("",quote=FALSE);print(paste(spacehere,"Section 32",sep=" "),quote=FALSE);
                               print("supfun");Hmisc::prn(df2);Hmisc::prn(initial.sample);Hmisc::prn(inner.rank);
                               Hmisc::prn(nopl);Hmisc::prn(formula)   }
          data <- df2                  # load factor subset
          dimdata.1 <- dim(data)[1]
          dimdata.2 <- dim(data)[2]                  # includes a column for Observation
          obsindex <- matrix(data[,1], nrow=dimdata.1, ncol=dimdata.2)             # column of original row number and randomized row number
          zlist <- vector("list", initial.sample)
          result <- matrix(0, nrow=dimdata.1, ncol=2)
          for(i in 1:initial.sample){
               zlist[[i]] <- result                                                        # this is the initialized 0,0 matrix
               bigvector <- 1:dimdata.1
               zlist[[i]][,1] <- sample(bigvector, dimdata.1, replace=FALSE)                  # sample permutation of observation numbers
           }  
           dataObs <- data[,1]
           medaugx <- matrix(1, nrow=initial.sample, ncol= 2)
           for(i in 1:initial.sample){
                zlist[[i]] <- result                                                        # this is the initialized 0,0 matrix
                zlist[[i]][,1] <- sample(dataObs,dimdata.1, replace=FALSE)                  # sample permutation of observation numbers
          }     # i
                           if(b.d <=33 ){ print("",quote=FALSE);print(paste(spacehere,"Section 33",sep=" "),quote=FALSE);
                                Hmisc::prn(zlist)   }

          for(i in 1:initial.sample){
               Potential <- matrix(0, nrow=initial.sample, ncol=inner.rank)
               dimpickm2 <- max(inner.rank, nopl)
               obsindex <- zlist[[i]][1:dimpickm2,1]
 
                           if(b.d <=35 ){ print("",quote=FALSE);print(paste(spacehere,"Section 35", i, sep=" "),quote=FALSE);
                                         Hmisc::prn(yesfactor);Hmisc::prn(inner.rank);Hmisc::prn(obsindex)   }

               Obs <- match(obsindex,data[,1])
               smalldata <- data[Obs,]
               #################################################################
               # Run linear model on inner group and predict to entire dataset #
               #################################################################
               this.form <- formula
               lmsmall <- stats::lm(formula=this.form, smalldata, singular.ok=TRUE)                             #    lm
               predsmall <- stats::predict(lmsmall, data, type="response", pred.var=1)                          # predict
               errorsmall <- data[, ycol] - predsmall
               sserrorsmall <- errorsmall^2
               rowz <- zlist[[i]]
               augx <- matrix(rowz[,1], nrow=dimdata.1, ncol=1)
               augx <- cbind(augx, sserrorsmall)
               augx <- augx[order(augx[,2]),]                            # order the squared errors
               medaugx[i,] <- augx[floor((dimdata.1 + inner.rank)/2),]   # place the median in a matrix

                           if(b.d <=37 ){print("",quote=FALSE); print(paste(spacehere,"Section 37 ", i, sep=" "),quote=FALSE);
                               Hmisc::prn(smalldata);Hmisc::prn(this.form);
                               Hmisc::prn(utils::head(augx));Hmisc::prn(utils::tail(augx));Hmisc::prn(medaugx[i,])    }
          }       #   i

          randset <- 1:initial.sample
          medaugx <- cbind(medaugx, randset)              # identify which permutation each one comes from
          minmed <- min(medaugx[,2])                      # identify the minimum of these
          locatemin <- medaugx[medaugx[,2] == minmed, ]   # this is the minimum of the medians; it could be more than one subset

                            if(b.d <=38 ){print("",quote=FALSE); print(paste(spacehere,"Section 38",sep=" "),quote=FALSE);
                                Hmisc::prn(utils::head(medaugx));
                                Hmisc::prn(utils::tail(medaugx));Hmisc::prn(minmed);Hmisc::prn(locatemin)    }

          if(is.matrix(locatemin))  locatemin <- locatemin[1,]     # resolve >1 equal minima
          lcm3 <- locatemin[3]       # pull out the randset to be used 
          zliststar <- zlist[[lcm3]]     # this is the permutation to use
         #
          rim <- c(zliststar[1:length(obsindex), 1])
         return(rim)
     }       # end supfun



#############################################    Main function starts here ##########################################
 
                           if(b.d <=31 ){ print("",quote=FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(inner.rank);
                                 Hmisc::prn(formula);Hmisc::prn(ycol)   }

     if(yesfactor){
          nlevels <- length(df1)
          rim.by.level <- vector("list", nlevels)
          rim <- NULL
          for(ii in 1:nlevels){
               holdrim <- supfun(df2=df1[[ii]], initial.sample=initial.sample, inner.rank=inner.rank,ycol,nopl=nopl,formula=formula)     # supfun
               rim.by.level[[ii]] <- holdrim
               rim <- c(rim, holdrim)
          }    #  ii
          rimout <- list(rim, rim.by.level)        # single vector of observation numbers and list of obs numbers by level
     }       # major function section: factor present
     else{
          rimout <- supfun(df2=df1, initial.sample=initial.sample, inner.rank=inner.rank, ycol, nopl=nopl, formula=formula)           #  supfun
     }       # major function section: no factor present

     return(rimout)
}
