#' @export
aStep1 <-
function (yesfactor, data, inner.rank, initial.sample, formula, ycol, nopl) 
{
     #                                    aStep1
     # 
     # VALUE      Produces rim for Step 1. If there are no factors, random choice 
     #            of inner.rank rows from entire dataset. If there are factors, use 
     #            variablelist and picksome to get a constrained choice of innersample.rows.
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            data           Data frame being analyzed by forward search. Presence of
     #                                 Observation column has no effect on output, but numbers 
     #                                 in Observation may not be 1 through n.
     #            inner.rank     Rank of X matrix of lm anlysis on entire dataset.
     #            initial.sample Number of random samples from which to take rim
     #            formula        Fixed effects formula of lm function
     #            ycol           Response column number
     #            nopl           n.obs.per.level
     #
     # NOTES
     #     Form inner groups available in data using variablelist.  rim will have nopl 
     #     observation from each of these inner groups plus enough randomly selected other 
     #     observations to make up a set of inner.rank observations, if any more are needed.
     #
     ##########################################################
     # Setup output and intermediate structures               #
     # pickm is a set of Observation numbers, not row numbers #
     ##########################################################
     dimdata.1 <- dim(data)[1]
     dimdata.2 <- dim(data)[2]                  # includes a column for Observation
     obsindex <- matrix(data[,1], nrow=dimdata.1, ncol=dimdata.1)             # column of original row number and randomized row number
     zlist <- vector("list", initial.sample)
     result <- matrix(0, nrow=dimdata.1, ncol=2)
     if(yesfactor){
          ss77 <- variablelist(datadf=data, verbose=FALSE)                   #   variablelist and picksome
          pickm <- picksome(subsetlist=ss77, nobs=dimdata.1, initial.sample, n.obs.per.level=nopl, rank=inner.rank)
          dimpickm <- dim(pickm)[2]
     }           # yesfactor TRUE
     #
     dataObs <- data[,1]
     medaugx <- matrix(1, nrow=initial.sample, ncol= 2)
     for(i in 1:initial.sample){
          zlist[[i]] <- result
          zlist[[i]][,1] <- sample(dataObs,dimdata.1)                                 # sample permutation of observation numbers
          if(yesfactor){
               rows.by.pickm <- zlist[[i]][1:dimpickm,1] <- pickm[i,]                 # sample formed by picksome
          }
     }     # i
     Potential <- matrix(0, nrow=initial.sample, ncol=inner.rank)
     if(yesfactor){
          obsindex <- zlist[[i]][1:dimpickm,1]
     }           # yesfactor TRUE
     else{
          obsindex <- zlist[[i]][1:inner.rank,1]
     }      
     for(i in 1:initial.sample){
          smalldata <- data[obsindex,]
          #################################################################
          # Run linear model on inner group and predict to entire dataset #
          #################################################################
          lmsmall <- stats::lm(formula, smalldata, singular.ok=TRUE)                                       #    lm
          predsmall <- stats::predict(lmsmall, data, type="response", pred.var=1)                          # predict
          errorsmall <- data[, ycol] - predsmall
          sserrorsmall <- errorsmall^2
          rowz <- zlist[[i]]
          augx <- matrix(rowz[,1], nrow=dimdata.1, ncol=1)
          augx <- cbind(augx, sserrorsmall)
          augx <- augx[order(augx[,2]),]
          medaugx[i,] <- augx[floor((dimdata.1 + inner.rank)/2),]
     }       #   i
     #
     randset <- 1:initial.sample
     medaugx <- cbind(medaugx, randset)
     minmed <- min(medaugx[,2])
     locatemin <- medaugx[medaugx[,2]==minmed,]
     if(is.matrix(locatemin)) locatemin <- locatemin[1,]
          zliststar <- zlist[[locatemin[3]]]
     #
     #######################################################################################
     # We have now created a matrix (Potential) of row numbers initial.sample x inner.rank #
     # Any row would meet the constraints regarding factors.                               #
     # Run lm() on each row and get predict() for each row against entire data base.       #
     # Calculate squared errors. Extract the median for each row of Potential. Add this    #
     # vector (cbind) to Potential. Identify the row of Potential whose first column       #
     # equals the minimm of the first column.  This is the rim we want.                    #
     #######################################################################################
     rim <- c(zliststar[1:inner.rank,1])
     return(rim)
}
