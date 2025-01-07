aStep1 <-
function (yesfactor, df1, df1.ls, inner.rank, initial.sample, formulaA, nofactform, ycol, b.d) 
{
     #                                    aStep1   
     # 
     # VALUE      Produces rim for Step 1. Uses candprep function to form matrix of candidate sets of observation numbers
     #            and runs lm and predictor to determine set with median sum of squared errors.
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            df1            Data frame being analyzed by forward search. Observation column has no effect on lm analysis
     #            df1.ls         List of df1 by factor subset or NULL
     #            inner.rank     Vector. Number of observations to pull from each factor subset, or from overall database otherwise
     #            initial.sample Number of random samples from which to take rim
     #            formulaA       Formula for all effects including factors and constructed variables    
     #            nofactform     Formula for all effects, omitting factors
     #            ycol           Response column number
     #            b.d            begin.diagnose Ranges from 25 to 45
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      aStep1           "    
# print("in aStep1") 

                           if(b.d <=31 ){ print("",quote=FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(inner.rank);
                                 Hmisc::prn(initial.sample);Hmisc::prn(dim(df1));Hmisc::prn(formulaA);Hmisc::prn(ycol)   }

     ############################################################################################################### 
     # Collect sufficient observation numbers from all factor subsets or the entire database if no factors present #
     # Assemble them into a matrix.  The number of rows is initial.sample and the number of columns depends on the # 
     # number of observations.                                                                            #
     # Then run lm() on each row and predict to the entire database. Determine the median of the squared errors of # 
     # each prediction. Select the row that produced the minimum of these mediums.                                 #
     # The formula for the prediction is complete, including constructed covariates and factors.                   #
     ############################################################################################################### 
     # Set up final holding variables
     nobs <- dim(df1)[1]
     Observation <- 1:initial.sample
     SSE <- rep(-99, initial.sample)
     pullN <- max(inner.rank)

                                if(b.d <= 31){print("", quote = FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                                      Hmisc::prn(inner.rank);Hmisc::prn(pullN)       }

     if(yesfactor){
          nfactsub <- length(inner.rank)
          hold.cands <- candprep(yf=yesfactor, dfa2=df1, fixd.ls=df1.ls, inner.rank=inner.rank, preprnk=pullN, 
                        in.sam=initial.sample, b.d=b.d)                                                       # candprep list
     }        #   yesfactor TRUE
     else{
          hold.cands <- candprep(yf=yesfactor, dfa2=df1, fixd.ls=NULL, inner.rank=inner.rank, preprnk=pullN, 
                       in.sam=initial.sample, b.d=b.d)                                                        # candprep matrix
          hold.cands <- hold.cands[[1]]
     }    # no factors present
     #
     ################################################################################################### 
     # Run the lm function on each row of hold.cands and predict to the entire database for each one   #
     # Identify the median of the squared errors of each predictio.                                    #
     # Calculate the sum of squares of the error for each one amd store it in SSE.                     #
     # Sort this matrix on the 2nd column and locate the median row.  Save this value. Sort the matrix #
     # on this value. The observations on the first row (lowest median SSE are the rim for Step 1      #                                                                  #
     ################################################################################################### 
     if(yesfactor){
          SSErim <- NULL
          rimout.ls <- vector("list", nfactsub)
          this.form <- nofactform
          for(j in 1:nfactsub){
               hold.cands.sub <- hold.cands[[j]]
               hold.cands.sub <- hold.cands.sub[,1:inner.rank[j]]
               SSE <- rep(-99, initial.sample)
               for(r in 1:initial.sample){
                    if(inner.rank[j]==1){
                         index <- hold.cands.sub[r]
                    }
                    else{
                         index <- hold.cands.sub[r,]
                    }
                    xmat <- match(index, df1[,1])
                    smalldata <- df1[xmat,]
                    indexuniv <- smalldata$holdISG[1]
                    universe <- df1[df1$holdISG==indexuniv,]
                    MED <- floor(dim(universe)[1]/2 + .99991)
                    MED[MED==0] <- 1

                    lmsmall <- stats::lm(formula=this.form, data=smalldata, singular.ok=TRUE)                             #    lm
                    predsmall <- stats::predict(lmsmall, data=universe)                                                    # predict

                    errorsmall <- df1[, ycol] - predsmall
                    sserrorsmall <- sort(errorsmall^2)  
                    SSE[r] <- sserrorsmall[MED]
               }     # r

                                if(b.d <= 33){print("", quote = FALSE);print(paste(spacehere,"Section 33",sep=" "),quote=FALSE);
                                      Hmisc::prn(predsmall);Hmisc::prn(SSE)       }

               hold.cands.sub <- cbind(hold.cands.sub, SSE)
               lastcol <- dim(hold.cands.sub)[2]
               hold.cands.sub <- hold.cands.sub[order(hold.cands.sub[,lastcol]),]
               uu <- hold.cands.sub[,-lastcol]

               if(inner.rank[j]==1) uu <- matrix(uu,ncol=1)

               SSErim <- c(SSErim,uu[1,])
               rimout.ls[[j]] <- sort(uu[1,])
          }    #  j 
          rimout <- sort(SSErim)
          nsubs <- length(df1.ls)
          temprims <- rimout
     }      #   yesfactor
     else{
          # hold.cands must be a matrix with last column of 0
          MED <- floor(nobs/2 + 0.00001)             # no factors, so must have at least inner.rank = 2
          for(i in 1:initial.sample){
               index <- hold.cands[i,]         # ith row is set of observation numbers, not row numbers
               dropone <- length(index) 
               index <- index[-dropone]
               smalldata <- df1[index, ]
               this.form <- formulaA

               lmsmall <- stats::lm(formula=this.form, data=smalldata, singular.ok=TRUE)                             #    lm
               predsmall <- stats::predict(lmsmall, data=df1)                                                    # predict 

               errorsmall <- df1[, ycol] - predsmall
               sserrorsmall <- sort(errorsmall^2)  
               SSE[i] <- sserrorsmall[MED]
          }     # i
          hold.cands <- cbind(hold.cands, SSE)
          lastcol <- dim(hold.cands)[2]
          hold.cands <- hold.cands[order(hold.cands[,lastcol]),]
          rimout2 <- hold.cands[1,]
          rimout <- rimout2[-lastcol]
          lastcol <- lastcol-1
          rimout <- rimout[-lastcol]
          nsubs <- length(df1.ls)
          temprims <- rimout
          rimout.ls <- NULL
     }     #   no factors present

                               if(b.d <= 35){print("", quote = FALSE);print(paste(spacehere,"Section 35",sep=" "),quote=FALSE);
                                      Hmisc::prn(rimout);Hmisc::prn(rimout.ls)       }

    allrims <- list(rimout, rimout.ls)
# print("leaving aStep1")

    return(allrims)
}
