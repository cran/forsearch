cStep1 <-
function (yesfactor, df1, df1.ls, inner.rank, initial.sample, formula, f.e, ycol, nopl, b.d) 
{
     #                                    cStep1   
     # REVISE ALL THIS
     # VALUE      Produces rim for Step 1. If there are factors, selects random sets from all factor subsets and runs lm
     #            and predictor to determine set with median sum of squared errors.
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            df1            Data frame being analyzed by forward search. Observation column has no effect on lm analysis
     #            df1.ls         List of df1 by factor subset or NULL
     #            inner.rank     Rank of lm analysis on dataset with or without factor variables, depending on yesfactor     
     #            initial.sample Number of random samples from which to take rim
     #            formula        Formula for all effects including factors and constructed variables    
     #            f.e            Right hand side of formula for Surv function
     #            ycol           Response column number
     #            nopl           n.obs.per.level 
     #            b.d            begin.diagnose Ranges from 25 to 45
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      cStep1           "    
  
                           if(b.d <=31 ){ print("",quote=FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(inner.rank);
                                 Hmisc::prn(initial.sample);Hmisc::prn(dim(df1));Hmisc::prn(formula);Hmisc::prn(ycol)   }


     ############################################################################################################### 
     # Collect sufficient observation numbers from all factor subsets or the entire database if no factors present #
     # Assemble them into a matrix.  The number of rows is initial.sample and the number of columns depends on the # 
     # number of observations and nopl.                                                                            #
     # Then run lm() on each row and predict to the entire database.  Select the Step 1 row which is the median of #
     # the sum of squared errors of the predictions. Formula contains all variables, including factors.            #
     ############################################################################################################### 

     # Set up final holding variables
     nobs <- dim(df1)[1]
     Observation <- 1:initial.sample
     Sampleno <- 1:initial.sample
     SSE <- -99

     sumSqError <- data.frame(Sampleno, SSE)     # from which to select median set based on SSE
     pullN <- max(inner.rank,nopl)

     if(yesfactor) nlevels <- length(df1.ls)
     else{
          nlevels <- 1
          df1.ls <- list(df1)
     }
     hold.cands <- matrix(-9L, nrow=initial.sample, ncol=pullN*nlevels)      # will be appended when obs identified by stratum 
     for(i in 1:initial.sample){
          hold.cands.this.subset <- NULL
          for(j in 1:nlevels){
               this.subset <- df1.ls[[j]]
               n.in.subset <- dim(this.subset)[1]    # may be more obs in some subsets than in others
               popu <- 1:n.in.subset
               hold.cands.temp <- sample(x=popu, size=pullN, replace=FALSE)                      # row numbers
               hold.cands.temp <- this.subset$Observation[hold.cands.temp]          # associated observation numbers
               hold.cands.this.subset <- c(hold.cands.this.subset, hold.cands.temp)      # append observation numbers to vector
          }    #   j
          len.this.subset <- length(hold.cands.this.subset)
          hold.cands[i,1:len.this.subset] <- sort(hold.cands.this.subset) 
     }     # i

                           if(b.d <=32 ){ print("",quote=FALSE);print(paste(spacehere,"Section 32",sep=" "),quote=FALSE);
                                 Hmisc::prn(inner.rank);Hmisc::prn(nopl);Hmisc::prn(nobs);Hmisc::prn(utils::head(sumSqError));Hmisc::prn(nlevels);
                                 Hmisc::prn(utils::head(hold.cands));Hmisc::prn(utils::tail(hold.cands))   }


     ################################################################################################# 
     # Run the coxph function on each row of hold.cands and predict to the entire database for each one #
     # Calculate the sum of squares of the error for each one amd store it in sumSqError[i,2]        #
     # Sort this matrix on the 2nd column and locate the median row.  The observations for this row  #
     # are the rim for Step 1                                                                        #
     ################################################################################################# 

     for(i in 1:initial.sample){
          index <- hold.cands[i,]
          smalldata <- df1[index,]

          thisform <- paste("survival::Surv(time=event.time, event=status)", f.e, sep=" ~ ")
          thisform <- stats::as.formula(thisform)

          lmsmall <- do.call(what=survival::coxph, args=list(formula=thisform, data=smalldata, 
               ties="efron", model=TRUE, singular.ok=FALSE, x=TRUE, y=TRUE))                                      # coxph

          predsmall <- stats::predict(object=lmsmall)                                                             #  predict

          errorsmall <- df1[, ycol] - predsmall
          predsmall2 <- predsmall^2
          sumSqError[i,2] <- sum(predsmall2)
     }     # i
     sumSqError <- sumSqError[order(sumSqError[,2]),]
     hold.cands <- hold.cands[sumSqError[,1],]
     MED <- round(initial.sample/2 + .000001)

                           if(b.d <=34 ){ print("",quote=FALSE);print(paste(spacehere,"Section 34",sep=" "),quote=FALSE);
                                 Hmisc::prn(smalldata);Hmisc::prn(thisform);Hmisc::prn(predsmall);Hmisc::prn(utils::head(sumSqError));
                                 Hmisc::prn(utils::tail(sumSqError))
                                 Hmisc::prn(MED);Hmisc::prn(utils::head(hold.cands));Hmisc::prn(utils::tail(hold.cands))    }

     rimout <- hold.cands[MED,]
     return(rimout)
}
