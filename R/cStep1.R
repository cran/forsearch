cStep1 <-
function (df1, df1.ls, inner.rank, initial.sample, f.e, cphties, ycol, b.d) 
{
     #                                    cStep1   
     # 
     # VALUE      Produces rim for Step 1. Uses candprep function to form matrix of candidate sets of observation numbers and runs coxph
     #            and predict to determine set with median sum of squared errors. 
     #
     # INPUT      df1              Data frame being analyzed by forward search. 
     #            df1.ls           List of df1 by factor subset or NULL
     #            inner.rank                Number of observations overall (yesfactor=F) or for each factor subset (yesfactor=T)
     #            initial.sample   Number of random samples from which to take rim
     #            f.e              Right hand side of formula for Surv function
     #            cphties          ties of coxph
     #            ycol             Response column number                                    
     #            b.d              begin.diagnose Ranges from 25 to 45
     #
     yesfactor=TRUE
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      cStep1           "    
# print(" ");print("in cStep1")
  
                           if(b.d <=39 ){ print("",quote=FALSE);print(paste(spacehere,"Section 39",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(dim(df1));
                                 Hmisc::prn(inner.rank);Hmisc::prn(initial.sample);Hmisc::prn(f.e);Hmisc::prn(cphties);Hmisc::prn(ycol)   }

     ########################################################################################################
     # coxph and/or predict.coxph won't work if only one predictor observation is used, so pull at least 2  #
     # THIS WILL LIKELY CREATE FALSE TRANSITION FROM STEP 1 TO STEP 2. IGNORE IT.                           #
     ########################################################################################################
     inner.rank[inner.rank==1] <- 2

     # Set up final holding variables
     nobs <- dim(df1)[1]
     Observation <- 1:initial.sample
     pullN <- max(inner.rank)     # random sample will have a common number of obs from each factor subset

                                if(b.d <= 40){print(" ", quote = FALSE);print(paste(spacehere,"Section 40",sep=" "),quote=FALSE);
                                      Hmisc::prn(inner.rank)       }

     ################################################
     # Ensure that hold.cands is a list of matrices #
     ################################################
     if(yesfactor){
          nfactsub <- length(df1.ls)
          hold.cands <- candprep(yf=yesfactor, fixd.ls=df1.ls, inner.rank=inner.rank, preprnk=pullN, in.sam=initial.sample, b.d=b.d)     # candprep list

                           if(b.d <=39 ){ print("",quote=FALSE);print(paste(spacehere,"Section 39",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(dim(df1))      }

     }        #   yesfactor TRUE
     #
     ################################################################################################### 
     # Run the coxph function on each row of hold.cands and predict to the entire database for each one   #
     # Identify the median of the squared errors of each predictio.                                    #
     # Calculate the sum of squares of the error for each one amd store it in SSE.                     #
     # Sort this matrix on the 2nd column and locate the median row.  Save this value. Sort the matrix #
     # on this value. The observations on the first row (lowest median SSE are the rim for Step 1      #                                                                  #
     ################################################################################################### 
     if(yesfactor){
          SSErim <- NULL
          rimout.ls <- vector("list", nfactsub)
          thisform <- paste("survival::Surv(time=event.time, event=status)", f.e, sep=" ~ ")
          thisform <- stats::as.formula(thisform)

                           if(b.d <=40 ){ print("",quote=FALSE);print(paste(spacehere,"Section 40",sep=" "),quote=FALSE);
                                 Hmisc::prn(hold.cands);Hmisc::prn(thisform)     }

          for(j in 1:nfactsub){
               hold.cands.sub <- hold.cands[[j]]                   # pick factor subset
               hold.cands.sub <- hold.cands.sub[ ,(1:inner.rank[j])]    #   no longer contains 0 column
               if(inner.rank[j]==1){
                    hold.cands.sub <- matrix(hold.cands.sub,ncol=1)
               }
               SSE <- rep(-99, initial.sample)
 

                           if(b.d <=42 ){ print("",quote=FALSE);print(paste(spacehere,"Section 42",sep=" "),quote=FALSE);
                                 Hmisc::prn(j);Hmisc::prn(hold.cands.sub)     }


              for(r in 1:initial.sample){
                    index <- hold.cands.sub[r,]              # index is set of observation numbers, not row numbers

                    xmat <- match(index, df1[,1])
                    smalldata <- df1[xmat,]
                    indexuniv <- smalldata$holdISG[1]
                    universe <- df1[df1$holdISG==indexuniv,]
                    MED <- floor(dim(universe)[1]/2 + .00001)

                    lmsmall <- do.call(what=survival::coxph, args=list(formula=thisform, data=smalldata,
                          ties=cphties, model=TRUE, singular.ok=TRUE, x=TRUE, y=TRUE) )                        #    coxph
                    predsmall <- stats::predict(object=lmsmall, newdata=universe)                       # predict.coxph

                    errorsmall <- universe[, ycol] - predsmall
                    sserrorsmall <- sort(errorsmall^2) + universe$wiggle 
                    SSE[r] <- sserrorsmall[MED]
               }     # r
                                if(b.d <= 49){print("", quote = FALSE);print(paste(spacehere,"Section 49",sep=" "),quote=FALSE);
                                      Hmisc::prn(j);Hmisc::prn(universe);Hmisc::prn(predsmall);Hmisc::prn(SSE)       }

               hold.cands.sub <- cbind(SSE, hold.cands.sub)
               hold.cands.sub <- hold.cands.sub[order(hold.cands.sub[,1]),]
               uu <- hold.cands.sub[,-1]
               SSErim <- cbind(SSErim, uu[1,])
          }    #  j 

          rimout <- sort(SSErim)
     }      #   yesfactor
 print("leaving cStep1")
#stop("cStep1")
     return(rimout)
}
