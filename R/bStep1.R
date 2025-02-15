bStep1 <-
function (yesfactor, df1, df1.ls, groups, inner.rank, source, initial.sample, nofactform=NULL, formulaA, randform, inc, ycol, b.d) 
{

     #                                    bStep1   
     # 
     # VALUE      Produces rim for Step 1. Uses candprep function to form matrix of candidate sets of observation numbers
     #            and runs lme and predictor to determine set of observation numbers with median sum of squared errors.
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            df1            Data frame being analyzed by forward search. 
     #            df1.ls         List of df1 by factor subset or NULL
     #            groups         groupISG
     #            inner.rank     Vector. Number of observations to pull from each factor subset, or from overall database otherwise
     #            source         Data frame of subsets. Only used for information, not calculations
     #            initial.sample Number of random samples from which to take rim
     #            nofactform     2-sided formula without factors                   NEEDED???
     #            formulaA       Formula for all effects including factors and constructed variables    
     #            randform       
     #            inc            Logical, TRUE causes relaxation of lmeControl
     #            ycol           Response column number
     #            b.d            begin.diagnose Ranges from 20 to 32
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      bStep1           "    
# print("in bStep1") 

                     if(b.d <=31 ){ print("",quote=FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                         Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(inner.rank);
                         Hmisc::prn(source);Hmisc::prn(initial.sample);Hmisc::prn(nofactform);Hmisc::prn(dim(df1));
                         Hmisc::prn(formulaA);Hmisc::prn(ycol);print("End of bStep1 argument listing")   }

     if(inc){
          ##################
          # Set lmeControl #
          ##################
          utils::str(lCtr <- nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, tolerance = 1e-2, 
               niterEM = 1000, msMaxEval = 1000, msTol = 1e-2, optimMethod = "L-BFGS-B",
               msVerbose = FALSE, returnObject = FALSE) )
               do.call(nlme::lmeControl, lCtr)
     }
     #
     Nlevels <- length(groups)
     nfacts <- length(inner.rank)
     ngroups <- length(df1.ls)
     listbyGroups <- vector("list", ngroups)
     hold.cands.pool <- vector("list", length(inner.rank))
     fixed <- unique(df1$fixedISG)
     nfixed <- length(fixed)
     newlist <- vector("list", nfixed)
     pullN <- max(inner.rank)

                 if(b.d <= 32){print("", quote = FALSE);print(paste(spacehere,"Section 32",sep=" "),quote=FALSE);
                                      Hmisc::prn(inner.rank);Hmisc::prn(pullN)       }

     ############################################################
     # Get a 2-level list each element of which is a matrix of  # 
     # sample observatons (initial.sample x factor/group subset #
     ############################################################
     if(yesfactor){
          for(mm in 1:nfixed){
               newlist[[mm]] <- df1[df1$fixedISG==fixed[mm],]
          }    #   mm
          hold.cands <- candprep(yf=yesfactor, dfa2=df1, fixd.ls=df1.ls, preprnk=pullN, 
                    inner.rank=inner.rank, in.sam=initial.sample, makearray=TRUE, b.d=b.d)           # candprep
     }     # yesfactor
     else{
               hold.cands <- candprep(yf=yesfactor, dfa2=df1, fixd.ls=df1.ls, preprnk=pullN, in.sam=initial.sample, 
                      inner.rank=inner.rank, makearray=TRUE, b.d=b.d)                                 # candprep
     }      # else, no factors

     ###########################################
     # now we have a set of candidate rims.  
     # Reduce the structure to a single matrix #
     ###########################################
     MED <- floor( dim(df1)[1]/2 + .00001 )
     SSE <- rep(-99, initial.sample)
     concatout <- matrix(-99, nrow=initial.sample, ncol= sum(inner.rank))
     for(r in 1:initial.sample){
          m <- 1
          subrim <- NULL
          for(k in 1:nfixed){
               for(j in 1:Nlevels){
                    piddly <- hold.cands[[j]][[k]][r,]
                    piddly <- piddly[1:inner.rank[m]  ]
                    subrim <- c(subrim, piddly)
                    m <- m + 1
               }    # in j
          }         # in k
          concatout[r,] <- sort(subrim)
     }              # in r

                 if(b.d <= 49){print("", quote = FALSE);print(paste(spacehere,"Section 49",sep=" "),quote=FALSE);
                                      Hmisc::prn(hold.cands);Hmisc::prn(concatout)       }

     for(r in 1:initial.sample){
          smalldata <- df1[concatout[r,],]
          lmesmall <- nlme::lme(fixed = formulaA, data = smalldata, random = randform)         # lme
          predsmall <- stats::predict(lmesmall, data=df1)                                # predict

          error2 <- (df1[, ycol] - predsmall)^2
          error2 <- sort(error2)
          mederror2 <- error2[MED]
          SSE[r] <- mederror2
     }   # r
     augmented <- cbind(SSE,concatout)
     augmented <- augmented[order(augmented[,1]),]
     smallestmed <- augmented[1,]
     smallestmed <- smallestmed[-1]

                 if(b.d <= 49){print("", quote = FALSE);print(paste(spacehere,"Section 49",sep=" "),quote=FALSE);
                                      Hmisc::prn(smalldata);Hmisc::prn(predsmall);Hmisc::prn(SSE);Hmisc::prn(smallestmed)       }

     return(smallestmed)

}
