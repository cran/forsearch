dStep1 <-
function (yesfactor, df1, inner.rank, initial.sample, formuladStep, fam, ycol, nopl, b.d) 
{
     #                                    dStep1   
     # 
     # VALUE      Produces rim for Step 1. If there are no factors, random choice 
     #            of inner.rank rows from entire dataset. If there are factors, use 
     #            variablelist and picksome to get a constrained choice of innersample.rows.
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            df1            Data frame being analyzed by forward search. Presence of
     #                                 Observation column has no effect on output, but numbers 
     #                                 in Observation may not be 1 through n. Data may be full dataset
     #                                 (no factors) or list of subsets (factors present)
     #            inner.rank     Rank of lm analysis on dataset with or without factor variables, depending on yesfactor     
     #            initial.sample Number of random samples from which to take rim
     #            formuladStep   Formula for fixed effects in lm always without factors
     #            fam            family
     #            ycol           Response column number
     #            nopl           n.obs.per.level 
     #            b.d            begin.diagnose Ranges from 25 to 45
     #
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      dStep1           "    
     
     ########################################################################################
     # Helper function for calculaton of deviance code for each observation based on family #
     ########################################################################################
     devianceCode <- function (obs, pred, ni=NULL, fam) 
     {
          if(obs/pred < 0) { out <- 0  }
          else{
               if(fam=="binomial"){
                    if(obs==0){ out <- ni*log(ni/(ni-pred))     }      # expects whole numbers, not proportions
                    else if(obs==ni){ out <- obs*log(obs/pred) }
                    else{ out <- obs*log(obs/pred) + (ni-obs)*log((ni-obs)/(ni-pred)) }
               }
               else if(fam=="Gamma"){ out <- -log(obs/pred) + (obs-pred)/pred  }
               else if(fam=="poisson"){ if(obs==0){ out <- pred }
                    else{ out <- obs*log(obs/pred) - obs + pred  }
               }
               else if(fam=="exponential"){ out <- -log(obs/pred) + obs*pred - 1 }
               else{
                    stop("family name not recognized")
               }
               if(is.na(out)){ out <- 0 }
               else{
                    if(abs(out)<= 10^(-12)){ out <- 0 }
                    else{ out <- 2 * out; vv <- obs - pred; vv <- vv/abs(vv)            #  1 or -1
                         out <- sqrt(out)*vv }
               }
          }
          return(out)
     }                                   # end of devianceCode function
     #
     ############################
     # Another support function #
     ############################
     supfun <- function(df2, initial.sample, inner.rank, ycol, nopl, formula ){

                          if(b.d <=32 ){ print("",quote=FALSE);print(paste(spacehere,"Section 32",sep=" "),quote=FALSE);
                               print("supfun");Hmisc::prn(utils::head(df2));Hmisc::prn(initial.sample);Hmisc::prn(inner.rank);
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
               #
               #########################################################################
               # Run general linear model on inner group and predict to entire dataset #
               #########################################################################
               this.form <- formula
               lmsmall <- stats::glm(formula=this.form, family=fam, data=smalldata, 
                        model=TRUE, x = TRUE, y=TRUE, singular.ok=TRUE)                                        #    glm
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

          rim <- c(zliststar[1:length(obsindex), 1])
          return(rim)
     }       # end supfun



#############################################    Main function starts here ##########################################

     formuladStep <- do.call("formula",list(formuladStep))

                           if(b.d <=31 ){ print("",quote=FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(df1);Hmisc::prn(inner.rank);
                                 Hmisc::prn(formuladStep);Hmisc::prn(fam);Hmisc::prn(ycol)   }

     if(yesfactor){
          nlevels <- length(df1)
          rim <- NULL
          for(ii in 1:nlevels){              # here is where we perform step 1 for each factor subset

               holdrim <- supfun(df2=df1[[ii]], initial.sample=initial.sample, inner.rank=inner.rank,ycol=ycol,nopl=nopl,formula=formuladStep)     # supfun

               rim <- c(rim, holdrim)
          }    #  ii
          rimout <- rim                      # single vector of observation numbers. Differs from aStep1 in no by subset list
     }       # major function section: factor present
     else{

          rimout <- supfun(df2=df1, initial.sample=initial.sample, inner.rank=inner.rank, ycol=ycol, nopl=nopl, formula=formuladStep)           #  supfun

     }       # major function section: no factor present

     return(rimout)
}
