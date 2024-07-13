eStep1 <-
function (df1, inner.rank, initial.sample, formula, start, contR, algo, ycol, b.d) 
{
     #                                    eStep1   
     # 
     # VALUE      Produces partial rim for Step 1. There will always be at least 1 factor subset.
     #
     # INPUT      df1            Subset data frame being analyzed. Presence of
     #                                 Observation, Phase, holdISG columns has no effect on output. 
     #            inner.rank     nobs     
     #            initial.sample Number of random samples from which to take rim
     #            formula        Formula for fixed effects in nls always without factors
     #            ycol           Response column number
     #            b.d            begin.diagnose Ranges from 25 to 45
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      eStep1        "    
     
                           if(b.d <=31 ){ print("",quote=FALSE);print(paste(spacehere,"Section 31",sep=" "),quote=FALSE);
                                 Hmisc::prn(df1);Hmisc::prn(inner.rank);
                                 Hmisc::prn(formula);Hmisc::prn(ycol);Hmisc::prn(start)   }

          data <- df1                  # load factor subset
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
           thisstart <- start
           thisform <- formula
           medaugx <- matrix(1, nrow=initial.sample, ncol= 2)
           for(i in 1:initial.sample){
                zlist[[i]] <- result                                                        # this is the initialized 0,0 matrix
                zlist[[i]][,1] <- sample(dataObs,dimdata.1, replace=FALSE)                  # sample permutation of observation numbers
          }     # i
#                           if(b.d <=33 ){ print("",quote=FALSE);print(paste(spacehere,"Section 33",sep=" "),quote=FALSE);
#                                Hmisc::prn(zlist)[[initial.sample]]   }

          for(i in 1:initial.sample){
               Potential <- matrix(0, nrow=initial.sample, ncol=inner.rank)
               obsindex <- sort(zlist[[i]][1:inner.rank,1])
 
                           if(b.d <=35 ){ print("",quote=FALSE);print(paste(spacehere,"Section 35,   Sample ", i, sep=" "),quote=FALSE);
                                         Hmisc::prn(inner.rank);Hmisc::prn(obsindex)   }

               Obs <- match(obsindex,data[,1])
               smalldata <- data[Obs,]
               ####################################################################
               # Run nonlinear model on inner group and predict to entire dataset #
               ####################################################################
               this.form <- formula
               lmsmall <- stats::nls(formula=this.form, data=smalldata, start = thisstart, 
                     control=contR, algorithm=algo, model=TRUE)                                                             # nls

               predsmall <- stats::predict(lmsmall, data, type="response", pred.var=1)                          # predict

               errorsmall <- data[, ycol] - predsmall
               sserrorsmall <- errorsmall^2
               rowz <- zlist[[i]]
               augx <- matrix(rowz[,1], nrow=dimdata.1, ncol=1)
               augx <- cbind(augx, sserrorsmall)
               augx <- augx[order(augx[,2]),]                            # order the squared errors
               medaugx[i,] <- augx[floor((dimdata.1 + inner.rank)/2),]   # place the median in a matrix

                           if(b.d <=37 ){print("",quote=FALSE); print(paste(spacehere,"Section 37    Sample: ", i, sep=" "),quote=FALSE);
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
}
