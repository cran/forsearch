dStep1 <-
function (yesfactor, df1, df1.ls, inner.rank, initial.sample, formuladStep, fam, ycol, b.d) 
{
     #                                    dStep1  
     # 
     # VALUE      Produces rim for Step 1. If there are no factors, random choice 
     #            of inner.rank rows from entire datadf2set. If there are factors, use 
     #            variablelist and picksome to get a constrained choice of innersample.rows.   REWRITE
     #
     # INPUT      yesfactor      Logical. TRUE if there are factors in the X matrix 
     #            df1            Data frame being analyzed by forward search. Presence of
     #                                 Observation column has no effect on output, but numbers 
     #                                 in Observation may not be 1 through n. Data may be full datadf2set
     #                                 (no factors) or list of subsets (factors present)
     #            df1.ls         List of df1 with elements the factor subsets, or NULL
     #            inner.rank     Rank of lm analysis on datadf2set with or without factor variables, depending on yesfactor     
     #            initial.sample Number of random samples from which to take rim
     #            formuladStep   Formula for fixed effects in lm always without factors
     #            fam            family
     #            ycol           Response column number
     #            b.d            b.d 
     #
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      dStep1           "    
# print(" ");print("in dStep1")
 
                           if(b.d <=20 ){ print("",quote=FALSE);print(paste(spacehere,"Section 20",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(dim(df1));
                                 Hmisc::prn(inner.rank);Hmisc::prn(initial.sample);
                                 Hmisc::prn(formuladStep);Hmisc::prn(ycol)   }

     # Set up final holding variables
     nbreak <- NULL
     nobs <- dim(df1)[1]
     Observation <- 1:initial.sample
     pullN <- max(inner.rank)     # random sample will have a common number of obs from each factor subset

                                if(b.d <= 21){print(" ", quote = FALSE);print(paste(spacehere,"Section 21",sep=" "),quote=FALSE);
                                      Hmisc::prn(inner.rank)       }

     ################################################
     # Ensure that hold.cands is a list of matrices #
     ################################################
     if(yesfactor){
          nfactsub <- length(df1.ls)
          hold.cands <- candprep(yf=yesfactor, dfa2 = df1, fixd.ls=df1.ls, preprnk=pullN, inner.rank=inner.rank, 
                                in.sam=initial.sample, b.d=b.d)                                               # candprep list
          for(i in 1:nfactsub){
               uu <- hold.cands[[i]]
               uu <- cbind(uu,0)
               hold.cands[[i]] <- uu
         }    # i

                           if(b.d <=39 ){ print("",quote=FALSE);print(paste(spacehere,"Section 39",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(dim(df1))      }

     }        #   yesfactor TRUE
     #
     ################################################################################################### 
     # Run the glm function on each row of hold.cands and predict to the entire database for each one   #
     # Identify the median of the squared errors of each predictio.                                    #
     # Calculate the sum of squares of the error for each one amd store it in SSE.                     #
     # Sort this matrix on the 2nd column and locate the median row.  Save this value. Sort the matrix #
     # on this value. The observations on the first row (lowest median SSE are the rim for Step 1      #                                                                  #
     ################################################################################################### 
     if(yesfactor){
          SSErim <- NULL
          nbreak <- rep(-99, nfactsub)        # for count of number of obs in each factor subset
          rimout.ls <- vector("list", nfactsub)

                           if(b.d <=40 ){ print("",quote=FALSE);print(paste(spacehere,"Section 40",sep=" "),quote=FALSE);
                                 Hmisc::prn(hold.cands)     }
          for(j in 1:nfactsub){
               hold.cands.sub <- hold.cands[[j]] 
               SSE <- rep(-99, initial.sample)
 
                           if(b.d <=42 ){ print("",quote=FALSE);print(paste(spacehere,"Section 42",sep=" "),quote=FALSE);
                                 Hmisc::prn(j);Hmisc::prn(hold.cands.sub)     }

               for(r in 1:initial.sample){
                    if(inner.rank[j]==1){
                         index <- hold.cands.sub[r]
                    }
                    else{
                         index <- hold.cands.sub[r,1:inner.rank[j]]              # index is set of observation numbers, not row numbers
                    }
                    xmat <- match(index, df1[,1])
                    smalldata <- df1[xmat,]
                    smalldata <- smalldata[1:inner.rank,]
                    indexuniv <- smalldata$holdISG[1]
                    universe <- df1[df1$holdISG==indexuniv,]
                    MED <- floor(dim(universe)[1]/2 + .00001)

                    lmsmall <- do.call(what=stats::glm, 
                                args=list(formula=formuladStep, family=fam, data=smalldata, 
                                 model=TRUE, x = TRUE, y=TRUE) )                                                       #    glm
                    predsmall <- stats::predict(lmsmall, newdata=universe, type="response")                            # predict
 
                    Devs <- rep(-99,dim(universe)[1])                   
                    holddeviance <- data.frame(universe, predsmall, Devs)
                    #
                    #############################################################
                    # Substitute a deviance for each pair of response/predsmall #
                    #############################################################         
                    nholddev <- length(predsmall)
                    if(fam[[1]]=="binomial"){
                         for(jj in 1:nholddev){
                              this.weight <- universe$bin.wts[jj]
                              this.y <- universe$proportion1[jj]*this.weight
                              this.pred <- predsmall[jj]*this.weight
                              if(this.pred==this.y){
                                   devsmall <- 0
                              }
                              else{
                                   devsmall <- devianceCode(this.y, this.pred, this.weight, fam=fam[[1]])                            # devianceCode
                              }
                              Devs[jj] <- devsmall^2 + universe$wiggle
                         }        # jj
                    }      # binomial
 
                    if(fam[[1]]=="Gamma"){
                         nholddev <- length(predsmall)
                         for(jj in 1:nholddev){
                              this.y <- universe[jj,ycol]
                              this.pred <- predsmall[jj]
                              devsmall <- devianceCode(obs=this.y, pred=this.pred, ni=NULL, fam=fam[[1]])       # devianceCode
                              Devs[jj] <- devsmall^2 + universe$wiggle
                         }        # jj 
                    }     # Gamma

                      if(fam[[1]]=="poisson"){
                         nholddev <- length(predsmall)
                         for(jj in 1:nholddev){
                              this.y <- universe[jj,ycol]
                              this.pred <- predsmall[jj]
                              devsmall <- devianceCode(obs=this.y, pred=this.pred, ni=NULL, fam=fam[[1]])       # devianceCode
                              Devs[jj] <- devsmall^2 + universe$wiggle
                         }   # jj
                         Devs <- sort(Devs)     # sort in order to determine median
                    }       # fam = poisson

                    SSE[r] <- Devs[MED]
               }     # r

                        if(b.d <= 49){print("", quote = FALSE);print(paste(spacehere,"Section 49",sep=" "),quote=FALSE);
                                      Hmisc::prn(j);Hmisc::prn(universe);Hmisc::prn(predsmall);Hmisc::prn(SSE)       }

              hold.cands.sub <- cbind(SSE, hold.cands.sub)
              hold.cands.sub <- hold.cands.sub[order(hold.cands.sub[,1]),]
              uu <- hold.cands.sub[,-1]
              uu <- uu[1,]
              uuv <- c(uu)
              uuv <- uuv[uuv > 0]
              SSErim <- c(SSErim, uuv)
              nbreak[j] <- length(uuv)
          }    #  j 
          SSErim <- SSErim[SSErim > 0]
          rimout <- sort(SSErim)
     }      #   yesfactor
     else{

          nfactsub <- 1
          hold.cands <- candprep(yf=yesfactor, dfa2=df1, preprnk=pullN, inner.rank=inner.rank, 
                            in.sam=initial.sample, b.d=b.d)                                                # candprep list
          hold.cands <- hold.cands[[1]]
                           if(b.d <=43 ){ print("",quote=FALSE);print(paste(spacehere,"Section 43",sep=" "),quote=FALSE);
                                 Hmisc::prn(yesfactor);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(dim(df1));
                                 Hmisc::prn(hold.cands)      }
     #
     ################################################################################################### 
     # Run the glm function on each row of hold.cands and predict to the entire database for each one   #
     # Identify the median of the squared errors of each predictio.                                    #
     # Calculate the sum of squares of the error for each one amd store it in SSE.                     #
     # Sort this matrix on the 2nd column and locate the median row.  Save this value. Sort the matrix #
     # on this value. The observations on the first row (lowest median SSE are the rim for Step 1      #                                                                  #
     ################################################################################################### 
          SSErim <- NULL
          rimout.ls <- vector("list", nfactsub)

                           if(b.d <=44 ){ print("",quote=FALSE);print(paste(spacehere,"Section 44",sep=" "),quote=FALSE);
                                 Hmisc::prn(hold.cands)     }

               hold.cands.sub <- hold.cands                   # pick factor subset
               SSE <- rep(-99, initial.sample)

                           if(b.d <=45 ){ print("",quote=FALSE);print(paste(spacehere,"Section 45",sep=" "),quote=FALSE);
                                 Hmisc::prn(hold.cands.sub)     }

              for(r in 1:initial.sample){
                    index <- hold.cands.sub[r,]              # index is set of observation numbers, not row numbers
                    xmat <- match(index, df1[,1])
                    smalldata <- df1[xmat,]
                    indexuniv <- smalldata$holdISG[1]
                    universe <- df1
                    thisform <- stats::formula(formuladStep)
                    MED <- floor(dim(universe)[1]/2 + .00001)

                    lmsmall <- stats::glm(formula=thisform, family=fam, data=smalldata, 
                                 model=TRUE, x = TRUE, y=TRUE)                                                        # glm
                    predsmall <- stats::predict(lmsmall, newdata=universe, type="response")                            # predict

                    Devs <- rep(-99,dim(universe)[1])                   
                    holddeviance <-data.frame(universe, predsmall, Devs)
                    # Substitute a deviance for each pair of response/predsmall #
                    #############################################################
                    # Substitute a deviance for each pair of response/predsmall #
                    #############################################################
                    nholddev <- length(predsmall)
                    if(fam[[1]]=="binomial"){
                         for(jj in 1:nholddev){
                              this.weight <- universe$bin.wts[jj]
                              this.y <- universe$proportion1[jj]*this.weight
                              this.pred <- predsmall[jj]*this.weight
                              if(this.pred==this.y){
                                   devsmall <- 0
                              }
                              else{
                                   devsmall <- devianceCode(this.y, this.pred, this.weight, fam=fam[[1]])                            # devianceCode
                              }
                              Devs[jj] <- devsmall^2 + universe$wiggle
                         }        # jj
                    }      # binomial
 
                   if(fam[[1]]=="Gamma"){
                         nholddev <- length(predsmall)
                         for(jj in 1:nholddev){
                              this.y <- universe[jj,ycol]
                              this.pred <- predsmall[jj]
                              devsmall <- devianceCode(obs=this.y, pred=this.pred, ni=NULL, fam=fam[[1]])       # devianceCode
                              Devs[jj] <- devsmall^2 + universe$wiggle
                         }        # jj 
                    }     # Gamma

                    if(fam[[1]]=="poisson"){
                         nholddev <- length(predsmall)
                         for(jj in 1:nholddev){
                              this.y <- universe[jj,ycol]
                              this.pred <- predsmall[jj]
                              devsmall <- devianceCode(obs=this.y, pred=this.pred, ni=NULL, fam=fam[[1]])       # devianceCode
                              Devs[jj] <- devsmall^2 + universe$wiggle
                         }   # jj
                         Devs <- sort(Devs)
                    }       # fam = poisson

                    SSE[r] <- Devs[MED]
               }     # r

                                if(b.d <= 47){print("", quote = FALSE);print(paste(spacehere,"Section 47",sep=" "),quote=FALSE);
                                      Hmisc::prn(universe);Hmisc::prn(predsmall);Hmisc::prn(SSE)       }

               hold.cands.sub <- cbind(SSE, hold.cands.sub)
               hold.cands.sub <- hold.cands.sub[order(hold.cands.sub[,1]),]
               uu <- hold.cands.sub[,-1]
               SSErim <- cbind(SSErim, uu[1,])

               SSErim <- SSErim[SSErim > 0]
               rimout <- sort(SSErim)
     }       # no factors
     out <- list(rimout, nbreak)


# print("leaving dStep1")
     return(out)
}
