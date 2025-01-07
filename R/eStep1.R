eStep1 <-
function (yf, df1=NULL, df1.ls=NULL, df1fact.list, inner.rank=NULL, initial.sample, formulaE=NULL, nofactform, 
                   start, nofactstart, inner.r.vector=NULL, contR=FALSE, algo, ycol, b.d) 
{
   #                                    eStep1   
     # 
     # VALUE      Produces rim for Step 1. Uses candprep function to form matrix of candidate sets of observation numbers
     #            and runs nls and predictor to determine set with median sum of squared errors. 
     #
     # INPUT      yf             yesfactor
     #            df1            Data frame being analyzed by forward search. 
     #            df1.ls         2-lovel list of df1 by sector and factor subset
     #            df1fact.list   1-level list by factor subset
     #            inner.rank     Number of observations to pull from entire database
     #            initial.sample Number of random samples from which to take rim
     #            formulaE       Formula for all effects including factors and constructed variables    
     #            nofactform     Formula for all effects, omitting factors
     #            start
     #            nofactstart    Vector of starting values, omitting factors
     #            inner.r.vector Vector of number to draw within each factor
     #            contR          Logical
     #            algo           nls algorithm
     #            ycol           Response column number
     #            b.d            begin.diagnose Ranges from 25 to 45
     #
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      eStep1           "    
# print("in eStep1") 

                           if(b.d <=20 ){ print("",quote=FALSE);print(paste(spacehere,"Section 20",sep=" "),quote=FALSE);
                                 Hmisc::prn(yf);Hmisc::prn(utils::head(df1));Hmisc::prn(df1.ls);Hmisc::prn(df1fact.list);Hmisc::prn(inner.rank);
                                 Hmisc::prn(initial.sample);Hmisc::prn(dim(df1));Hmisc::prn(formulaE);
                                 Hmisc::prn(start);Hmisc::prn(inner.r.vector);Hmisc::prn(ycol);
                                 print("end of initial b.d in eStep1")   }
#stop("initial eStep1")
     ###################################
     # Set up final holding variables  #
     ###################################
     nobs <- dim(df1)[1]

     usection <- unique(df1$Section)     # will pull a number of obs to cover the sections completely 
     nusection <- length(usection)

     pullN <- inner.rank     # default. +1 to avoid nls crash
     if(yf){
          pullN <- max(inner.r.vector)            # candprep will pull this many for each factor/section subset
     }
     hold.cands <- candprep(yf=yf, dfa2=df1, fixd.ls=df1.ls, inner.rank=inner.r.vector, preprnk=pullN, 
                   in.sam=initial.sample, makearray=TRUE, b.d=b.d)                                       # candprep as 2-level list

     nfactorlevels <- length(hold.cands)
     SSE.outer <- rep(-99, initial.sample)
     SSE.inner <- matrix(-99, nrow=initial.sample, ncol=nfactorlevels)

     if(yf){
          ###############################################
          # Collapse all the sections of each subfactor #
          ###############################################
          hold.cands.sub <- vector("list", nfactorlevels)
          for(i in 1:nfactorlevels){
               tempmatrix <- NULL
               for(j in 1:nusection){
                    index <- 1:inner.r.vector[i]
                    tempmatrix <- cbind(tempmatrix, hold.cands[[i]][[j]][,index]    )
                    hold.cands.sub[[i]] <- tempmatrix
               }    #   j
          }    #    i
          #########################################################################
          # Calculate the SSE for each sample subset WITHIN its own factor subset #
          #########################################################################
          for(j in 1:nfactorlevels){
               hold.cands.c <- hold.cands.sub[[j]]   # candidate matrix within 1 factor subset
               markerhold <- df1[hold.cands.c[1,],][1,]
               markerhold <- markerhold$holdISG
               universe <- df1[df1$holdISG==markerhold,]
               ##############################################
               # If nls controls have been weakened and nls #
               # still doesn't converge, give notice and    #
               # ignore those candidates.                   #
               ##############################################
               for(r in 1:initial.sample){
                    xmat <- hold.cands.c[r,]       # 1 row out of initial.sample
                    #############################################
                    # Add a wild card to each candidate set     #
                    # if the candidate set is too small         #      SHOULD WE REALLY BE DOING THIS?
                    # Seems to need 1 more than starting number #
                    #############################################
                    if(length(xmat) <= length(nofactstart)){
                         uniObs <- universe$Observation
                         locat <- pmatch(xmat,uniObs)
                         wild <- uniObs[-locat]
                         wild <- sample(wild,1)
                         xmat <- c(wild, xmat)
                    }
                    smalldata <- df1[xmat,]
                    MED <- floor(dim(universe)[1]/2 + .00001)
                    MED[MED==0] <- 1
                    if(contR){
                         a.warn.occurred <- FALSE
                         tryCatch(
                           expr=(lmsmall <- stats::nls(formula=nofactform, data=smalldata, start=nofactstart, 
                           algorithm = algo, control=list(maxiter=1000, tol=.001, warnOnly=FALSE)))
                           , warning=function(w) {a.warn.occurred <- TRUE}
                         )                                                                                     #  nls
                         if(a.warn.occurred){                            
                              messprop <- "nls failed to converge"
                              print(messprop, quote=FALSE)
                              SSE.inner[r,j] <- -99
                         }
                         else{
                              predsmall <- stats::predict(lmsmall, newdata=universe)                            # predict
                              errorsmall <- df1[, ycol] - predsmall
                              sserrorsmall <- sort(errorsmall^2)  
                              SSE.inner[r,j] <- sserrorsmall[MED]
                         }
                    }    #  contR TRUE
                    else{
                         lmsmall <- stats::nls(formula=nofactform, data=smalldata, start=start, algorithm=algo)                                                 # nls
                         predsmall <  - stats::predict(lmsmall, newdata=universe)                                # predict
 
                         errorsmall <- df1[, ycol] - predsmall
                         sserrorsmall <- sort(errorsmall^2)  
                         SSE.inner[r,j] <- sserrorsmall[MED]
                    }# contR FALSE
               }     #   r
          }          #   j
          sumup <- SSE.inner == -99
          sumsumup <- sum(sumup)
          propsum <- 100 * sumsumup/(initial.sample * nfactorlevels)
          if(propsum > 0){
               print(" ", quote=FALSE)
               print("Testing squared errors within factors,", quote=FALSE)
               summess <- paste(sumsumup, "(", propsum, "%)", "candidate Step 1 observation sets failed to converge", sep=" ")
               print(summess, quote = FALSE)
               print(" ", quote = FALSE)
               maxSSEinner <- max(SSE.inner)
               index <- SSE.inner < 0               
               SSE.inner[index] <- maxSSEinner + 10
          }        # propsum > 0
          else{
               print(" ", quote=FALSE)
               print("Testing squared errors within factors, nls calculations for all candidate sets converged.", quote=FALSE)
               print(" ", quote=FALSE)
          }
     }     #    end of yf TRUE    EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
     #
     ###################################################
     # Compute outer SSE for both yf TRUE and yf FALSE #
     ###################################################
     hold.cands.b <- NULL
     for(jj in 1:nfactorlevels){
          for(kk in 1:nusection){
               uu <- hold.cands[[jj]][[kk]]
               dimuu2 <- dim(uu)[2]
               uu <- uu[,-dimuu2]
               hold.cands.b <- cbind(hold.cands.b, uu)
          }   #  kk
     }        #  jj  
     #
     ##############################################
     # Calculate the outer median squared error   #
     # If nls controls have been weakened and nls #
     # still doesn't converge, give notice and    #
     # ignore those candidates.                   #
     ##############################################
     universe <- df1
     MED <- floor(dim(universe)[1]/2 + .00001)
     MED[MED==0] <- 1
     for(r in 1:initial.sample){
          xmat <- hold.cands.b[r,]
          smalldata <- df1[xmat,]
          if(contR){
               a.warn.occurred <- FALSE
               an.error.occurred <- FALSE
               tryCatch(
                    expr=(lmsmall <- stats::nls(formula=formulaE, data=smalldata, start=start, 
                    algorithm = algo, control=list(maxiter=1000, tol=.001, warnOnly=FALSE)))
                    , warning=function(w) {a.warn.occurred <- TRUE}
                    , error=function(e) {an.error.occurred <- TRUE}
               )                                                                              #  nls
               if(a.warn.occurred | an.error.occurred){                            
                    messprop <- "nls failed to converge"
                    print(messprop, quote=FALSE)
                    SSE.outer[r] <- -999
                    if(a.warn.occurred){print("It was a warning") }
                    if(an.error.occurred) {print("It was an error")}
                    stop("failed")
               }
               else{
                    predsmall <- stats::predict(lmsmall, newdata=universe)                        # predict
                    errorsmall <- df1[, ycol] - predsmall
                    sserrorsmall <- sort(errorsmall^2)  
                    SSE.outer[r] <- sserrorsmall[MED]
              }
          }      #   contR
          else{
               lmsmall <- stats::nls(formula=formulaE, data=smalldata, start=start, algorithm=algo)  # nls
               predsmall <- stats::predict(lmsmall, newdata=universe)                                # predict

               errorsmall <- df1[, ycol] - predsmall
               sserrorsmall <- sort(errorsmall^2)  
               SSE.outer[r] <- sserrorsmall[MED]
          }    # not contR
     }     # r
     ##########################################################################################
     # Combine into a single SSE.inner and penalize any candidate set that failed to converge #
     ##########################################################################################
     SSE.inner <- apply(SSE.inner,1,mean)
     SSE.inner[SSE.inner < 0]<- max(SSE.inner) + 10
     #
     sumup <- SSE.outer == -99
     sumsumup <- sum(sumup)
     propsum <- 100 * sumsumup/initial.sample
     if(propsum > 0){
          print("Testing outer squared errors,", quote=FALSE)
          summess <- paste(sumsumup, "(", propsum, "%)", "candidate Step 1 observation sets failed to converge", sep=" ")
          print(summess, quote = FALSE)
          print(" ", quote = FALSE)
          maxSSEouter <- max(SSE.outer)
          index <- SSE.outer < 0               
          SSE.outer[index] <- maxSSEouter + 10
     }   # propsum > 0
     else{
          print("Testing outer squared errors, nls calculations converged for all candidate sets.", quote=FALSE)
          print(" ", quote=FALSE)
     }
          #
          #################################################################
          # Penalize any candidate set whose SSE.outer failed to converge #
          # Combine SSE.inner to SSE.outer if there are factors           #
          #################################################################
          SSE.outer[SSE.outer < 0] <- max(SSE.outer) + 10
          if(yf){
               SSE <- (SSE.inner + SSE.outer)/2
          }
          else{
               SSE <- SSE.outer
          }
          ##############################################################
          # Attach SSE to hold.cands.c and sort to identify lowest MED #
          ##############################################################
          rimoutx <- cbind(SSE, hold.cands.b)
          rimoutx <- rimoutx[order(rimoutx[,1]),]
          rimoutx <- rimoutx[,-1]
          rimoutx <- rimoutx[1,]
# print("leaving eStep1")

     return(sort(rimoutx))
}
