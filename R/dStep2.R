dStep2 <-
function(f2, dfa2, randm2, ms, finalm, fbg, b.d, rnk2, ycol, fam) 
{
     #                                            dStep2
     #
     # VALUE        A list of 2 levels. The first of these is a complete list of the rim during the search.
     #                 Second level of the primary list is the saved glm output for subsequent extraction of statistics. 
     #                 For each level  m  of the factor subsets, select the  rnk  observations with the smallest squared errors
     #                 Then pool the results and add enough of the remaining observations to bring the total to  m+1.
     #
     # INPUT 
     #     fam     Family with link
     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        dStep2               "

     nobs <- dim(dfa2)[1]

                            if(b.d <=60 ){print("",quote=FALSE); print(paste(spacer,"Section 60",sep=" "),quote=FALSE);
                                Hmisc::prn(f2);Hmisc::prn(utils::head(dfa2));Hmisc::prn(utils::tail(dfa2));Hmisc::prn(dim(dfa2));
                                Hmisc::prn(ms);Hmisc::prn(rnk2)    }


     fooResult <- vector("list", nobs)
     for(i in ms:(nobs-1)){
          remainder <- NULL
          diff2 <- -999
          fixdat.mod <- data.frame(dfa2,diff2)
          sbsts <- unique(fixdat.mod$holdISG)
          nsubs <- length(sbsts) 
          rim <- finalm[[i-1]]
          thisdata <- fixdat.mod[rim,]

                            if(b.d <=61 ){print("",quote=FALSE); print(paste(spacer,"Section 61",sep=" "),quote=FALSE);
                                Hmisc::prn(thisdata)    }

          f3 <- f2
          thislme <- stats::glm(formula=f3, family=fam, data=thisdata, y=TRUE, x=TRUE)                                      # glm
          fooResult[[i]] <- thislme

          thispredict <- stats::predict(thislme, fixdat.mod)
          fixdat.mod$diff2 <- (fixdat.mod[, ycol] - thispredict)^2
          fixdat.mod <- fixdat.mod[order(fixdat.mod$diff2),]
          firstobs <- NULL

                          if(b.d <=67 ){ print("",quote=FALSE);print(paste(spacer,"Section 67",sep=" "),quote=FALSE);
                                   Hmisc::prn(thislme);Hmisc::prn(thispredict);Hmisc::prn(fixdat.mod$diff2)   }    
 
          for(j in 1:nsubs){
               candidates <- fixdat.mod[fixdat.mod$holdISG==sbsts[j],]
               firstobs <- rbind(firstobs, candidates[1:rnk2,])
               candidates <- candidates[-(1:rnk2),]
               remainder <- rbind(remainder, candidates)
          }     #   j
          remainder <- remainder[order(remainder$diff2),]
          needed <- i - rnk2*nsubs
          needed.obs <- remainder$Observation[1:needed] 
          finalm[[i]] <- c(firstobs$Observation, needed.obs)
     }    #   i

                          if(b.d <=80 ){ print("",quote=FALSE);print(paste(spacer,"Section 80",sep=" "),quote=FALSE);
                                   Hmisc::prn(finalm);Hmisc::prn(fooResult)   }    
 
    outlist <- list(finalm, fooResult)
    return(outlist)
}
