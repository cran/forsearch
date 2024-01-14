cStep2 <-
function (f.e, finalm, dfa2, ms, rnk2, ss, b.d) 
{
     #                                            cStep2
     #
     # VALUE        An updated list of the rim during Step 2 NOT HAVING A LIST OF LISTS ANY MORE
     #                 Second element of the primary list is saved lm output for subsequent extraction of statistics. 
     #                 For each level of the factor subsets, select the rnk observations with the smallest squared errors
     #                 Then pool the results and add enough of the remaining observations to bring the total to m+1.
     #
     # INPUT 
     #       f.e               cont.form.rhs
     #       fbl               fixdat.by.level  list of observations by factor level all status levels
     #       finalm            See VALUE above. finalm argument is the same but only for Step 1 values
     #       dfa2              Complete data frame being analyzed by forward search. Presence of Observation column has 
     #                             no effect on output
     #       ms                First subset to be defined
     #       rnk2              Rank of X matrix. For factors, this is rank with factors removed; ie, rank for each factor subset.
     #       ss                skip.step1
     #       b.d               Number at whidh to begin diagnostic listings
     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        cStep2               "

     nobs <- dim(dfa2)[1]
     # Delete observations by factor subset from initial set of observations (from Step 1).
     if(is.null(ss)){finalm[[ms-1]] <- finalm[[ms-1]][[1]]}

                            if(b.d <=60 ){print("",quote=FALSE); print(paste(spacer,"Section 60",sep=" "),quote=FALSE);
                                Hmisc::prn(f.e);Hmisc::prn(finalm);Hmisc::prn(utils::head(dfa2));Hmisc::prn(utils::tail(dfa2));Hmisc::prn(dim(dfa2));
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
          td.et <- thisdata$event.time
          td.st <- thisdata$status
          xform <- paste("survival::Surv(time=td.et, event=td.st)", f.e, sep=" ~ ")
          xform <- stats::as.formula(xform)

                            if(b.d <=61 ){print("",quote=FALSE); print(paste(spacer,"Section 61",sep=" "),quote=FALSE);
                                Hmisc::prn(xform);Hmisc::prn(thisdata)    }

          thiscph <- survival::coxph(formula=xform, data=thisdata, ties="efron", model=TRUE, singular.ok=TRUE, x=TRUE, y=TRUE)     # coxph
          fooResult[[i]] <- thiscph

          thispredict <- stats::predict(thiscph, fixdat.mod)
          fixdat.mod$diff2 <- (fixdat.mod$event.time - thispredict)^2
          fixdat.mod <- fixdat.mod[order(fixdat.mod$diff2),]
          firstobs <- NULL

                          if(b.d <=67 ){ print("",quote=FALSE);print(paste(spacer,"Section 67",sep=" "),quote=FALSE);
                                   Hmisc::prn(thiscph);Hmisc::prn(thispredict);Hmisc::prn(fixdat.mod$diff2)   }    
 
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
