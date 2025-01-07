dStep2 <-
function(yf, f2, dfa2, onlyfactor=FALSE, ms, finalm, fbg, rnk2, ycol, fam, b.d) 
{
     #                                            dStep2   
     #
     # VALUE        A list of 2 levels. The first of these is a complete list of the rim during the search.
     #                 Second level of the primary list is the saved glm output for subsequent extraction of statistics. 
     #                 For each level  m  of the factor subsets, select the  rnk  observations with the smallest squared errors
     #                 within each factor subset, if yf=TRuE.  Otherwise, only overall 
     #                 Then pool the results and add enough of the remaining observations to bring the total to  m+1.
     #
     # INPUT   yf          yesfactor
     #         f2          formula
     #         dfa2        data frame with a holdISG column that is '_' if no factor subsets
     #         onlyfactor  TRUE indicates that there are no continuous independent variables. Causes use of agony variable     
     #         ms
     #         finalm
     #         fbg         dfa2 as list by factor subset or NULL
     #         rnk2        if yf=TRUE, vector of observations per factor subset in primary stage; otherwise
     #                            number from entire dataset
     #         ycol        number of response column
     #         fam         Family with link
     #         b.d
     #
     spacer <- "XXXXXXXXXXXXXXXXXXXXXXXXXXX        dStep2               "
    
     nobs <- dim(dfa2)[1]

                            if(b.d <=60 ){print("",quote=FALSE); print(paste(spacer,"Section 60",sep=" "),quote=FALSE);
                                Hmisc::prn(f2);Hmisc::prn(ms);Hmisc::prn(rnk2);Hmisc::prn(dfa2)    }

     ################################################
     # Fill in last values for finalm and fooResult #
     ################################################
     finalm[[nobs]] <- 1:nobs
     fooResult <- vector("list", nobs)
     fooResult[[nobs]] <- stats::glm(formula=f2, family=fam, data=dfa2, y=TRUE, x=TRUE)              # glm
     #
     ######################################E#####################################
     # Create all other entries to finalm and fooResult within each value of ms #
     ############################################################################
     rimtotal <- NULL

     for(i in ms:(nobs-1)){
          ############################################################
          # Generate current fooResult and predict to entire dataset #
          ############################################################
          thisrim <- finalm[[i-1]]   
          galaxy <- dfa2
          smalldata <- dfa2[thisrim,]

          thisglm <- fooResult[[i]] <- stats::glm(formula=f2, family=fam, data=smalldata, y=TRUE, x=TRUE)             # glm
          predsmally <- stats::predict.glm(thisglm, newdata=galaxy, type="response")                                # predict

          nholddev <- length(predsmally)
          Devs <- rep(-99,nholddev)
                    ##############################################################
                    # Substitute a deviance for each pair of response/predsmally #
                    ##############################################################         
                    if(fam[[1]]=="binomial"){
                         for(jj in 1:nholddev){
                              this.weight <- galaxy$bin.wts[jj]
                              this.y <- galaxy$proportion1[jj]*this.weight
                              this.pred <- predsmally[jj]*this.weight
                              if(this.pred==this.y){
                                   devsmall <- 0
                              }
                              else{
                                   devsmall <- devianceCode(this.y, this.pred, this.weight, fam=fam[[1]])       # devianceCode
                              }
                              Devs[jj] <- devsmall^2 + galaxy$wiggles[jj]     # * onlyfactor
                         }        # jj
                    }      # binomial
                    if(fam[[1]]=="Gamma"){
                         for(jj in 1:nholddev){
                              this.y <- galaxy[jj,ycol]
                              this.pred <- predsmally[jj]
                              devsmall <- devianceCode(obs=this.y, pred=this.pred, ni=NULL, fam=fam[[1]])       # devianceCode
                              Devs[jj] <- devsmall^2
                         }        # jj 
                    }     # Gamma

                      if(fam[[1]]=="poisson"){
                         for(jj in 1:nholddev){
                              this.y <- galaxy[jj,ycol]
                              this.pred <- predsmally[jj]
                              devsmall <- devianceCode(obs=this.y, pred=this.pred, ni=NULL, fam=fam[[1]])       # devianceCode
                              Devs[jj] <- devsmall^2
                         }   # jj
                    }       # fam = poisson 
          ##################################################################################
          # Assign the Devs to the observations in galaxy, sort and select this subgroup #
          ##################################################################################
          galaxy <- cbind(Devs, dfa2)
          galaxy <- galaxy[order(galaxy$Devs),]
          #
          if(yf){
               ###################################################################################
               # Pull primary within each factor subset and then concatenate them for next stage #
               ###################################################################################
               nfactsub <- length(fbg)
               cumrim <- NULL
               namesISG <- unique(dfa2$holdISG)
               for(j in 1:nfactsub){
                    keeprim <- NULL 
                    thisuniv <- galaxy[galaxy$holdISG == namesISG[j],]
                    pullobs <- thisuniv[1:rnk2[j],2]
                    cumrim <- c(cumrim,pullobs)
               }    #  j
          }    #   yf
          else{
               ###############################################
               # Pull the primary from the entire data frame #
               ###############################################

               galaxy2 <- galaxy[1:rnk2,]
               cumrim <- galaxy2[,2]
          } # no factors   

                           if(b.d <= 69 ){ print("",quote=FALSE);print(paste(spacer,"Section 69",sep=" "),quote=FALSE);
                                Hmisc::prn(galaxy);Hmisc::prn(rnk2);print(" ");print("The following is rim entering dStep2:");
                                Hmisc::prn(thisrim); print(" ");print("The following is primary in Step2:");Hmisc::prn(cumrim);     stop("in dStep2 b.d 69")   }

               ##################################################
               # Append enough observations to bring total to i #    
               ##################################################
               needed <- i - length(cumrim)
               rimout <- cumrim
               if(needed > 0){
                    galaxy <- galaxy[order(galaxy[,2]),]     # reorder
                    galaxy <- galaxy[-cumrim,]
                    galaxy <- galaxy[order(galaxy[,1]),]
                    additional <- galaxy[1:needed,2]
                    rimout <- c(rimout, additional)
               }     #   needed > 0
               finalm[[i]] <- sort(rimout)
     }    #  i in ms: nobs-1

     outlist <- list(finalm, fooResult)

     return(outlist)
}
