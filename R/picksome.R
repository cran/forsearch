#' @export
picksome <-
function(subsetlist, nobs, initial.sample, n.obs.per.level, rank, verbose=FALSE)
{
     #                          picksome
     # 
     # VALUE              Matrix of observation numbers (dimensions:  initial.sample x (rank * n.obs.per.level) that 
     #                          could satisfy Step 1         
     #
     # INPUT        subsetlist      List, the names of which are codes for levels of character variables (all possible 
     #                                  crosses). Created by variablelist()
     #              nobs            Number of observations in dataset
     #              initial.sample  Number of random samples of observation numbers desired 
     #              n.obs.per.level Number of observations to take from each level
     #              rank            Rank of the X-matrix of the analysis
     #              verbose
     #
     MC <- match.call()
#print("In picksome")
     if(verbose) {
          print("", quote = FALSE)
          print("Running picksome", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     namesSL <- names(subsetlist)
     lennamesSL <- length(namesSL)
     ################################################################################
     # Set up and populate matrix of observation numbers from categorical variables #
     ################################################################################
     out <- array(rep(0,initial.sample*lennamesSL*n.obs.per.level), dim=c(initial.sample, lennamesSL, n.obs.per.level))
     for(i in 1:initial.sample){
          for(j in 1:lennamesSL){   
               ufromSL <- subsetlist[[j]][,1]                                                # observations in this level

               nufromSL <- length(ufromSL)                                                   # maximum number of available observations in level
               if(nufromSL < n.obs.per.level) stop("Too many observations requested per level")
               out[i,j,] <- sample(ufromSL, n.obs.per.level, replace=FALSE)                   # sample
          }            #  j
     }     # i
     out1 <- matrix(c(out),nrow=initial.sample)
     #
     ##############################################################################################
     # Add observations from total set of observations not yet included in each row of matrix out #
     ##############################################################################################
     p <- rank - lennamesSL * n.obs.per.level     # difference between rank of lm H matrix and number of observations obtained by product
     if(p > 0){              # needs augmentation 
          out2 <- matrix(0, nrow=initial.sample, ncol = p)
          for(i in 1:initial.sample){
               allobs <- 1:nobs
               possibles <- allobs[-1*out1[i,]]
               poss <- sample(possibles, p, replace=FALSE)
               out2[i,] <- poss[1:p]
          }
          out <- cbind(out1, out2)
     }
     else out <- out1     
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running picksome", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(out)
}
