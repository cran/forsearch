picksome <-
function(subsetlist, nobs, initial.sample, n.obs.per.level, rank)
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
     #
     namesSL <- names(subsetlist)
     nopl <- n.obs.per.level
     lennamesSL <- length(namesSL)
     ################################################################################
     # Set up and populate matrix of observation numbers from categorical variables #
     ################################################################################
     nsample <- rank*nopl
     out <- array(0, dim=c(lennamesSL, nsample, initial.sample))
     for(k in 1:initial.sample){
          for(i in 1:lennamesSL){   
               ufromSL <- subsetlist[[i]][,1]
               nufromSL <- length(ufromSL)                                                   # maximum number of available observations in level
               if(nufromSL < nopl) stop("Too many observations requested per level")
               out[i,,k] <- sample(ufromSL, nsample, replace=FALSE)                         # sample put into array
          }            #  j
     }     # i
     out1 <- matrix(0,nrow=initial.sample, ncol=lennamesSL*nsample)                       # collapse out into matrix out1
     for(i in 1:initial.sample){
          outtemp <- NULL
          for(j in 1:lennamesSL){
               outtemp <- c(outtemp, out[j,,i])
          }     # j
          out1[i,] <- outtemp
     }    # i

#print("leaving picksome")
    return(out1)
}
