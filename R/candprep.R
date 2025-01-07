candprep <-
function (yf, dfa2=NULL, fixd.ls=NULL, preprnk, inner.rank, in.sam, makearray=FALSE, b.d) 
{
     #                                      candprep
     #
     # VALUE    Samples of observations in the form of an array or list, SEE BELOW
     #               This is an oversample that will later be culled in xStep1.
     #
     # INPUT     yf          Logical.  TRUE if there are factors in the dataset
     #           dfa2        Data frame. Used only by non makearray. See below
     #           fixd.ls     List , See below
     #           preprnk     Single constant. Rank that accounts for construction variables.
     #           inner.rank  Vector. Number of observations to pull from each source.
     #           in.sam      initial.sample
     #           makearray   Logical. TRUE causes output to be an array for forsearch_lme
     #           b.d         begin diagnose

     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      candprep       "    
# print("in candprep")

                     if(b.d <=34 ){ print("",quote=FALSE);print(paste(spacehere,"Section 34",sep=" "),quote=FALSE);
                                 Hmisc::prn(yf);Hmisc::prn(dfa2);Hmisc::prn(fixd.ls);Hmisc::prn(preprnk);Hmisc::prn(inner.rank);
                                 Hmisc::prn(in.sam);Hmisc::prn(makearray)   }

     ########################################################
     # Restructure the input to a commonly structured list  #
     #     one with 4 dimensions. nA, nB, nC and nD are     #
     # in.sam, preprnk, nfactsub, ngroups                   #
     ########################################################
     nA <- in.sam
     nB <- preprnk

     if(makearray){
              samplist <- fixd.ls

     }             #  makearray = TRUE
     else{
         ###############
         # No makelist #
         ###############
         if(yf){ nD <- 1; nC <-length(fixd.ls)
               outlist <- vector("list", nD)
               Xlist  <- vector("list", nC)
               outlist[[1]] <- Xlist 
               for(j in 1:nC){ 
                    outlist[[1]][[j]] <- fixd.ls[[j]]    # yes.see input below
               }
               samplist <- outlist
          }        # factors present
          else{ 
              nD <- 1 
              nC <- 1
              outlist <- vector("list", nD)
              Xlist  <- vector("list", nC)
              Xlist[[1]] <- dfa2                 #yes, see input below
              outlist[[1]] <- Xlist
              samplist <- outlist
          }        # no factors
     }             #  makearray = FALSE

                     if(b.d <=35 ){ print("",quote=FALSE);print(paste(spacehere,"Section 35",sep=" "),quote=FALSE);
                                 ;Hmisc::prn(samplist)     }

     #
     ##############################################
     # Sample the database per the specifications #
     #   samplist[[2]][[5]] is for group 2        #
     ##############################################
     nD <- length(samplist)
     nC <- length(samplist[[1]])
     candarray <- vector("list",nD)               # outer list
     for(i in 1:nD){
          candarray[[i]] <- vector("list", nC)    # inner list  
     }
     for(j in 1:nD){
          for(k in 1:nC){
               thisset <- samplist[[j]][[k]]    # OK
               tempsamp <- matrix(-99, nrow=in.sam, ncol=preprnk + 1)
               for(m in 1:in.sam){
                    tempsamp[m,] <- c(sample(thisset[,1], preprnk, replace=FALSE), 0)
               }   # m
               candarray[[j]][[k]] <- tempsamp
          }        # k
     }             # j
 
                     if(b.d <=36 ){ print("",quote=FALSE);print(paste(spacehere,"Section 36",sep=" "),quote=FALSE);
                                 Hmisc::prn(candarray)     }

     #
     ####################################
     # Convert to list if not makearray #
     ####################################
     if(!makearray){
          outf <- vector("list", nC)
          for(i in 1:nC){
               outf[[i]] <- candarray[[1]][[i]]
          }
          candarray <- outf

                      if(b.d <=38 ){ print("",quote=FALSE);print(paste(spacehere,"Section 38",sep=" "),quote=FALSE);
                                 Hmisc::prn(candarray)   }
     }   # if ! makearray
     #
#######################################################################################
#######################################################################################
########### INPUT

# If makearray, input is a list, regardless of yesfactor
# Example with 2 groups, 3 factor subsets
# [[3]][[2]]
#    Observation G1  F1       C1        y       agony groupISG fixedISG
# 51          51  8 300 1.826639 330.5434 0.001157884       _8     _300
# 52          52  8 300 5.528472 366.2486 0.004975073       _8     _300
# 53          53  8 300 2.079564 331.8921 0.003539868       _8     _300

# If not makearray, 
#   if yf, input is a list by factor subset (fixd.ls)
#   if not yf, input is a data frame        (dfa2)
#############################################
########### OUTPUT

# if makearray, output is 2-layer list  
#   if yf   7 groups, 5 factor subsets
#   [[7]][[5]]
#         [,1] [,2] [,3]
#   [1,]  140   70    0
#   [2,]   70  140    0

# if not yf,  3 groups, 0 factors (note the 1 default)
#[[3]]
#[[3]][[1]]
#   Observation        C1        C2  G1         y       wiggle fixedISG groupISG      comboISG
#17          17 0.7297850 0.2752169 200  8.050512 0.0078237682    _None     _200 _None_F,G_200
#18          18 0.3655336 0.6287630 200 11.890049 0.0014210136    _None     _200 _None_F,G_200

1# if not makearray, output is a list
#   if yf, the list is a set of matrices by factor subset
#   if not yf, the list is a 1 level matrix
############################################
# STRATEGY
# Make a common 4-level matrix out of each input and perform the same sampling loops
########################################################################################

# print("leaving candprep")

     return(candarray)
}
