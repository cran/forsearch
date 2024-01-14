variablelist <-
function(datadf, prank)
{
     #                          variablelist
     #
     # VALUE     List of 2-column matrices, each of which is a level of the level1 * level2 * level3 * . . . 
     #               possible levels of the character variables
     #
     # INPUT          datadf     Data frame of independent variables in analysis
     #                prank       Rank of X matrix continuous variables 
     #
# prn(datadf)
    dimdata <- dim(datadf)
     nnrows <- dimdata[1]
     nncols <- dimdata[2]
     SubsetCode <- rep("_", nnrows)
     for(j in 2:nncols){
          nlevs <- levels(datadf[,j])
          if(!is.null(nlevs)){
               for(i in 1:nnrows){
                    SubsetCode[i] <- paste(SubsetCode[i], datadf[i,j], sep="_") 
               }
          }
     }
     Subsetsdf <- data.frame(datadf[,1], SubsetCode) 

     ###########################################################################################
     # Message regarding removal of observations when it is the only example of a factor level #
     ###########################################################################################
     tableSub <- table(SubsetCode)
     dntableSub <- dimnames(tableSub)$SubsetCode
     ncells <- length(tableSub)
     anyonlyone <- any(tableSub < prank)

#     if(anyonlyone){
#          print("There is as least one instance where a level of the possibly crossed factors contains too few observations.")
#          print("It may be necessary to remove any such observations before forsearch analysis.")
#          print("The factor codes and their number of observations are:")
#          for(i in 1:ncells){
#               if(tableSub[[i]]==1){
#                    singleobs <- Subsetsdf[Subsetsdf$SubsetCode==dntableSub[i],]
#                    Hmisc::prn(singleobs)
#               }     # if
#          }     # for
#     }     # anyonlyone
     uSubsetCode <- sort(unique(SubsetCode))
     nSC <- length(uSubsetCode)
     Subsetlist <- vector("list", nSC)
     names(Subsetlist) <- uSubsetCode
     for(i in 1:nSC){
          Subsetlist[[i]] <- Subsetsdf[Subsetsdf[,2] == uSubsetCode[i],]
     }
     return(Subsetlist)
}
