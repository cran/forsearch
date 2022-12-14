#' @export
variablelist <-
function(datadf, verbose=TRUE)
{
     #                          variablelist
     #
     # VALUE     List of 2-column matrices, each of which is a level of the level1 * level2 * level3 * . . . 
     #               possible levels of the character variables
     #
     # INPUT          datadf     Data frame of independent variables in analysis
     #                verbose
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running variablelist", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     dimdata <- dim(datadf)
     nrows <- dimdata[1]
     ncols <- dimdata[2]
     SubsetCode <- rep("_", nrows)
     for(j in 2:ncols){
          nlevs <- levels(datadf[,j])
          if(!is.null(nlevs)){
               for(i in 1:nrows){
                    SubsetCode[i] <- paste(SubsetCode[i], datadf[i,j], sep="_") 
               }
          }
     }
     Subsetsdf <- data.frame(datadf[,1], SubsetCode) 
     uSubsetCode <- sort(unique(SubsetCode))
     nSC <- length(uSubsetCode)
     Subsetlist <- vector("list", nSC)
     names(Subsetlist) <- uSubsetCode
     for(i in 1:nSC){
          Subsetlist[[i]] <- Subsetsdf[Subsetsdf[,2] == uSubsetCode[i],]
     }
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running variablelist", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     return(Subsetlist)
}
