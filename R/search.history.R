#' @export
search.history <-
function(list1, verbose=TRUE)
{
     #                                 search.history
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running search.history", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     ##################################################
     # Obtain 1st element of list1 and set up outputs #
     ##################################################
     x1 <- list1$"Rows in stage"
     lenx1 <- length(x1)
     IN <- OUT <- vector("list",lenx1) 
     lenIN <- rep(0,lenx1)
     #
     #########################################################
     # Determine changes in observation sets as m progresses #
     #########################################################
     uu <- x1[[1]]
     if(!is.null(uu)) IN[[1]] <- uu
     for(i in 2:lenx1){
                        # uu is the current element; vv is the previous one
          uu <- x1[[i]]
          if(!is.null(uu)){                         #  keep going until find some IN entries
               vv <- x1[[i-1]]
               if(is.null(vv)){
                    IN[[i]] <- uu                   # vv can be NULL only at first nontrivial addition to set
               }
               else{
                    rr <- match(vv, uu, nomatch = 0)# rr and tt are vectors of POSITIONS of matches; rr is positions of elements of vv in uu
                    if(sum(rr)>0){
                         IN[[i]] <- uu[-1*rr]            # IN will be the later one(s) that don't match the earlier ones
                    }
                    else{
                         IN[[i]] <- uu                   #  all new observations replace all old ones
                    }
                    tt <- match(uu, vv, nomatch = 0)
                    if(sum(tt)>0){
                         OUT[[i]] <- vv[-1*tt]           # OUT wll be the early one(s) that don't match the later ones 
                    }
                    else{
                         OUT[[i]]<- vv                   #  all old ones eliminated
                    }
               }       #    vv not null
          }       #   uu  not null 
     }      #  i    
     #
     ################################################
     # Put IN and OUT into convenient output format #
     ################################################
     for(i in 1:lenx1){
          lenIN[i] <- length(IN[[i]])
     }
     maxIN <- max(lenIN)
     histIN <- matrix(0,nrow=lenx1,ncol=maxIN)
     histOUT <- matrix(0,nrow=lenx1,ncol=(maxIN-1))
     for(i in 1:lenx1){
          uu <- c(IN[[i]],rep(0,maxIN))
          uu <- uu[1:maxIN]
          histIN[i,] <- uu
          vv <- c(OUT[[i]],rep(0,lenx1))
          vv <- vv[1:(maxIN-1)]
          histOUT[i,] <- vv
     }      #   i
     lasthist <- dim(histOUT)[2]
     # histIN is already in the minimalist columnar format

     INname <- paste("IN",1:maxIN,sep="")
     histIN <- as.data.frame(histIN)
     names(histIN) <- INname
     history <- histIN

     if(sum(histOUT)>0){
          orijOUT <- dim(histOUT)[2]
          for(j in (maxIN-1):2){
               if(sum(histOUT[,j])==0) {histOUT <- histOUT[,-j]}
               else {break}
          }    # j
          histOUT <- as.matrix(histOUT,nrow=maxIN)
          dimOUT2 <- dim(histOUT)[2]
          if(dimOUT2 > 0){
               XX <- rep(" //",lenx1)
               history <- cbind(histIN,XX,histOUT)
               history <- as.data.frame(history)
               OUTname <- paste("OUT",1:dimOUT2,sep="")
               names(history) <- c(INname,"XX",OUTname)
          }     # dimOUT > 0
     }   # sum histOUT > 0  
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running search.history", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
     outout <- list(history=history, Call=MC)
     return(outout)
}
