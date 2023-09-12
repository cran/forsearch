aStep2 <-
function (thislm, data, ycol, thisi) 
{
#                                            aStep2
#
# VALUE        List. first element is updated rim during Step 2. No constraints with respect to factors. 
#              Avoid going beyone rows.in.model file size by restriction in calling routine.
#
# INPUT thislm            lm object
#       data              Data frame being analyzed by forward search. Presence of Observation column has 
#                             no effect on output
#       ycol              Response column number, including 1 for Observation
#       thisi             Current iteration in Step 2; also, number of integers in output
#
# NOTE               Called by bStep2
#
     #####################################################################
     # Calculate the error for each observation using the current subset #
     # Select next subset.                                               #
     #####################################################################
#prn(data)
     dimx1 <- dim(data)[1]
#prn(dimx1)
     preds <- stats::predict(thislm, data, pred.var=1)                                 #   predict
     errors <- data[, ycol] - preds
     errors2 <- errors * errors
#prn(errors2)
     medaugx <- matrix(1:dimx1, nrow=dimx1, ncol=2, byrow=FALSE)
     medaugx[,2] <- errors2             
#prn(medaugx)
     medaugx <- medaugx[order(medaugx[,2]),]
#prn(medaugx)
     rim <- sort(medaugx[1:thisi, 1])
#prn(medaugx[1:thisi,2])
     score <- sum(medaugx[1:thisi,2])
#prn(score)
     outlist <- list(rim=rim, score=score)
#print("leaving aStep2")
     return(outlist)
}
