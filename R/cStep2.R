#' @export
cStep2 <-
function(df1, rim, formula.elements, r2){
     #
     #                                      cStep2
     #
     # INPUT
     #       df1         Complete data frame. This data frame is never changed in order that the row number and Observation number coincide
     #       rim         Vector of row numbers already in subset
     #               df1.rim is the subset of df1 that is already defined 
     #               df1.id is a 3-column matrix of N rows, containing the Observation number, an indicator of whether the element in the 1st column is in df1.id
     #               (1=Yes, 0=No), and the sum of squared residuals that results from
     #               running coxph on each subset that is made up of df1.rim and one of the other rows of df1. There is a vector containing a logical indicator 
     #               of whether this set of observations has no redundancies. This matrix is initialized by setting the 3rd
     #               column to -9. The vector is initialized to TRUE. 
     #      r2       Argument for Hmisc::redun. Higher values more difficult to declare redundancy.

     # OUTLINE OF ROUTINE TO DEFINE STEP 2

     # df1 is predefined
     # subset df1 to df1.rim using rim
     # define df1.id
     # Iterate on rows of df1.id:
     #    Define a test data frame as cbind(df1.rim, a row of df1). Skip rows with df1.id col2=1.
     #    Calculate coxph and extract residuals using method
     #    Square and sum. Enter result into appropriate row of df1.id, column 3
     #    Eliminate rows of df1.id with col 3 < 0.
     #    Sort df1.id by col 3.
     # Add col 1 of sorted df1.id to rim.
     #
     beg.diag.3 <- 100                                                       # Set diagnostic start here
     spacehere <- "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      cStep2           "    

     formula.rhs.rim <- paste(formula.elements, collapse=" + ")
     dimx1 <- dim(df1)[1]
     df1.rim <- df1[rim,]   # subset to rim

     df1.id <- matrix(0,nrow=dimx1,ncol=3)
     df1.id[,1] <- df1$Observation
     df1.id[rim,2] <- 1                    # set defaults
     df1.id[,3] <- -9

                                               if(beg.diag.3 <=10 ){ print(paste(spacehere,"Section 10",sep=" "),quote=FALSE);Hmisc::prn(df1.rim)    }
     Redundant.free <- rep(TRUE,dimx1)     # initialize

     coxph.out05 <- NULL
     for(uk in 1:(dimx1-1)){
          if(df1.id[uk,2] == 0){                               # just changed from !=
               tempdata <- rbind(df1.rim, df1[uk,])


                                               if(beg.diag.3 <=15 ){ print(paste(spacehere,"Section 15",sep=" "),quote=FALSE);Hmisc::prn(uk);
                                                     Hmisc::prn(df1.rim);Hmisc::prn(df1[uk,]);Hmisc::prn(utils::head(tempdata));Hmisc::prn(utils::tail(tempdata))   }
               rim.time <- df1.rim$event.time
               rim.status <- df1.rim$status

               xformc <- paste("survival::Surv(time=rim.time, event=rim.status, type='right')", formula.rhs.rim, sep=" ~  ")    # Surv
               formula.increment <- stats::as.formula(xformc)
                                               if(beg.diag.3 <=20 ){ print(paste(spacehere,"Section 20",sep=" "),quote=FALSE);Hmisc::prn(formula.increment);
                                                  Hmisc::prn(df1.rim);Hmisc::prn(tempdata[,1])   }
               #
               #####################################################
               # test for redundant variables before running coxph #
               #####################################################
               thisdata <- df1[tempdata[,1],]
               checkform <- paste("~ ",formula.rhs.rim, sep=" ")
               checkform <- stats::as.formula(checkform)
               redun.out <- Hmisc::redun(checkform, tempdata, r2=r2)
               redun.out4 <- redun.out[[4]][1]
               gotone <- match(redun.out4, formula.elements, nomatch = -99)
               if(gotone > -99) 
               {
                    Redundant.free[uk] <- FALSE
               }     # if gotone 

               if(Redundant.free[uk]){
                    xformc <- paste("survival::Surv(time=event.time, event=status, type='right')", formula.rhs.rim, sep=" ~  ")    # Surv
                    formula.increment <- stats::as.formula(xformc)
                    coxph.out05 <- survival::coxph(formula=formula.increment, data=tempdata, 
                        ties = "efron", singular.ok = TRUE, model = TRUE, x = TRUE, y = TRUE )                               # coxph
                    resids <- stats::residuals(coxph.out05, type="martingale")
                    df1.id[uk,3] <- sum(resids^2)          
               }    # Redundant.free
               else stop("not redundant free in cStep2")

          }       # df1.id[uk,2] !=1
                                               if(beg.diag.3 <=30 ){ print(paste(spacehere,"Section 30",sep=" "),quote=FALSE);Hmisc::prn(resids);
                                                     Hmisc::prn(df1.id)   }
     }    # uk
     index <- df1.id[,3] > 0 & Redundant.free
                                               if(beg.diag.3 <=40 ){ print(paste(spacehere,"Section 40",sep=" "),quote=FALSE);Hmisc::prn(index)   }
     remaining <- df1.id[,2]
     remaining0 <- remaining[remaining==0]
     n0 <- length(remaining0)
     if(n0==1){
          lastone <- df1.id[df1.id[,2]==0,]
          lastone <- as.data.frame(lastone)
          rim <- c(rim, lastone[1,1])
     }
     else{
          if(any(index)){
               df1.id <- df1.id[index,]
               # return rim, not in numerical order
               if(is.matrix(df1.id)){
                    df1.id <- df1.id[order(df1.id[,3]),]
                    rim <- c(rim, df1.id[1,1])
               }
               else{
                    rim <- c(rim, df1.id[1])
               }
          }    # any index
     }    #  n0 ==1
                                                 if(beg.diag.3 <=90 ){ print(paste(spacehere,"Section 90",sep=" "),quote=FALSE);Hmisc::prn(df1.id);
                                                           Hmisc::prn(rim);print("leaving cStep2")   }

     return(rim)
}
