#' @export
plotdiag.ANOX2 <-
function(forn, 
     anova.rows=NULL,
     ylab.extend=c("proportionality", "variance"),
     maintitle="Put main title here", 
     subtitle="Put subtitle here", 
     caption="Put caption here",
     wmf="Put_stored_name_here", 
     Cairo=TRUE,
     printgraph=TRUE,
     legend="Dummy legend name",
     verbose=TRUE)
{
     #                          plotdiag.ANOX2
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.ANOX2", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }    # verbose
 
     #################
     # Plot function #
     #################
#                                                           print("OK to plotdiag.ANOX2    10")
     plotB1 <- function(data, xcol, ycol, cov1col,
                     facetcol=NULL, facetdir, 
                     mtitle,
                     stitle,
                     horlabel,
                     vertlabel,
                     cap,
                     filewidth=5, fileheight=5,
                     legendname,
                     verbose)
     {

          XVAR <- data[,xcol]
          YVAR <- data[,ycol]
          COV1 <- data[,cov1col]
          dfplot <- data.frame(XVAR,YVAR,COV1)
          FACET <- data[,facetcol]
          dfplot <- data.frame(dfplot,FACET)
          out <- ggplot2::ggplot(data=dfplot,ggplot2::aes(XVAR,YVAR,COV1)) + ggplot2::geom_point(ggplot2::aes(shape = COV1))
          out$labels$shape <- legendname

          if(!is.null(facetcol)){
               if(facetdir=="h") out <- out + ggplot2::facet_grid(.~FACET)
               if(facetdir=="v") out <- out + ggplot2::facet_grid(FACET~.)
               if(facetdir=="w") out <- out + ggplot2::facet_wrap(~FACET)
          }      
  
          ############################################
          #    Add titles, axis labels, and caption. #
          ############################################
          out <- out + ggplot2::ggtitle(mtitle,subtitle=stitle) + ggplot2::xlab(horlabel) + ggplot2::ylab(vertlabel) + ggplot2::labs(caption=cap,legend=legendname)
          out <- out + ggplot2::labs(color=legendname)

          ############################
          # Print and save the graph #
          ############################
          if(Cairo){
               Cairo::CairoWin(width = 7, height = 7, pointsize = 12, record = getOption("graphics.record"),
                 rescale = c("R", "fit", "fixed"), bg = "transparent", canvas = "white", gamma = getOption("gamma"),
                 xpos = NA, ypos = NA, buffered = getOption("windowsBuffered"), restoreConsole = FALSE)
          }      # Cairo

          print(out)             # this line plots the graph

          if(printgraph){
                filename <- paste(wmf,".wmf",sep="")
                ggplot2::ggsave(filename,width=filewidth, height=fileheight)
                grDevices::dev.off()
          }    # printgraph
          #
          if(verbose) {
               print("", quote = F)
               print("Finished running plotB1", quote = F)
               print("", quote = F)
               print(date(), quote = F)
               print("", quote = F)
          }    #  verbose
          return()
     }      #  plotB1
     #
     ############################                                         #   prepstuff
     # Preparation for plotting #
     # Remove any na's          #
     # Remove initial zero rows #
     ############################
     prepstuff <- function(rightforn, yle){
#                                                                        print("OK to plotdiag.ANOX2    30")
          df2 <- rightforn
          indexna <- is.na(df2)
          indexna <- apply(indexna,1,any)
          df2 <- df2[!indexna,]
          ndf2 <- dim(df2)[1]

          for(ir in 1:ndf2){
               uu <- df2[1,-1] 
               if(sum(uu)==0){df2 <- df2[-1,]} else break
          }      #   ir


          nrows <- dim(df2)[1]
          df3 <- df2[,-1]
          ncols <- dim(df3)[2]

          if(is.null(ncols)){
               col1 <- df2[,1]
               col2 <- df3
               col3 <- "1"                       # default for only 1 anova row
               pvals <- as.data.frame(tibble::tibble(col1,col2,col3))
               wmf2 <- paste(wmf,".wmf",sep="")    
#                                                                      print("OK to plotdiag.ANOX2    40")
               ##############################################
               # Call for plot using support function above #
               ##############################################
               if(yle=="p"){
                    plotB1(data=pvals, xcol=1, ycol=2, cov1col=3,
                        mtitle=maintitle,
                        stitle=subtitle,
                        horlabel="Subset size m",
                        vertlabel="p-values for test of proportionality",
                        cap=caption,
                        legendname=legend,
                        verbose=verbose)

               }    #   yle=="p"
               if(yle=="v"){
                    plotB1(data=pvals, xcol=1, ycol=2, cov1col=3,
                        mtitle=maintitle,
                        stitle=subtitle,
                        horlabel="Subset size m",
                        vertlabel="p-values for ANOVA test of null hypothesis",
                        cap=caption,
                        legendname=legend,
                        verbose=verbose)
               }   # if yle == "v"
          }  # is.null ncols
          else{                                 # now for case in which anova rows have been subset
               if(yle=="p"){
#                                                                    print("OK to plotdiag.ANOX2    50")
                    cola <- df2[,1] 
                    col1 <- rep(cola, times=ncols)
                    namesdf3 <- names(df3)
                    col2 <- c(unlist(df3),use.names = FALSE)
                    col3 <- NULL
                    nnames <- length(namesdf3)

                    for(i in 1:nnames){
                         col3 <- c(col3,rep(namesdf3[i],times=nrows))
                    }    #  i

                    betacoeffs <- as.data.frame(tibble::tibble(col1,col2,col3))
                    wmf2 <- paste(wmf,".wmf",sep="")    
#                                                                print("OK to plotdiag.ANOX2    60")
   
                    ##############################################
                    # Call for plot using support function above #
                    ##############################################

                    plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                        mtitle=maintitle,
                        stitle=subtitle,
                        horlabel="Subset size m",
                        vertlabel="p-values for test of proportionality",
                        cap=caption,
                        legendname=legend,
                        verbose=verbose)
               }   # if yle = "p"
               #
               if(yle=="v"){
#                                                                    print("OK to plotdiag.ANOX2    61")
                    cola <- df2[,1] 
                    col1 <- rep(cola, times=ncols)
                    namesdf3 <- names(df3)
                    col2 <- c(unlist(df3),use.names = FALSE)
                    col3 <- NULL
                    nnames <- length(namesdf3)

                    for(i in 1:nnames){
                         col3 <- c(col3,rep(namesdf3[i],times=nrows))
                    }    #  i

                    betacoeffs <- as.data.frame(tibble::tibble(col1,col2,col3))
                    wmf2 <- paste(wmf,".wmf",sep="")    
#                                                                print("OK to plotdiag.ANOX2    62")
   
                    ##############################################
                    # Call for plot using support function above #
                    ##############################################

                    plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                        mtitle=maintitle,
                        stitle=subtitle,
                        horlabel="Subset size m",
                        vertlabel="p-values for ANOVA test of null hypothesis",
                        cap=caption,
                        legendname=legend,
                        verbose=verbose)
               }   # if yle = "p"
          }     #   if   else
     return()
     }     # End of preparation function #




##################################################   Main function    ##############################################

     grouped <- FALSE       #   possible future development
     #############################################################
     # Extract each set of coefficients and form into data frame #
     #############################################################
     ylab.extend <- substr(ylab.extend,1,1)                        # first letter 
     if(!(ylab.extend=="p" | ylab.extend =="v" ))stop("ylab.extend is not recognized")

     if(ylab.extend=="p"){    
          df1 <- forn$"Proportionality Test"                  # df1 already has a column of m
          if(is.character(df1))stop("Proportionality not tested.")
     }
     if(ylab.extend=="v"){
          df1 <- forn$"ANOVA"                            # df1 already has a column of m
     }
     mcolumn <- df1[,1]
     ndf1 <- dim(df1)[2]
     df1 <- df1[,-1]
     if(ndf1==2){
          df1 <- matrix(df1,nrow=length(df1), ncol=1)
          df1 <- as.data.frame(df1)
     }
     namesdf1 <- names(df1)
     if(!is.null(anova.rows)){ namesdf1 <- names(df1)[anova.rows]  }
     df8 <- NULL
     if(is.null(anova.rows)){
          df8 <- df1
     } 
     else{
          ccj <- anova.rows
          lencc <- length(ccj)
          df8 <- data.frame(df1[,ccj[1]])     # initialize df8
          if(lencc > 1){
               for(j in 2:lencc){
                    ccj <- anova.rows[j]
                    df8 <- data.frame(df8, df1[,ccj])
               }
          }    # lencc > 1
     }  #  if else
#                                                                              print("OK to plotdiag.ANOX2    20")
     df8 <- cbind(mcolumn,df8)
     df8 <- as.data.frame(df8)
     names(df8)<- c("m",namesdf1)
     prepstuff(rightforn=df8, yle=ylab.extend)                     #    call to prepstuff with anova rows set
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.ANOX2", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
