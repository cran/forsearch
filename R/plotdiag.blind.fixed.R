#' @export
plotdiag.blind.fixed <-
function(forn, 
     coeff.codenums=NULL,
     maintitle="Put main title here", 
     subtitle="Put subtitle here", 
     caption="Put caption here",
     wmf="Put_stored_name_here", 
     Cairo=TRUE,
     printgraph=TRUE,
     legend="Dummy legend name",
     verbose=TRUE)
{
     #                          plotdiag.blind.fixed
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.blind.fixed", quote = FALSE)
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
     ############################
     # Preparation for plotting #
     # Remove initial zero rows #
     ############################
     prepstuff <- function(rightforn,gg){
          df2 <- rightforn
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
               col3 <- "1"                       # default for only 1 coefficient
               betacoeffs <- as.data.frame(tibble::tibble(col1,col2,col3))
               wmf2 <- paste(wmf,".wmf",sep="")    
               ##############################################
               # Call for plot using support function above #
               ##############################################
               plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Estimated beta coefficient, blinded",
                     cap=caption,
                     legendname=legend,
                     verbose=verbose)

          }
          else{
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
               ##############################################
               # Call for plot using support function above #
               ##############################################
               plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Estimated beta coefficients, blinded",
                     cap=caption,
                     legendname=legend,
                     verbose=verbose)

          }    # else for is null cols
     return()
     }
     # End of preparation function #

##################################################   Main function    ##############################################

     grouped <- FALSE       #   possible future development

     #############################################################
     # Extract each set of coefficients and form into data frame #
     # Note location of zeros                                    #
     #############################################################
     df1 <- forn$"Fixed parameter estimates"                    # df1 already has a column of m
     Observation <- df1[,1]
     #############################
     # Remove observation number #
     # Replace NAs with 0        #
     #############################
     dfab <- df1[,-1]
     dimdfab <-dim(dfab)
     dimdfab1 <- dimdfab[1]
     dimdfab2 <- dimdfab[2]

     index0 <- matrix(FALSE, nrow=dimdfab1, ncol=dimdfab2)
     index0[dfab==0] <- TRUE
     index0 <- apply(index0,1,all)
     for(nn in 1:dimdfab2){
          index <- is.na(dfab[,nn])
          dfab[index,nn]<- 0
     }

     rowmean <- apply(dfab,2,mean)
     rowmean <- matrix(rowmean,nrow=dimdfab1,ncol=dimdfab2,byrow=TRUE)
     df8 <- dfab - rowmean
     df8 <- data.frame(Observation, df8)
     df8 <- df8[!index0,]
     prepstuff(rightforn=df8)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.blind.fixed", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     
     }
     return()
}
