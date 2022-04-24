#' @export
plotdiag.params.random <-
function(forn, coeff.codenums=NULL,
     asfacets=FALSE, facetdir=c("h","v"), 
     maintitle="Put maintitle here", 
     subtitle="Put subtitle here", 
     caption="Put caption here",
     wmf="Put_stored_name_here", 
     Cairo=TRUE,
     printgraph=TRUE,
     legend="Dummy legend name",
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.params.random
     #
     # VALUE      Plot of the diagnostic statistics resulting from a forward search of a database.  Shows the influence of each observation on each independent variable coefficient.
     #                 For databases with more than 6 independent variables (including intercept), must subset the parameters in the plot.
     #                 Handles linear models and mixed effects (grouped data) models.  The same subset of independent variables will be evaluated in each subgroup.
     #
     # INPUT    forn            File resulting from run of forsearch_lme( ).
     #          coeff.codenums  Vector of code numbers, subset of those resulting from identifyCoeffs( ).
     #          maintitle       Graph main title
     #          subtitle        Graph subtitle
     #          caption         Graph caption
     #          wmf             Graph title in storage space for ungrouped plots; omit ".wmf"; ".wmf" and subgroup appendix (if needed) will be added in function
     #          Cairo           TRUE causes use of Cairo graphics
     #          printgraph      TRUE causes graph to be printed in a Windows metafile and closes the device
     #          legend          Legend name
     #
     #          diagnose        Logical. TRUE causes printing of diagnostic content
     #          verbose         Logical. TRUE causes printing of program ID before and after running.
     #
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.params.random", quote = FALSE)
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
     plotB1 <- function(data, xcol, ycol, cov1col=NULL, 
                     facetcol=NULL, facetdir=NULL, 
                     mtitle,
                     stitle,
                     horlabel,
                     vertlabel,
                     cap,
                     filewidth=5, fileheight=5,
                     legendname)
     {

          XVAR <- data[,xcol]
          YVAR <- data[,ycol]
          COV1 <- data[,cov1col]
          dfplot <- data.frame(XVAR,YVAR,COV1)
                      if(diagnose)Hmisc::prn(dfplot)
          FACET <- data[,facetcol]
          dfplot <- data.frame(dfplot,FACET)
                     if(diagnose)Hmisc::prn(dfplot)
          if(is.null(facetcol)){
               out <- ggplot2::ggplot(data=dfplot,ggplot2::aes(XVAR,YVAR,COV1)) + ggplot2::geom_point(ggplot2::aes(shape = COV1))
               out$labels$shape <- legendname
          }
          else{
               if(!(facetdir=="v" | facetdir=="h")) stop("facetdir must be 'v' or 'h'")
               ncells <- length(unique(FACET))
               out <- ggplot2::ggplot(data=dfplot, ggplot2::aes(XVAR,YVAR) ) + ggplot2::geom_point(ggplot2::aes(shape = FACET))
               if(facetdir=="h") out <- out + ggplot2::facet_wrap(~FACET,ncol=ncells,scales="free_y")
               if(facetdir=="v") out <- out + ggplot2::facet_wrap(~FACET,nrow=ncells,scales="free_y")
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
     }      #  plotB1
     #
     ############################
     # Preparation for plotting #
     # Remove initial zeros     #
     ############################
     prepstuff <- function(df2, usefacets){            #   df2 is a data frame
          nrows <- dim(df2)[1]
          ncols <- dim(df2)[2] - 1         # columns without observation number
          cola <- df2[,1] 
          if(ncols > 1){
               col1 <- rep(cola, times=ncols)

               df3 <- df2[,-1]
               namesdf3 <- names(df3)
               col2 <- c(unlist(df3), use.names = FALSE)
               col3 <- NULL
               nnames <- length(namesdf3)
               for(i in 1:nnames){
                    col3 <- c(col3,rep(namesdf3[i],nrows))
               }    #  i
          }     # ncols > 1
          else{
               col1 <- cola
               col2 <- c(df2[,2])
               col3 <- rep(names(df2)[2],nrows)
               usefacets <- FALSE
          }
          betacoeffs <- as.data.frame(tibble::tibble(col1,col2,col3))
                                  if(diagnose) Hmisc::prn(betacoeffs)
          index <- betacoeffs[,2] > 0
          betacoeffs <- betacoeffs[index,]
          wmf2 <- paste(wmf,".wmf",sep="")    
          #
          ##############################################
          # Call for plot using support function above #
          ##############################################
          if(usefacets){
               plotB1(data=betacoeffs, xcol=1, ycol=2, facetcol=3, facetdir=facetdir,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Root mean square of random coefficients",
                     cap=caption,
                     legendname=legend)

          }
          else
               plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Root mean square of random coefficients",
                     cap=caption,
                     legendname=legend)

     }
     # End of preparation function #




##################################################   Main function    ##############################################


     #############################################################
     # Extract each set of coefficients and form into data frame #
     #############################################################
     df2 <- forn$"Random parameter estimates" 
     if(!is.null(coeff.codenums)){
          df2 <- df2[,coeff.codenums]
     }
     m <- 1:(dim(df2)[1])                      #   add observation numbers
     df2 <- cbind(m,df2)
     prepstuff(df2,usefacets=asfacets)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.params.random", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
