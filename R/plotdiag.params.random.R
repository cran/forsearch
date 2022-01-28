#' @export
plotdiag.params.random <-
function(forn, coeff.codenums, 
     maintitle="Put maintitle here", 
     subtitle="Put subtitle here", 
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
     plotB1 <- function(data, xcol, ycol, cov1col,
                     facetcol=NULL, facetdir, 
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
          Hmisc::prn(out)             # this line plots the graph

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
     prepstuff <- function(df2){
          ndf2 <- dim(df2)[1]
          for(ir in 1:ndf2){
               uu <- df2[1,-1] 
               if(sum(uu)==0){df2 <- df2[-1,]} else break
          }      #   ir
                                  if(diagnose) Hmisc::prn(df2)

          nrows <- dim(df2)[1]
          df3 <- df2[,-1]
          ncols <- dim(df3)[2]
          cola <- df2[,1] 
          col1 <- rep(cola, times=ncols)
          namesdf3 <- names(df3)
          col2 <- c(unlist(df3),use.names = FALSE)
          col3 <- NULL
          nnames <- length(namesdf3)
          for(i in 1:nnames){
               col3 <- c(col3,rep(namesdf3[i],nrows))
          }    #  i
          betacoeffs <- as.data.frame(tibble::tibble(col1,col2,col3))
                                  if(diagnose) Hmisc::prn(betacoeffs)

          cap2 <- "Result of forward search"
          wmf2 <- paste(wmf,".wmf",sep="")    
          #
          ##############################################
          # Call for plot using support function above #
          ##############################################
          plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Estimated beta coefficients",
                     cap=cap2,
                     legendname=legend)

     }
     # End of preparation function #




##################################################   Main function    ##############################################


     #############################################################
     # Extract each set of coefficients and form into data frame #
     #############################################################
     lencc <- length(coeff.codenums)
     if(lencc<2) stop("Must call for plot of at least 2 coefficients (coeff.codes)")

     df1 <- forn$"Random parameter estimates" 
     ccj <- coeff.codenums[1]
     df2 <- data.frame(df1[[ccj]])
     for(j in 2:lencc){
          ccj <- coeff.codenums[j]
          df2 <- cbind(df2, df1[[ccj]])
     }
     m <- 1:(dim(df2)[1])
     df2 <- cbind(m,df2)
     prepstuff(df2)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.params.random", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
