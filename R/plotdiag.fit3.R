#' @export
plotdiag.fit3 <-
function(forn, 
     maintitle="Put main title here", 
     subtitle="Put subtitle here", 
     caption="Put caption here",
     wmf="Put_stored_name_here", 
     Cairo=TRUE,
     printgraph=TRUE,
     legend="Dummy legend name",
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.fit3
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.fit3", quote = FALSE)
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
                     diagnose,verbose)
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
                 rescale = c("R", "fit3", "fixed"), bg = "transparent", canvas = "white", gamma = getOption("gamma"),
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
     }      #  plotB1
     #
     ############################
     # Preparation for plotting #
     ############################
     prepstuff <- function(rightforn){
          df2 <- rightforn
          wmf2 <- paste(wmf,".wmf",sep="")    
          ##############################################
          # Call for plot using support function above #
          ##############################################
          plotB1(data=df2, xcol=1, ycol=3, cov1col=2,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Statistic value",
                     cap=caption,
                     legendname=legend,
                     diagnose=diagnose,verbose=verbose)

     }
     # End of preparation function #




##################################################   Main function    ##############################################
     #############################################################
     # Extract each set of coefficients and form into data frame #
     #############################################################
     df1 <- forn$"Fit statistics"                    # df1 already has a column of m
     mcolumn <- df1[,1]
     ndf1 <- dim(df1)[1]
     df1 <- df1[,-1]
     df2 <- c(as.matrix(df1))
     covs <- c(rep("AIC", times=ndf1), rep("BIC", times=ndf1), rep("logLik", times=ndf1))
     df2 <- data.frame(mcolumn,covs,df2)  
     index <- df2[,3] != 0
     df2 <- df2[index,]
     prepstuff(rightforn=df2)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.fit3", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
