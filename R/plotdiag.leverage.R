#' @export
plotdiag.leverage <-
function(forn,  
     hilos=c(1,0),
     maintitle = "Put main title here", 
     subtitle = "Put subtitle here", 
     caption="Put caption here",
     wmf = "Put_graph_title_here", 
     Cairo=TRUE,
     printgraph = TRUE,
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.leverage
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.leverage", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
############################################## Main function begins after functions plotD1( ) and prepstuff( ) ########################
     #################
     # Plot function #
     #################
     plotD1 <- function(data, xcol, ycol, cov2col, df3, 
          highslows,
          mtitle, stitle, cap, 
          horlabel, vertlabel, 
          filewidth, fileheight,  
          diagnose, verbose)
     {

     XVAR <- data[,xcol]
     YVAR <- data[,ycol]
     COV2 <- data[,cov2col]
     SD <- data$SD
     N <- data$N
     dfplot <- data.frame(COV2,XVAR,YVAR,SD,N)
                if(diagnose)Hmisc::prn(dfplot)
     upper <- dfplot$YVAR + dfplot$SD
     lower <- dfplot$YVAR - dfplot$SD
     dfplot <- data.frame(dfplot,upper,lower)

      out <- ggplot2::ggplot(data=dfplot, ggplot2::aes(x=XVAR, y=YVAR, group=COV2, color=COV2))
      out <- out + ggplot2::geom_line() + ggplot2::geom_errorbar(ggplot2::aes(x=XVAR,ymin=lower,ymax=upper),width=0.2)+ ggplot2::theme(legend.position="none")

      highs <- highslows[1]
                                      if(diagnose) Hmisc::prn(df3)
      if(highs>0){        
           dim2 <- max(XVAR)
           for(ihigh in 1:highs){
                out <- out + ggplot2::annotate("text", x=dim2+3,y=df3[ihigh,3],label=as.character(df3[ihigh,2]))
           }
      }
      #
      ############################################
      #    Add titles, axis labels, and caption. #
      ############################################
      out <- out + ggplot2::ggtitle(mtitle,subtitle=stitle) + ggplot2::xlab(horlabel) + ggplot2::ylab(vertlabel) + ggplot2::labs(caption=cap)
                   if(diagnose)    Hmisc::prn(as.character(out))
      #
      #############################
      # Print and save the graph. #
      #############################
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
        print("Finished running plotD1", quote = F)
        print("", quote = F)
        print(date(), quote = F)
        print("", quote = F)
      }
     }
# End plot function #

     ############################
     # Preparation for plotting #
     ############################
     prepstuff <- function(rightforn,gg){
          df1 <- rightforn$Leverage
          SD <- 0
          N <- 1
          df2 <- data.frame(df1, SD, N)
                            if(diagnose) Hmisc::prn(df2)
                            if(diagnose) temphist <- search.history(rightforn)
          df2order <- df2[order(df2$m),]
          maxm <- max(df2[,1])
          df3 <- df2[df2[,1]==maxm,,]
          df3 <- df3[   order(df3[,3],decreasing=TRUE)   ,]

          print("", quote = F)
          print("Observation leverages in order of final value", quote=FALSE)
          print(utils::head(df3[,2:3],n=10L), quote=FALSE)     # for identifing outliers on graph
          print("", quote = F)

          wmf2 <- paste(wmf,".wmf",sep="")    
          if(grouped){
               wmf2 <- paste(wmf," Subgroup ",gg,".wmf",sep="")    
          }
               plotD1(data=df2, xcol=1, ycol=3, cov2col=2, df3=df3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     highslows=hilos,
                     horlabel="Subset size m",
                     vertlabel="Leverage",
                     cap=caption, filewidth=5,fileheight=5,
                     diagnose=diagnose,verbose=verbose)
     }
     # End of preparation function #



##################################################   Main function    ##############################################

     grouped <- FALSE      #  possible future development

     ################################################################
     # Extract each subgroup for plotting if forn is a grouped list #
     ################################################################
     if(grouped){
          nnames_forn <- length(forn)
          for(gg in 1:nnames_forn){
               prepstuff(rightforn=forn[[gg]], gg)
          }         #    gg
     }              #    grouped
     else{
          prepstuff(rightforn=forn, gg="")
     }              # not grouped
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.leverage", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
