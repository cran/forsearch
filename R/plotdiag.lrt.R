#' @export
plotdiag.lrt <-
function(forn,  
     maintitle= "Put main title here", 
     subtitle= "Put subtitle here" , 
     caption="Put caption here",
     wmf = "Put_graph_filename_here", 
     Cairo=TRUE,
     printgraph = TRUE,
     addline=c("none","loess","straight"),
     verbose=TRUE)
{
     #                          plotdiag.lrt
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.lrt", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     if(length(addline)>1)  addline <- "none"
     #################
     # Plot function #
     #################
     plotB2 <- function(data,
               mtitle2,
               stitle2,
               cap,
               filewidth=5, fileheight=5,
               wmfname){

          XVAR <- data[,1]
          YVAR <- data[,2]
          dfplot <- as.data.frame(tibble::tibble(XVAR,YVAR))
#                     if(diagnose)Hmisc::prn(dfplot)
          temp <- stats::lm(formula=YVAR ~ XVAR,data=dfplot)
          FIT <- temp$fitted.values
          out <- ggplot2::ggplot(data=dfplot,ggplot2::aes(XVAR,YVAR)) + ggplot2::geom_point()
          if(substring(addline,1)=="l"){
                 out <- out + ggplot2::geom_smooth()
          }
          if(substring(addline,1)=="s"){
                dfplot <- data.frame(dfplot,FIT)
                out <- out + ggplot2::geom_line(ggplot2::aes(XVAR,FIT)) 
          }
          #    Add titles, axis labels, and caption. 
               horlabel <- "Subset size m"
               vertlabel <- "Likelihood ratio Test"
          out <- out + ggplot2::ggtitle(mtitle2,subtitle=stitle2) + ggplot2::xlab(horlabel) + ggplot2::ylab(vertlabel) + ggplot2::labs(caption=cap)
#                 if(diagnose)Hmisc::prn(as.character(out))   
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
             print("", quote = FALSE)
             print("Finished running plotB1", quote = FALSE)
             print("", quote = FALSE)
             print(date(), quote = FALSE)
             print("", quote = FALSE)
          }

     }    #    End of plot function

     ############################
     # Preparation for plotting #
     ############################
     prepstuff <- function(rightforn,gg){
          s2 <- rightforn$"Likelihood ratio test"     # a data frame
          s2 <- s2[-1*c(1,2),]
          wmf2 <- paste(wmf,".wmf",sep="")    
          ##############################################
          # Call for plot using support function above #
          ##############################################
          plotB2(data=s2,
               mtitle2=maintitle,
               stitle2=subtitle,
               cap=caption)
     }
     # End of preparation function #

##################################################   Main function    ##############################################

     grouped <- FALSE       #  possible future development

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
          print("Finished running plotdiag.lrt", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
