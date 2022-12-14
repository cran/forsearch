#' @export
plotdiag.deviances <-
function(forn,  
     devtype,
     maintitle="Put main title here", 
     subtitle= "Put subtitle here", 
     caption="Put caption here",
     wmf="Put_plot_file_title_here", 
     Cairo=TRUE,
     printgraph = TRUE,
     addline=c("none","loess","straight"),
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.deviances
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.deviances", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     #################
     # Plot function #
     #################
     plotB2 <- function(data,
               mtitle2,
               stitle2,
               filewidth=5, fileheight=5,
               cap){

          XVAR <- data[,1]
          YVAR <- data[,2]
          dfplot <- as.data.frame(tibble::tibble(XVAR,YVAR))
                     if(diagnose)Hmisc::prn(dfplot)
          temp <- stats::lm(formula=YVAR ~ XVAR,data=dfplot)
          FIT <- temp$fitted.values
          out <- ggplot2::ggplot(data=dfplot,ggplot2::aes(XVAR,YVAR)) + ggplot2::geom_point()
          addline <- substr(addline,1,1)
          if(addline=="l"){
                 out <- out + ggplot2::geom_smooth()
          }
          if(addline=="s"){
                dfplot <- data.frame(dfplot,FIT)
                out <- out + ggplot2::geom_line(ggplot2::aes(XVAR,FIT)) 
          }
          #    Add titles, axis labels, and caption. 
               horlabel <- "Subset size m"
               vertlabel <- "Null deviance"
               if(devtype == "R") vertlabel <- "Residual deviance"
          out <- out + ggplot2::ggtitle(mtitle2, subtitle=stitle2) + ggplot2::xlab(horlabel) + ggplot2::ylab(vertlabel) + ggplot2::labs(caption=cap)
                 if(diagnose)Hmisc::prn(as.character(out))   
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
          df1 <- rightforn
          df2 <- df1
          ndf2 <- length(df2)
          column1 <- 1:ndf2
          for(ir in 1:ndf2){
               if(df2[1]==0){
                    df2 <- df2[-1]
                    column1 <- column1[-1]
               }
               else break
          }    #  ir
                        if(diagnose) Hmisc::prn(df2)
          nrows <- length(df1)
          s2 <- data.frame(column1,df2)
          s2 <- s2[-1*c(1,2),]
          wmf2 <- paste(wmf,".wmf",sep="")    
          if(grouped){
               wmf2 <- paste(wmf,"Subgroup ", gg, ".wmf", sep="")
          }
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

     grouped <- FALSE              #   possible future development

     ################################################################
     # Extract each subgroup for plotting if forn is a grouped list #
     ################################################################
     if(devtype != "R" & devtype != "N") stop("Please indicate devtype='R' or devtype='N'")
     if(devtype=="R")forn2 <-forn$"Residual deviance"
     else forn2 <- forn$"Null deviance"
 
     if(grouped){
          nnames_forn <- length(forn)
          for(gg in 1:nnames_forn){
               prepstuff(rightforn=forn[[gg]], gg)
          }         #    gg
     }              #    grouped
     else{
          prepstuff(rightforn=forn2, gg="")
     }              # not grouped
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.deviances", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
