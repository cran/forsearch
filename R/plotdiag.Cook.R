#' @export
plotdiag.Cook <-
function(forn,  
     maintitle="Put main title here", 
     subtitle= "Put subtitle here", 
     wmf="Put_plot_file_title_here", 
     Cairo=TRUE,
     printgraph = TRUE,
     loess=FALSE,
     subdiag=FALSE, subverb=FALSE,
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.Cook
     #
     # VALUE      Plot of the diagnostic statistics resulting from a forward search of a database.  Shows the influence of each observation on the estmate of variance.
     #                 Handles linear models and mixed effects (grouped data) models.  The same subset of independent variables will be evaluated in each subgroup.
     #
     # INPUT    forn         File resulting from run of forsearch_lm( ) or forsearch_lme( ), the latter for mixed effects models.
     #          maintitle    Graph main title
     #          subtitle     Graph subtitle
     #          wmf          Graph title in storage space for ungrouped plots; omit ".wmf"; ".wmf" and subgroup appendix (if needed) will be added in function
     #          Cairo        TRUE causes use of Cairo graphics
     #          printgraph   TRUE causes graph to be printed in a Windows metafile and closes the device
     #          loess        Logical T calls for loess fit to points
     #          subdiag      Logical. TRUE causes printing of diagnostic content of called subfunctions
     #          subverb      Logical. TRUE causes printing of subfunction ID before and after running.
     #
     #          diagnose     Logical. TRUE causes printing of diagnostic content
     #          verbose      Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.Cook", quote = FALSE)
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
                     if(subdiag)Hmisc::prn(dfplot)
          temp <- stats::lm(formula=YVAR ~ XVAR,data=dfplot)
          FIT <- temp$fitted.values
          out <- ggplot2::ggplot(data=dfplot,ggplot2::aes(XVAR,YVAR)) + ggplot2::geom_point()
          if(loess){
                 out <- out + ggplot2::geom_smooth()
          }
          else{
                dfplot <- data.frame(dfplot,FIT)
                out <- out + ggplot2::geom_line(ggplot2::aes(XVAR,FIT)) 
          }
          #    Add titles, axis labels, and caption. 
               horlabel <- "Subset size m"
               vertlabel <- "Modified Cook distance"
          out <- out + ggplot2::ggtitle(mtitle2, subtitle=stitle2) + ggplot2::xlab(horlabel) + ggplot2::ylab(vertlabel) + ggplot2::labs(caption=cap)
                 if(subdiag)Hmisc::prn(as.character(out))   
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
          #
          if(subverb) {
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
          df1 <- rightforn$"Modified Cook distance"
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
          cap2 <- "Result of forward search"
          wmf2 <- paste(wmf,".wmf",sep="")    
          if(grouped){
               cap2 <- 0
               wmf2 <- paste(wmf,"Subgroup ", gg, ".wmf", sep="")
          }
          ##############################################
          # Call for plot using support function above #
          ##############################################
          plotB2(data=s2,
               mtitle2=maintitle,
               stitle2=subtitle,
               cap=cap2)
     }

     # End of preparation function #


##################################################   Main function    ##############################################

     grouped <- FALSE              #   possible future development

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
          print("Finished running plotdiag.Cook", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
