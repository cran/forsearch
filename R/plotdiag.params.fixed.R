#' @export
plotdiag.params.fixed <-
function(forn, coeff.codenums=NULL,
     maintitle="Put main title here", 
     subtitle="Put subtitle here", 
     wmf="Put_stored_name_here", 
     Cairo=TRUE,
     printgraph=TRUE,
     legend="Dummy legend name",
     subdiag=FALSE, subverb=FALSE,
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.params.fixed
     #
     # VALUE      Plot of the diagnostic statistics resulting from a forward search of a database.  Shows the influence of each observation on each independent variable coefficient.
     #                 For databases with more than 6 independent variables (including intercept), must subset the parameters in the plot.
     #                 Handles linear models and mixed effects (grouped data) models.  The same subset of independent variables will be evaluated in each subgroup.
     #
     # INPUT    forn            File resulting from run of forsearch_lm( ) or forsearch_lme( ), the latter for mixed effects models.
     #          coeff.codenums  Vector of index numbers (codes) of fixed effects to print, or NULL to print them all
     #          maintitle       Graph main title
     #          subtitle        Graph subtitle
     #          wmf             Graph title in storage space for ungrouped plots; omit ".wmf"; ".wmf" and subgroup appendix (if needed) will be added in function
     #          Cairo           TRUE causes use of Cairo graphics
     #          printgraph      TRUE causes graph to be printed in a Windows metafile and closes the device
     #          legend          Legend name
     #          subdiag         Logical. TRUE causes printing of diagnostic content of called subfunctions
     #          subverb         Logical. TRUE causes printing of subfunction ID before and after running.
     #
     #          diagnose        Logical. TRUE causes printing of diagnostic content
     #          verbose         Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.params.fixed", quote = FALSE)
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
               print("", quote = F)
               print("Finished running plotB1", quote = F)
               print("", quote = F)
               print(date(), quote = F)
               print("", quote = F)
          }    #  subverb
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
                                  if(diagnose) Hmisc::prn(betacoeffs)

               cap2 <- "Result of forward search"
               wmf2 <- paste(wmf,".wmf",sep="")    
               ##############################################
               # Call for plot using support function above #
               ##############################################
               plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Estimated beta coefficient",
                     cap=cap2,
                     legendname=legend,
                     diagnose=subdiag,verbose=subverb)

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
                                  if(diagnose) Hmisc::prn(betacoeffs)
               cap2 <- "Result of forward search"
               wmf2 <- paste(wmf,".wmf",sep="")    
               ##############################################
               # Call for plot using support function above #
               ##############################################

               plotB1(data=betacoeffs, xcol=1, ycol=2, cov1col=3,
                     mtitle=maintitle,
                     stitle=subtitle,
                     horlabel="Subset size m",
                     vertlabel="Estimated beta coefficients",
                     cap=cap2,
                     legendname=legend,
                     diagnose=subdiag,verbose=subverb)

          }    # else for is null cols
     }
     # End of preparation function #




##################################################   Main function    ##############################################

     grouped <- FALSE       #   possible future development

     #############################################################
     # Extract each set of coefficients and form into data frame #
     #############################################################
     df1 <- forn$"Fixed parameter estimates"                    # df1 already has a column of m
     mcolumn <- df1[,1]
     ndf1 <- dim(df1)[2]
     df1 <- df1[,-1]
     if(ndf1==2){
          df1 <- matrix(df1,nrow=length(df1), ncol=1)
          df1 <- as.data.frame(df1)
     }
     namesdf1 <- names(df1)
     if(!is.null(coeff.codenums)) namesdf1 <- names(df1)[coeff.codenums]
     df8 <- NULL
     if(is.null(coeff.codenums)){
          df8 <- df1
     } 
     else{
          ccj <- coeff.codenums
          lencc <- length(ccj)
          df8 <- data.frame(df1[,ccj[1]])     # initialize df8
          if(lencc > 1){
               for(j in 2:lencc){
                    ccj <- coeff.codenums[j]
                    df8 <- data.frame(df8, df1[,ccj])
               }
          }    # lencc > 1
     }  #  if else
     df8 <- cbind(mcolumn,df8)
     df8 <- as.data.frame(df8)
     names(df8)<- c("m",namesdf1)
     prepstuff(rightforn=df8)
     #
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running plotdiag.params.fixed", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
