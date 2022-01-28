#' @export
plotdiag.residuals <-
function(forn, squared=FALSE,
     hilos=c(1,0),
     maintitle="Put main title here", 
     subtitle="Put subtitle here", 
     wmf = "Put_graph_title_here", 
     Cairo=TRUE,
     printgraph = TRUE,
     legend="Dummy legend name",
     subdiag=FALSE, subverb=FALSE,
     diagnose=FALSE, verbose=TRUE)
{
     #                          plotdiag.residuals
     #
     # VALUE      Plot of the diagnostic statistics resulting from a forward search of a database.  Shows the influence of each observation on the scaled residuals.
     #                 For databases with more than 6 independent variables (including intercept), must subset the parameters in the plot. This doesn't affect
     #                 fitting of the model; this has already been done in the forward search procedure that formed the input to this function.
     #                 Handles linear models and mixed effects (grouped data) models.  The same subset of independent variables will be evaluated in each subgroup.
     #
     # INPUT    forn         File (list) resulting from run of forsearch1( ) or myforsearch2( ), the latter for mixed effects models.
     #          squared      Logical. TRUE causes residuals to be squared.
     #          hilos        Vector with number of high and number of low responses to flag in graph 
     #          maintitle    Graph main title
     #          subtitle     Graph subtitle
     #          wmf          Graph title in storage space for each plot; omit ".wmf"; ".wmf" and subgroup appendix (if needed) will be added in function
     #          Cairo        TRUE causes use of Cairo graphics
     #          legend       Legend name, if needed
     #          subdiag      Logical. TRUE causes printing of diagnostic content of called subfunctions
     #          subverb      Logical. TRUE causes printing of subfunction ID before and after running.
     #
     #          diagnose     Logical. TRUE causes printing of diagnostic content
     #          verbose      Logical. TRUE causes printing of program ID before and after running.
     #
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running plotdiag.residuals", quote = FALSE)
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
     plotD1 <- function(data, df3, xcol, ycol, dim1, dim2,
          cov2col, SE=FALSE, loess=T, 
          facetcol=NULL, facetdir, 
          mtitle, stitle, cap, horlabel, vertlabel, 
          filewidth=5, fileheight=5, legendname, 
          highslows, subdiag, subverb)
     {
     XVAR <- data[,xcol]
     YVAR <- data[,ycol]
     if(squared) {YVAR <- YVAR^2}
     COV2 <- data[,cov2col]
     SD <- data$SD
     N <- data$N
     dfplot <- data.frame(COV2,XVAR,YVAR,SD,N)
                if(diagnose)Hmisc::prn(dfplot)
     if(!is.null(facetcol)){
         FACET <- data[,facetcol]
         dfplot <- data.frame(dfplot,FACET)
                    if(diagnose)Hmisc::prn(dfplot)
     }
                if(diagnose)Hmisc::prn(dfplot)
     if(SE){
            upper <- dfplot$YVAR + dfplot$SD/sqrt(dfplot$N)
            lower <- dfplot$YVAR - dfplot$SD/sqrt(dfplot$N)
      }
      else    {
            upper <- dfplot$YVAR + dfplot$SD
            lower <- dfplot$YVAR - dfplot$SD
      }
      dfplot <- data.frame(dfplot,upper,lower)

      out <- ggplot2::ggplot(data=dfplot, ggplot2::aes(x=XVAR, y=YVAR, group=COV2, color=COV2))
      out <- out + ggplot2::geom_line() + ggplot2::geom_errorbar(ggplot2::aes(x=XVAR,ymin=lower,ymax=upper),width=0.2)+ ggplot2::theme(legend.position="none")

      highs <- highslows[1]
      lows <- highslows[2]
      if(highs>0){        
           for(ihigh in 1:highs){
                out <- out + ggplot2::annotate("text", x=dim2+3,y=df3[ihigh,2],label=as.character(df3[ihigh,1]))
           }
      }
      dimdf3 <- dim(df3)[1]
      if(lows>0){        
          for(ilow in 1:lows){
                out <- out + ggplot2::annotate("text", x=dim2+3,y=df3[dimdf3+1-ilow,2],label=as.character(df3[dimdf3+1-ilow,1]))
          }
      }
      ############################
      # Produce graph in facets? #
      ############################
      if(!is.null(facetcol)){
          if(facetdir=="h") out <- out + ggplot2::facet_grid(.~FACET)
          if(facetdir=="v") out <- out + ggplot2::facet_grid(FACET~.)
      }
      ############################################
      #    Add titles, axis labels, and caption. #
      ############################################
      out <- out + ggplot2::ggtitle(mtitle,subtitle=stitle) + ggplot2::xlab(horlabel) + ggplot2::ylab(vertlabel) + ggplot2::labs(caption=cap,legend=legendname)
      out <- out + ggplot2::labs(color=legendname)
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

      Hmisc::prn(out)             # this line plots the graph

      if(printgraph){
           filename <- paste(wmf,".wmf",sep="")
           ggplot2::ggsave(filename,width=filewidth, height=fileheight)
           grDevices::dev.off()
      }    # printgraph
    #
    if(subverb) {
        print("", quote = F)
        print("Finished running plotD1", quote = F)
        print("", quote = F)
        print(date(), quote = F)
        print("", quote = F)
    }
}
# End plot function #

     #####################################
     # Preparation function for plotting #
     #####################################
     prepstuff <- function(rightforn,gg){
          df1 <- rightforn$"Standardized residuals"
          ndf1 <- dim(df1)[2]
          mm <- 1:ndf1
          for(ir in 1:ndf1){
               uu <- df1[,1] 
               if(sum(uu)==0){
                    df1 <- df1[,-1]
                    mm <- mm[-1]
               } else break
          }      #   ir
          if(squared){
               df1 <- df1*df1
          }
          dimdf1 <- dim(df1)
          dim1 <- dimdf1[1]
          dim2 <- dimdf1[2]
 
          df1order <- order(df1[,dim2], decreasing=TRUE)
          lastresids <- df1[,dim2]
          lastresids <- sort(lastresids,decreasing=TRUE)
          residsdf3 <- data.frame(df1order,lastresids)
          names(residsdf3)[1] <- "Observation"
          print("", quote = F)
          print("Highest residuals when all included", quote=FALSE)
          print(utils::head(residsdf3), quote=FALSE)     # for identifing outliers on graph
          print("", quote = F)
          print("Lowest residuals when all included", quote=FALSE)
          print(utils::tail(residsdf3), quote=FALSE)     # for identifing outliers on graph
          print("", quote = F)

          COV2 <- rep(1:dim1, each=dim2)
          XVAR <- rep(mm, times=dim1)
          resids <- c(t(df1))
          SD <- 0
          N <- 1
          df2 <- as.data.frame(tibble::tibble(COV2, XVAR, resids, SD, N))
                                     if(subdiag)Hmisc::prn(df2)
          cap2 <- ""
          wmf2 <- paste(wmf,".wmf",sep="")    
          if(grouped){
               cap2 <- paste("Subgroup",gg) 
               wmf2 <- paste(wmf," Subgroup ",gg,".wmf",sep="")    
          }
               if(squared){vl <- "Scaled squared residuals"}
               else{vl <- "Scaled residuals"}
               plotD1(data=df2, df3=residsdf3, xcol=2, ycol=3, 
                     cov2col=1,dim1=dim1, dim2=dim2,
                     mtitle=maintitle,
                     stitle=subtitle,
                     highslows=hilos,
                     horlabel="Subset size m",
                     vertlabel=vl,
                     cap=cap2,
                     legendname=legend,
                     subdiag=subdiag,subverb=subverb)
     }


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
          print("Finished running plotdiag.residuals", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
     }
}
