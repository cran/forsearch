bStep1 <-
function(fixed, yf, mnf, nOuter, yobs, s.o, nopl, nobs, i.s, fbg, b.d, verbose=TRUE)
{
#                                                             bStep1
#
# VALUE		Initial set of observation numbers (Step 1) for use in forsearch_lme
#
# INPUT
#
#         begin.diagnose	Numeric indicator of first diagnostic to print for this function. 0 causes all diagnostics to print.
#         verbose             Logical. TRUE causes print of function identifiers at start and end of function
#
# NOTE: Calls aStep1
#
     MC <- match.call()
     if(verbose) {
          print("", quote = FALSE)
          print("Running bStep1", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote = FALSE)
          print("", quote = FALSE)
          print("Call:", quote = FALSE)
          print(MC)
          print("", quote = FALSE)
     }
     yesfactor <- yf
     maxnfixdat <- mnf
     nufixdatOuter<- nOuter
     saved.obsnums <- s.o
     n.obs.per.level <- nopl
     initial.sample <- i.s
     fixdat.by.group <- fbg
     begin.diagnose <- b.d 
     spacer <- rep(" ", 20)
     #
     ####################################################################
     # Set up variables to hold partial results. This is augmenting set #
     ####################################################################
     rim.all.subgroups <- vector("list", nufixdatOuter)                        # holds vectors of final determination of each rim in each subgroup kk runs over columns within a row
     current.obs.by.group <- vector("list", nobs)                              # contains rim for each stage within each subgroup          dd runs down column, over rows
     for(kk in 1:nufixdatOuter){
          rim.all.subgroups[[kk]] <- current.obs.by.group                      # this only sets up a list within a list structure, not data
     }
     newrims <- rim.all.subgroups                                              # for translation back to original observation numbers
     #
     ####################################################
     # Replace initial.sample with inner.initial.sample #
     # Define mstart                                    #
     ####################################################
     inner.initial.sample <- max(100, round(initial.sample/nufixdatOuter))
     for(kk in 1:nufixdatOuter){
          fixdatkk <- fixdat.by.group[[kk]]
          lm.inner <- stats::lm(fixed, data=fixdatkk, singular.ok=TRUE)              # first level of rim.all.subgroups to be filled        lm
          rim <- aStep1(yesfactor, data=fixdatkk, inner.rank=lm.inner$rank, initial.sample=inner.initial.sample, 
                  formula=fixed, ycol=yobs, nopl=n.obs.per.level)                                                               #    aStep 1
          current.obs.by.group[[kk]] <- rim
     }     #   kk
     mstart <- lm.inner$rank
                                                              if(begin.diagnose <= 20){print(c(spacer,"Section 20"), quote=FALSE);Hmisc::prn(c(mstart,maxnfixdat))}
     #
     ###############################################
     # Install step 1 results in rim.all.subgroups #
     ###############################################
     compltd <- rep(0, nufixdatOuter)              # 0 = not completed, 1 = completed 
     for(tt in 1:nufixdatOuter){
          cobg <- current.obs.by.group[[tt]]
          rim.all.subgroups[[tt]][[mstart]] <- cobg
          if(length(cobg) == dim(fixdat.by.group[[tt]])[1]){
               compltd[tt] <- 1
          } 
     }             #   tt
                                                              if(begin.diagnose <= 25){print(c(spacer,"Section 25"), quote=FALSE);
                                                                Hmisc::prn(rim.all.subgroups[[1]][[mstart]]);Hmisc::prn(rim.all.subgroups[[2]][[mstart]])}
     # 

     out <- list(mstart=mstart, rim=rim.all.subgroups)
     if(verbose) {
          print("", quote = FALSE)
          print("Finished running bStep1", quote = FALSE)
          print("", quote = FALSE)
          print(date(), quote=FALSE)
          print("", quote = FALSE)
     }
     return(out)
}
