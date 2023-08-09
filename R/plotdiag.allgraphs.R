#' @export
plotdiag.allgraphs <-
function (object, mt = " ", st = " ", cpt = " ", cc = NULL, ccrand = NULL, Cairo=TRUE) 
{
#                                  plotdiag.allgraphs
#
# VALUE            Complete run of all plotdiag.xxx functions 
#
# INPUT
#      object           Name of object file contaiining data for plotting
# 

     MC <- object$Call
     if(is.null(MC))stop("object is not the name of a forsearch object file")

     print("The Call from this database is")
     print(MC)
     analysis <- substr(MC,start=11, stop=13)[1]
     if(!any(analysis==c("lm","lme","glm","cph"))){
          Hmisc::prn(analysis)
          stop("Call does not indicate a recognized underlying analysis function")
     }
     if(analysis=="cph"){
#          plotdiag.AICX(object, maintitle = mt, subtitle = st, caption = cpt, wmf="AICX", Cairo=Cairo)
#          plotdiag.Cook(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Cook", Cairo=Cairo)  
#          plotdiag.deviance.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviance residuals", Cairo=Cairo)
#          plotdiag.deviances(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviances", Cairo=Cairo)
#          plotdiag.fit3(object, maintitle = mt, subtitle = st, caption = cpt, wmf="fit3", Cairo=Cairo)
          plotdiag.leverage(object, maintitle = mt, subtitle = st, caption = cpt, wmf="leverage", Cairo=Cairo)         
          plotdiag.params.fixed(object, coeff.codenums=cc, maintitle = mt, subtitle = st, caption = cpt, wmf="params fixed", Cairo=Cairo)
#          plotdiag.params.random(object, coeff.codenums=ccrand, maintitle = mt, subtitle = st, caption = cpt, wmf="params random"
#          plotdiag.phihatx(object, maintitle = mt, subtitle = st, caption = cpt, wmf="phihatx", Cairo=Cairo)
#          plotdiag.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="residuals", Cairo=Cairo)         
#          plotdiag.s2(object, maintitle = mt, subtitle = st, caption = cpt, wmf="s2", Cairo=Cairo)
#          plotdiag.tstats(object, maintitle = mt, subtitle = st, caption = cpt, wmf="tstats", Cairo=Cairo)  

          plotdiag.loglik(object, maintitle = mt, subtitle = st, caption = cpt, wmf="loglik", Cairo=Cairo)  
          plotdiag.Wald(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Wald", Cairo=Cairo)  
          plotdiag.lrt(object, maintitle = mt, subtitle = st, caption = cpt, wmf="lrt", Cairo=Cairo) 
 
          plotdiag.loglik(object, maintitle = mt, subtitle = st, caption = cpt, wmf="loglik", Cairo=Cairo)  


          Hmisc::prn(names(object))
          print(search.history(object))
     }

     if(analysis=="lm"){
#          plotdiag.AICX(object, maintitle = mt, subtitle = st, caption = cpt, wmf="AICX", Cairo=Cairo)
          plotdiag.Cook(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Cook", Cairo=Cairo)  
#          plotdiag.deviance.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviance residuals", Cairo=Cairo)
#          plotdiag.deviances(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviances", Cairo=Cairo)
#          plotdiag.fit3(object, maintitle = mt, subtitle = st, caption = cpt, wmf="fit3", Cairo=Cairo)
          plotdiag.leverage(object, maintitle = mt, subtitle = st, caption = cpt, wmf="leverage", Cairo=Cairo)         
          plotdiag.params.fixed(object, coeff.codenums=cc, maintitle = mt, subtitle = st, caption = cpt, wmf="params fixed", Cairo=Cairo)
#          plotdiag.params.random(object, coeff.codenums=ccrand, maintitle = mt, subtitle = st, caption = cpt, wmf="params random"
#          plotdiag.phihatx(object, maintitle = mt, subtitle = st, caption = cpt, wmf="phihatx", Cairo=Cairo)
          plotdiag.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="residuals", Cairo=Cairo)         
          plotdiag.s2(object, maintitle = mt, subtitle = st, caption = cpt, wmf="s2", Cairo=Cairo)
          plotdiag.tstats(object, maintitle = mt, subtitle = st, caption = cpt, wmf="tstats", Cairo=Cairo)  

          Hmisc::prn(names(object))
          print(search.history(object))
     }

     if(analysis=="lme"){
#          plotdiag.AICX(object, maintitle = mt, subtitle = st, caption = cpt, wmf="AICX", Cairo=Cairo)
          plotdiag.Cook(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Cook", Cairo=Cairo)  
#          plotdiag.deviance.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviance residuals", Cairo=Cairo)
#          plotdiag.deviances(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviances", Cairo=Cairo)
          plotdiag.fit3(object, maintitle = mt, subtitle = st, caption = cpt, wmf="fit3", Cairo=Cairo)
          plotdiag.leverage(object, maintitle = mt, subtitle = st, caption = cpt, wmf="leverage", Cairo=Cairo)          
          plotdiag.params.fixed(object, coeff.codenums=cc, maintitle = mt, subtitle = st, caption = cpt, wmf="params fixed", Cairo=Cairo)
          plotdiag.params.random(object, coeff.codenums=ccrand, maintitle = mt, subtitle = st, caption = cpt, wmf="params random", Cairo=Cairo)
#          plotdiag.phihatx(object, maintitle = mt, subtitle = st, caption = cpt, wmf="phihatx", Cairo=Cairo)
          plotdiag.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="residuals", Cairo=Cairo)         
#          plotdiag.s2(object, maintitle = mt, subtitle = st, caption = cpt, wmf="s2", Cairo=Cairo)
          plotdiag.tstats(object, maintitle = mt, subtitle = st, caption = cpt, wmf="tstats", Cairo=Cairo)  

          Hmisc::prn(names(object))
          print(search.history(object))
     }

     if(analysis=="glm"){
          plotdiag.AICX(object, maintitle = mt, subtitle = st, caption = cpt, wmf="AICX", Cairo=Cairo)
#          plotdiag.Cook(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Cook", Cairo=Cairo)                              Should be doing this?  
          plotdiag.deviance.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="Deviance residuals", Cairo=Cairo)
          plotdiag.deviances(object, devtype="R", maintitle = mt, subtitle = st, caption = cpt, wmf="Deviances type R", Cairo=Cairo)
          plotdiag.deviances(object, devtype="N", maintitle = mt, subtitle = st, caption = cpt, wmf="Deviances type N", Cairo=Cairo)
#          plotdiag.fit3(object, maintitle = mt, subtitle = st, caption = cpt, wmf="fit3", Cairo=Cairo)
          plotdiag.leverage(object, maintitle = mt, subtitle = st, caption = cpt, wmf="leverage", Cairo=Cairo)          
          plotdiag.params.fixed(object, coeff.codenums=cc, maintitle = mt, subtitle = st, caption = cpt, wmf="params fixed", Cairo=Cairo)
#          plotdiag.params.random(object, coeff.codenums=ccrand, maintitle = mt, subtitle = st, caption = cpt, wmf="params random", Cairo=Cairo)
#          plotdiag.residuals(object, maintitle = mt, subtitle = st, caption = cpt, wmf="residuals", Cairo=Cairo)         
#          plotdiag.s2(object, maintitle = mt, subtitle = st, caption = cpt, wmf="s2", Cairo=Cairo)
          plotdiag.tstats(object, maintitle = mt, subtitle = st, caption = cpt, wmf="tstats", Cairo=Cairo)  

          Hmisc::prn(names(object))
          print(search.history(object))

          plotdiag.phihatx(object, maintitle = mt, subtitle = st, caption = cpt, wmf="phihatx")
     }
     return()
}
