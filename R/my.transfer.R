my.transfer <-
function (app=FALSE) 
{

 dump(list="brMask",                      file="c:/temp/blindreviewHold/brMask.R", append=app) 
 dump(list="unmask",                      file="c:/temp/blindreviewHold/unmask.R", append=app) 

 dump(list="forsearch_lm",                file="c:/temp/forsearchHold/forsearch_lm.R", append=app) 
 dump(list="forsearch_lme",               file="c:/temp/forsearchHold/forsearch_lme.R", append=app)
 dump(list="forsearch_glm",               file="c:/temp/forsearchHold/forsearch_glm.R", append=app)   
 dump(list="forsearch_cph",               file="c:/temp/forsearchHold/forsearch_cph.R", append=app) 
 dump(list="forsearch_nls",               file="c:/temp/forsearchHold/forsearch_nls.R", append=app) 
# dump(list="forsearch_nme",               file="c:/temp/forsearchHold/forsearch_nme.R", append=app)

 dump(list="candprep",                    file="c:/temp/forsearchHold/candprep.R")                   
 dump(list="devianceCode",                file="c:/temp/forsearchHold/devianceCode.R")                   
 dump(list="aStep1",                      file="c:/temp/forsearchHold/aStep1.R")                   
 dump(list="aStep2",                      file="c:/temp/forsearchHold/aStep2.R")                   
 dump(list="bStep1",                      file="c:/temp/forsearchHold/bStep1.R")                   
 dump(list="bStep2",                      file="c:/temp/forsearchHold/bStep2.R")                   
 dump(list="cStep1",                      file="c:/temp/forsearchHold/cStep1.R")                   
 dump(list="cStep2",                      file="c:/temp/forsearchHold/cStep2.R") 
 dump(list="dStep1",                      file="c:/temp/forsearchHold/dStep1.R") 
 dump(list="dStep2",                      file="c:/temp/forsearchHold/dStep2.R") 
 dump(list="eStep1",                      file="c:/temp/forsearchHold/eStep1.R") 
 dump(list="eStep2",                      file="c:/temp/forsearchHold/eStep2.R") 
# dump(list="fStep1",                      file="c:/temp/forsearchHold/fStep1.R") 
# dump(list="fStep2",                      file="c:/temp/forsearchHold/fStep2.R") 

 dump(list="plotdiag.AICX",               file="c:/temp/forsearchHold/plotdiag.AICX.R", append=app)
 dump(list="plotdiag.ANOX2",              file="c:/temp/forsearchHold/plotdiag.ANOX2.R", append=app)
 dump(list="plotdiag.blind.fixed",        file="c:/temp/forsearchHold/plotdiag.blind.fixed.R", append=app)
 dump(list="plotdiag.Cook",               file="c:/temp/forsearchHold/plotdiag.Cook.R", append=app)
 dump(list="plotdiag.deviance.residuals", file="c:/temp/forsearchHold/plotdiag.deviance.residuals.R", append=app)
 dump(list="plotdiag.deviances",          file="c:/temp/forsearchHold/plotdiag.deviances.R", append=app)
 dump(list="plotdiag.fit3",               file="c:/temp/forsearchHold/plotdiag.fit3.R", append=app)
 dump(list="plotdiag.leverage",           file="c:/temp/forsearchHold/plotdiag.leverage.R", append=app)
 dump(list="plotdiag.loglik",             file="c:/temp/forsearchHold/plotdiag.loglik.R", append=app)
 dump(list="plotdiag.lrt",                file="c:/temp/forsearchHold/plotdiag.lrt.R", append=app)
 dump(list="plotdiag.params.fixed",       file="c:/temp/forsearchHold/plotdiag.params.fixed.R", append=app)
 dump(list="plotdiag.params.random",      file="c:/temp/forsearchHold/plotdiag.params.random.R", append=app)
 dump(list="plotdiag.phihatx",            file="c:/temp/forsearchHold/plotdiag.phihatx.R", append=app)
 dump(list="plotdiag.residuals",          file="c:/temp/forsearchHold/plotdiag.residuals.R", append=app)
 dump(list="plotdiag.s2",                 file="c:/temp/forsearchHold/plotdiag.s2.R", append=app)
 dump(list="plotdiag.tstats",             file="c:/temp/forsearchHold/plotdiag.tstats.R", append=app)                   
 dump(list="plotdiag.Wald",               file="c:/temp/forsearchHold/plotdiag.Wald.R", append=app)

 dump(list="plotdiag.allgraphs",          file="c:/temp/forsearchHold/plotdiag.allgraphs.R", append=app)
 
 dump(list="showme",                      file="c:/temp/forsearchHold/showme.R", append=app)       

 dump(list="identifyCoeffs",              file="c:/temp/forsearchHold/identifyCoeffs.R", append=app)
 dump(list="identifyFixedCoeffs",         file="c:/temp/forsearchHold/identifyFixedCoeffs.R", append=app)
 dump(list="logist3",                     file="c:/temp/forsearchHold/logist3.R", append=app)                   
 dump(list="search.history",              file="c:/temp/forsearchHold/search.history.R", append=app)                   

 dump(list="variablelist",                file="c:/temp/forsearchHold/variablelist.R")                   

 dump(list="my.transfer",                 file="c:/temp/forsearchHold/my.transfer.R")

# dump(list="Alfalfa.O.forlme",            file="c:/temp/forsearchHold/Alfalfa.O.forlme.R")
# dump(list="Machines.O.forlme2",          file="c:/temp/forsearchHold/Machines.O.forlme2.R")
# dump(list="crossdata.for1",              file="c:/temp/forsearchHold/crossdata.for1.R")                   
# dump(list="Alfalfa",                     file="c:/temp/forsearchHold/Alfalfa.R")
# dump(list="Alfalfa.O.forlme2",           file="c:/temp/forsearchHold/Alfalfa.O.forlme2.R")
# dump(list="Machines.O",                  file="c:/temp/forsearchHold/Machines.O.R")
# dump(list="crossdata",                   file="c:/temp/forsearchHold/crossdata.R")
# dump(list="Alfalfa.O.A",                 file="c:/temp/forsearchHold/Alfalfa.O.A.R")
# dump(list="micem1",                      file="c:/temp/forsearchHold/micem1.R")
# dump(list="micem1.for",                  file="c:/temp/forsearchHold/micem1.for.R")
# dump(list="train.for3",                  file="c:/temp/forsearchHold/train.for3.R")

# dump(list="paradigm.A.lm",               file="c:/temp/forsearchHold/paradigm.A.lm.R") 

}
