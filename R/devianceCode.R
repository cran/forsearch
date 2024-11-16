devianceCode <-
function (obs, pred, ni=NULL, fam) 
{
          if(is.na(obs)) {
print("obs was NA")
               out <- -88  }
          else{
               if(fam=="binomial"){
                    if(obs==0){ out <- ni*log(ni/(ni-pred))     }      # expects whole numbers, not proportions
                    else if(obs==ni){ out <- obs*log(obs/pred) }
                    else{ out <- obs*log(obs/pred) + (ni-obs)*log((ni-obs)/(ni-pred)) }
               }
               if(fam=="Gamma"){ out <- -log(obs/pred) + (obs-pred)/pred  }

               if(fam=="poisson"){ if(obs==0){ out <- pred }
                    else{ out <- obs*log(obs/pred) - obs + pred  }
               }
 #              if(fam=="exponential"){ out <- -log(obs/pred) + obs*pred - 1 }

               if(is.na(out)){ out <- 0 }
               else{
                    if(abs(out)<= 10^(-12)){ out <- 0 }
                    else{ out <- 2 * out; vv <- obs - pred; vv <- vv/abs(vv)            #  1 or -1
                         out <- sqrt(out)*vv }
               }
          }
#prn(obs)
#prn(pred)
#prn(out)
#stop("dev fn")
          return(out)
}
