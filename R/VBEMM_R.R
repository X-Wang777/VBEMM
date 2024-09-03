#' @export     
VBEMM <- function(data,K,Q,EMlambda0=2,Vlambda0=1,EMeta0=-1,delta0=1,ueta0=-2,ulambda0=0,veta0=10,
                  vlambda0=10,Mlambdainitial=2,Vlambdainitial=1,Metainitial=-1,rhoinitial=runif(2^K),deltainitial=1,
                  alphainitial=1,piinitial=runif(2^K),lambdaindex=0,T=5000,trunc_num=0,e0=0.0001,model){
  if(model=="DINA"){
    time1 =as.POSIXlt(Sys.time())
    r <- VBEMM_DINA(data,K,Q,EMlambda0=EMlambda0,Vlambda0_value=Vlambda0, EMeta0=EMeta0,delta0_value=delta0,
               ueta0=ueta0,ulambda0=ulambda0,veta0=veta0,vlambda0=vlambda0,Mlambdainitial=Mlambdainitial,
               Vlambdainitial_value=Vlambdainitial,Metainitial=Metainitial,deltainitial=deltainitial,
               trunc_num=trunc_num,T=T,e0=e0) 
    time2 =as.POSIXlt(Sys.time())
    time = difftime(time2, time1, units="secs")
  }else if(model=="DINO"){
    time1 =as.POSIXlt(Sys.time())
    r <- VBEMM_DINO(data,K,Q,EMlambda0=EMlambda0,Vlambda0_value=Vlambda0, EMeta0=EMeta0,delta0_value=delta0,
                    ueta0=ueta0,ulambda0=ulambda0,veta0=veta0,vlambda0=vlambda0,Mlambdainitial=Mlambdainitial,
                    Vlambdainitial_value=Vlambdainitial,Metainitial=Metainitial,deltainitial=deltainitial,
                    trunc_num=trunc_num,T=T,e0=e0) 
    time2 =as.POSIXlt(Sys.time())
    time = difftime(time2, time1, units="secs")
    
  }else{
    time1 =as.POSIXlt(Sys.time())
    D=2^K-1
    comb <- matrix(0,D,K)##types of combnation
    comb[1:K,1]<-c(1:K)
    tt=1;comlen <- c()
    for(k in 1:K){
      index <- utils::combn(K,k)
      len <- dim(index)[2]
      comb[tt:(tt+len-1),1:k] <- t(index)
      comlen <- c(comlen,rep(k,len))
      tt=tt+len
    }
    if(sum(lambdaindex)==0){
     lambdaindex <- Q_saturatedLCDM(Q)
    }
    r <- VBEMM_LCDM(data,K,Q,lambdaindex,comb,comlen,EMlambda0_value=EMlambda0,Vlambda0_value=Vlambda0,
                    EMeta0=EMeta0,delta0_value=delta0,ueta0=ueta0,ulambda0=ulambda0,veta0=veta0,vlambda0=vlambda0,
                    Mlambdainitial=Mlambdainitial,Vlambdainitial_value=Vlambdainitial,Metainitial=Metainitial, 
                    deltainitial=deltainitial,T=T,e0=e0) 
    time2 =as.POSIXlt(Sys.time())
    time = difftime(time2, time1, units="secs")
  }
  r$time=time
  SD_lambda <- t(sqrt(apply(r$Vlambda,1,diag)))
  SD_pi <- sqrt(r$delta*(sum(r$delta)-r$delta)/(sum(r$delta)^2*(sum(r$delta+1))))
  r$SD=list(lambda_star=SD_lambda,pi=SD_pi)
  return(r)
}

