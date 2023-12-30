#' @export     
VBEMM <- function(data,K,Q,EMlambda0=2,Vlambda0=1,EMeta0=-1,delta0=1,ueta0=-2,ulambda0=0,veta0=10,vlambda0=10,
               Mlambdainitial=2,Vlambdainitial=1,Metainitial=-1,rhoinitial=runif(2^K),deltainitial=1,
               alphainitial=1,piinitial=runif(2^K),lambdaindex=0,T=5000,trunc_num=0,e0=0.0001,model,output=TRUE){
  if(model=="DINA"){
    r <- VI_DINA(data,K,Q,EMlambda0,Vlambda0,EMeta0,delta0,ueta0,ulambda0,veta0,vlambda0,Mlambdainitial,
                 Vlambdainitial,Metainitial,rhoinitial,deltainitial,alphainitial,piinitial,trunc_num,T,e0,output)
  }else if(model=="DINO"){
    r <- VI_DINO(data,K,Q,EMlambda0,Vlambda0,EMeta0,delta0,ueta0,ulambda0,veta0,vlambda0,Mlambdainitial,
                 Vlambdainitial,Metainitial,rhoinitial,deltainitial,alphainitial,piinitial,trunc_num,T,e0,output)
  }else{
    r <- VI_LCDM(data,K,Q,EMlambda0,Vlambda0,EMeta0,delta0,ueta0,ulambda0,veta0,vlambda0,Mlambdainitial,
                 Vlambdainitial,Metainitial,rhoinitial,deltainitial,alphainitial,piinitial,lambdaindex,
                 trunc_num,T,e0,output)
  }
  return(r)
}

