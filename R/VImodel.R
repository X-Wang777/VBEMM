VI_DINA <- function(data,K,Q,EMlambda0=2,Vlambda0=1,EMeta0=-1,delta0=1,uyi0=-2,ula0=0,vyi0=10,vla0=10,
                    Mlambdainitial=2,Vlambdainitial=1,Metainitial=-1,rhoinitial=runif(2^K),deltainitial=1,
                    alphainitial=1,piinitial=runif(2^K),trunc_num=0,T=5000,e0=0.0001,output=TRUE){
  time1 =as.POSIXlt(Sys.time())
  ######################
  N <- dim(data)[1]##num of subjects
  J <- dim(data)[2]##num of items
  L <- 2^K
  Vlambda0<-diag(2)*Vlambda0
  Vlambdainitial=diag(2)*Vlambdainitial
  delta0=rep(delta0,L)
  ystar <- 2*data-1###0,1 response data y transform to 1,-1 ystar
  alpha_all <- attribute_profile(K)##all attribute profile
  #####some definations
  alpha <- array(NA,c(T,N,K))
  eta <- matrix(0,N,J)
  pi <- array(NA,c(T,L))
  Vlambda <- array(0,c(T,J,2,2))
  Mlambda <- array(0,c(T,J))
  Meta <- array(0,c(T,J))
  xi <- array(0,c(T,L,J))
  Nalpha <- c()
  rho <- array(NA,c(T,N,L))
  delta <- array(NA,c(T,L))
  h <-array(NA,c(J,L,2))
  phi1 <- c()
  phi2 <- c()
  Mlambda0 <- c()
  Meta0 <- c()
  VMlambda0 <- c()
  VMeta0 <- c()
  #initial value
  for(j in 1:J){
    Vlambda[1,j,,] <- Vlambdainitial
  }
  Mlambda[1,] <- Mlambdainitial
  Meta[1,]<- Metainitial
  
  rhoinitial=rhoinitial/sum(rhoinitial)
  rho[1,,] <- matrix(1,N,1)%*%t(rhoinitial)
  delta[1,] <- deltainitial
  pi[1,] <- piinitial/sum(piinitial)
  alpha[1,,]=alphainitial
  t=1###loop index
  PH0 <- (1-pnorm(trunc_num,ula0,sqrt(vla0))) 
  
  H <- H_allprofile_dina(alpha_all,Q)
  for(j in 1:J){
    h[j,,] <- cbind(rep(1,L),H[,j]) 
  }
  xi[t,,]=ke_1dim(h,cbind(Meta[t,],Mlambda[t,]),Vlambda[t,,,])
  e=1
  elbo <-c()
  elbo[1]=0
  ####################
  #######main loop
  while((t<T) & (e >= e0)){
    t <- t+1
    ###1,update rho(pars for z_i)
    phi1 <- digamma(delta[t-1,])
    phi2 <- digamma(sum(delta[t-1,])) 
    result=update(ystar,h,cbind(Meta[t-1,],Mlambda[t-1,]),Vlambda[t-1,,,],xi[t-1,,],phi1,rep(phi2,L),delta0,c(EMeta0,EMlambda0),Vlambda0)
    rho[t,,] <- result$rho
    delta[t,] <- result$delta
    Vlambda[t,,,] <- result$Vlambda
    Meta[t,] <- result$lambda[,1]
    Mlambda[t,] <- result$lambda[,2]
    xi[t,,] <- result$xi
    ###update prior
    VMeta0[t] <- (J+1/(vyi0))^(-1)
    Meta0[t]= VMeta0[t]*(uyi0/vyi0+sum(Meta[t,]))
    EMeta0 <- Meta0[t]
    
    VMlambda0[t]<- (J+1/(vla0))^(-1)
    Mlambda0[t]= VMlambda0[t]*(ula0/vla0+sum(Mlambda[t,]))
    PHt <- (1-pnorm(trunc_num,Mlambda0[t],sqrt(VMlambda0[t])))
    df <- dnorm(trunc_num,Mlambda0[t],sqrt(VMlambda0[t]))
    EMlambda0 <- Mlambda0[t]+VMlambda0[t]*df/PHt
    
    elbo[t]=ELBO_1dim(ystar,Vlambda[t,,,],xi[t,,],rho[t,,],Vlambda0,delta[t,],delta0,c(uyi0, ula0),c(vyi0,vla0),c(Meta0[t], Mlambda0[t]),c(VMeta0[t],VMlambda0[t]),PH0,PHt,df)
    e=abs(elbo[t]-elbo[t-1])
    
    if(output){
      cat('\rIteration = ',t-1,sprintf(', ELBO change = %0.5f',e))
    }
  }
  time2 =as.POSIXlt(Sys.time())
  time=difftime(time2, time1, units="secs")
  pi_est=Est_pi(delta[t,])
  alpha_all=attribute_profile(K)
  alpha_est= Est_alpha(rho[t,,],alpha_all)
  
  para_est =list(eta=Meta[t,],lambda=Mlambda[t,],pi=pi_est,alpha=alpha_est)
  info=list(iter=t-1,timeused=time,change_diff=e,elbo=elbo[t])
  extra=list(Meta=Meta[t,],Mlambda=Mlambda[t,],Vlambda=Vlambda[t,,,],rho=rho[t,,],delta=delta[t,],Meta0=Meta0[t],VMeta0=VMeta0,Mlambda0=Mlambda0[t],VMlambda0=VMlambda0)
  whole_step <- list(rho=rho[1:t,,],delta=delta[1:t,],Mlambda=Mlambda[1:t,],Vlambda=Vlambda[1:t,,,],Meta=Meta[1:t,],Meta0=Meta0[1:t],Mlambda0=Mlambda0[1:t],elbo=elbo[1:t])
  r <- list(info=info,para_est=para_est,extra=extra,iter_data=whole_step)
  
  return(r)
}
VI_DINO <- function(data,K,Q,EMlambda0=2,Vlambda0=1,EMeta0=-1,delta0=1,uyi0=-2,ula0=0,vyi0=10,vla0=10,
                    Mlambdainitial=2,Vlambdainitial=1,Metainitial=-1,rhoinitial=runif(2^K),deltainitial=1,
                    alphainitial=1,piinitial=runif(2^K),trunc_num=0,T=5000,e0=0.0001,output=TRUE){
  time1 =as.POSIXlt(Sys.time())
  ######################
  N <- dim(data)[1]##num of subjects
  J <- dim(data)[2]##num of items
  L <- 2^K
  Vlambda0<-diag(2)*Vlambda0
  Vlambdainitial=diag(2)*Vlambdainitial
  delta0=rep(delta0,L)
  ystar <- 2*data-1###0,1 response data y transform to 1,-1 ystar
  alpha_all <- attribute_profile(K)##all attribute profile
  #####some definations
  alpha <- array(NA,c(T,N,K))
  eta <- matrix(0,N,J)
  pi <- array(NA,c(T,L))
  Vlambda <- array(0,c(T,J,2,2))
  Mlambda <- array(0,c(T,J))
  Meta <- array(0,c(T,J))
  xi <- array(0,c(T,L,J))
  Nalpha <- c()
  rho <- array(NA,c(T,N,L))
  delta <- array(NA,c(T,L))
  h <-array(NA,c(J,L,2))
  phi1 <- c()
  phi2 <- c()
  Mlambda0 <- c()
  Meta0 <- c()
  VMlambda0 <- c()
  VMeta0 <- c()
  #initial value
  for(j in 1:J){
    Vlambda[1,j,,] <- Vlambdainitial
  }
  Mlambda[1,] <- Mlambdainitial
  Meta[1,]<- Metainitial
  
  rhoinitial=rhoinitial/sum(rhoinitial)
  rho[1,,] <- matrix(1,N,1)%*%t(rhoinitial)
  delta[1,] <- deltainitial
  pi[1,] <- piinitial/sum(piinitial)
  alpha[1,,]=alphainitial
  t=1###loop index
  PH0 <- (1-pnorm(trunc_num,ula0,sqrt(vla0))) 
  
  H <- H_allprofile_dino(alpha_all,Q)
  for(j in 1:J){
    h[j,,] <- cbind(rep(1,L),H[,j]) 
  }
  xi[t,,]=ke_1dim(h,cbind(Meta[t,],Mlambda[t,]),Vlambda[t,,,])
  e=1
  elbo <-c()
  elbo[1]=0
  ####################
  #######main loop
  while((t<T) &(e>=e0)){
    t <- t+1
    ###1,update rho(pars for z_i)
    phi1 <- digamma(delta[t-1,])
    phi2 <- digamma(sum(delta[t-1,])) 
    result=update(ystar,h,cbind(Meta[t-1,],Mlambda[t-1,]),Vlambda[t-1,,,],xi[t-1,,],phi1,rep(phi2,L),delta0,c(EMeta0,EMlambda0),Vlambda0)
    rho[t,,] <- result$rho
    delta[t,] <- result$delta
    Vlambda[t,,,] <- result$Vlambda
    Meta[t,] <- result$lambda[,1]
    Mlambda[t,] <- result$lambda[,2]
    xi[t,,] <- result$xi
    ###update prior
    VMeta0[t] <- (J+1/(vyi0))^(-1)
    Meta0[t]= VMeta0[t]*(uyi0/vyi0+sum(Meta[t,]))
    EMeta0 <- Meta0[t]
    
    VMlambda0[t]<- (J+1/(vla0))^(-1)
    Mlambda0[t]= VMlambda0[t]*(ula0/vla0+sum(Mlambda[t,]))
    PHt <- (1-pnorm(trunc_num,Mlambda0[t],sqrt(VMlambda0[t])))
    df <- dnorm(trunc_num,Mlambda0[t],sqrt(VMlambda0[t]))
    EMlambda0 <- Mlambda0[t]+VMlambda0[t]*df/PHt
    ###update elbo
    elbo[t]=ELBO_1dim(ystar,Vlambda[t,,,],xi[t,,],rho[t,,],Vlambda0,delta[t,],delta0,c(uyi0, ula0),c(vyi0,vla0),c(Meta0[t], Mlambda0[t]),c(VMeta0[t],VMlambda0[t]),PH0,PHt,df)
    e=abs(elbo[t]-elbo[t-1])
    if(output){
      cat('\rIteration = ',t-1,sprintf(', ELBO change = %0.5f',e))
    }
  }
  time2 =as.POSIXlt(Sys.time())
  time=difftime(time2, time1, units="secs")
  pi_est=Est_pi(delta[t,])
  alpha_all=attribute_profile(K)
  alpha_est= Est_alpha(rho[t,,],alpha_all)
  
  para_est =list(eta=Meta[t,],lambda=Mlambda[t,],pi=pi_est,alpha=alpha_est)
  info=list(iter=t,timeused=time,change_diff=e,elbo=elbo[t])
  extra=list(Meta=Meta[t,],Mlambda=Mlambda[t,],Vlambda=Vlambda[t,,,],rho=rho[t,,],delta=delta[t,],Meta0=Meta0[t],VMeta0=VMeta0,Mlambda0=Mlambda0[t],VMlambda0=VMlambda0)
  whole_step <- list(rho=rho[1:t,,],delta=delta[1:t,],Mlambda=Mlambda[1:t,],Vlambda=Vlambda[1:t,,,],Meta=Meta[1:t,],Meta0=Meta0[1:t],Mlambda0=Mlambda0[1:t],elbo=elbo[1:t])
  r <- list(info=info,para_est=para_est,extra=extra,iter_data=whole_step)
  return(r)
}
VI_LCDM <- function(data,K,Q,EMlambda0=2,Vlambda0=1,EMeta0=-1,delta0=1,uyi0=-2,ula0=0,vyi0=10,vla0=10,
                    Mlambdainitial=2,Vlambdainitial=1,Metainitial=-1,rhoinitial=runif(2^K),deltainitial=1,
                    alphainitial=1,piinitial=runif(2^K),lambdaindex,trunc_num=0,T=5000,e0=0.0001,output=TRUE){
  time1 =as.POSIXlt(Sys.time())
  ######################
  N <- dim(data)[1]##num of subjects
  J <- dim(data)[2]##num of items
  L <- 2^K
  D <- L-1
  EMlambda0 <- rep(EMlambda0,D)
  Vlambda0=diag(L)*Vlambda0
  Mlambdainitial=rep(Mlambdainitial,D)
  Vlambdainitial=diag(L)*Vlambdainitial
  delta0=rep(delta0,L)
  comb <- matrix(0,D,K)##types of combnation
  comb[1:K,1]<-c(1:K)
  tt=1
  comlen <- c()
  for(k in 1:K){
    index <- utils::combn(K,k)
    len <- dim(index)[2]
    comb[tt:(tt+len-1),1:k] <- t(index)
    comlen <- c(comlen,rep(k,len))
    tt=tt+len
  }
  ystar <- 2*data-1###0,1 response data y transform to 1,-1 ystar
  alpha_all <- attribute_profile(K)##all attribute profile
  #####some definations
  alpha <- array(NA,c(T,N,K))
  eta <- matrix(0,N,J)
  pi <- array(NA,c(T,L))
  Vlambda <- array(0,c(T,J,L,L))
  Mlambda <- array(0,c(T,J,D))
  Meta <- array(0,c(T,J))
  xi <- array(0,c(T,L,J))
  Nalpha <- c()
  rho <- array(NA,c(T,N,L))
  delta <- array(NA,c(T,L))
  h <-array(NA,c(J,L,L))
  phi1 <- c()
  phi2 <- c()
  Meta0 <- c()
  #VMeta0 <- c()
  Mlambda0 <- array(NA,c(T,D))
  VMlambda0 <- c()
  PHt <- c()
  df <- c()
  lambdaindd <- matrix(0,J,D)
  lambdalen <- apply(lambdaindex,1,sum)
  for(j in 1:J){
    lambdaindd[j,1:lambdalen[j]] <- which(lambdaindex[j,]==1)
  }
  
  #initial value
  
  for(j in 1:J){
    indd <- lambdaindd[j,1:lambdalen[j]]
    Vlambda[1,j,,] <- Vlambdainitial
  }
  Mlambda[1,,] <- Mlambdainitial
  Meta[1,]<- Metainitial
  rhoinitial=rhoinitial/sum(rhoinitial)
  rho[1,,] <- matrix(1,N,1)%*%t(rhoinitial)
  delta[1,] <- deltainitial
  pi[1,] <- piinitial/sum(piinitial)
  alpha[1,,]=alphainitial
  t=1###loop index
  H <- H_allprofile_lcdm(alpha_all,Q,comb,comlen)
  for(j in 1:J){
    h[j,,] <- cbind(rep(1,L),H[j,,]) 
  }
  xi[t,,]=ke_1dim(h,cbind(Meta[t,],Mlambda[t,,]),Vlambda[t,,,])
  e=1
  elbo <-c()
  elbo[1] <- 0
  lengthd <- apply(lambdaindex,2,sum)
  for(d in 1:D){
    VMlambda0[d]<- (lengthd[d]+1/(vla0))^(-1)
  }
  VMeta0 <- (J+1/(vyi0))^(-1)
  PH0 <- (1-pnorm(trunc_num,ula0,sqrt(vla0))) 
  
  ####################
  #######main loop
  while((t<T) &(e>e0)){
    t <- t+1
    ###1,update rho(pars for z_i)
    phi1 <- digamma(delta[t-1,])
    phi2 <- digamma(sum(delta[t-1,])) 
    
    result=update_lcdm(ystar,h,cbind(Meta[t-1,],Mlambda[t-1,,]),Vlambda[t-1,,,],xi[t-1,,],phi1,rep(phi2,L),delta0,c(EMeta0,EMlambda0),Vlambda0,cbind(rep(0,J),lambdaindd),(lambdalen+1))
    rho[t,,] <- result$rho
    delta[t,] <- result$delta
    Vlambda[t,,,] <- result$Vlambda
    Meta[t,] <- result$lambda[,1]
    Mlambda[t,,] <- result$lambda[,2:L]
    xi[t,,] <- result$xi
    ####update prior
    Meta0[t]= VMeta0*(uyi0/vyi0+sum(Meta[t,]))
    EMeta0 <- Meta0[t]
    for(d in 1:D){
      inl <- which(Mlambda[t,,d]!=0)
      Mlambda0[t,d] <- VMlambda0[d]*(ula0/(vla0)+sum(Mlambda[t,inl,d]))
      PHt[d] <- (1-pnorm(trunc_num,Mlambda0[t,d],sqrt(VMlambda0[d])))
      df[d]  <- dnorm(trunc_num,Mlambda0[t,d],sqrt(VMlambda0[d]))
      EMlambda0[d] <- Mlambda0[t,d]+VMlambda0[d]*df[d]/PHt[d]
    }
    
    elbo[t]=ELBO_lcdm(ystar,Vlambda[t,,,],xi[t,,],rho[t,,],Vlambda0,delta[t,],delta0,cbind(rep(0,J),lambdaindd),(lambdalen+1),c(uyi0, rep(ula0,D)),c(vyi0,rep(vla0,D)),c(Meta0[t], Mlambda0[t,]),c(VMeta0,VMlambda0),rep(PH0,D),PHt,df)
    e=abs(elbo[t]-elbo[t-1])
    if(output){
      cat('\rIteration = ',t-1,sprintf(', ELBO change = %0.5f',e))
    }
  }
  time2 =as.POSIXlt(Sys.time())
  time=difftime(time2, time1, units="secs")
  
  pi_est=Est_pi(delta[t,])
  alpha_all=attribute_profile(K)
  alpha_est= Est_alpha(rho[t,,],alpha_all)
  
  para_est =list(eta=Meta[t,],lambda=Mlambda[t,,],pi=pi_est,alpha=alpha_est)
  info=list(iter=t,timeused=time,change_diff=e,elbo=elbo[t])
  extra=list(Meta=Meta[t,],Mlambda=Mlambda[t,,],Vlambda=Vlambda[t,,,],rho=rho[t,,],delta=delta[t,],Meta0=Meta0[t],VMeta0=VMeta0,Mlambda0=Mlambda0[t,],VMlambda0=VMlambda0)
  whole_step <- list(rho=rho[1:t,,],delta=delta[1:t,],lambda=Mlambda[1:t,,],Vlambda=Vlambda[1:t,,,],eta=Meta[1:t,],Meta0=Meta0[1:t],Mlambda0=Mlambda0[1:t,],elbo=elbo[1:t])
  r <- list(info=info,para_est=para_est,extra=extra,iter_data=whole_step)
  return(r)
}
