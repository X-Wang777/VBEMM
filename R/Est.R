
Est_alpha <- function(post,alpha_all){
  N <- dim(post)[1]
  K <-dim(alpha_all)[2]
  alpha <- matrix(NA,N,K)
  for(i in 1:N){
    ind <- which(post[i,]==max(post[i,]))[1]
    alpha[i,]<- alpha_all[ind,]
  }
  return(alpha)
}
Est_pi <- function(delta){
  pi=delta/sum(delta)
  return(pi)
}
AAR <- function(alpha0,alpha){
  K <-dim(alpha)[2]
  N <- dim(alpha)[1]
  aar <- c()
  for(k in 1:K){
    aar[k] <- 1-sum(abs(alpha[,k]-alpha0[,k]))/N
  }
  return(aar)
}
PAR <- function(alpha0,alpha){
  N <- dim(alpha)[1]
  num <- 0
  for(i in 1:N){
    if(all(alpha[i,]==alpha0[i,])){num=num+1}
  }
  return(num/N)
}

PARn <- function(alpha0,alpha,n){
  N <- dim(alpha)[1]
  num <- 0
  for(i in 1:N){
    if(sum(abs(alpha[i,]-alpha0[i,]))<=n){num=num+1}
  }
  return(num/N)
}

