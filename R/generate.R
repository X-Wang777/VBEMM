attribute_profile <- function(K){
  alphaprofile <- array(0,c(2^K,K))
  in1 <- 1
  for(k in 1:K){
    index <- utils::combn(K,k)
    in2 <- dim(index)[2]
    for(i in 1:in2){
      alphaprofile[in1+i,index[,i]] <- 1
    }
    in1=in1+in2 
  }
  return(alphaprofile)
}

Qgenerate <- function(J,K,prob){
  Q1 <- diag(K)
  Q2 <- matrix(stats::rbinom((J-2*K)*K,1,prob),J-2*K,K)
  for(i in 1:(J-2*K)){
    if(sum(Q2[i,])==0){
      Q2[i,1] <- 1
    }
  }
  Q <- rbind(Q1,Q1,Q2)
  return(Q)
}
###find an attribute for one person
fun=function(x){r <- which(x==1);return(r)}

###generate attribute for N person
attribute_allsubject <- function(N,K,pai=rep(1/2^K,2^K)){
  allprofile <-attribute_profile(K) 
  index <- t(stats::rmultinom(N,1,pai))
  index <- apply(index,1,fun)
  alpha <- matrix(0,N,K)
  for(i in 1:N){
    alpha[i,] <- allprofile[index[i],]
  }
  return(alpha)
}##N*K matrix
h_onesubject_lcdm_R <- function(alpha,Q){
  J <- nrow(Q)
  K <- ncol(Q)
  alphamatrix <- matrix(rep(alpha,J),J,K,2)
  L <- 2^K-1
  QL <- matrix(NA,J,L)
  alphaL <- matrix(NA,J,L)
  t <- 1
  for(k in 1:K){
    QL[,t] <- Q[,k]
    alphaL[,t] <- alphamatrix[,k]
    t <- t+1
  }
  for(i in 2:K){
    ind <- utils::combn(K,i)
    num <- dim(ind)[2]
    for(l in 1:num){
      QL[,t] <- apply(Q[,ind[,l]],1,prod)
      alphaL[,t] <- apply(alphamatrix[,ind[,l]],1,prod)
      t <- t+1
    }
  }
  alphaL=alphaL*QL
  return(alphaL)
}

LCDMdatagenerate <- function(N,J,K,lambda,eta,Q,sigma=0){
    V <- matrix(sigma,K,K)
    diag(V) <- rep(1,K)
    alpha <- MASS::mvrnorm(N,rep(0,K),V)
    alpha[alpha<0]=0
    alpha[alpha>0]=1
  r <- array(NA,c(N,J))
  for(i in 1:N){
    h <- h_onesubject_lcdm_R(alpha[i,],Q)
    p <- 1/(1+exp(-rowprod(lambda,h)-eta))
    for(j in 1:J){
      r[i,j] <- stats::rbinom(1,1,p[j])
    }
  }
  re <- list(y=r,alpha=alpha,h=h)
  return(re)
}
eta_onesubject_dina <- function(alpha,Q){
  J <- nrow(Q)
  K <- ncol(Q)
  alphamatrix <- matrix(rep(alpha,J),J,K,2)^Q
  eta <- apply(alphamatrix,1,prod)
  return(eta)
}
###DINA model generate response data
##s and g are vectors(J*1)
DINAdatagenerate <- function(N,J,K,s,g,Q,sigma=0){
    V <- matrix(sigma,K,K)
    diag(V) <- rep(1,K)
    alpha <- MASS::mvrnorm(N,rep(0,K),V)
    alpha[alpha<0]=0
    alpha[alpha>0]=1
  r <- array(NA,c(N,J))
  for(i in 1:N){
    eta <- eta_onesubject_dina(alpha[i,],Q)
    p <- (1-s)^eta*g^(1-eta)
    for(j in 1:J){
      r[i,j] <- stats::rbinom(1,1,p[j])
    }
  }
  re <- list(y=r,alpha=alpha)
  return(re)
}
eta_onesubject_dino <- function(alpha,Q){
  J <- nrow(Q)
  K <- ncol(Q)
  alphamatrix <- (1-matrix(rep(alpha,J),J,K,2))^Q
  eta <- 1-apply(alphamatrix,1,prod)
  return(eta)
}
###DINO model generate response data
##s and g are vectors(J*1)
DINOdatagenerate <- function(N,J,K,s,g,Q,sigma=0){
    V <- matrix(sigma,K,K)
    diag(V) <- rep(1,K)
    alpha <- MASS::mvrnorm(N,rep(0,K),V)
    alpha[alpha<0]=0
    alpha[alpha>0]=1
  r <- array(NA,c(N,J))
  for(i in 1:N){
    eta <- eta_onesubject_dino(alpha[i,],Q)
    p <- (1-s)^eta*g^(1-eta)
    for(j in 1:J){
      r[i,j] <- stats::rbinom(1,1,p[j])
    }
    
  }
  re <- list(y=r,alpha=alpha)
  return(re)
}

Q_saturatedLCDM <- function(Q){
  J <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K-1
  QL <- matrix(NA,J,L)
  t <- 1
  for(k in 1:K){
    QL[,t] <- Q[,k]
    t <- t+1
  }
  for(i in 2:K){
    ind <- combn(K,i)
    num <- dim(ind)[2]
    for(l in 1:num){
      QL[,t] <- apply(Q[,ind[,l]],1,prod)
      t <- t+1
    }
  }
  return(QL)
}
logit <- function(a){
  r <- log(a/(1-a))
  return(r)
}
