
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
