// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat maipow(arma::mat a,arma::mat b){
  int m=a.n_rows;
  int n=a.n_cols;
  arma::mat re(m,n);
  for(int i=0; i<m ;i++){
    for(int j=0;j<n;j++){
      re(i,j) = pow(a(i,j),b(i,j));
    }
  }
  return re;
}
//[[Rcpp::export]]
arma::vec rowprod(arma::mat a,arma::mat b){
  int J=a.n_rows;
  arma::vec re(J);
  for(int j=0;j<J;j++){
    re(j)=sum((a.row(j)%b.row(j)));
  }
  return re;
}
//[[Rcpp::export]]
double sigmod(double x){
  double re;
  re=1/(1+exp(-x));
  return re;
}
//[[Rcpp::export]]
arma::mat Sigmod(arma::mat x){
  arma::mat re;
  re=1/(1+exp(-x));
  return re;
}
//[[Rcpp::export]]
double LA1dima(double x){
  double re;
  if(x==0) re=(double)1/8;
  else re=(1/(1+exp(-x))-0.5)/(2*x);
  return re;
}
//[[Rcpp::export]]
arma::mat LAmatrix(arma::mat x){
  int N=x.n_rows;
  int M=x.n_cols;
  arma::mat re(N,M);
  for(int i=0;i<N;i++){
    for(int j=0;j<M;j++){
      if(x(i,j)==0)
        re(i,j)=(double)1/8;
      else
        re(i,j)=(sigmod(x(i,j))-0.5)/(2*x(i,j));
    }
    
  }
  return re;
}
// [[Rcpp::export]]
arma::mat H_allprofile_dina(arma::mat alfa,arma::mat Q){
  int J=Q.n_rows;
  int K=Q.n_cols;
  int L=alfa.n_rows;
  arma::mat QN(L,K);
  arma::mat h(L,J);
  arma::mat re;
  for(int j=0;j<J;j++){
    QN=ones(L,1)*Q.row(j);
    re=maipow(alfa,QN);
    re=cumprod(re,1);
    h.col(j)=re.col(K-1);
  }
  return h;
}
//[[Rcpp::export]]
arma::mat H_allprofile_dino(arma::mat alfa,arma::mat Q){
  int J=Q.n_rows;
  int K=Q.n_cols;
  int L=alfa.n_rows;
  arma::mat QN(L,K);
  arma::mat h(L,J);
  arma::mat re;
  for(int j=0;j<J;j++){
    QN=ones(L,1)*Q.row(j);
    re=maipow(1-alfa,QN);
    re=cumprod(re,1);
    h.col(j)=re.col(K-1);
  }
  h=1-h;
  return h; 
}
//[[Rcpp::export]]
arma::cube H_allprofile_lcdm(arma::mat alfa,arma::mat Q,arma::umat comb,arma::ivec leng){
  int J=Q.n_rows;
  int K=Q.n_cols;
  int L=alfa.n_rows;
  int D=comb.n_rows;
  arma::uvec index;
  arma::mat QD(J,D);
  arma::mat alfaD(L,D);
  arma::cube h(J,L,D);
  h.zeros(J,L,D);
  for(int k=0;k<K;k++){
    alfaD.col(k)=alfa.col(k);
    QD.col(k)=Q.col(k);
  }
  for(int d=K;d<D;d++){
    arma::uvec ind=comb.row(d).t();
    arma::uvec inde;
    inde=ind.head(leng(d))-1;
    arma::mat ad=alfa.cols(inde);
    arma::mat qd=Q.cols(inde);
    ad=cumprod(ad,1);
    qd=cumprod(qd,1);
    alfaD.col(d)=ad.col(leng(d)-1);
    QD.col(d)=qd.col(leng(d)-1);
  }
  for(int j=0;j<J;j++){
    arma::mat qj;
    qj=ones(L,1)*QD.row(j);
    h.row(j)=alfaD%qj;
  }
  return h; 
}

//[[Rcpp::export]]
arma::mat ke_1dim(arma::cube h,arma::mat lambda,arma::cube Vlambda){
  int L=h.n_cols;
  int J=lambda.n_rows;
  arma::mat kesainew(L,J);
  arma::mat Vlambdanewj;
  arma::vec lambdanewj;
  arma::mat hj;
  arma::mat kk;
  for(int j=0;j<J;j++){
    Vlambdanewj=Vlambda.row(j);
    lambdanewj=(lambda.row(j)).t();
    hj=h.row(j);
    for(int i=0;i<L;i++){
      kk=hj.row(i)*(Vlambdanewj+lambdanewj*lambdanewj.t())*hj.row(i).t();
      kesainew(i,j)=sqrt(kk(0,0));
    }
  }
  return kesainew;
}

// [[Rcpp::export]]
Rcpp::List update(arma::mat y,arma::cube h,arma::mat lambda,arma::cube Vlambda,arma::mat xi,arma::vec phi1,arma::vec phi2,arma::vec delta0,arma::vec lambda0,arma::mat Vlambda0){
  int N=y.n_rows;
  int J=y.n_cols;
  int L=xi.n_rows;
  arma::mat Vlambdanewj;
  arma::vec lambdanewj;
  arma::mat yy0;
  arma::vec yy1(J);
  arma::vec yy2(J);
  arma::vec yy3;
  arma::mat yy;
  arma::mat rho(N,L);
  arma::mat hl;
  arma::mat yy4;
  arma::vec delta(L);
  for(int l=0;l<L;l++){
    hl=h.col(l);
    for(int j=0;j<J;j++){
      Vlambdanewj=Vlambda.row(j);
      lambdanewj=(lambda.row(j)).t();
      yy0=0.5*lambda.row(j)*hl.row(j).t();
      yy1(j)=yy0(0,0);
      yy0=log(sigmod(xi(l,j)))-0.5*xi(l,j)-LA1dima(xi(l,j))*(hl.row(j)*(Vlambdanewj+lambdanewj*lambdanewj.t())*hl.row(j).t()-pow(xi(l,j),2));
      yy2(j)=yy0(0,0);
    }
    yy=ones(N,1)*yy1.t();
    yy=yy%y;
    yy4=cumsum(yy,1);
    yy3=yy4.col(J-1);
    rho.col(l)=yy3+sum(yy2)+phi1(l)-phi2(l); 
  }
  rho=exp(rho);
  yy0=cumsum(rho,1);
  yy0=yy0.col(L-1)*ones(1,L);
  rho=rho/yy0;
  //update delta
  for(int d=0;d<L;d++){
    delta[d]=sum(rho.col(d))+delta0[d];
  }
  //update lambda,xi
  arma::mat lambdanew;
  lambdanew.zeros(J,2);
  arma::cube Vlambdanew;
  Vlambdanew.zeros(J,2,2);
  arma::mat xinew(L,J);
  arma::mat kk;
  for(int j=0;j<J;j++){
    arma::mat Vlambdanewj=Vlambdanew.row(j);
    arma::vec lambdanewj=lambdanew.row(j).t();
    arma::mat hstar;
    hstar.zeros(2,2);
    arma::vec yh;
    yh.zeros(2);
    arma::vec lamh;
    lamh.zeros(2);
    arma::mat hj;
    hj=h.row(j);
    arma::mat hl(2,2);
    arma::mat vv;
    for(int l=0;l<L;l++){
      hl=LA1dima(xi(l,j))*hj.row(l).t()*hj.row(l);
      for(int i=0;i<N;i++){
        hstar=hstar+rho(i,l)*hl;
        yh=yh+rho(i,l)*y(i,j)*hj.row(l).t();
      }
      
    }
    
    Vlambdanewj=inv(Vlambda0)+2*hstar;
    Vlambdanewj=inv(Vlambdanewj);
    lambdanewj=Vlambdanewj*(0.5*yh+inv(Vlambda0)*lambda0);
    Vlambdanew.row(j)=Vlambdanewj;
    lambdanew.row(j)=lambdanewj.t();
    for(int d=0;d<L;d++){
      kk=hj.row(d)*(Vlambdanewj+lambdanewj*lambdanewj.t())*hj.row(d).t();
      xinew(d,j)=sqrt(kk(0,0));
    }
  }
  arma::mat ke;
  ke=log(Sigmod(xinew))+LAmatrix(xinew)%pow(xinew,2);
  ke=cumsum(ke,1);
  ke=ones(N*1)*ke.col(J-1).t();
  
  return List::create(Named("rho") = rho, Named("delta") = delta,
                      Named("Vlambda") = Vlambdanew, Named("lambda") = lambdanew,
                      Named("xi") = xinew);
  
}

// [[Rcpp::export]]
Rcpp::List update_lcdm(arma::mat y,arma::cube h,arma::mat lambda,arma::cube Vlambda,arma::mat xi,arma::vec phi1,arma::vec phi2,arma::vec delta0,arma::vec lambda0,arma::mat Vlambda0,arma::umat lambdaindd,arma::ivec lambdalen){
  int N=y.n_rows;
  int J=y.n_cols;
  int L=xi.n_rows;
  arma::mat Vlambdanewj;
  arma::vec lambdanewj;
  arma::mat yy0;
  arma::vec yy1(J);
  arma::vec yy2(J);
  arma::vec yy3;
  arma::mat yy;
  arma::mat rho(N,L);
  arma::mat hl;
  arma::mat yy4;
  arma::vec delta(L);
  for(int l=0;l<L;l++){
    
    for(int j=0;j<J;j++){
      hl=h.col(l);
      arma::uvec index=lambdaindd.row(j).t();
      index=index.head(lambdalen(j));
      arma::vec la(lambdalen(j));
      arma::mat va(lambdalen(j),lambdalen(j));
      hl=hl.cols(index);
      Vlambdanewj=Vlambda.row(j);
      lambdanewj=(lambda.row(j)).t();
      yy0=0.5*lambdanewj(index).t()*hl.row(j).t();
      yy1(j)=yy0(0,0);
      yy0=log(sigmod(xi(l,j)))-0.5*xi(l,j)-LA1dima(xi(l,j))*(hl.row(j)*(Vlambdanewj(index,index)+lambdanewj(index)*lambdanewj(index).t())*hl.row(j).t()-pow(xi(l,j),2));
      yy2(j)=yy0(0,0);
    }
    yy=ones(N,1)*yy1.t();
    yy=yy%y;
    yy4=cumsum(yy,1);
    yy3=yy4.col(J-1);
    rho.col(l)=yy3+sum(yy2)+phi1(l)-phi2(l); 
  }
  rho=exp(rho);
  yy0=cumsum(rho,1);
  yy0=yy0.col(L-1)*ones(1,L);
  rho=rho/yy0;
  //update delta
  for(int d=0;d<L;d++){
    delta[d]=sum(rho.col(d))+delta0[d];
  }
  //update lambda,xi
  arma::mat lambdanew;
  lambdanew.zeros(J,L);
  arma::cube Vlambdanew;
  Vlambdanew.zeros(J,L,L);
  arma::mat xinew(L,J);
  arma::mat kk;
  for(int j=0;j<J;j++){
    arma::mat Vlambdanewj=Vlambdanew.row(j);
    arma::vec lambdanewj=lambdanew.row(j).t();
    arma::uvec index=lambdaindd.row(j).t();
    index=index.head(lambdalen(j));
    arma::vec la(lambdalen(j));
    arma::mat va(lambdalen(j),lambdalen(j));
    
    arma::mat hstar;
    hstar.zeros(lambdalen(j),lambdalen(j));
    arma::vec yh;
    yh.zeros(lambdalen(j));
    arma::vec lamh;
    lamh.zeros(lambdalen(j));
    arma::mat hj;
    hj=h.row(j);
    hj=hj.cols(index);
    arma::mat hl(lambdalen(j),lambdalen(j));
    arma::mat vv;
    for(int l=0;l<L;l++){
      hl=LA1dima(xi(l,j))*hj.row(l).t()*hj.row(l);
      for(int i=0;i<N;i++){
        hstar=hstar+rho(i,l)*hl;
        yh=yh+rho(i,l)*y(i,j)*hj.row(l).t();
      }
      
    }
    
    Vlambdanewj(index,index)=inv(Vlambda0(index,index))+2*hstar;
    Vlambdanewj(index,index)=inv(Vlambdanewj(index,index));
    lambdanewj(index)=Vlambdanewj(index,index)*(0.5*yh+inv(Vlambda0(index,index))*lambda0(index));
    Vlambdanew.row(j)=Vlambdanewj;
    lambdanew.row(j)=lambdanewj.t();
    for(int d=0;d<L;d++){
      kk=hj.row(d)*(Vlambdanewj(index,index)+lambdanewj(index)*lambdanewj(index).t())*hj.row(d).t();
      xinew(d,j)=sqrt(kk(0,0));
    }
  }
  arma::mat ke;
  ke=log(Sigmod(xinew))+LAmatrix(xinew)%pow(xinew,2)-0.5*xinew;
  ke=cumsum(ke,1);
  ke=ones(N*1)*ke.col(J-1).t();
  
  return List::create(Named("rho") = rho, Named("delta") = delta,
                      Named("Vlambda") = Vlambdanew, Named("lambda") = lambdanew,
                      Named("xi") = xinew);
  
}

// [[Rcpp::export]]
double ELBO_1dim(arma::mat y,arma::cube Vlambda,arma::mat xi,arma::mat rho,arma::mat Vlambda0,arma::vec delta,arma::vec delta0,arma::vec Mlambda0,arma::vec VMlambda0,arma::vec Mlambda0t,arma::vec VMlambda0t,double ph0,double pht,double df){
  int N=y.n_rows;
  int J=y.n_cols;
  arma::mat ke;
  double re;
  
  arma::vec Vj(J);
  for(int j=0;j<J;j++){
    arma::mat Vlambdanewj=Vlambda.row(j);
    Vj(j)=log(det(Vlambdanewj))-log(det(Vlambda0));
  }
  
  ke=log(Sigmod(xi))+LAmatrix(xi)%pow(xi,2)-0.5*xi;
  ke=cumsum(ke,1);
  ke=ones(N*1)*ke.col(J-1).t();
  re=sum(sum(rho%(ke-log(rho))))+0.5*sum(Vj)+sum(lgamma(delta)-lgamma(delta0))+lgamma(sum(delta0))-lgamma(sum(delta));
  re=re+0.5*sum(log(VMlambda0t)-log(VMlambda0))+log(pht)-log(ph0)+
    0.5*(VMlambda0t(0)*(1/VMlambda0t(0)-1/VMlambda0(0))-pow(Mlambda0t(0)-Mlambda0(0),2)/VMlambda0(0)+
    (VMlambda0t(1)-Mlambda0t(1)*VMlambda0t(1)*df/pht)*(1/VMlambda0t(1)-1/VMlambda0(1))-pow(Mlambda0t(1)-Mlambda0(1),2)/VMlambda0(1));
  return re;
}

// [[Rcpp::export]]
double ELBO_lcdm(arma::mat y,arma::cube Vlambda,arma::mat xi,arma::mat rho,arma::mat Vlambda0,arma::vec delta,arma::vec delta0,arma::umat lambdaindd,arma::ivec lambdalen,arma::vec Mlambda0,arma::vec VMlambda0,arma::vec Mlambda0t,arma::vec VMlambda0t,arma::vec ph0,arma::vec pht,arma::vec df){
  int N=y.n_rows;
  int J=y.n_cols;
  arma::mat ke;
  int L=Vlambda0.n_cols;
  double re;
  
  arma::vec Vj(J);
  for(int j=0;j<J;j++){
    arma::mat Vlambdanewj=Vlambda.row(j);
    arma::uvec index=lambdaindd.row(j).t();
    index=index.head(lambdalen(j));
    Vj(j)=log(det(Vlambdanewj(index,index)))-log(det(Vlambda0(index,index)));
  }
  
  ke=log(Sigmod(xi))+LAmatrix(xi)%pow(xi,2)-0.5*xi;
  ke=cumsum(ke,1);
  ke=ones(N*1)*ke.col(J-1).t();
  re=sum(sum(rho%(ke-log(rho))))+0.5*sum(Vj)+sum(lgamma(delta)-lgamma(delta0))+lgamma(sum(delta0))-lgamma(sum(delta));
  re=re+0.5*sum(log(VMlambda0t)-log(VMlambda0))+sum(log(pht)-log(ph0))+
    0.5*(VMlambda0t(0)*(1/VMlambda0t(0)-1/VMlambda0(0))-pow(Mlambda0t(0)-Mlambda0(0),2)/VMlambda0(0));
  arma::vec re1(L-1);
  for(int d=1;d<(L);d++){
    re1(d-1)=0.5*((VMlambda0t(d)-Mlambda0t(d)*VMlambda0t(d)*df(d-1)/pht(d-1))*(1/VMlambda0t(d)-1/VMlambda0(d))-pow(Mlambda0t(d)-Mlambda0(d),2)/VMlambda0(d));
  }
  re=re+sum(re1);
  return re;
}

