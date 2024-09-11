// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
double MAX=1e300;
double MIN=1e-300;
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
arma::mat Trans_10to2_mat(int K,const arma::ivec& CL) {
  int Col=CL.n_elem;
  arma::mat alpha(Col,K);
  for(int i=0;i<Col;i++){
    arma::colvec alphai(K);double cl=CL(i);
    for(int k=0;k<K;k++){
      double twopow = pow(2,K-k-1);
      alphai(k) = (twopow<=cl);
      cl = cl - twopow*alphai(k);
    }
    alpha.row(i)=alphai.t();
  }
  return alpha;
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
// [[Rcpp::export]]
arma::vec H_allprofile_dina_oneitem(arma::mat alfa,arma::vec Q){
  int K=alfa.n_cols;
  int L=alfa.n_rows;
  arma::mat QN(L,K);
  arma::vec h(L);
  arma::mat re;
  QN=ones(L,1)*Q.t();
  re=maipow(alfa,QN);
  re=cumprod(re,1);
  h=re.col(K-1);
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
arma::mat H_allprofile_oneitem_lcdm(arma::mat alfa,arma::vec Q,arma::umat comb,arma::ivec leng){
  int K=alfa.n_cols;
  int L=alfa.n_rows;
  int D=comb.n_rows;
  arma::uvec index;
  arma::vec QD(D);
  arma::mat alfaD(L,D);
  arma::mat h(L,D);
  h.zeros(L,D);
  for(int k=0;k<K;k++){
    alfaD.col(k)=alfa.col(k);
    QD(k)=Q(k);
  }
  for(int d=K;d<D;d++){
    arma::uvec ind=comb.row(d).t();
    arma::uvec inde;
    inde=ind.head(leng(d))-1;
    arma::mat ad=alfa.cols(inde);
    arma::vec qd=Q(inde);
    ad=cumprod(ad,1);
    QD(d)=prod(qd);
    alfaD.col(d)=ad.col(leng(d)-1);
  }
  arma::mat qj;
  qj=ones(L,1)*QD.t();
  h=alfaD%qj;
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


void update(arma::mat y,arma::cube h,arma::mat& lambda,arma::cube& Vlambda,arma::mat& xi,arma::vec& phi1,arma::vec& phi2,
             arma::vec& delta,arma::vec& delta0,arma::vec& Elambda0t,arma::mat& Vlambda0,arma::mat& rho){
  int N=y.n_rows;int J=y.n_cols;int L=xi.n_rows;
  arma::mat Vlambdanewj;
  arma::vec lambdanewj;
  arma::mat yy0;
  arma::vec yy1(J);
  arma::vec yy2(J);
  arma::vec yy3;
  arma::mat yy;
  arma::mat hl;
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
    rho.col(l)=sum((ones(N,1)*yy1.t())%y,1)+sum(yy2)+phi1(l)-phi2(l); 
  }
  rho=exp(rho);
  yy0=sum(rho,1)*ones(1,L);
  rho=rho/yy0;
  //update delta
  for(int d=0;d<L;d++){
    delta(d)=sum(rho.col(d))+delta0(d);
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
    arma::vec yh;yh.zeros(2);
    arma::vec lamh;lamh.zeros(2);
    arma::mat hj;
    hj=h.row(j);
    arma::mat hl(2,2);
    arma::mat vv;
    for(int l=0;l<L;l++){
      hl=LA1dima(xi(l,j))*hj.row(l).t()*hj.row(l);
      hstar=hstar+sum(rho.col(l))*hl;
      yh=yh+sum(rho.col(l)%y.col(j))*hj.row(l).t();
    }
    
    Vlambdanewj=inv(Vlambda0)+2*hstar;
    Vlambdanewj=inv(Vlambdanewj);
    lambdanewj=Vlambdanewj*(0.5*yh+inv(Vlambda0)*Elambda0t);
    Vlambdanew.row(j)=Vlambdanewj;
    lambdanew.row(j)=lambdanewj.t();
    for(int d=0;d<L;d++){
      kk=hj.row(d)*(Vlambdanewj+lambdanewj*lambdanewj.t())*hj.row(d).t();
      xinew(d,j)=sqrt(kk(0,0));
    }
  }
  arma::mat ke;
  ke=log(Sigmod(xinew))+LAmatrix(xinew)%pow(xinew,2);
  ke=ones(N*1)*sum(ke,1).t();
  Vlambda = Vlambdanew;lambda = lambdanew;xi=xinew;
}

// [[Rcpp::export]]
void update_lcdm(arma::mat y,arma::cube h,arma::mat& lambda,arma::cube& Vlambda,arma::mat& xi,arma::vec& phi1,arma::vec& phi2,
                  arma::vec& delta0,arma::vec& delta,arma::mat& rho,arma::vec& lambda0,arma::mat& Vlambda0,arma::umat& lambdaindd,arma::ivec& lambdalen){
  int N=y.n_rows;int J=y.n_cols;int L=xi.n_rows;
  arma::mat Vlambdanewj;
  arma::vec lambdanewj;
  arma::mat yy0;
  arma::vec yy1(J);
  arma::vec yy2(J);
  arma::vec yy3;
  arma::mat yy;
  //arma::mat rho(N,L);
  arma::mat hl;
  arma::mat yy4;
  //arma::vec delta(L);
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
    rho.col(l)=sum((ones(N,1)*yy1.t())%y,1)+sum(yy2)+phi1(l)-phi2(l);
  }
  rho=exp(rho);
  yy0=sum(rho,1)*ones(1,L);
  rho=rho/yy0;
  //update delta
  for(int d=0;d<L;d++){
    delta(d)=sum(rho.col(d))+delta0(d);
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
      hstar=hstar+sum(rho.col(l))*hl;
      yh=yh+sum(rho.col(l)%y.col(j))*hj.row(l).t();
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
  ke=sum(ke,1);
  ke=ones(N*1)*ke.t();
  Vlambda=Vlambdanew;lambda=lambdanew;xi=xinew;
  
}


// [[Rcpp::export]]
double ELBO_1dim(arma::mat y,arma::cube Vlambda,arma::mat xi,arma::mat rho,arma::mat Vlambda0,arma::vec delta,arma::vec delta0,
                 arma::vec lambda00,arma::vec Vlambda00,arma::vec lambda0t,arma::vec Vlambda0t,double ph0,double pht,double df){
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
  re=re+0.5*sum(log(Vlambda0t)-log(Vlambda00))+log(pht)-log(ph0)+
    0.5*(Vlambda0t(0)*(1/Vlambda0t(0)-1/Vlambda00(0))-pow(lambda0t(0)-lambda00(0),2)/Vlambda00(0)+
    (Vlambda0t(1)-lambda0t(1)*Vlambda0t(1)*df/pht)*(1/Vlambda0t(1)-1/Vlambda00(1))-pow(lambda0t(1)-lambda00(1),2)/Vlambda00(1));
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
// [[Rcpp::export]]
double ELBO_lcdm_main_inter(arma::mat y,arma::cube Vlambda,arma::mat xi,arma::mat rho,arma::mat Vlambda0,arma::vec delta,arma::vec delta0,
                            arma::umat lambdaindd,arma::ivec lambdalen,arma::vec Mlambda0,arma::vec VMlambda0,arma::vec Mlambda0t,
                            arma::vec VMlambda0t,arma::vec ph0,arma::vec pht,arma::vec df){
  int N=y.n_rows;
  int J=y.n_cols;
  arma::mat ke;
  //int L=Vlambda0.n_cols;
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
  arma::vec re1(2);
  for(int d=1;d<3;d++){
    re1(d-1)=0.5*((VMlambda0t(d)-Mlambda0t(d)*VMlambda0t(d)*df(d-1)/pht(d-1))*(1/VMlambda0t(d)-1/VMlambda0(d))-pow(Mlambda0t(d)-Mlambda0(d),2)/VMlambda0(d));
  }
  re=re+sum(re1);
  return re;
}

// [[Rcpp::export]]
arma::vec digamma_vec(arma::vec a){
  int L=a.n_elem;
  arma::vec b(L);
  for(int i=0;i<L;i++){
    b(i)= R::digamma(a(i));
  }
  return b;
}
// [[Rcpp::export]]
arma::mat Estimate_alpha(arma::mat rho,arma::mat alpha_all,int N,int K){
  arma::mat alpha(N,K);
  for(int i=0;i<N;i++){
    arma::uvec ind = arma::find(rho.row(i) == max(rho.row(i)));
    if((ind).n_elem==1){
     alpha.row(i) = alpha_all.rows(arma::find(rho.row(i) == max(rho.row(i)))); 
    }
    if((ind).n_elem>1){
      arma::mat alpha1 = alpha_all.rows(arma::find(rho.row(i) == max(rho.row(i)))); 
      alpha.row(i)=alpha1.row(0);
    }
  }
  return(alpha);
}

// [[Rcpp::export]]
Rcpp::List VBEMM_DINA(arma::mat data,int K,arma::mat Q,double EMlambda0=2,double Vlambda0_value=1,double EMeta0=-1,double delta0_value=1,
                      double ueta0=-2,double ulambda0=0,double veta0=10,double vlambda0=10, double Mlambdainitial=2,double Vlambdainitial_value=1,
                      double Metainitial=-1,double deltainitial=1,double trunc_num=0,int T=5000,double e0=0.0001) 
{
  int N = data.n_rows;int J = data.n_cols;int L = pow(2,K);
  arma::mat ystar = 2*data-1;//0,1 response data y transform to 1,-1 ystar
  arma::mat I2;I2.eye(2,2);
  arma::mat Vlambda0 = I2*Vlambda0_value;
  arma::vec delta0 = rep(delta0_value,L);
  arma::ivec CL=seq(0,L-1);
  arma::mat alpha_all = Trans_10to2_mat(K, CL);//all attribute profile
  //some definations
  arma::cube h(J,L,2);
  double Mlambda0,Meta0,VMlambda0,VMeta0;
  arma::vec lambda00(2);lambda00(0)=ueta0;lambda00(1)=ulambda0;
  arma::vec Vlambda00(2);Vlambda00(0)=veta0;Vlambda00(1)=vlambda0;
  arma::vec lambda0t(2);
  arma::vec Elambda0t(2);
  arma::vec Vlambda0t(2);
  //initial value
  Elambda0t(0)=EMeta0;Elambda0t(1)=EMlambda0;
  arma::mat Vlambdainitial = I2*Vlambdainitial_value;
  arma::cube Vlambda(J,2,2);
  for(int j=0;j<J;j++){
    Vlambda.row(j) = Vlambdainitial;
  }
  arma::vec Mlambda = arma::ones<arma::vec>(J)*Mlambdainitial;
  arma::vec Meta = arma::ones<arma::vec>(J)*Metainitial;
  arma::mat lambda(J,2);
  lambda.col(0) = Meta;lambda.col(1)=Mlambda;
  arma::vec rhoinitial = arma::randu<arma::vec>(L);
  rhoinitial = rhoinitial/sum(rhoinitial);
  //rhoinitial=arma::ones<arma::vec>(L)/32;
  arma::mat rho = arma::ones(N,1)*rhoinitial.t();
  arma::vec delta = arma::ones<arma::vec>(L)*deltainitial;
  arma::vec pi = arma::randu<arma::vec>(L);
  pi = pi/sum(pi);//pi=arma::ones<arma::vec>(L)/32;
  arma::mat alpha=arma::ones(N,K);
  double PH0 = R::pnorm(trunc_num,  ulambda0, sqrt(vlambda0),1,0);//lower = true, log = false
  PH0 = 1-PH0;
  arma::mat H = H_allprofile_dina(alpha_all,Q);
  arma::mat hj=h.row(1);
  hj.col(0)=arma::ones<arma::vec>(L);
  for(int j=0;j<J;j++ ){
    hj.col(1)=H.col(j);
    h.row(j) = hj; 
  }
  arma::mat xi(L,J);
  xi=ke_1dim(h,lambda,Vlambda);
  double PHt,df,e;
  arma::vec phi1;arma::vec phi2;double phi2_uni;
  arma::vec elbo(T);elbo(0)=0;
  for(int t=1;t<T;t++){
    //1,update rho(pars for z_i)
    phi1 = digamma_vec(delta);
    phi2_uni = R::digamma(sum(delta));
    phi2 = arma::ones<arma::vec>(L)*phi2_uni;
    update(ystar,h,lambda,Vlambda,xi,phi1,phi2,delta,delta0,Elambda0t,Vlambda0,rho);
    Meta =lambda.col(0);
    Mlambda = lambda.col(1);
    //update prior
    VMeta0 = 1/(J+1/(veta0));
    Meta0 = VMeta0*(ueta0/veta0+sum(lambda.col(0)));
    EMeta0 = Meta0;
    
    VMlambda0 = 1/(J+1/(vlambda0));
    Mlambda0 = VMlambda0*(ulambda0/vlambda0+sum(lambda.col(1)));
    PHt = R::pnorm(trunc_num, Mlambda0,sqrt(VMlambda0),  1,  0 );PHt=1-PHt;
    df = R::dnorm(trunc_num,Mlambda0,sqrt(VMlambda0),0);
    EMlambda0 = Mlambda0+VMlambda0*df/PHt;
    lambda0t(0)=Meta0;lambda0t(1)=Mlambda0;
    Elambda0t(0)=EMeta0;Elambda0t(1)=EMlambda0;
    Vlambda0t(0)=VMeta0;Vlambda0t(1)=VMlambda0;
    elbo(t)=ELBO_1dim(ystar,Vlambda,xi,rho,Vlambda0,delta,delta0,lambda00,Vlambda00,lambda0t,Vlambda0t,PH0,PHt,df);
    e=abs(elbo(t)-elbo(t-1));
    Rcout << "Iteration = "<< t ;
    Rcout <<";ELBO change = "<< e<<'\r';
    if(e<e0){
      break;
    }
  }
  pi=delta/sum(delta);
  alpha=Estimate_alpha(rho,alpha_all, N, K);
  return Rcpp::List::create(Rcpp::Named("Est_alpha",alpha),
                            Rcpp::Named("Est_eta",Meta),
                            Rcpp::Named("Est_lambda",Mlambda),
                            Rcpp::Named("Est_pi",pi),
                            Rcpp::Named("delta",delta),
                            Rcpp::Named("rho",rho),
                            Rcpp::Named("Vlambda",Vlambda)
  );
  
}

// [[Rcpp::export]]
Rcpp::List VBEMM_DINO(arma::mat data,int K,arma::mat Q,double EMlambda0=2,double Vlambda0_value=1,double EMeta0=-1,double delta0_value=1,
                      double ueta0=-2,double ulambda0=0,double veta0=10,double vlambda0=10, double Mlambdainitial=2,double Vlambdainitial_value=1,
                      double Metainitial=-1,double deltainitial=1,double trunc_num=0,int T=5000,double e0=0.0001) 
{
  int N = data.n_rows;int J = data.n_cols;int L = pow(2,K);
  arma::mat ystar = 2*data-1;//0,1 response data y transform to 1,-1 ystar
  arma::mat I2;I2.eye(2,2);
  arma::mat Vlambda0 = I2*Vlambda0_value;
  arma::vec delta0 = rep(delta0_value,L);
  arma::ivec CL=seq(0,L-1);
  arma::mat alpha_all = Trans_10to2_mat(K, CL);//all attribute profile
  //some definations
  arma::cube h(J,L,2);
  double Mlambda0,Meta0,VMlambda0,VMeta0;
  arma::vec lambda00(2);lambda00(0)=ueta0;lambda00(1)=ulambda0;
  arma::vec Vlambda00(2);Vlambda00(0)=veta0;Vlambda00(1)=vlambda0;
  arma::vec lambda0t(2);
  arma::vec Elambda0t(2);
  arma::vec Vlambda0t(2);
  //initial value
  Elambda0t(0)=EMeta0;Elambda0t(1)=EMlambda0;
  arma::mat Vlambdainitial = I2*Vlambdainitial_value;
  arma::cube Vlambda(J,2,2);
  for(int j=0;j<J;j++){
    Vlambda.row(j) = Vlambdainitial;
  }
  arma::vec Mlambda = arma::ones<arma::vec>(J)*Mlambdainitial;
  arma::vec Meta = arma::ones<arma::vec>(J)*Metainitial;
  arma::mat lambda(J,2);
  lambda.col(0) = Meta;lambda.col(1)=Mlambda;
  arma::vec rhoinitial = arma::randu<arma::vec>(L);
  rhoinitial = rhoinitial/sum(rhoinitial);
  //rhoinitial=arma::ones<arma::vec>(L)/32;
  arma::mat rho = arma::ones(N,1)*rhoinitial.t();
  arma::vec delta = arma::ones<arma::vec>(L)*deltainitial;
  arma::vec pi = arma::randu<arma::vec>(L);
  pi = pi/sum(pi);//pi=arma::ones<arma::vec>(L)/32;
  arma::mat alpha=arma::ones(N,K);
  double PH0 = R::pnorm(trunc_num,  ulambda0, sqrt(vlambda0),1,0);//lower = true, log = false
  PH0 = 1-PH0;
  arma::mat H = H_allprofile_dino(alpha_all,Q);
  arma::mat hj=h.row(1);
  hj.col(0)=arma::ones<arma::vec>(L);
  for(int j=0;j<J;j++ ){
    hj.col(1)=H.col(j);
    h.row(j) = hj; 
  }
  arma::mat xi(L,J);
  xi=ke_1dim(h,lambda,Vlambda);
  double PHt,df,e;
  arma::vec phi1;arma::vec phi2;double phi2_uni;
  arma::vec elbo(T);elbo(0)=0;
  for(int t=1;t<T;t++){
    //1,update rho(pars for z_i)
    phi1 = digamma_vec(delta);
    phi2_uni = R::digamma(sum(delta));
    phi2 = arma::ones<arma::vec>(L)*phi2_uni;
    update(ystar,h,lambda,Vlambda,xi,phi1,phi2,delta,delta0,Elambda0t,Vlambda0,rho);
    Meta =lambda.col(0);
    Mlambda = lambda.col(1);
    //update prior
    VMeta0 = 1/(J+1/(veta0));
    Meta0 = VMeta0*(ueta0/veta0+sum(lambda.col(0)));
    EMeta0 = Meta0;
    
    VMlambda0 = 1/(J+1/(vlambda0));
    Mlambda0 = VMlambda0*(ulambda0/vlambda0+sum(lambda.col(1)));
    PHt = R::pnorm(trunc_num, Mlambda0,sqrt(VMlambda0),  1,  0 );PHt=1-PHt;
    df = R::dnorm(trunc_num,Mlambda0,sqrt(VMlambda0),0);
    EMlambda0 = Mlambda0+VMlambda0*df/PHt;
    lambda0t(0)=Meta0;lambda0t(1)=Mlambda0;
    Elambda0t(0)=EMeta0;Elambda0t(1)=EMlambda0;
    Vlambda0t(0)=VMeta0;Vlambda0t(1)=VMlambda0;
    elbo(t)=ELBO_1dim(ystar,Vlambda,xi,rho,Vlambda0,delta,delta0,lambda00,Vlambda00,lambda0t,Vlambda0t,PH0,PHt,df);
    e=abs(elbo(t)-elbo(t-1));
    Rcout << "Iteration = "<< t ;
    Rcout <<";ELBO change = "<< e<<'\r';
    if(e<e0){
      break;
    }
  }
  pi=delta/sum(delta);
  alpha=Estimate_alpha(rho,alpha_all, N, K);
  return Rcpp::List::create(Rcpp::Named("Est_alpha",alpha),
                            Rcpp::Named("Est_eta",Meta),
                            Rcpp::Named("Est_lambda",Mlambda),
                            Rcpp::Named("Est_pi",pi),
                            Rcpp::Named("delta",delta),
                            Rcpp::Named("rho",rho),
                            Rcpp::Named("Vlambda",Vlambda)
  );
  
}
// [[Rcpp::export]]
Rcpp::List VBEMM_LCDM(arma::mat data,int K,arma::mat Q,arma::umat lambdaindex,arma::umat comb,arma::ivec comlen,double EMlambda0_value=2,double Vlambda0_value=1,
                      double EMeta0=-1,double delta0_value=1,double ueta0=-2,double ulambda0=0,double veta0=10,double vlambda0=10, double Mlambdainitial=2,
                      double Vlambdainitial_value=1,double Metainitial=-1,double deltainitial=1,int T=5000,double e0=0.0001) 
{
  int N = data.n_rows;int J = data.n_cols;int L = pow(2,K);int D=L-1;
  arma::mat ystar = 2*data-1;//0,1 response data y transform to 1,-1 ystar
  arma::mat IL;IL.eye(L,L);
  arma::mat Vlambda0 = IL*Vlambda0_value;
  arma::vec delta0 = rep(delta0_value,L);
  arma::ivec CL=seq(0,L-1);
  arma::mat alpha_all = Trans_10to2_mat(K, CL);//all attribute profile
  arma::vec EMlambda0 = rep(EMlambda0_value,D);
  //arma::vec Mlambdainitial=rep(Mlambdainitial_value,D);
  arma::mat Vlambdainitial=IL*Vlambdainitial_value;
  //some definations
  arma::cube h(J,L,L);
  arma::mat alpha(N,K);
  arma::umat lambdaindd=arma::zeros<arma::umat>(J,D);
  arma::ivec lambdalen(J); //= (sum(lambdaindex,0)).t();
  arma::uvec lambdaindd_j;arma::uvec  kk;
  for(int j=0;j<J;j++){
    kk=arma::find(lambdaindex.row(j)>0);
    lambdalen(j)=kk.n_elem;
    lambdaindd_j=(lambdaindd.row(j)).t();
    lambdaindd_j.head(lambdalen(j)) = arma::find(lambdaindex.row(j)==1)+1;
    lambdaindd.row(j)=lambdaindd_j.t();
  }
  arma::umat lamind_all=arma::zeros<arma::umat>(J,L);
  lamind_all.cols(1,L-1)=lambdaindd;
  arma::ivec lambdalen_c=lambdalen+1;
  //initial value
  double Meta0,VMeta0;
  arma::vec Mlambda0(2),VMlambda0(2);
  arma::vec lambda00=rep(ulambda0,3);lambda00(0)=ueta0;//lambda00.tail(D)=rep(ulambda0,D);
  arma::vec Vlambda00=rep(vlambda0,3);Vlambda00(0)=veta0;//Vlambda00.tail(D)=rep(vlambda0,D);
  arma::vec lambda0t(3);
  arma::vec Elambda0t(L);Elambda0t(0)=EMeta0;Elambda0t.tail(D)=EMlambda0;
  arma::vec Vlambda0t(3);
  arma::cube Vlambda(J,L,L);
  for(int j=0;j<J;j++ ){
    Vlambda.row(j) = Vlambdainitial;
  }
  arma::mat Mlambda = arma::ones(J,D)*Mlambdainitial;
  arma::vec Meta = arma::ones<arma::vec>(J)*Metainitial;
  arma::mat lambda(J,L);
  lambda.col(0) = Meta;lambda.cols(1,L-1)=Mlambda;
  arma::vec rhoinitial = arma::randu<arma::vec>(L);
  rhoinitial = rhoinitial/sum(rhoinitial);
  //rhoinitial=arma::ones<arma::vec>(L)/32;
  arma::mat rho = arma::ones(N,1)*rhoinitial.t();
  arma::vec delta = arma::ones<arma::vec>(L)*deltainitial;
  arma::vec pi = arma::randu<arma::vec>(L);
  pi = pi/sum(pi);//pi=arma::ones<arma::vec>(L)/32;
  arma::cube H = H_allprofile_lcdm(alpha_all,Q,comb,comlen);
  arma::mat hj(L,L);arma::mat H_j(L,D);
  hj.col(0)=arma::ones<arma::vec>(L);
  for(int j=0;j<J;j++ ){
    H_j=H.row(j);
    hj.cols(1,L-1)=H_j;
    h.row(j) = hj; 
  }
  arma::mat xi=ke_1dim(h,lambda,Vlambda);
  int lengthd=0;
  for(int d=0;d<K;d++){
    kk=arma::find(lambdaindex.col(d)>0);
    lengthd=lengthd+(kk).n_elem;
  }
  VMlambda0(0)= 1/(lengthd+1/(vlambda0));
  lengthd=0;
  for(int d=K;d<D;d++){
    kk=arma::find(lambdaindex.col(d)>0);
    lengthd=lengthd+(kk).n_elem;
  }
  VMlambda0(1)= 1/(lengthd+1/(vlambda0));
  VMeta0 = 1/(J+1/(veta0));
  double PH0 = R::pnorm(0,  ulambda0, sqrt(vlambda0),1,0);//lower = true, log = false
  PH0 = 1-PH0;
  arma::vec PH0_vec=rep(PH0,2);
  PH0_vec(1)=1-R::pnorm(R_NegInf,  ulambda0, sqrt(vlambda0),1,0);
  arma::vec elbo(T);elbo(0)=0;
  double e=1.0;  
  arma::vec PHt(2),df(2);
  arma::vec phi1;arma::vec phi2;double phi2_uni;
  for(int t=1;t<T;t++){
    //1,update rho(pars for z_i)
    phi1 = digamma_vec(delta);
    phi2_uni = R::digamma(sum(delta));
    phi2 = arma::ones<arma::vec>(L)*phi2_uni;
    
    update_lcdm(ystar,h,lambda,Vlambda,xi,phi1,phi2,delta0,delta,rho,Elambda0t,Vlambda0,lamind_all,lambdalen_c);
    Meta =lambda.col(0);
    Mlambda = lambda.cols(1,L-1);
    //update prior
    Meta0= VMeta0*(ueta0/veta0+sum(Meta));
    EMeta0 = Meta0;
    arma::uvec inl;arma::vec Mlambda_d;
    double lam_sum=0;
    for(int d=0;d<K;d++){
      inl = arma::find(Mlambda.col(d)!=0);
      Mlambda_d=Mlambda.col(d);
      lam_sum = lam_sum+sum(Mlambda_d.elem(inl));
    }
    Mlambda0(0) = VMlambda0(0)*(ulambda0/(vlambda0)+lam_sum);
    PHt(0) =1- R::pnorm(0, Mlambda0(0),sqrt(VMlambda0(0)),  1,  0 );
    df(0) = R::dnorm(0,Mlambda0(0),sqrt(VMlambda0(0)),0);
    double Elam_main = Mlambda0(0)+VMlambda0(0)*df(0)/PHt(0);
    EMlambda0.head(K) = arma::ones<arma::vec>(K)*Elam_main;
    lam_sum=0;
    for(int d=K;d<D;d++){
      inl = arma::find(Mlambda.col(d)!=0);
      Mlambda_d=Mlambda.col(d);
      lam_sum = lam_sum+sum(Mlambda_d.elem(inl));
    }
    Mlambda0(1) = VMlambda0(1)*(ulambda0/(vlambda0)+lam_sum);
    PHt(1) =1- R::pnorm(R_NegInf, Mlambda0(0),sqrt(VMlambda0(0)),  1,  0 );
    df(1) = R::dnorm(R_NegInf,Mlambda0(1),sqrt(VMlambda0(1)),0);
    double Elam_inter = Mlambda0(1)+VMlambda0(1)*df(1)/PHt(1);
    EMlambda0.tail(D-K) = arma::ones<arma::vec>(D-K)*Elam_inter;
    
    lambda0t(0)=Meta0;lambda0t.tail(2)=Mlambda0;
    Elambda0t(0)=EMeta0;Elambda0t.tail(D)=EMlambda0;
    Vlambda0t(0)=VMeta0;Vlambda0t.tail(2)=VMlambda0;
    
    elbo(t)=ELBO_lcdm_main_inter(ystar,Vlambda,xi,rho,Vlambda0,delta,delta0,lamind_all,lambdalen_c,lambda00,
         Vlambda00,lambda0t,Vlambda0t,PH0_vec,PHt,df);
    e=abs(elbo(t)-elbo(t-1));
    Rcout << "Iteration = "<< t ;
    Rcout <<";ELBO change = "<< e<<'\r';
    if(e<e0){
      break;
    }                  
  }
  
  pi=delta/sum(delta);
  alpha=Estimate_alpha(rho,alpha_all, N, K);
  return Rcpp::List::create(Rcpp::Named("Est_alpha",alpha),
                            Rcpp::Named("Est_eta",Meta),
                            Rcpp::Named("Est_lambda",Mlambda),
                            Rcpp::Named("Est_pi",pi),
                            Rcpp::Named("delta",delta),
                            Rcpp::Named("rho",rho),
                            Rcpp::Named("Vlambda",Vlambda)
  );
  
}


// [[Rcpp::export]]
Rcpp::List VBEMM_LCDM_by_dim(arma::mat data,int K,arma::mat Q,arma::umat lambdaindex,arma::umat comb,arma::ivec comlen,double EMlambda0_value=2,double Vlambda0_value=1,
                       double EMeta0=-1,double delta0_value=1,double ueta0=-2,double ulambda0=0,double veta0=10,double vlambda0=10, double Mlambdainitial=2,
                       double Vlambdainitial_value=1,double Metainitial=-1,double deltainitial=1,double trunc_num=0,int T=5000,double e0=0.0001) 
{
  int N = data.n_rows;int J = data.n_cols;int L = pow(2,K);int D=L-1;
  arma::mat ystar = 2*data-1;//0,1 response data y transform to 1,-1 ystar
  arma::mat IL;IL.eye(L,L);
  arma::mat Vlambda0 = IL*Vlambda0_value;
  arma::vec delta0 = rep(delta0_value,L);
  arma::ivec CL=seq(0,L-1);
  arma::mat alpha_all = Trans_10to2_mat(K, CL);//all attribute profile
  arma::vec EMlambda0 = rep(EMlambda0_value,D);
  //arma::vec Mlambdainitial=rep(Mlambdainitial_value,D);
  arma::mat Vlambdainitial=IL*Vlambdainitial_value;
  //some definations
  arma::cube h(J,L,L);
  arma::mat alpha(N,K);
  arma::umat lambdaindd=arma::zeros<arma::umat>(J,D);
  arma::ivec lambdalen(J); //= (sum(lambdaindex,0)).t();
  arma::uvec lambdaindd_j;arma::uvec  kk;
  for(int j=0;j<J;j++){
    kk=arma::find(lambdaindex.row(j)>0);
    lambdalen(j)=kk.n_elem;
    lambdaindd_j=(lambdaindd.row(j)).t();
    lambdaindd_j.head(lambdalen(j)) = arma::find(lambdaindex.row(j)==1)+1;
    lambdaindd.row(j)=lambdaindd_j.t();
  }
  arma::umat lamind_all=arma::zeros<arma::umat>(J,L);
  lamind_all.cols(1,L-1)=lambdaindd;
  arma::ivec lambdalen_c=lambdalen+1;
  //initial value
  double Meta0,VMeta0;
  arma::vec Mlambda0(D),VMlambda0(D);
  arma::vec lambda00=rep(ulambda0,L);lambda00(0)=ueta0;//lambda00.tail(D)=rep(ulambda0,D);
  arma::vec Vlambda00=rep(vlambda0,L);Vlambda00(0)=veta0;//Vlambda00.tail(D)=rep(vlambda0,D);
  arma::vec lambda0t(L);
  arma::vec Elambda0t(L);Elambda0t(0)=EMeta0;Elambda0t.tail(D)=EMlambda0;
  arma::vec Vlambda0t(L);
  arma::cube Vlambda(J,L,L);
  for(int j=0;j<J;j++ ){
    Vlambda.row(j) = Vlambdainitial;
  }
  arma::mat Mlambda = arma::ones(J,D)*Mlambdainitial;
  arma::vec Meta = arma::ones<arma::vec>(J)*Metainitial;
  arma::mat lambda(J,L);
  lambda.col(0) = Meta;lambda.cols(1,L-1)=Mlambda;
  arma::vec rhoinitial = arma::randu<arma::vec>(L);
  rhoinitial = rhoinitial/sum(rhoinitial);
  //rhoinitial=arma::ones<arma::vec>(L)/32;
  arma::mat rho = arma::ones(N,1)*rhoinitial.t();
  arma::vec delta = arma::ones<arma::vec>(L)*deltainitial;
  arma::vec pi = arma::randu<arma::vec>(L);
  pi = pi/sum(pi);//pi=arma::ones<arma::vec>(L)/32;
  arma::cube H = H_allprofile_lcdm(alpha_all,Q,comb,comlen);
  arma::mat hj(L,L);arma::mat H_j(L,D);
  hj.col(0)=arma::ones<arma::vec>(L);
  for(int j=0;j<J;j++ ){
    H_j=H.row(j);
    hj.cols(1,L-1)=H_j;
    h.row(j) = hj; 
  }
  arma::mat xi=ke_1dim(h,lambda,Vlambda);
  arma::ivec lengthd(D);
  for(int d=0;d<D;d++){
    kk=arma::find(lambdaindex.col(d)>0);
    lengthd(d)=(kk).n_elem;
    VMlambda0(d)= 1/(lengthd(d)+1/(vlambda0));
  }
  VMeta0 = 1/(J+1/(veta0));
  double PH0 = R::pnorm(trunc_num,  ulambda0, sqrt(vlambda0),1,0);//lower = true, log = false
  PH0 = 1-PH0;
  arma::vec PH0_vec=rep(PH0,D);
  arma::vec elbo(T);elbo(0)=0;
  double e=1.0;  
  arma::vec PHt(D),df(D);
  arma::vec phi1;arma::vec phi2;double phi2_uni;
  for(int t=1;t<T;t++){
    //1,update rho(pars for z_i)
    phi1 = digamma_vec(delta);
    phi2_uni = R::digamma(sum(delta));
    phi2 = arma::ones<arma::vec>(L)*phi2_uni;
    
    update_lcdm(ystar,h,lambda,Vlambda,xi,phi1,phi2,delta0,delta,rho,Elambda0t,Vlambda0,lamind_all,lambdalen_c);
    Meta =lambda.col(0);
    Mlambda = lambda.cols(1,L-1);
    //update prior
    Meta0= VMeta0*(ueta0/veta0+sum(Meta));
    EMeta0 = Meta0;
    arma::uvec inl;arma::vec Mlambda_d;
    for(int d=0;d<D;d++){
      inl = arma::find(Mlambda.col(d)!=0);
      Mlambda_d=Mlambda.col(d);
      Mlambda0(d) = VMlambda0(d)*(ulambda0/(vlambda0)+sum(Mlambda_d.elem(inl)));
      PHt(d) =1- R::pnorm(trunc_num, Mlambda0(d),sqrt(VMlambda0(d)),  1,  0 );
      df(d) = R::dnorm(trunc_num,Mlambda0(d),sqrt(VMlambda0(d)),0);
      EMlambda0(d) = Mlambda0(d)+VMlambda0(d)*df(d)/PHt(d);
    }
    lambda0t(0)=Meta0;lambda0t.tail(D)=Mlambda0;
    Elambda0t(0)=EMeta0;Elambda0t.tail(D)=EMlambda0;
    Vlambda0t(0)=VMeta0;Vlambda0t.tail(D)=VMlambda0;
    
    elbo(t)=ELBO_lcdm(ystar,Vlambda,xi,rho,Vlambda0,delta,delta0,lamind_all,lambdalen_c,lambda00,
         Vlambda00,lambda0t,Vlambda0t,PH0_vec,PHt,df);
    e=abs(elbo(t)-elbo(t-1));
    Rcout << "Iteration = "<< t ;
    Rcout <<";ELBO change = "<< e<<'\r';
    if(e<e0){
      break;
    }                  
  }
  
  pi=delta/sum(delta);
  alpha=Estimate_alpha(rho,alpha_all, N, K);
  return Rcpp::List::create(Rcpp::Named("Est_alpha",alpha),
                            Rcpp::Named("Est_eta",Meta),
                            Rcpp::Named("Est_lambda",Mlambda),
                            Rcpp::Named("Est_pi",pi),
                            Rcpp::Named("delta",delta),
                            Rcpp::Named("rho",rho),
                            Rcpp::Named("Vlambda",Vlambda)
  );
  
}


