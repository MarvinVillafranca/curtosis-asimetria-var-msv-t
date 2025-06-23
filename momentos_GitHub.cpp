#include <RcppArmadillo.h>
#include<RcppArmadilloExtensions/sample.h>
#include<Rcpp.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
using namespace std;
using namespace Rcpp;
using namespace arma;

/*#====================================================vec====================================================*/
// [[Rcpp::export]] 
vec vec_Armadillo(mat A)
{
  int rows = A.n_rows;
  int cols = A.n_cols;
  vec d(rows*cols);
  for(int i=0;i<cols;i++)
  {
    d.subvec(i*rows,(i+1)*rows-1) = A.col(i);
  }
  return d;
}

/*#====================================================in_vec====================================================*/
// [[Rcpp::export]] 
mat inv_vec_Armadillo(vec A,int rows,int cols)
{
  mat d(rows,cols);
  for(int i=0;i<cols;i++)
  {
    d.col(i) = A.subvec(i*rows,(i+1)*rows-1);
  }
  return d;
}

/*================================Funcion Vt con RcppArmadillo================================*/
// [[Rcpp::export]]
mat Vt_Armadillo(mat alpha, int t)
{
  int p = alpha.n_rows;
  mat Vt = diagmat(exp(alpha.col(t-1) / 2)); // Directly create diagonal matrix
  return Vt;
}

/*================================Funcion Vt inversa con RcppArmadillo================================*/
// [[Rcpp::export]]
mat solve_Vt_Armadillo(mat alpha, int t)
{
  int p = alpha.n_rows;
  mat solve_Vt = diagmat(exp(-alpha.col(t-1) / 2)); // Directly create diagonal matrix
  return solve_Vt;
}

/*================================Funcion mt con RcppArmadillo================================*/
// [[Rcpp::export]]
vec mt_Armadillo(mat alpha, mat phi,vec diag_phi, mat solve_SenSnn, int n, int t)
{
  int p = phi.n_rows;
  vec mt;
  if (t < n)
  {
    mt = solve_SenSnn * (alpha.col(t) - diag_phi%alpha.col(t - 1));
  }
  else
  {
    mt = arma::zeros<arma::vec>(p);  // Correcta inicialización de mt
  }
  return mt;
}

/*#========================================Funcion Mt con Rcpp========================================*/
// [[Rcpp::export]]
vec Mt_Armadillo(mat alpha, mat phi,vec diag_phi, mat solve_SenSnn, int n, int t)
{
  return Vt_Armadillo(alpha, t) * mt_Armadillo(alpha, phi,diag_phi, solve_SenSnn, n, t);
}

//==================================St======================================*/
mat St_Armadillo(mat See,mat SeeSenSnnSne,int n,int t)
{
  if(t<n)
  {
    return SeeSenSnnSne;
  }
  else
  {
    return See;
  }
}

//==================================Priori Minesota==================================
mat PrioriMinesota(double lambda,double teta,int n,int p, int k, mat See)
{
  mat V(p*p*k+p,p*p*k+p);
  V.zeros();
  int contador = 0;
  
  for(int l=1;l<=k;l++)
  {
    for(int j=0;j<p;j++)
    {
      for(int i=0;i<p;i++)
      {
        if(i==j)
        {
          V(p + contador,p + contador) =  ((l*l)/(lambda*lambda));
        }
        else
        {
          V(p + contador,p + contador) = (l*l*See(j,j))/(lambda*teta*lambda*teta*See(i,i));
        }
        contador = contador + 1;
      } 
    } 
  }
  return V;
}

//==================================Estimacion de A1,...,Ak======================================*/
// [[Rcpp::export]]
mat estimacionA0_Ak_Armadillo(mat alpha,mat Yn,mat lambda,mat phi,vec diag_phi,mat sigma,mat COV_prior,int n, int p,int k, mat Yk,mat A1_temp,double lambda_conts,double teta_conts)
{
  mat See = sigma.submat(0,0,(p-1),(p-1));
  mat Sen = sigma.submat(0,p,(p-1),(2*p-1));
  mat Sne = sigma.submat(p,0,(2*p-1),(p-1));
  mat Snn = sigma.submat(p,p,(2*p-1),(2*p-1));
  mat solve_Snn = inv_sympd(Snn);
  mat solve_SnnSne = solve_Snn*Sne;
  mat solve_SenSnn = solve_SnnSne.t();
  mat SeeSenSnnSne = See - Sen*solve_SnnSne;
  
  mat Y = join_rows(Yk,Yn);
  vec Yt(k*p+1);
  Yt.ones();
  int T = n-k;
  vec Mb;
  Mb.zeros(k*p*p+p);
  mat Vb;
  mat inv_Vb;
  Vb.zeros(k*p*p+p,k*p*p+p);
  mat Vt;
  vec Mt;
  mat Ct;
  mat inv_Ct;
  mat diag;
  diag.eye(p,p);
  mat temp;
  for(int i=0;i<n;i++)
  {
    Vt = Vt_Armadillo(alpha,i+1);
    Mt = pow(lambda(0,i),-0.5)*Mt_Armadillo(alpha,phi,diag_phi,solve_SenSnn,n,i+1);
    Ct = pow(lambda(0,i),-1)*Vt*St_Armadillo(See,SeeSenSnnSne,n,i+1)*Vt;
    for(int i1=0;i1<p;i1++)
    {
      for(int i2=(i1+1);i2<p;i2++)
      {
        Ct(i2,i1) = Ct(i1,i2);
      } 
    }
    inv_Ct = inv_sympd(Ct);
    for(int j=0;j<k;j++)
    {
      Yt.subvec(j*p+1,(j+1)*p) = Y.col(i+k-1-j);
    }
    temp = kron(Yt.t(),diag);
    Mb = Mb + temp.t()*inv_Ct*(Yn.col(i)-Mt);
    Vb = Vb + temp.t()*inv_Ct*temp;
  }
  //Esto asegura la simetria para que no se pierdan datos
  for(int i=0;i<p;i++)
  {
    for(int j=(i+1);j<p;j++)
    {
      Vb(j,i) = Vb(i,j);
    } 
  }
  
  mat PM = PrioriMinesota(lambda_conts,teta_conts,n,p,k,COV_prior);
  mat diag_PM = PM.diag();
  inv_Vb = inv_sympd(Vb+PM);
  Mb = inv_Vb*Mb;
  vec vec_A = mvnrnd(Mb,inv_Vb);
  mat A = inv_vec_Armadillo(vec_A,p,k*p+1);
  return A;
}

/*=================================yt-v-A1yt-1,....Akyt-k==========================*/
// [[Rcpp::export]]
mat Y_var_Armadillo(mat Yn,mat Yk,mat VAR,int p,int n,int k)
{
  mat Y(p,n+k);
  Y.submat(0,0,p-1,k-1) = Yk;
  Y.submat(0,k,p-1,n+k-1) = Yn;
  vec y = vec_Armadillo(Y);
  mat temp(p,n);
  vec Yt;
  Yt.ones(k*p+1);
  
  for(int i=0;i<n;i++)
  {
    for(int j=(k-1);j>=0;j--)
    {
      Yt.subvec((1+(k-j-1)*p),(k-j)*p) = y.subvec((j*p+i*p),((j+1)*p-1+i*p));
    }
    temp.col(i) = Yn.col(i) - VAR*Yt;
  }
  return temp;
}

/*=================================Sigma w==========================*/
// [[Rcpp::export]]
mat Sigma_w_Armadillo(mat See,mat sigma0,mat phi,double v,int p)
{
  mat temp(p,p);
  for(int i=0;i<p;i++)
  {
    temp(i,i) = (v/(v-2))*exp(0.5*sigma0(i,i))*See(i,i);
    for(int j=i+1;j<p;j++)
    {
      temp(i,j) = (v/(v-2))*exp(0.125*(sigma0(i,i)+sigma0(j,j)+2*sigma0(i,j)))*See(i,j);
      temp(j,i) = temp(i,j);
    } 
  }
  return temp;
}

/*=================================Propiedad 4==========================*/
// [[Rcpp::export]]
double p4_Armadillo(mat See,mat sigma0,double v,int i,int j,int k,int l)
{
  double temp;
  temp = (pow(v,2)/((v-2)*(v-4)))*(See(i,j)*See(k,l)+See(i,k)*See(j,l)+See(i,l)*See(j,k))*exp(0.125*(sigma0(i,i)+sigma0(j,j)+sigma0(k,k)+sigma0(l,l)+2*(sigma0(i,j)+sigma0(i,k)+sigma0(i,l)+sigma0(j,k)+sigma0(j,l)+sigma0(k,l))));
  return temp;
}

/*=================================Propiedad 5==========================*/
// [[Rcpp::export]]
double p5_Armadillo(mat See,mat sigma0,double v,int i,int j)
{
  double temp;
  temp = ((3*pow(v,2)*See(i,i)*See(i,j))/((v-2)*(v-4)))*exp(0.125*(9*sigma0(i,i)+sigma0(j,j)+6*sigma0(i,j)));
  return temp;
}

/*=================================Propiedad 6==========================*/
// [[Rcpp::export]]
double p6_Armadillo(mat See,mat sigma0,double v,int i,int j)
{
  double temp;
  temp = ((pow(v,2)*(See(i,i)*See(j,j)+2*pow(See(i,j),2)))/((v-2)*(v-4)))*exp(0.5*(sigma0(i,i)+sigma0(j,j)+2*sigma0(i,j)));
  return temp;
}

/*=================================Propiedad 7==========================*/
// [[Rcpp::export]]
double p7_Armadillo(mat See,mat sigma0,double v,int i,int j,int k)
{
  double temp;
  temp = ((pow(v,2)*(See(i,i)*See(j,k)+2*See(i,j)*See(i,k)))/((v-2)*(v-4)))*exp(0.125*(4*sigma0(i,i)+sigma0(j,j)+sigma0(k,k)+2*(2*sigma0(i,j)+2*sigma0(i,k)+sigma0(j,k))));
  return temp;
}

/*=================================Propiedad 8==========================*/
// [[Rcpp::export]]
double p8_Armadillo(mat See,mat sigma0,double v,int i)
{
  double temp;
  temp = 3*(pow(v,2)/((v-2)*(v-4)))*exp(2*sigma0(i,i))*pow(See(i,i),2);
  return temp;
}

/*============================selecciona la propiedad==========================*/
// [[Rcpp::export]]
double Ew_Armadillo(mat See,mat sigma0,double v,int p,int i, int j, int k, int l)
{
  double temp;
  if(i!=j & i!=k & i!=l & j!=k & j!=l & k!=l) /*1 i != j != k != l*/ 
  {
    temp = p4_Armadillo(See,sigma0,v,i,j,k,l);
  }
  else if(i==j & i==k & i!=l)                   /*2 tres y uno*/
  {
    temp = p5_Armadillo(See,sigma0,v,i,l);
  }
  else if(i==k & i==l & i!=j)                   /*3 tres y uno*/
  {
    temp = p5_Armadillo(See,sigma0,v,i,j);
  }
  else if(i==l & i==j & i!=k)                   /*4 tres y uno*/
  {
    temp = p5_Armadillo(See,sigma0,v,i,k);
  }
  else if(j==k & j==l & j!=i)                   /*5 tres y uno*/
  {
    temp = p5_Armadillo(See,sigma0,v,j,i);
  }
  else if(i==j & k==l & i!=k)                 /*6 i=j, k=l, i!=k*/
  {
    temp = p6_Armadillo(See,sigma0,v,i,k);
  }
  else if(i==l & k==j & i!=k)                 /*7 i=l, k=j, i!=k*/
  {
    temp = p6_Armadillo(See,sigma0,v,i,k);
  }
  else if(i==k & j==l & i!=j)                   /*8 i=k, j=l, i!=j**/
  {
    temp = p6_Armadillo(See,sigma0,v,i,j);
  }
  else if(i==j & i!=k & i!=l & k!=l)                 /*9 dos iguales y dos distintos*/
  {
    temp = p7_Armadillo(See,sigma0,v,i,k,l);
  }
  else if(i==k & i!=l & i!=j & l!=j)                        /*10 dos iguales y dos distintos*/
  {
    temp = p7_Armadillo(See,sigma0,v,i,j,l); 
  }
  else if(i==l & i!=j & i!=k & j!=k)                        /*11 dos iguales y dos distintos*/
  {
    temp = p7_Armadillo(See,sigma0,v,i,j,k);
  }
  else if(j==k & j!=i & j!=l & i!=l)                        /*12 dos iguales y dos distintos*/
  {
    temp = p7_Armadillo(See,sigma0,v,j,i,l);
  }
  else if(j==l & j!=i & j!=k & i!=k)                        /*13 dos iguales y dos distintos*/
  {
    temp = p7_Armadillo(See,sigma0,v,j,i,k);
  }
  else if(k==l & k!=i & k!=j & i!=j)                        /*14 dos iguales y dos distintos*/
  {
    temp = p7_Armadillo(See,sigma0,v,k,i,j);
  }
  else if(i==j & i==k & i==l & j==k & j==l & k==l)   /*15 todos iguales*/
  {
    temp = p8_Armadillo(See,sigma0,v,i);
  }
  
  
  return temp;
}

/*=================================Kurtosis Total==========================*/
// [[Rcpp::export]]
double Kurtosis_total_Armadillo(mat See,mat sigma0,mat phi,double v,int p)
{
  mat L_w = chol(Sigma_w_Armadillo(See,sigma0,phi,v,p));
  L_w = inv(L_w.t());
  mat L_w_K = kron(kron(L_w,L_w),kron(L_w,L_w));
  int P = pow(p,4);
  mat indices(P,p);
  vec E_w(P);
  int i1 = 0;
  int i2 = 0;
  int i3 = 0;
  int i4 = 0;
  for(int i = 0;i<P;i++)
  {
    if(i%(P/p)==0 & i!=0) 
    {
      i1 = i1 + 1;
      if(i1==p)
      {
        i1=0;
      }
    }
    if(i%(P/(p*p))==0 & i!=0)
    {
      i2 = i2 + 1;
      if(i2==p)
      {
        i2=0;
      }
    }
    if(i%(P/(p*p*p))==0 & i!=0)
    {
      i3 = i3 + 1;
      if(i3==p)
      {
        i3=0;
      }
    }
    if(i4==p)
    {
      i4 = 0;
    }
    indices(i,0) = i1;
    indices(i,1) = i2;
    indices(i,2) = i3;
    indices(i,3) = i4;
    E_w(i) = Ew_Armadillo(See,sigma0,v,p,i1,i2,i3,i4);
    i4 = i4 + 1;
  }
  vec temp;
  double Kurtosis_total = 0;
  for(int i=0;i<P;i++)
  {
    temp = L_w_K.row(i)*E_w;
    Kurtosis_total = Kurtosis_total + pow(temp(0),2);
  }
  return Kurtosis_total;
}

/*=================================Covarianza de yt ==========================*/
// [[Rcpp::export]]
mat COV_yt_Armadillo(mat COV_w,mat A,int k,int p)
{
  mat I1;
  I1.eye(p,p);
  mat J(p,k*p);
  J.zeros(p,k*p);
  J.submat(0,0,p-1,p-1) = I1;
  mat I2;
  I2.eye(pow(k*p,2),pow(k*p,2));
  mat Sw;
  Sw.zeros(k*p,k*p);
  Sw.submat(0,0,p-1,p-1) = COV_w;
  vec vec_COV_yt = kron(J,J)*inv(I2-kron(A,A))*vec_Armadillo(Sw);
  mat COV_yt = inv_vec_Armadillo(vec_COV_yt,p,p);
  return COV_yt;
}
