// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  
using namespace Rcpp;
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <vector>
#include "Classes_C.h"
// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#else 
#define omp_get_thread_num() 0 
#endif 

/*  LSMLE - Functions C++  */

arma::mat Psi_LSMLE_Beta_Cpp(NumericVector mu_hat, NumericVector phi_hat, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  double q = 1-alpha;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  NumericVector phi_q=phi_hat/q;
  
  NumericVector aq=mu_hat * phi_q;
  NumericVector bq=(1-mu_hat) * phi_q;
  
  NumericVector y_star = Rcpp::log(y)-Rcpp::log((1-y));
  NumericVector mu_star = Rcpp::digamma(aq)-Rcpp::digamma(bq);
  NumericVector Tb = pow(link.d_linkmu(mu_hat),-1);
  
  //PONTO ALTERACAO
  //NumericVector f_q_star = Rcpp::pow(degbeta_C(y_star,mu_hat,phi_q),alpha);
  NumericVector f_q_star = Rcpp::pow(degbeta_CR(y,mu_hat,phi_q),alpha);
  
  NumericMatrix Phi_q_Tb_fqstar = Rcpp::diag(phi_q*Tb*f_q_star);
  NumericVector diff = y_star-mu_star;
  arma::mat diff_arma =as<arma::vec>(diff);
  arma::mat Phi_q_Tb_fqstar_arma =as<arma::mat>(Phi_q_Tb_fqstar);
  arma::mat result = X.t()*Phi_q_Tb_fqstar_arma*diff_arma;
  return(result);
}

arma::mat Psi_LSMLE_Gamma_Cpp(NumericVector mu_hat, NumericVector phi_hat, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  double q = 1-alpha;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  NumericVector phi_q=phi_hat/q;
  NumericVector aq=mu_hat * phi_q;
  NumericVector bq=(1-mu_hat) * phi_q;
  NumericVector y_dagger = Rcpp::log((1-y));
  NumericVector y_star = Rcpp::log(y)-y_dagger;
  NumericVector mu_star = Rcpp::digamma(aq)-Rcpp::digamma(bq);
  NumericVector mu_dagger = Rcpp::digamma(bq)-Rcpp::digamma(phi_q);
  
  NumericVector Tg = pow(link.d_linkphi(phi_hat),-1);
  
  //PONTO ALTERACAO
  //NumericVector f_q_star = Rcpp::pow(degbeta_C(y_star,mu_hat,phi_q),alpha);
  NumericVector f_q_star = Rcpp::pow(degbeta_CR(y,mu_hat,phi_q),alpha);
  
  NumericVector eta =  mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger);
  NumericMatrix Tg_fqstar = Rcpp::diag(Tg*f_q_star/q);
  arma::mat eta_arma =as<arma::vec>(eta);
  arma::mat Tg_fqstar_arma =as<arma::mat>(Tg_fqstar);
  arma::mat result = Z.t()*Tg_fqstar_arma*eta_arma;
  return(result);
}

//' @useDynLib RobustBetaReg, .registration=TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
arma::mat Psi_LSMLE_Cpp(arma::vec Theta, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  int k = X.n_cols;
  int m = Z.n_cols;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  arma::vec Beta=Theta.subvec(0,k-1);
  arma::vec Gamma=Theta.subvec(k,k+m-1);
  NumericVector mu_hat = link.inv_linkmu(eta(X,Beta));
  NumericVector phi_hat = link.inv_linkphi(eta(Z,Gamma));

  arma::mat psi_beta = Psi_LSMLE_Beta_Cpp(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi);
  arma::mat psi_gamma = Psi_LSMLE_Gamma_Cpp(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi);
  arma::mat psi = join_cols(psi_beta,psi_gamma);
  return(psi.t());
}

//' @useDynLib RobustBetaReg, .registration=TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
double L_alpha(arma::vec Theta, NumericVector y, arma::mat X, arma::mat Z, double alpha, StringVector link_mu, StringVector link_phi){
  double L_q;
  double q = 1-alpha;
  int k = X.n_cols;
  int m = Z.n_cols;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  arma::vec Beta=Theta.subvec(0,k-1);
  arma::vec Gamma=Theta.subvec(k,k+m-1);
  NumericVector mu_hat = link.inv_linkmu(eta(X,Beta));
  NumericVector phi_hat = link.inv_linkphi(eta(Z,Gamma));
  NumericVector phi_q=phi_hat/q;
  NumericVector f_q_star = degbeta_CR(y,mu_hat,phi_q);
  
  if(alpha==0){
    L_q = sum(Rcpp::log(f_q_star));
  }else{
    L_q = sum((Rcpp::pow(f_q_star,alpha)-1)/alpha);
  }
  return(L_q);
}

// [[Rcpp::export]]
arma::mat Psi_LSMLE_Jacobian_C(arma::vec Theta, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  int k = X.n_cols;
  int m = Z.n_cols;
  double q = 1-alpha;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  arma::vec Beta=Theta.subvec(0,k-1);
  arma::vec Gamma=Theta.subvec(k,k+m-1);
  NumericVector mu_hat, phi_hat, phi_q, aq, bq, y_star, y_dagger, mu_star, mu_dagger;
  mu_hat = link.inv_linkmu(eta(X,Beta));
  phi_hat = link.inv_linkphi(eta(Z,Gamma));
  phi_q = phi_hat/q;
  aq=mu_hat*phi_q;
  bq=(1-mu_hat)*phi_q;
  y_dagger=log(1-y);
  y_star=log(y)-y_dagger;
  mu_star = Rcpp::digamma(aq)-Rcpp::digamma(bq);
  mu_dagger = Rcpp::digamma(bq)-Rcpp::digamma(phi_q);
  
  NumericVector d_linkmu, d2_linkmu, d_linkphi, d2_linkphi; 
  d_linkmu = link.d_linkmu(mu_hat);
  d2_linkmu = link.d2_linkmu(mu_hat);
  d_linkphi = link.d_linkphi(phi_hat);
  d2_linkphi = link.d2_linkphi(phi_hat);
  
  NumericMatrix Tb, Tg;
  Tb = Rcpp::diag(pow(d_linkmu,-1));
  Tg = Rcpp::diag(pow(d_linkphi,-1));
  
  NumericVector u_mu, u_phi, u_mumu, u_muphi,  u_phiphi,  fqstar;
  u_mu = phi_q*(y_star-mu_star);
  u_phi=(mu_hat*(y_star-mu_star)+y_dagger-mu_dagger)/q;
  u_mumu=-pow(phi_q,2)*(trigamma(aq)+trigamma(bq));
  u_muphi=((y_star-mu_star)-phi_q*(mu_hat*trigamma(aq)-(1-mu_hat)*trigamma(bq)))/q;
  u_phiphi=(trigamma(phi_q)-trigamma(aq)*pow(mu_hat,2)-trigamma(bq)*pow(1-mu_hat,2))/pow(q,2);
  fqstar = pow(degbeta_C(y_star,mu_hat,phi_q),alpha);
  //arma::vec u_mu_test =as<arma::vec>(u_phi);
  
  NumericVector core1, core2, core3;
  core1 = fqstar*(-(d2_linkmu/d_linkmu)*u_mu+u_mumu+alpha*pow(u_mu,2));
  core2 = fqstar*(u_muphi+alpha*u_mu*u_phi);
  core3 = fqstar*(-(d2_linkphi/d_linkphi)*u_phi+u_phiphi+alpha*pow(u_phi,2));
  
  arma::vec Core1 =as<arma::vec>(core1);
  arma::vec Core2 =as<arma::vec>(core2);
  arma::vec Core3 =as<arma::vec>(core3);
  arma::mat Tbb = as<arma::mat>(Tb);
  arma::mat Tgg = as<arma::mat>(Tg);
  
  mat J, J11, J12, J21, J22;
  J11 = X.t()*Tbb*diagmat(Core1)*Tbb*X;
  J12 = X.t()*Tbb*diagmat(Core2)*Tgg*Z;
  J22 = Z.t()*Tgg*diagmat(Core3)*Tgg*Z;
  
  J = join_cols(join_rows(J11,J12),join_rows(J12.t(),J22));
  return(J);
}


/*    Saddlepoint        */

// [[Rcpp::export]]
double p_LSMLE_Cpp(double y_star, double mu, double phi, double mu_0, double phi_0, double alpha, double Kx, double Kz)
{
  double q, phi_q, a, b, a0, b0, u1, u2, muA, muB;
  q=1-alpha;
  phi_q=phi/q;
  a=mu*phi_q;
  b=phi_q-a;
  a0=mu_0*phi_0;
  b0=(1-mu_0)*phi_0;  

  muA = R::digamma(a)-R::digamma(b);
  muB = R::digamma(b)-R::digamma(phi_q);
  u1=Kx*phi_q*(y_star-muA);
  u2=Kz*(mu*(y_star-muA)-R::log1pexp(y_star)-muB)/q;
  double f = exp(-(R::lbeta(a,b)+b*y_star+phi_q*R::log1pexp(-y_star))); 
  double f0 = exp(-(R::lbeta(a0,b0)+b0*y_star+phi_0*R::log1pexp(-y_star)));
  double result = (u1+u2)*pow(f,alpha)+log(f0);
  
  return(result);
}

// [[Rcpp::export]]
NumericVector dp_LSMLE_Cpp(NumericVector y_star, double mu, double phi, double mu_0, double phi_0, double alpha, double Kx, double Kz)
{
  double q, phi_q, a, b, a0, b0, u1, u2, du1, du2, muA, muB, f, f0, df, dlnf0;
  q=1-alpha;
  phi_q=phi/q;
  a=mu*phi_q;
  b=phi_q-a;
  a0=mu_0*phi_0;
  b0=(1-mu_0)*phi_0;  
  muA = R::digamma(a)-R::digamma(b);
  muB = R::digamma(b)-R::digamma(phi_q);
  int n = y_star.size();
  NumericVector result (n);
  for(int i=0;i<n;i++)
  {
    u1=phi_q*(y_star[i]-muA);
    u2=(mu*(y_star[i]-muA)-R::log1pexp(y_star[i])-muB)/q;
    du1=phi_q;
    du2=(mu-exp(y_star[i])/(1+exp(y_star[i])))/q;
    f = exp(-(R::lbeta(a,b)+b*y_star[i]+phi_q*R::log1pexp(-y_star[i])));
    df = alpha*phi_q*(mu+(mu-1)*exp(y_star[i]))*pow(f,alpha)/(1+exp(y_star[i]));
    dlnf0 = phi_0*(mu_0+(mu_0-1)*exp(y_star[i]))/(1+exp(y_star[i]));
    result[i] = (Kx*du1+Kz*du2)*pow(f,alpha)+(Kx*u1+Kz*u2)*df+dlnf0;
  }
  return(result);
}

double dp_LSMLE(double y, double mu, double phi, double mu_0, double phi_0, double alpha, double Kx, double Kz)
{
  double q, phi_q, a, b, a0, b0, u1, u2, du1, du2, muA, muB, f, f0, df, dlnf0;
  double result, EXPY;
  q=1-alpha;
  phi_q=phi/q;
  a=mu*phi_q;
  b=phi_q-a;
  a0=mu_0*phi_0;
  b0=(1-mu_0)*phi_0;  
  muA = R::digamma(a)-R::digamma(b);
  muB = R::digamma(b)-R::digamma(phi_q);
  EXPY=exp(y);
  
  u1=phi_q*(y-muA);
  u2=(mu*(y-muA)-R::log1pexp(y)-muB)/q;
  du1=phi_q;
  du2=(mu-EXPY/(1+EXPY))/q;
  f = exp(-(R::lbeta(a,b)+b*y+phi_q*R::log1pexp(-y)));
  df = alpha*phi_q*(mu+(mu-1)*EXPY)*pow(f,alpha)/(1+EXPY);
  dlnf0 = phi_0*(mu_0+(mu_0-1)*EXPY)/(1+EXPY);
  result = (Kx*du1+Kz*du2)*pow(f,alpha)+(Kx*u1+Kz*u2)*df+dlnf0;
  
  return(result);
}

double * DXfVector_Cpp(double y_star, double mu, double phi, double alpha) {
  
  double EGB, AP, AMP, AMPm1, EXPY;
  EGB=exp(-alpha*(R::lbeta(mu*phi,(1-mu)*phi)+(1-mu)*phi*y_star+phi*log(1+exp(-y_star))));
  AP=alpha*phi;
  AMP=AP*mu;
  AMPm1=AP*(mu-1);
  EXPY=exp(y_star);
  
  double df = (AMP+AMPm1*EXPY)*EGB/(1+EXPY);
  double df2 = (EGB/pow(1+EXPY,2))*(pow(AMP,2)+pow(EXPY,2)*pow(AMPm1,2)+AP*EXPY*(2*AMPm1*mu-1));
  double df3 = (EGB/pow(1+EXPY,3))*(pow(AMP,2)*(AMP+EXPY*(AMP-AP-2))+(pow(AMPm1,2))*pow(EXPY,2)*(AMP+AMPm1*EXPY+2)+AP*(2*AMPm1*mu-1)*EXPY*(AMP+EXPY*(AMPm1-1)+1));
  
  double fa = EGB/pow(1+EXPY,4);
  double a1 = pow(AMP,2)*(AMP-AP-2)*EXPY*(AMP+EXPY*(AMP-AP-2)+1);
  double a12 = pow(AMP,3)*(AMP+EXPY*(AMPm1-3));
  double a21 = pow(AMPm1,2)*(AMP+2)*exp(2*y_star)*(AMP+EXPY*(AMPm1-1)+2);
  double a22 = pow(AMPm1*EXPY,3)*(AMP+AMPm1*EXPY+3);
  double a31 = AP*(2*AMPm1*mu-1)*(AMP+1)*EXPY*(AMP+EXPY*(AMPm1-2)+1);
  double a32 = AP*(2*AMPm1*mu-1)*(AMPm1-1)*pow(EXPY,2)*(AMP+EXPY*(AMPm1-1)+2);
  
  double df4 = fa*(a1+a12+a21+a22+a31+a32);
  
  double *result = (double *) calloc(4,sizeof(double));
  result[0]=df;
  result[1]=df2;
  result[2]=df3;
  result[3]=df4;
  return result ;
}

double* DXlnfVector_Cpp(double y_star, double mu, double phi){
  
  double EXPY, AP;
  EXPY = exp(y_star);
  AP = phi;
  
  double dlnf = AP*(mu+(mu-1)*EXPY)/(1+EXPY);
  double d2lnf = -AP*EXPY/pow(1+EXPY,2);
  double d3lnf = AP*EXPY*(EXPY-1)/pow(1+EXPY,3);
  double d4lnf = -AP*EXPY*(-4*EXPY+(pow(EXPY,2))+1)/pow(1+EXPY,4);
  
  //List result = List::create(dlnf,d2lnf, d3lnf, d4lnf);
  double *result = (double *) calloc(4,sizeof(double));
  result[0]=dlnf;
  result[1]=d2lnf;
  result[2]=d3lnf;
  result[3]=d4lnf;
  return result;
}

double DerivativeVector_Cpp(double y_star, double mu, double phi, double mu_0,double phi_0, double alpha, double Kx,double Kz, double p0) {
  
  double phi_q,a,b, muA, muB, EXPY;
  double q=1-alpha;
  double *Df = (double *) calloc(4,sizeof(double));
  double *Dlnf = (double *) calloc(4,sizeof(double));
  
  double ep0=exp(p0);
  
  phi_q=phi/q;
  a=mu*phi_q;
  b=phi_q-a;
  muA = R::digamma(a)-R::digamma(b);
  muB = R::digamma(b)-R::digamma(phi_q);
  EXPY=exp(y_star);
  
  double u1= phi_q*(y_star-muA);
  double u2= (mu*(y_star-muA)-std::log1p(exp(y_star))-muB)/q;
  double du1= phi_q;
  double du2= (mu-EXPY/(1+EXPY))/q;
  double d2u2=-EXPY/(q*pow((1+EXPY),2));
  double d3u2=(1/q)*EXPY*(EXPY-1)/pow(1+EXPY,3);
  double d4u2=(-1/q)*EXPY*(-4*EXPY+pow(EXPY,2)+1)/pow(1+EXPY,4);
  
  double f =  exp(-alpha*(R::lbeta(a,b)+b*y_star+phi_q*log(1+exp(-y_star))));
  Df = DXfVector_Cpp(y_star, mu, phi_q, alpha);
  double df = Df[0];
  double d2f = Df[1];
  double d3f = Df[2];
  double d4f = Df[3];
  
  Dlnf = DXlnfVector_Cpp(y_star, mu_0, phi_0);
  double dlnf = Dlnf[0];
  double d2lnf = Dlnf[1];
  double d3lnf = Dlnf[2];
  double d4lnf = Dlnf[3];
  
  double D2 = (Kz*d2u2)*f+2*(Kx*du1+Kz*du2)*df+(Kx*u1+Kz*u2)*d2f+d2lnf;
  double D3 = (Kz*d3u2)*f+3*(Kz*d2u2)*df+3*(Kx*du1+Kz*du2)*d2f+(Kx*u1+Kz*u2)*d3f+d3lnf;
  double D4 = (Kz*d4u2)*f+4*(Kz*d3u2)*df+6*(Kz*d2u2)*d2f+4*(Kx*du1+Kz*du2)*d3f+(Kx*u1+Kz*u2)*d4f+d4lnf;
  
  double a1=D4/(8*pow(D2,2))-5*(pow(D3,2))/(24*pow(D2,3));
  double result = pow(2*pi/std::abs(D2),0.5)*ep0*(1+a1);
  
  free(Df);
  free(Dlnf);
  
  return result;
}

double NQ_LSMLE(double y, double mu, double phi, double mu_0, double phi_0, double alpha, double Kx, double Kz)
{
  double result, D1, D2, q, phi_q, a, b, a0, b0, u1, u2, du1, du2, d2u2, muA, muB, f, f0, df, d2f, dlnf0, d2lnf0;
  double AP, AMP, AMPm1, EXPY, EGB;
  q=1-alpha;
  phi_q=phi/q;
  a=mu*phi_q;
  b=phi_q-a;
  a0=mu_0*phi_0;
  b0=(1-mu_0)*phi_0;  
  muA = R::digamma(a)-R::digamma(b);
  muB = R::digamma(b)-R::digamma(phi_q);
  
  AP=alpha*phi_q;
  AMP=AP*mu;
  AMPm1=AP*(mu-1);
  EXPY=exp(y);
  EGB=exp(-alpha*(R::lbeta(a,b)+b*y+phi_q*log(1+exp(-y))));
  
  u1=phi_q*(y-muA);
  u2=(mu*(y-muA)-R::log1pexp(y)-muB)/q;
  du1=phi_q;
  du2=(mu-EXPY/(1+EXPY))/q;
  d2u2=-EXPY/(q*pow((1+EXPY),2));
  
  df = AP*(mu+(mu-1)*EXPY)*EGB/(1+EXPY);
  d2f = (EGB/pow(1+EXPY,2))*(pow(AMP,2)+pow(EXPY,2)*pow(AMPm1,2)+AP*EXPY*(2*AMPm1*mu-1));
  dlnf0 = phi_0*(mu_0+(mu_0-1)*EXPY)/(1+EXPY);
  d2lnf0 = -phi_0*EXPY/pow(1+EXPY,2);
  
  D1 = (Kx*du1+Kz*du2)*EGB+(Kx*u1+Kz*u2)*df+dlnf0;
  D2 = (Kz*d2u2)*EGB+2*(Kx*du1+Kz*du2)*df+(Kx*u1+Kz*u2)*d2f+d2lnf0;
  
  result=D1/D2;
  
  return(result);
}

double uniroot_C(double lower, double upper, double mu, double phi, double mu_0, double phi_0, double alpha, double Kx, double Kz)
{
  double y0, y1, dy, NQ;
  int intMax = 200;
  int k=1;
  double tol = 1e-5;//tolerance
  bool run = true;
  
  //IntegerVector inty = seq((int)lower, (int)upper);
  //NumericVector yy = as<NumericVector>(inty);
  
  y0=(upper+lower)/2; //initial guess
  while(run)
    {
      NQ=NQ_LSMLE(y0,mu,phi,mu_0,phi_0,alpha,Kx,Kz);
      y1=y0-NQ;
      if(std::abs(y1-y0)<tol)
        {
          dy=dp_LSMLE(y1,mu,phi,mu_0, phi_0, alpha, Kx,Kz);
          if(std::abs(dy)<tol)
            {
              run=false;
            }
        }
      if(k>=intMax){run=false;}
      y0=y1;
      k++;
    }
  return(y1);
}

// [[Rcpp::export]]
NumericVector La_Cpp(NumericVector mu, NumericVector phi, NumericVector mu_0,NumericVector phi_0, double alpha, NumericVector Kx,NumericVector Kz, int thrd)
{
  int i;
  int size = mu.size();
  double y0, p0;
  NumericVector result (size);
  double *v = (double *) calloc(size,sizeof(double));
  //Inputs
  double *Vmu = (double *) calloc(size,sizeof(double));
  double *Vphi = (double *) calloc(size,sizeof(double));
  double *Vmu_0 = (double *) calloc(size,sizeof(double));
  double *Vphi_0 = (double *) calloc(size,sizeof(double));
  double *VKx = (double *) calloc(size,sizeof(double));
  double *VKz = (double *) calloc(size,sizeof(double));
  //Vector assembling
  for(i=0;i<size;i++)
  {
    Vmu[i]=mu[i];
    Vphi[i]=phi[i];
    Vmu_0[i]=mu_0[i];
    Vphi_0[i]=phi_0[i];
    VKx[i]=Kx[i];
    VKz[i]=Kz[i];
  }
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(size,v,Vmu,Vphi,Vmu_0,Vphi_0,VKx,VKz,alpha) private(i,y0,p0)
{
#pragma omp for
  for(i=0;i<size;i++)
    {
      y0=uniroot_C(-25, 25, Vmu[i], Vphi[i], Vmu_0[i], Vphi_0[i], alpha, VKx[i], VKz[i]);
      p0=p_LSMLE_Cpp(y0, Vmu[i], Vphi[i], Vmu_0[i], Vphi_0[i], alpha, VKx[i], VKz[i]);     
      v[i]=DerivativeVector_Cpp(y0, Vmu[i], Vphi[i],Vmu_0[i],Vphi_0[i],alpha,VKx[i],VKz[i],p0);
    }
} 
for(i=0; i<size;i++)
  {
    result[i]=v[i];
  }
  free (v);
  free (Vmu);
  free (Vphi);
  free (Vmu_0);
  free (Vphi_0);
  free (VKx);
  free (VKz);

  return(result);
}


// [[Rcpp::export]]
NumericVector La_CppB(NumericVector p0,NumericVector y0, NumericVector mu, NumericVector phi, NumericVector mu_0,NumericVector phi_0, double alpha, NumericVector Kx,NumericVector Kz, int thrd)
{
  int i;
  int size = y0.size();
  NumericVector result (size);
  double *v = (double *) calloc(size,sizeof(double));
  //Inputs
  double *Vp0 = (double *) calloc(size,sizeof(double));
  double *Vy0 = (double *) calloc(size,sizeof(double));
  double *Vmu = (double *) calloc(size,sizeof(double));
  double *Vphi = (double *) calloc(size,sizeof(double));
  double *Vmu_0 = (double *) calloc(size,sizeof(double));
  double *Vphi_0 = (double *) calloc(size,sizeof(double));
  double *VKx = (double *) calloc(size,sizeof(double));
  double *VKz = (double *) calloc(size,sizeof(double));
  //Inicializa vetor
  for(i=0;i<40;i++)
  {
    Vp0[i]=p0[i];
    Vy0[i]=y0[i];
    Vmu[i]=mu[i];
    Vphi[i]=phi[i];
    Vmu_0[i]=mu_0[i];
    Vphi_0[i]=phi_0[i];
    VKx[i]=Kx[i];
    VKz[i]=Kz[i];
  }
#if _OPENMP
  omp_set_num_threads(thrd);//define number of threads 
#endif  
#pragma omp parallel shared(size,v,Vy0,Vmu,Vphi,Vmu_0,Vphi_0,VKx,VKz,Vp0,alpha) private(i)
{
#pragma omp for
  for(i=0;i<40;i++)
  {
    v[i]=DerivativeVector_Cpp(Vy0[i], Vmu[i], Vphi[i],Vmu_0[i],Vphi_0[i],alpha,VKx[i],VKz[i],Vp0[i]);
    //v[i]=teste2(Vy0[i], Vmu[i], Vphi[i],Vmu_0[i],Vphi_0[i],alpha,VKx[i],VKz[i],Vp0[i]);
  }
} 
for(i=0; i<40;i++)
{
  result[i]=v[i];
}
free (v);
free (Vp0);
free (Vy0);
free (Vmu);
free (Vphi);
free (Vmu_0);
free (Vphi_0);
free (VKx);
free (VKz);

return(result);
}


/*  LMDPDE - Functions C++  */

arma::mat Psi_LMDPDE_Beta_Cpp(NumericVector mu_hat, NumericVector phi_hat, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  double q = 1-alpha;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  //Vectors
  NumericVector a0=mu_hat*phi_hat;
  NumericVector b0=(1-mu_hat)*phi_hat;
  NumericVector a_alpha=a0*(1+alpha);
  NumericVector b_alpha=b0*(1+alpha);
  
  NumericVector y_star = Rcpp::log(y)-Rcpp::log((1-y));
  NumericVector mu_star = Rcpp::digamma(a0)-Rcpp::digamma(b0);
  NumericVector mu_star_alpha = Rcpp::digamma(a_alpha)-Rcpp::digamma(b_alpha);
  NumericVector Phi_Tb = phi_hat*pow(link.d_linkmu(mu_hat),-1);
  NumericVector diff1 = y_star-mu_star;
  NumericVector diff2 = mu_star_alpha-mu_star;
  
  //Matrixes
  //NumericMatrix Phi_Tb_Walpha = Rcpp::diag(Phi_Tb*Rcpp::pow(degbeta_C(y_star,mu_hat,phi_hat),alpha));
  NumericMatrix Phi_Tb_Walpha = Rcpp::diag(Phi_Tb*Rcpp::pow(degbeta_CR(y,mu_hat,phi_hat),alpha));
  NumericMatrix Phi_Tb_Calpha = Rcpp::diag(Phi_Tb*Rcpp::exp(Rcpp::lbeta(a_alpha,b_alpha)-(1+alpha)*Rcpp::lbeta(a0,b0)));
  //arma Matrixes
  arma::mat Phi_Tb_Walpha_arma =as<arma::mat>(Phi_Tb_Walpha);
  arma::mat Phi_Tb_Calpha_arma =as<arma::mat>(Phi_Tb_Calpha);
  arma::mat D1 =as<arma::vec>(diff1);
  arma::mat D2 =as<arma::vec>(diff2);
  
  arma::mat result = X.t()*(Phi_Tb_Walpha_arma*D1-Phi_Tb_Calpha_arma*D2);
  return(result);
}


arma::mat Psi_LMDPDE_Gamma_Cpp(NumericVector mu_hat, NumericVector phi_hat, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  double q = 1-alpha;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  //Vectors
  NumericVector a0, b0, a_alpha, b_alpha, y_dagger, y_star, mu_star, mu_star_alpha, mu_dagger, mu_dagger_alpha, Tg, eta, E_eta;
  
  a0=mu_hat*phi_hat;
  b0=(1-mu_hat)*phi_hat;
  a_alpha=a0*(1+alpha);
  b_alpha=b0*(1+alpha);
  y_dagger = Rcpp::log((1-y));
  y_star = Rcpp::log(y)-y_dagger;
  mu_star = Rcpp::digamma(a0)-Rcpp::digamma(b0);
  mu_star_alpha = Rcpp::digamma(a_alpha)-Rcpp::digamma(b_alpha);
  mu_dagger = digamma(b0)-digamma(phi_hat);
  mu_dagger_alpha = digamma(b_alpha)-digamma(phi_hat*(1+alpha));
  Tg = pow(link.d_linkphi(phi_hat),-1);
  eta =  mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger);
  E_eta =  mu_hat*(mu_star_alpha-mu_star)+(mu_dagger_alpha-mu_dagger);
  //Matrixes
  //NumericMatrix Tg_Walpha = Rcpp::diag(Tg*Rcpp::pow(degbeta_C(y_star,mu_hat,phi_hat),alpha));
  NumericMatrix Tg_Walpha = Rcpp::diag(Tg*Rcpp::pow(degbeta_CR(y,mu_hat,phi_hat),alpha));
  NumericMatrix Tg_Calpha = Rcpp::diag(Tg*Rcpp::exp(Rcpp::lbeta(a_alpha,b_alpha)-(1+alpha)*Rcpp::lbeta(a0,b0)));
  //arma Matrixes
  arma::mat Tg_Walpha_arma =as<arma::mat>(Tg_Walpha);
  arma::mat Tg_Calpha_arma =as<arma::mat>(Tg_Calpha);
  arma::mat eta_arma =as<arma::vec>(eta);
  arma::mat E_eta_arma =as<arma::vec>(E_eta);
  
  arma::mat result = Z.t()*(Tg_Walpha_arma*eta_arma-Tg_Calpha_arma*E_eta_arma);
  return(result);
}

// [[Rcpp::export]]
arma::mat Psi_LMDPDE_Cpp(arma::vec Theta, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  int k = X.n_cols;
  int m = Z.n_cols;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  arma::vec Beta=Theta.subvec(0,k-1);
  arma::vec Gamma=Theta.subvec(k,k+m-1);
  NumericVector mu_hat = link.inv_linkmu(eta(X,Beta));
  NumericVector phi_hat = link.inv_linkphi(eta(Z,Gamma));
  
  arma::vec psi_beta = Psi_LMDPDE_Beta_Cpp(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi);
  arma::vec psi_gamma = Psi_LMDPDE_Gamma_Cpp(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi);
  arma::mat psi = join_cols(psi_beta,psi_gamma);
  return(psi.t());
}


//' @useDynLib RobustBetaReg, .registration=TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
double D_alpha(arma::vec Theta, NumericVector y, arma::mat X, arma::mat Z, double alpha, StringVector link_mu, StringVector link_phi){
  double D_q;
  double q = 1-alpha;
  int k = X.n_cols;
  int m = Z.n_cols;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  arma::vec Beta=Theta.subvec(0,k-1);
  arma::vec Gamma=Theta.subvec(k,k+m-1);
  NumericVector mu_hat = link.inv_linkmu(eta(X,Beta));
  NumericVector phi_hat = link.inv_linkphi(eta(Z,Gamma));
  
  if(alpha==0){
    D_q = sum(Rcpp::log(degbeta_CR(y,mu_hat,phi_hat)));
  }else{
    NumericVector a0 = mu_hat*phi_hat;
    NumericVector b0 = (1-mu_hat)*phi_hat;
    NumericVector a_alpha = a0*(1+alpha);
    NumericVector b_alpha = b0*(1+alpha);
    NumericVector E_alpha = Rcpp::exp(Rcpp::lbeta(a_alpha,b_alpha)-(1+alpha)*Rcpp::lbeta(a0,b0));
    D_q=sum((1+alpha)/(alpha)*Rcpp::pow(degbeta_CR(y,mu_hat,phi_hat),alpha)-E_alpha);
  }

  return(D_q);
}



// [[Rcpp::export]]
arma::mat Psi_LMDPDE_Jacobian_C(arma::vec Theta, NumericVector y, arma::mat X,arma::mat Z, double alpha, StringVector link_mu,StringVector link_phi)
{
  int k = X.n_cols;
  int m = Z.n_cols;
  Link link;
  link.setLinkMu(link_mu);
  link.setLinkPhi(link_phi);
  arma::vec Beta=Theta.subvec(0,k-1);
  arma::vec Gamma=Theta.subvec(k,k+m-1);
  NumericVector mu_hat, phi_hat, a0, b0, a_alpha, b_alpha, y_star, y_dagger, mu_star, mu_star_alpha, mu_dagger, mu_dagger_alpha;
  mu_hat = link.inv_linkmu(eta(X,Beta));
  phi_hat = link.inv_linkphi(eta(Z,Gamma));

  a0=mu_hat*phi_hat;
  b0=(1-mu_hat)*phi_hat;
  a_alpha=a0*(1+alpha);
  b_alpha=b0*(1+alpha);
  y_dagger=log(1-y);
  y_star=log(y)-y_dagger;
  mu_star = Rcpp::digamma(a0)-Rcpp::digamma(b0);
  mu_star_alpha = Rcpp::digamma(a_alpha)-Rcpp::digamma(b_alpha);
  mu_dagger = Rcpp::digamma(b0)-Rcpp::digamma(phi_hat);
  mu_dagger_alpha = Rcpp::digamma(b_alpha)-Rcpp::digamma(phi_hat*(1+alpha));
  NumericVector d_linkmu, d2_linkmu, d_linkphi, d2_linkphi; 
  d_linkmu = link.d_linkmu(mu_hat);
  d2_linkmu = link.d2_linkmu(mu_hat);
  d_linkphi = link.d_linkphi(phi_hat);
  d2_linkphi = link.d2_linkphi(phi_hat);
  
  NumericMatrix Tb, Tg;
  Tb = Rcpp::diag(pow(d_linkmu,-1));
  Tg = Rcpp::diag(pow(d_linkphi,-1));
  
  NumericVector core1, c_star,k_mumu, k_mu, u_mu, u_phi, u_mumu,f_q_star;
  c_star=exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0));
  k_mumu=pow(phi_hat,2)*c_star*((1+alpha)*(pow(mu_star_alpha-mu_star,2)+trigamma(a_alpha)+trigamma(b_alpha))-trigamma(a0)-trigamma(b0));
  k_mu = phi_hat*c_star*(mu_star_alpha-mu_star);
  u_mu = phi_hat*(y_star-mu_star); 
  u_phi = mu_hat*(y_star-mu_star)+y_dagger-mu_dagger;
  u_mumu = -pow(phi_hat,2)*(trigamma(a0)+trigamma(b0));
  f_q_star=pow(degbeta_C(y_star,mu_hat,phi_hat),alpha);
    
  core1 = f_q_star*(u_mumu+alpha*pow(u_mu,2))-(k_mumu+(d2_linkmu/d_linkmu)*(f_q_star*u_mu-k_mu));  
    
  NumericVector core2, D1, dd, k_muphi, u_muphi; 
  D1 = mu_hat*(mu_star_alpha-mu_star)+(mu_dagger_alpha-mu_dagger);
  dd = (1+alpha)*D1*(mu_star_alpha-mu_star)+((1+alpha)*(mu_hat*trigamma(a_alpha)-(1-mu_hat)*trigamma(b_alpha))-(mu_hat*trigamma(a0)-(1-mu_hat)*trigamma(b0)));
  k_muphi = c_star*(mu_star_alpha-mu_star+phi_hat*dd);
  u_muphi = (y_star-mu_star)-phi_hat*(mu_hat*trigamma(a0)-(1-mu_hat)*trigamma(b0));
  
  core2=f_q_star*(u_muphi+alpha*u_mu*u_phi)-k_muphi;
  
  NumericVector core3, nu_alpha, nu_0, u_phiphi, k_phi, k_phiphi;
  
  nu_alpha=trigamma(a_alpha)*pow(mu_hat,2)+trigamma(b_alpha)*pow(1-mu_hat,2)-trigamma(phi_hat*(1+alpha));
  nu_0=trigamma(a0)*pow(mu_hat,2)+trigamma(b0)*pow(1-mu_hat,2)-trigamma(phi_hat);
  u_phiphi=trigamma(phi_hat)-(trigamma(a0)*pow(mu_hat,2)+trigamma(b0)*pow(1-mu_hat,2));
  k_phi=c_star*D1;
  k_phiphi = c_star*((1+alpha)*(pow(mu_hat*(mu_star_alpha-mu_star)+mu_dagger_alpha-mu_dagger,2)+nu_alpha)-nu_0);
  
  core3=f_q_star*(u_phiphi+alpha*pow(u_phi,2))-(d2_linkphi/d_linkphi)*(u_phi*f_q_star-k_phi)-k_phiphi;
  
  arma::vec Core1 =as<arma::vec>(core1);
  arma::vec Core2 =as<arma::vec>(core2);
  arma::vec Core3 =as<arma::vec>(core3);
  arma::mat Tbb = as<arma::mat>(Tb);
  arma::mat Tgg = as<arma::mat>(Tg);
  
  mat J, J11, J12, J21, J22;
  J11 = X.t()*Tbb*diagmat(Core1)*Tbb*X;
  J12 = X.t()*Tbb*diagmat(Core2)*Tgg*Z;
  J22 = Z.t()*Tgg*diagmat(Core3)*Tgg*Z;
  
  J = join_cols(join_rows(J11,J12),join_rows(J12.t(),J22));
  
  return(J);
}

// [[Rcpp::export]]
NumericVector SQV_Cpp(arma::mat zq, double n, double p)
{
    arma::mat temp =  sqrt(sum(pow(diff(zq,1,0),2),1))/(sqrt(n)*p);
    NumericVector result = wrap(temp.t()); 
    return(result);
}


/*    Test Functions    */

//' @useDynLib RobustBetaReg, .registration=TRUE
//' @importFrom Rcpp evalCpp
// [[Rcpp::export]]
int OpenMPTest(int numThrd)
{
  int j=0;
#if _OPENMP
  omp_set_num_threads(numThrd);//define number of threads 
#endif   
  #pragma omp parallel firstprivate(j)
  {
    j=omp_get_thread_num();
    Rprintf("The id of this Thread is %i\n", j);
  }
//Rcout << "O id da Thread j eh: " << j << "\n";
return(j);
}