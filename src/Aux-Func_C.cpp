// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  
using namespace Rcpp;
#include <math.h>
#include <vector>
#include "Classes_C.h"


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
  NumericVector f_q_star = Rcpp::pow(degbeta_C(y_star,mu_hat,phi_q),alpha);
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
  NumericVector f_q_star = Rcpp::pow(degbeta_C(y_star,mu_hat,phi_q),alpha);
  NumericVector eta =  mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger);
  NumericMatrix Tg_fqstar = Rcpp::diag(Tg*f_q_star/q);
  arma::mat eta_arma =as<arma::vec>(eta);
  arma::mat Tg_fqstar_arma =as<arma::mat>(Tg_fqstar);
  arma::mat result = Z.t()*Tg_fqstar_arma*eta_arma;
  return(result);
}

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
  
  arma::vec psi_beta = Psi_LSMLE_Beta_Cpp(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi);
  arma::vec psi_gamma = Psi_LSMLE_Gamma_Cpp(mu_hat,phi_hat,y,X,Z,alpha,link_mu,link_phi);
  arma::mat psi = join_cols(psi_beta,psi_gamma);
  return(psi.t());
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

  //delete[] Tb; 
  
  return(J);
}


// [[Rcpp::export]]
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
  NumericMatrix Phi_Tb_Walpha = Rcpp::diag(Phi_Tb*Rcpp::pow(degbeta_C(y_star,mu_hat,phi_hat),alpha));
  NumericMatrix Phi_Tb_Calpha = Rcpp::diag(Phi_Tb*Rcpp::exp(Rcpp::lbeta(a_alpha,b_alpha)-(1+alpha)*Rcpp::lbeta(a0,b0)));
  //arma Matrixes
  arma::mat Phi_Tb_Walpha_arma =as<arma::mat>(Phi_Tb_Walpha);
  arma::mat Phi_Tb_Calpha_arma =as<arma::mat>(Phi_Tb_Calpha);
  arma::mat D1 =as<arma::vec>(diff1);
  arma::mat D2 =as<arma::vec>(diff2);
  
  arma::mat result = X.t()*(Phi_Tb_Walpha_arma*D1-Phi_Tb_Calpha_arma*D2);
  return(result);
}

// [[Rcpp::export]]
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
  NumericMatrix Tg_Walpha = Rcpp::diag(Tg*Rcpp::pow(degbeta_C(y_star,mu_hat,phi_hat),alpha));
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

