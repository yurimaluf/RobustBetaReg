// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;
#include <math.h>
#include <vector>

const double pi = 3.141593;

NumericVector ddnorm(NumericVector x)
{
  return(-x*exp(-pow(x,2)/2)/sqrt(2*pi));
}

NumericVector eta(arma::mat x,arma::mat theta) 
{
  NumericVector result =wrap(x*theta);
  return(result);
}

NumericVector degbeta_C(NumericVector y_star, NumericVector mu, NumericVector phi)
{
  NumericVector a0=mu*phi;
  NumericVector b0=(1-mu)*phi;
  NumericVector result = Rcpp::exp(-(Rcpp::lbeta(a0,b0)+b0*y_star+phi*Rcpp::log(1+Rcpp::exp(-y_star))));
  return(result);
}


class Link{
public:
  Rcpp::StringVector LinkMu; // Name of mean link function
  Rcpp::StringVector LinkPhi;// Name of precision link function
  NumericVector linkmu(NumericVector mu);// Returns mean link function
  NumericVector d_linkmu(NumericVector mu);// Returns derivative of mean link function
  NumericVector d2_linkmu(NumericVector mu);// Returns second derivative of mean link function
  NumericVector inv_linkmu(NumericVector eta);// Returns inverse of mean link function
  NumericVector linkphi(NumericVector phi);// Returns precision link function
  NumericVector d_linkphi(NumericVector phi);// Returns derivative of precision link function
  NumericVector d2_linkphi(NumericVector phi);// Returns second derivative of precision link function
  NumericVector inv_linkphi(NumericVector eta);// Returns inverse of precision link function
  
  void setLinkMu( Rcpp::StringVector Nome );
  void setLinkPhi( Rcpp::StringVector Nome );
};

void Link::setLinkMu( Rcpp::StringVector Nome)
{
  LinkMu=Nome[0];
}

void Link::setLinkPhi( Rcpp::StringVector Nome)
{
  LinkPhi=Nome[0];
}

NumericVector Link::linkmu(NumericVector mu)
{
  mu=pmax(pmin(mu,1-2.22e-16),2.22e-16);
  if(LinkMu[0]=="logit")
  {
    return(log(mu)-log(1-mu));
  }
  if(LinkMu[0]=="probit")
  {
    return(Rcpp::qnorm(mu,0,1,true,false));
  }
  if(LinkMu[0]=="cloglog")
  {
    return(log(-log(1-mu)));
  }
  if(LinkMu[0]=="cauchit")
  {
    return(tan(pi*(mu-0.5)));
  }
  if(LinkMu[0]=="loglog")
  {
    return(-log(-log(mu)));
  }
  return(0);
}

NumericVector Link::d_linkmu(NumericVector mu)
{
  mu=pmax(pmin(mu,1-2.22e-16),2.22e-16);
  if(LinkMu[0]=="logit")
  {
    return(pow(mu*(1-mu),-1));
  }
  if(LinkMu[0]=="probit")
  {
    return(Rcpp::dnorm(Rcpp::qnorm(mu,0,1,true,false),0,1,false));
  }
  if(LinkMu[0]=="cloglog")
  {
    return(pow(-(1-mu)*log(1-mu),-1));
  }
  if(LinkMu[0]=="cauchit")
  {
    return(pi*pow(cos(pi*(mu-0.5)),-2));
  }
  if(LinkMu[0]=="loglog")
  {
    return(pow(mu*log(mu),-1));
  }
  return(0);  
}

NumericVector Link::d2_linkmu(NumericVector mu)
{
  mu=pmax(pmin(mu,1-2.22e-16),2.22e-16);
  if(LinkMu[0]=="logit")
  {
    return((2*mu-1)/(pow(mu*(1-mu),2)));
  }
  if(LinkMu[0]=="probit")
  {
    return(-ddnorm(Rcpp::qnorm(mu,0,1,true,false))/pow(Rcpp::dnorm(Rcpp::qnorm(mu,0,1,true,false),0,1,false),3));
  }
  if(LinkMu[0]=="cloglog")
  {
    return(-(log(1-mu)+1)/pow((1-mu)*log(1-mu),2));
  }
  if(LinkMu[0]=="cauchit")
  {
    return(2*pi*tan(pi*(mu-0.5))*pow(cos(pi*(mu-0.5)),-2));
  }
  if(LinkMu[0]=="loglog")
  {
    return((log(mu)+1)/pow(mu*log(mu),2));
  }
  return(0);  
}

NumericVector Link::inv_linkmu(NumericVector eta)
{
  if(LinkMu[0]=="logit")
  {
    return(pmax(pmin(exp(eta-log1p(exp(eta))),1-2.22e-16),2.22e-16));
    //return(exp(eta-log1p(exp(eta))));
  }
  if(LinkMu[0]=="probit")
  {
    return(pmax(pmin(Rcpp::pnorm(eta,0,1,true,false),1-2.22e-16),2.22e-16));
  }
  if(LinkMu[0]=="cloglog")
  {
    return(pmax(pmin(1-exp(-exp(-eta)),1-2.22e-16),2.22e-16));
  }
  if(LinkMu[0]=="cauchit")
  {
    return(pmax(pmin(0.5+atan(eta)/pi,1-2.22e-16),2.22e-16));
  }
  if(LinkMu[0]=="loglog")
  {
    return(pmax(pmin(exp(-exp(-eta)),1-2.22e-16),2.22e-16));
  }
  return(0);
}

NumericVector Link::linkphi(NumericVector phi)
{
  if(LinkPhi[0]=="log")
  {
    return(log(phi));
  }
  if(LinkPhi[0]=="identity")
  {
    return(phi);
  }
  if(LinkPhi[0]=="sqrt")
  {
    return(sqrt(phi));
  }
  return(0);
}

NumericVector Link::d_linkphi(NumericVector phi)
{
  if(LinkPhi[0]=="log")
  {
    return(pow(phi,-1));
  }
  if(LinkPhi[0]=="identity")
  {
    NumericVector v(phi.size(),1);
    return(v);
  }
  if(LinkPhi[0]=="sqrt")
  {
    return(pow(2*sqrt(phi),-1));
  }
  return(0);
}

NumericVector Link::d2_linkphi(NumericVector phi)
{
  if(LinkPhi[0]=="log")
  {
    return(-pow(phi,-2));
  }
  if(LinkPhi[0]=="identity")
  {
    NumericVector v(phi.size(),0);
    return(v);
  }
  if(LinkPhi[0]=="sqrt")
  {
    return(-0.25*pow(phi,3/2));
  }
  return(0);
}

NumericVector Link::inv_linkphi(NumericVector eta)
{
  if(LinkPhi[0]=="log")
  {
    return(exp(eta));
  }
  if(LinkPhi[0]=="identity")
  {
    return(eta);
  }
  if(LinkPhi[0]=="sqrt")
  {
    return(pow(eta,2));
  }
  return(0);
}

