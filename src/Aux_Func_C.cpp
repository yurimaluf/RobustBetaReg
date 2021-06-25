// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;
#include <math.h>
#include <vector>

const double pi = 3.141593;

// [[Rcpp::export]]
List rcpp_hello_world() {
  
  CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
  NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
  List z            = List::create( x, y ) ;
  
  return z ;
}

// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}

// [[Rcpp::export]]
int times4(int x) {
  return x * 4;
}

// [[Rcpp::export]]
arma::mat a19(arma::mat x) {
  return( x.t() * x ) ;
}

// [[Rcpp::export]]
arma::mat a20(arma::mat x) {
  return( inv(x.t() * x) ) ;
}

// [[Rcpp::export]]
NumericVector eta(arma::mat x,arma::mat theta) {
  NumericVector result =wrap(x*theta);
  return(result);
}

// [[Rcpp::export]]
NumericVector degbeta_C(NumericVector y_star, NumericVector mu, NumericVector phi)
{
  NumericVector a0=mu*phi;
  NumericVector b0=(1-mu)*phi;
  NumericVector result = Rcpp::exp(-(Rcpp::lbeta(a0,b0)+b0*y_star+phi*Rcpp::log(1+Rcpp::exp(-y_star))));
  return(result);
}
