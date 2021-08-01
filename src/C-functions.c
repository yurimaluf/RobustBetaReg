#include <R.h>
#include <Rinternals.h>
#include <omp.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else 
#define omp_get_thread_num() 0 
#endif 

const double pi = 3.1415926535;


/* Mu:  logit link function  */
SEXP C_logit_link(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] = log(mu)-log(1-mu);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
  {
  REAL(out)[j]=v[j];
  }
free(v);
UNPROTECT(1);
return out;
}

SEXP C_logit_linkinv(SEXP eta_ , SEXP thrd_) {
  
  int n = length(eta_);
  double eta, temp;
  int i; 
  double *v = (double *) calloc(n,sizeof(double));
  int thrd=asReal(thrd_);
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,eta,temp)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    eta = REAL(eta_)[i];
    temp =exp(eta-log1p(exp(eta)));
    temp=temp>1-2.22e-16 ? 1-2.22e-16 : temp;
    temp=temp<2.22e-16 ? 2.22e-16 : temp;
   v[i]=temp;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_logit_deta_dmu(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for(i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =pow(mu*(1-mu),-1);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_logit_d2eta_dmu2(SEXP mu_, SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =(2*mu-1)/(pow(mu*(1-mu),2));
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}


/* Mu: Complement log-log link function  */

SEXP C_cloglog_link(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] = log(-log(1-mu)) ;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_cloglog_linkinv(SEXP eta_ , SEXP thrd_) {
  int n = length(eta_);
  double eta, temp;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,eta,temp)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    eta = REAL(eta_)[i];
    temp =1-exp(-exp(-eta));
    temp=temp>1-2.22e-16 ? 1-2.22e-16 : temp;
    temp=temp<2.22e-16 ? 2.22e-16 : temp;
    v[i] =temp;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_cloglog_deta_dmu(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =pow(-(1-mu)*log(1-mu),-1);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_cloglog_d2eta_dmu2(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =-(log(1-mu)+1)/pow((1-mu)*log(1-mu),2);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

/* Mu: cauchit link function */

SEXP C_cauchit_link(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] = tan(pi*(mu-0.5)) ;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
} 

SEXP C_cauchit_linkinv(SEXP eta_ , SEXP thrd_) {
  int n = length(eta_);
  double eta, temp;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,eta,temp)
{
#pragma omp for
  for(i = 0; i < n; i++) {
    eta = REAL(eta_)[i];
    temp =0.5+atan(eta)/pi;
    temp=temp>1-2.22e-16 ? 1-2.22e-16 : temp;
    temp=temp<2.22e-16 ? 2.22e-16 : temp;
    v[i] =temp;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_cauchit_deta_dmu(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for(i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =pi*pow(cos(pi*(mu-0.5)),-2);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_cauchit_d2eta_dmu2(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =2*pi*tan(pi*(mu-0.5))*pow(cos(pi*(mu-0.5)),-2);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
} 

/* Mu: loglog link function */

SEXP C_loglog_link(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] = -log(-log(mu)) ;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_loglog_linkinv(SEXP eta_ , SEXP thrd_) {
  int n = length(eta_);
  double eta;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,eta)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    eta=REAL(eta_)[i];
    v[i] =exp(-exp(-eta));
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_loglog_deta_dmu(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =-pow(mu*log(mu),-1);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
}

SEXP C_loglog_d2eta_dmu2(SEXP mu_ , SEXP thrd_) {
  int n = length(mu_);
  double mu;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,mu)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    mu = REAL(mu_)[i];
    mu=mu>1-2.22e-16 ? 1-2.22e-16 : mu;
    mu=mu<2.22e-16 ? 2.22e-16 : mu;
    v[i] =(log(mu)+1)/pow(mu*log(mu),2);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
free(v);
UNPROTECT(1);
return out;
} 

/* Phi: log link function  */
SEXP C_log_link(SEXP phi_ , SEXP thrd_) {
  int n = length(phi_);
  int i; 
  int thrd=asReal(thrd_);
  double phi;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,phi)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    phi = REAL(phi_)[i];
    v[i] = log(phi);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
  {
  REAL(out)[j]=v[j];
  }
UNPROTECT(1);
free(v);
return out;
}

SEXP C_log_linkinv(SEXP eta_ , SEXP thrd_) {
  int n = length(eta_);
  int i; 
  int thrd=asReal(thrd_);
  double eta;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(n) private(i,eta)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    eta =REAL(eta_)[i];
   v[i] =exp(eta);
  } 
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
}

SEXP C_log_deta_dphi(SEXP phi_ , SEXP thrd_) {
  int n = length(phi_);
  int i; 
  int thrd=asReal(thrd_);
  double phi;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(n) private(i,phi)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    phi = REAL(phi_)[i];
    v[i] =pow(phi,-1);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
}

SEXP C_log_d2eta_dphi2(SEXP phi_ , SEXP thrd_) {
  int n = length(phi_);
  int i; 
  int thrd=asReal(thrd_);
  double phi;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(n) private(i,phi)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    phi = REAL(phi_)[i];
    v[i] =-pow(phi,-2);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
} 

/* Phi: sqrt link function  */
SEXP C_sqrt_link(SEXP phi_ , SEXP thrd_) {
  int n = length(phi_);
  int i; 
  int thrd=asReal(thrd_);
  double phi;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,phi)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    phi = REAL(phi_)[i];
    v[i] = sqrt(phi) ;
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
}

SEXP C_sqrt_linkinv(SEXP eta_ , SEXP thrd_) {
  int n = length(eta_);
  int i; 
  int thrd=asReal(thrd_);
  double eta;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,eta)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    eta = REAL(eta_)[i];
    v[i] =pow(eta,2);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
}

SEXP C_sqrt_deta_dphi(SEXP phi_ , SEXP thrd_) {
  int n = length(phi_);
  int i; 
  int thrd=asReal(thrd_);
  double phi;
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,phi)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    phi = REAL(phi_)[i];
    v[i] =pow(2*sqrt(phi),-1);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
}

SEXP C_sqrt_d2eta_dphi2(SEXP phi_ , SEXP thrd_) {
  int n = length(phi_);
  double phi;
  int i; 
  int thrd=asReal(thrd_);
  double *v = (double *) calloc(n,sizeof(double));
#if _OPENMP
  omp_set_num_threads(thrd); 
#endif  
#pragma omp parallel shared(v,n) private(i,phi)
{
#pragma omp for
  for (i = 0; i < n; i++) {
    phi=REAL(phi_)[i];
    v[i] =-0.25*pow(sqrt(phi),-3);
  }
}
SEXP out = PROTECT(allocVector(REALSXP, n));
for(int j=0;j<n;j++)
{
  REAL(out)[j]=v[j];
}
UNPROTECT(1);
free(v);
return out;
} 


/* Test Function */
SEXP testeC(SEXP x_, SEXP y_ , SEXP thrd_) {
  
  double x = asReal(x_);
  double y = asReal(y_);
  
  double sum = x + y;
  
  return ScalarReal(sum);
}
