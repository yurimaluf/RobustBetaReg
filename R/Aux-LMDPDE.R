#' Score Function LMDPDE - Beta
#' 
#' Modified Score Function LMDPDE - Beta
#' 
#' For more details see: 
#' 
#' @param Beta ...
#' @param Gamma ...
#' @param y ...
#' @param X ...
#' @param Z ...
#' @param alpha ...
#' 
#' @return Return value ...
#'
Psi_Beta=function(Beta,Gamma,y,X,Z,alpha)
{
  limit=1e-10
  n=length(y)
  mu_hat=h1(X,Beta)
  phi_hat=h2(Z,Gamma)
  #options(digits=7)
  
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)
  b_alpha=b0*(1+alpha)
  
  y_star=log(y/(1-y))
  mu_star=suppressWarnings(digamma(a0)-digamma(b0))
  mu_star_alpha=suppressWarnings(digamma(a_alpha)-digamma(b_alpha))
  C_alpha=diag(exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0)))
  Tb=diag(x=sapply(mu_hat,function(k) (k-k^2)))
  W_alpha=diag(degbeta(y_star,mu_hat,phi_hat)^(alpha))
  Phi=diag(phi_hat)
  
  return(t(X)%*%Phi%*%Tb%*%(W_alpha%*%(y_star-mu_star)-C_alpha%*%(mu_star_alpha-mu_star)))
}

#' Score Function LMDPDE - Gamma
#' 
#' Modified Score Function LMDPDE - Gamma
#' 
#' For more details see: 
#' 
#' @param Beta ...
#' @param Gamma ...
#' @param y ...
#' @param X ...
#' @param Z ...
#' @param alpha ...
#' 
#' @return Return value ...
#'
Psi_Gamma=function(Beta,Gamma,y,X,Z,alpha)
{
  n=length(y)
  mu_hat=h1(X,Beta)
  phi_hat=h2(Z,Gamma)
  
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)
  b_alpha=b0*(1+alpha)
  
  y_star=log(y/(1-y))
  y_dagger=log(1-y)
  mu_dagger=suppressWarnings(digamma(b0)-digamma(phi_hat))
  mu_star=suppressWarnings(digamma(a0)-digamma(b0))
  mu_star_alpha=suppressWarnings(digamma(a_alpha)-digamma(b_alpha))
  kappa_alpha=suppressWarnings(mu_hat*(mu_star_alpha-mu_star)+digamma(b_alpha)-digamma(b0)+digamma(phi_hat)-digamma(phi_hat*(1+alpha)))
  eta=mu_hat*(y_star-mu_star)+y_dagger-mu_dagger
  C_alpha=diag(exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0)))
  Tg=diag(phi_hat)
  W_alpha=diag(degbeta(y_star,mu_hat,phi_hat)^(alpha))
  
  return(t(Z)%*%Tg%*%(W_alpha%*%eta-C_alpha%*%kappa_alpha))
}


#' Modified Score Vector - LMDPDE
#' 
#' Modified Score Function LMDPDE - Gamma
#' 
#' For more details see: 
#' 
#' @param Theta ...
#' @param y ...
#' @param X ...
#' @param Z ...
#' @param alpha ...
#' 
#' @return Return value ...
#'
Psi_LMDPDE=function(Theta,y,X,Z,alpha)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  Beta=Theta[1:ncol(X)]
  Gamma=Theta[1:ncol(Z)+ncol(X)]
  
  psi_beta=Psi_Beta(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha)
  psi_gamma=Psi_Gamma(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha)
  
  return(c(psi_beta,psi_gamma))
}

#' Robust Point Estimation - LMDPDE
#' 
#' Robust point estimation by estimating equation via LSMLE
#' 
#' For more details see: 
#' 
#' @param Theta ...
#' @param y ...
#' @param x ...
#' @param z ...
#' @param alpha ...
#' 
#' @return Return value ...
#'
Robst.LMDPDE.Beta.Reg=function(y,x,z,start_theta,alpha,tolerance,maxit)
{
  theta=list()
  if(missing(tolerance)){tolerance=1e-3}
  if(missing(maxit)){maxit=150}
  if(missing(start_theta)){
    if(dim(z)[2]==1)
    {
      mle=suppressWarnings(betareg(y~x[,-1]|1))  
    }
    else{mle=suppressWarnings(betareg(y~x[,-1]|z[,-1]))}
    start_theta=as.numeric(c(mle$coefficients$mean,mle$coefficients$precision))
  }
  theta$x=rep(0,length(start_theta))
  theta$fvec=10
  theta$msg=theta$error=NULL
  #browser()
  
  #options(warn = 2) #Convert warnings in errors 
  # theta=tryCatch(Newton.Raphson(start_theta,FUN = Psi_LMDPDE,alpha=alpha,y=y,X=x,Z=z,details = T,tol=tolerance,M=maxit),error=function(e){
  #   theta$msg<-e$message
  #   theta$error=T
  #   return(theta)})
  # theta$x=theta$sol
  #theta=within(theta,rm(sol))
  
  theta=tryCatch(nleqslv(start_theta,Psi_LMDPDE,y=y,X=x,Z=z,alpha=alpha,control=list(ftol=tolerance,maxit=maxit),jacobian=TRUE,method="Newton"),error=function(e){
    theta$msg<-e$message
    return(theta)})
  theta$converged=F
  if(all(abs(theta$fvec)<tolerance) & !all(theta$fvec==0)){theta$converged=T}
  return(theta)
}


#' Auto Selecting tuning parameter algorithm
#' 
#' The algorithm of auto-selecting tuning.
#' 
#' For more details see: 
#' 
#' @param y The numeric response vector (with values in (0,1)).
#' @param x The numeric regressor matrix for mean model.
#' @param z The numeric regressor matrix for precision model, defaulting to an intercept only.
#' @param L A parameter of auto selecting algorithm of tuning parameter (default L=0.02).
#' @param M A integer parameter value of auto selecting algorithm of tuning parameter (default M=3).
#' @param tolerance The function value tolerance.
#' 
#' @return Return a tuning given by data driven algorithm selection
#'
Opt.Tuning.LMDPDE=function(y,x,z,control)
{
  #browser()
  if(missing(control)){control=robustbetareg.control()}
  LMDPDE.list=LMDPDE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.5,0.02)
  K=length(alpha_tuning)
  M=control$M
  L=control$L
  n=length(y)
  unstable=F
  sqv.unstable=T
  ponto.inicial.robst=ponto.inicial.temp=Initial.points(y,x,z)
  p=length(ponto.inicial.robst)
  if(dim(z)[2]!=1)
  {
    est.log.lik=tryCatch(suppressWarnings(betareg(y~x[,-1]|z[,-1])),error=function(e) NULL)
  }else
  {
    est.log.lik=tryCatch(suppressWarnings(betareg(y~x[,-1]|1)),error=function(e) NULL)
  }
  Est.param=as.numeric(c(est.log.lik$coefficients$mean,est.log.lik$coefficients$precision))
  for(k in 1:(M+1))
  {
    LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=Est.param,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(!LMDPDE.par$converged)
    {
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.temp,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    }
    if(!LMDPDE.par$converged)
    {
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.robst,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    }
    if(LMDPDE.par$converged)
    {
      ponto.inicial.temp=Par.q(LMDPDE.par)
    }
    if(any(is.na(Z.q(LMDPDE.par))) || is.null(SE.q(LMDPDE.par)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LMDPDE.list[[k]]<-LMDPDE.par
    zq.t<-unname(rbind(zq.t,Z.q(LMDPDE.par)))
  }
  sqv=SQV(zq.t,n,p)
  if(all(sqv<=L))
  {
    LMDPDE.par.star<-LMDPDE.list[[1]]
    LMDPDE.par.star$sqv=sqv
    return(LMDPDE.par.star)
  }
  k=k+1
  while(sqv.unstable)
  {
    LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=Est.param,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(!LMDPDE.par$converged)
    {
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.temp,control = control),error=function(e) {LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    }
    if(!LMDPDE.par$converged)
    {
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.robst,control = control),error=function(e) {LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    }
    if(LMDPDE.par$converged)
    {
      ponto.inicial.temp=Par.q(LMDPDE.par)
    }
    if(any(is.na(Z.q(LMDPDE.par))) || is.null(SE.q(LMDPDE.par)))
    {
      unstable=T
      break
    }
    LMDPDE.list[[k]]=LMDPDE.par
    zq.t=unname(rbind(zq.t,Z.q(LMDPDE.par)))
    sqv=SQV(zq.t,n,p)
    sqv.test=sqv[(k-M):(k-1)]
    if(all(sqv.test<=L) || k==K )
    {
      sqv.unstable=F
      #break
    }
    k=k+1
  }
  if(k>=K || unstable)
  {
    LMDPDE.par.star=LMDPDE.list[[1]]
    LMDPDE.par.star$sqv=sqv
    LMDPDE.par.star$message="Lack of stability"
  }else{
    LMDPDE.par.star=LMDPDE.list[[(k-M)]]
    LMDPDE.par.star$sqv=sqv
  }
  return(LMDPDE.par.star)
}


#' Sandwich Matrix - LMDPDE
#' 
#' Covariance matrix of LSMDPDE.
#' 
#' For more details see: 
#' 
#' @param mu The numeric mean parameter.
#' @param phi The numeric precision parameter.
#' @param X The numeric regressor matrix for mean model.
#' @param Z The numeric regressor matrix for precision model, defaulting to an intercept only.
#' @param alpha The tuning with values (0,1), for robust estimation. When alpha is equal to zero, it is equivalent to MLE. 
#' 
#' @return Return a tuning given by data driven algorithm selection
#'
LMDPDE_Cov_Matrix=function(mu,phi,X,Z,alpha)
{
  n=length(mu)
  a0=mu*phi
  b0=(1-mu)*phi
  a_alpha=(1+alpha)*a0
  b_alpha=(1+alpha)*b0
  a_2alpha=(1+2*alpha)*a0
  b_2alpha=(1+2*alpha)*b0
  
  K_alpha=exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0))
  K_2alpha=exp(lbeta(a_2alpha,b_2alpha)-(1+2*alpha)*lbeta(a0,b0))
  
  Tb=diag((mu-mu^2))
  Tg=diag(phi)
  
  mu_alpha_star=digamma(a_alpha)-digamma(b_alpha)
  mu_alpha_dagger=digamma(b_alpha)-digamma(a_alpha+b_alpha)
  mu_star=digamma(a0)-digamma(b0)
  mu_dagger=digamma(b0)-digamma(a0+b0)
  nu_alpha=mu^2*trigamma(a_alpha)+(1-mu)^2*trigamma(b_alpha)-trigamma(a_alpha+b_alpha)
  
  mu_2alpha_star=digamma(a_2alpha)-digamma(b_2alpha)
  mu_2alpha_dagger=digamma(b_2alpha)-digamma(a_2alpha+b_2alpha)
  kappa_mu=phi*K_alpha*(mu_alpha_star-mu_star)
  kappa_phi=K_alpha*(mu*(mu_alpha_star-mu_star)+mu_alpha_dagger-mu_dagger)
  #E.u2.phi_2alpha=mu^2*(trigamma(a_2alpha)+trigamma(b_2alpha)+(mu_2alpha_star-mu_star)^2)+2*mu*((mu_2alpha_star-mu_star)*(mu_2alpha_dagger-mu_dagger)-trigamma(b_2alpha))+trigamma(b_2alpha)-trigamma(a_2alpha+b_2alpha)+(mu_2alpha_dagger-mu_dagger)^2
  E.u2.phi_2alpha=(mu*(mu_2alpha_star-mu_star)+(mu_2alpha_dagger-mu_dagger))^2+mu^2*trigamma(a_2alpha)+(1-mu)^2*trigamma(b_2alpha)-trigamma(a_2alpha+b_2alpha)
  
  Lambda_mu_mu=diag(phi^2*K_alpha*((mu_alpha_star-mu_star)^2+trigamma(a_alpha)+trigamma(b_alpha)))
  Lambda_mu_phi=diag(phi*K_alpha*(mu*trigamma(a_alpha)-(1-mu)*trigamma(b_alpha)+mu*(mu_alpha_star-mu_star)^2+(mu_alpha_star-mu_star)*(mu_alpha_dagger-mu_dagger)))
  Lambda_phi_phi=diag(K_alpha*((mu*(mu_alpha_star-mu_star)+mu_alpha_dagger-mu_dagger)^2+nu_alpha))
  
  Sigma_mu_mu=diag(phi^2*K_2alpha*(trigamma(a_2alpha)+trigamma(b_2alpha)+(mu_2alpha_star-mu_star)^2)-kappa_mu^2)
  Sigma_mu_phi=diag(K_2alpha*phi*(mu*trigamma(a_2alpha)-(1-mu)*trigamma(b_2alpha)+mu*(mu_2alpha_star-mu_star)^2+(mu_2alpha_star-mu_star)*(mu_2alpha_dagger-mu_dagger))-kappa_mu*kappa_phi)
  Sigma_phi_phi=diag(K_2alpha*E.u2.phi_2alpha-kappa_phi^2)
  
  Lambda_beta_beta=t(X)%*%Tb%*%Lambda_mu_mu%*%Tb%*%X
  Lambda_beta_gamma=t(X)%*%Tb%*%Lambda_mu_phi%*%Tg%*%Z
  Lambda_gamma_gamma=t(Z)%*%Tg%*%Lambda_phi_phi%*%Tg%*%Z
  
  Lambda=rbind(cbind(Lambda_beta_beta,Lambda_beta_gamma),cbind(t(Lambda_beta_gamma),Lambda_gamma_gamma))
  
  Sigma_beta_beta=t(X)%*%Tb%*%Sigma_mu_mu%*%Tb%*%X
  Sigma_beta_gamma=t(X)%*%Tb%*%Sigma_mu_phi%*%Tg%*%Z
  Sigma_gamma_gamma=t(Z)%*%Tg%*%Sigma_phi_phi%*%Tg%*%Z
  
  Sigma=rbind(cbind(Sigma_beta_beta,Sigma_beta_gamma),cbind(t(Sigma_beta_gamma),Sigma_gamma_gamma))
  
  V=n*solve(Lambda)%*%Sigma%*%t(solve(Lambda))
  
  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma
  
  result$Cov=V/n
  result$Std.Error=c(sqrt(diag(V/n)))
  
  return(result)
}



