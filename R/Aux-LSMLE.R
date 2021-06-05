#' Score Function LSMLE - Beta
#' 
#' Modified Score Function LSMLE - Beta
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
Psi_Beta_LSMLE=function(Beta,Gama,y,X,Z,alpha)
{
  q=1-alpha
  n=length(y)
  X=as.matrix(X)
  Z=as.matrix(Z)
  mu_hat=h1(X,Beta)
  phi_hat=h2(Z,Gama)
  phi_q=phi_hat/q
  
  a.q=mu_hat*phi_q
  b.q=(1-mu_hat)*phi_q
  
  y_star=log(y/(1-y))
  Phi_q=diag(phi_q)
  mu_star=suppressWarnings(digamma(a.q)-digamma(b.q)) 
  Tb=diag(x=sapply(mu_hat,function(k) (k-k^2)))
  f_q_star=diag((degbeta(y_star=y_star,mu=mu_hat,phi=phi_q))^(alpha))
  
  return(t(X)%*%Phi_q%*%Tb%*%f_q_star%*%(y_star-mu_star))
}

#' Score Function LSMLE - Gamma
#' 
#' Modified Score Function LSMLE - Gamma
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
Psi_Gamma_LSMLE=function(Beta,Gamma,y,X,Z,alpha)
{
  q=1-alpha
  mu_hat=h1(X,Beta)
  phi_hat=h2(Z,Gamma)
  phi_q=phi_hat/q
  
  a.q=mu_hat*phi_q
  b.q=(1-mu_hat)*phi_q
  
  y_star=log(y/(1-y))
  y_dagger=log(1-y)
  mu_star=suppressWarnings(digamma(a.q)-digamma(b.q)) 
  mu_dagger=suppressWarnings(digamma(b.q)-digamma(phi_q))
  
  Tg=diag(phi_q)
  eta=mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger)
  f_q_star=diag((degbeta(y_star,mu_hat,phi_q)^(alpha)))
  
  return(t(Z)%*%Tg%*%f_q_star%*%eta)
}

#' Modified Score Vector - LSMLE
#' 
#' Modified Score Function LSMLE - Gamma
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
Psi_LSMLE=function(Theta,y,X,Z,alpha)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  Beta=Theta[1:ncol(X)]
  Gamma=Theta[1:ncol(Z)+ncol(X)]
  
  psi_beta=Psi_Beta_LSMLE(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha)
  psi_gamma=Psi_Gamma_LSMLE(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha)
  
  return(c(psi_beta,psi_gamma))
}


#' Robust Point Estimation - LSMLE
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
Robst.LSMLE.Beta.Reg=function(y,x,z,start_theta,alpha,tolerance,maxit)
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
  theta$msg=NULL
  #options(warn = 2) #Converte warnings em erros 
  # theta=tryCatch(Newton.Raphson(start_theta,FUN = Psi_LSMLE,alpha=alpha,y=y,X=x,Z=z,details = T,tol=tolerance,M=maxit),error=function(e){
  #   theta$msg<-e$message
  #   return(theta)})
  # theta$x=theta$sol
  # theta=within(theta,rm(sol))
  theta=tryCatch(nleqslv(start_theta,Psi_LSMLE,y=y,X=x,Z=z,alpha=alpha,control=list(ftol=tolerance,maxit=maxit),jacobian=TRUE,method="Newton"),error=function(e){
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
Opt.Tuning.LSMLE=function(y,x,z,control)
{
  if(missing(control)){control=robustbetareg.control()}
  LSMLE.list=list()
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
    LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=Est.param,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.Par)})
    if(!LSMLE.par$converged)
    {
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.temp,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.Par)})
    }
    if(!LSMLE.par$converged)
    {
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.temp,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.Par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=Par.q(LSMLE.par)
    }
    if(is.null(LSMLE.par) || any(is.na(Z.q(LSMLE.par))) || is.null(SE.q(LSMLE.par)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LSMLE.list[[k]]<-LSMLE.par
    zq.t<-unname(rbind(zq.t,Z.q(LSMLE.par)))
  }
  sqv=SQV(zq.t,n,p)
  if(all(sqv<=L))
  {
    LSMLE.par.star<-LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  k=k+1
  while(sqv.unstable)
  {
    LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=Est.param,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.temp,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],alpha.optimal=F,start_theta=ponto.inicial.robst,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=Par.q(LSMLE.par)  
    }
    if(any(is.na(Z.q(LSMLE.par))) || is.null(SE.q(LSMLE.par)))
    {
      unstable=T
      break
    }
    LSMLE.list[[k]]=LSMLE.par
    zq.t=unname(rbind(zq.t,Z.q(LSMLE.par)))
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
    LSMLE.par.star=LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$message="Lack of stability"
  }else{
    LSMLE.par.star=LSMLE.list[[(k-M)]]
    LSMLE.par.star$sqv=sqv
  }
  return(LSMLE.par.star)
}


#' Sandwich Matrix - LSMLE
#' 
#' Covariance matrix of LSMLE.
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
LSMLE_Cov_Matrix=function(mu,phi,X,Z,alpha)
{
  n=length(mu)
  q=1-alpha
  
  a.0=mu*phi
  b.0=(1-mu)*phi
  #
  phi.k0=phi/(1-alpha)
  a.k0=mu*phi.k0
  b.k0=(1-mu)*phi.k0
  #
  phi.k1=phi*(1+alpha)/(1-alpha)
  a.k1=mu*phi.k1
  b.k1=(1-mu)*phi.k1
  
  C.k0=diag(exp(q*lbeta(a.k0,b.k0)-lbeta(a.0,b.0)))
  C.k1=diag(exp(lbeta(a.k1,b.k1)-lbeta(a.0,b.0)-2*alpha*lbeta(a.k0,b.k0)))
  
  mu_star.k0=digamma(a.k0)-digamma(b.k0)
  mu_dagger.k0=digamma(b.k0)-digamma(phi.k0)
  mu_star.k1=digamma(a.k1)-digamma(b.k1)
  mu_dagger.k1=digamma(b.k1)-digamma(phi.k1)
  
  #Logit mean link function 
  g_beta.mu.dev=1/(mu-mu^2)
  #Log precision link funtion
  g_gamma.phi.dev=1/phi
  
  Tb=diag(1/g_beta.mu.dev)
  Tg=diag(1/(g_gamma.phi.dev))
  Phi=diag(phi)
  Q.inv=diag(n)/q
  
  Lambda_mu_mu=diag(trigamma(a.k0)+trigamma(b.k0))
  Lambda_mu_phi=diag(mu*trigamma(a.k0)-(1-mu)*trigamma(b.k0))
  Lambda_phi_phi=diag(mu^2*trigamma(a.k0)+(1-mu)^2*trigamma(b.k0)-trigamma(phi.k0))
  
  Lambda_beta_beta=t(X)%*%Tb%*%Q.inv%*%(Phi^2)%*%C.k0%*%Lambda_mu_mu%*%Tb%*%X
  Lambda_beta_gamma=t(X)%*%Tb%*%Q.inv%*%Phi%*%C.k0%*%Lambda_mu_phi%*%Tg%*%Z
  Lambda_gamma_gamma=t(Z)%*%Tg%*%Q.inv%*%C.k0%*%Lambda_phi_phi%*%Tg%*%Z
  
  Lambda=-1*rbind(cbind(Lambda_beta_beta,Lambda_beta_gamma),cbind(t(Lambda_beta_gamma),Lambda_gamma_gamma))
  
  Sigma_mu_mu=diag(trigamma(a.k1)+trigamma(b.k1)+(mu_star.k1-mu_star.k0)^2)
  Sigma_mu_phi=diag(mu*trigamma(a.k1)-(1-mu)*trigamma(b.k1)+mu*(mu_star.k1-mu_star.k0)^2+(mu_star.k1-mu_star.k0)*(mu_dagger.k1-mu_dagger.k0))
  #Sigma_phi_phi=diag(mu^2*(trigamma(a.k1)+trigamma(b.k1)+(mu_star.k1-mu_star.k0)^2)+2*mu*((mu_star.k1-mu_star.k0)*(mu_dagger.k1-mu_dagger.k0)-trigamma(b.k1))+trigamma(b.k1)-trigamma(a.k1+b.k1)+(mu_dagger.k1-mu_dagger.k0)^2)
  Sigma_phi_phi=diag(((mu*(mu_star.k1-mu_star.k0)+(mu_dagger.k1-mu_dagger.k0))^2)+(mu^2)*trigamma(a.k1)+((1-mu)^2)*trigamma(b.k1)-trigamma(phi.k1))
  
  Sigma_beta_beta=t(X)%*%Tb%*%(Q.inv^2)%*%(Phi^2)%*%C.k1%*%Sigma_mu_mu%*%Tb%*%X
  Sigma_beta_gamma=t(X)%*%Tb%*%(Q.inv^2)%*%Phi%*%C.k1%*%Sigma_mu_phi%*%Tg%*%Z
  Sigma_gamma_gamma=t(Z)%*%Tg%*%(Q.inv^2)%*%C.k1%*%Sigma_phi_phi%*%Tg%*%Z
  
  Sigma=rbind(cbind(Sigma_beta_beta,Sigma_beta_gamma),cbind(t(Sigma_beta_gamma),Sigma_gamma_gamma))
  
  V=n*solve(Lambda)%*%Sigma%*%t(solve(Lambda))
  
  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma
  result$Cov=V/n
  result$Std.Error=c(sqrt(diag(V/n)))
  
  return(result)
}

