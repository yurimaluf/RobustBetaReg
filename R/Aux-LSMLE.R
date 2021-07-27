# Auto Selecting tuning parameter algorithm
Opt.Tuning.LSMLE=function(y,x,z,link,link.phi,control)
{
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LSMLE.list=LSMLE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.5,0.02)
  K=length(alpha_tuning)
  M=control$M
  L=control$L
  n=length(y)
  unstable=F
  sqv.unstable=T
  est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
  #Est.param=c(est.log.lik$coefficients$mean,est.log.lik$coefficients$precision)
  Est.param=do.call("c",est.log.lik$coefficients)
  ponto.inicial.robst=ponto.inicial.temp=Initial.points(y,x,z)
  names(ponto.inicial.robst)=names(ponto.inicial.temp)=names(Est.param)=c(colnames(x),colnames(z))
  p=length(Est.param)
  for(k in 1:(M+1))
  {
    control$start=Est.param
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
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
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  #browser()
  k=k+1
  while(sqv.unstable)
  {
    control$start=Est.param
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
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
    }
    k=k+1
  }
  if(k>=K || unstable)
  {
    LSMLE.par.star=LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
  }else{
    LSMLE.par.star=LSMLE.list[[(k-M)]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
  }
  return(LSMLE.par.star)
}

#Sandwich Matrix - LSMLE
LSMLE_Cov_Matrix=function(mu,phi,X,Z,alpha,linkobj)
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
  Tb=diag(inverse(linkobj$linkfun.mu$d.linkfun(mu)))
  #Log precision link funtion
  Tg=diag(inverse(linkobj$linkfun.phi$d.linkfun(phi)))
  
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
  
  #V=n*solve(Lambda)%*%Sigma%*%t(solve(Lambda))
  V=n*MASS::ginv(Lambda)%*%Sigma%*%t(MASS::ginv(Lambda))
  
  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma
  result$Cov=V/n
  result$Std.Error=c(sqrt(diag(V/n)))
  
  return(result)
}

#Hat matrix
hatvalues.LSMLE=function(object)
{
  mu_hat=object$fitted.values$mu.predict
  phi_hat=object$fitted.values$phi.predict
  y=object$y
  X=object$model$mean
  linkobj=make.link(link.mu=object$link,link.phi=object$link.phi)
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  
  y_star=log(y)-log(1-y)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  W.PHI=diag(x=phi_hat*V_star*(inverse(d.link.mu))^2)
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  return(diag(H))
}

sup_K_psi_C=function(Theta_0,ind_free,ind_fix,eta_0,Beta,Gamma,X,Z,alpha,linkobj,thrd)
{
  m=ncol(X)
  k=ncol(Z)
  #browser()
  #Vector assembling 
  T_0=rep(0,(m+k))
  T_0[ind_free]=Theta_0
  T_0[ind_fix]=eta_0
  B_0=T_0[1:m]
  G_0=T_0[(m+1):(m+k)]
  
  l0=rep(0,(m+k))
  return(suppressWarnings(-(stats::nlminb(l0,K2_psi_C,Beta=Beta,Gamma=Gamma,Beta_0=B_0,Gamma_0=G_0,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd))$objective))
}

#Kumulant of Psi function
K2_psi_C=function(l,Beta,Gamma,Beta_0,Gamma_0,X,Z,alpha,linkobj,thrd)
{
  k=ncol(X)
  m=ncol(Z)
  l1=l[1:k]
  l2=l[(k+1):(k+m)]
  Psi.moment=LaplaceApx2_C(Beta,Gamma,Beta_0,Gamma_0,l1,l2,X,Z,alpha,linkobj,thrd=thrd)
  #if(any(Psi.moment<0)){browser()}
  return(mean(log(Psi.moment)))
}

LaplaceApx2_C=function(Beta,Gamma,Beta_0,Gamma_0,l1,l2,X,Z,alpha,linkobj,thrd)
{
  mu=mu_0=linkobj$linkfun.mu$inv.link(X%*%Beta)
  phi=phi_0=linkobj$linkfun.phi$inv.link(Z%*%Gamma)
  if(!all(Beta==Beta_0))
  {
    mu_0=linkobj$linkfun.mu$inv.link(X%*%Beta_0)  
  }
  if(!all(Gamma==Gamma_0))
  {
    phi_0=linkobj$linkfun.phi$inv.link(Z%*%Gamma_0)  
  }
  Kx=(X%*%l1)/linkobj$linkfun.mu$d.linkfun(mu)
  Kz=(Z%*%l2)/linkobj$linkfun.phi$d.linkfun(phi)
  n=length(mu)
  La=NULL
  #browser()
  La=La_Cpp(mu, phi, mu_0, phi_0, alpha, Kx, Kz,thrd=thrd) 
  return(La)
}

