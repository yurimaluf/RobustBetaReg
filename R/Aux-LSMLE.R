# Score Function LSMLE - Beta
Psi_Beta_LSMLE=function(Beta,Gamma,y,X,Z,alpha,linkobj)
{
  q=1-alpha
  mu_hat=linkobj$linkfun.mu$inv.link(X%*%Beta)
  phi_hat=linkobj$linkfun.phi$inv.link(Z%*%Gamma)
  phi_q=phi_hat/q
  
  a.q=mu_hat*phi_q
  b.q=(1-mu_hat)*phi_q
  
  y_star=log(y)-log(1-y)
  #Phi_q=diag(phi_q)
  mu_star=suppressWarnings(digamma(a.q)-digamma(b.q)) 
  Tb=inverse(linkobj$linkfun.mu$d.linkfun(mu_hat))
  f_q_star=(degbeta(y_star=y_star,mu=mu_hat,phi=phi_q))^(alpha)
  Phi_q.Tb.f_q_star=diag(phi_q*Tb*f_q_star)
  
  return(t(X)%*%Phi_q.Tb.f_q_star%*%(y_star-mu_star))
}

# Score Function LSMLE - Gamma
Psi_Gamma_LSMLE=function(Beta,Gamma,y,X,Z,alpha,linkobj)
{
  q=1-alpha
  mu_hat=linkobj$linkfun.mu$inv.link(X%*%Beta)
  phi_hat=linkobj$linkfun.phi$inv.link(Z%*%Gamma)
  phi_q=phi_hat/q
  
  a.q=mu_hat*phi_q
  b.q=(1-mu_hat)*phi_q
  
  y_dagger=log(1-y)
  y_star=log(y)-y_dagger
  mu_star=suppressWarnings(digamma(a.q)-digamma(b.q)) 
  mu_dagger=suppressWarnings(digamma(b.q)-digamma(phi_q))
  
  Tg=inverse(linkobj$linkfun.phi$d.linkfun(phi_hat))/q
  eta=mu_hat*(y_star-mu_star)+(y_dagger-mu_dagger)
  f_q_star=degbeta(y_star,mu_hat,phi_q)^(alpha)
  Tg.f_q_star=diag(Tg*f_q_star)
  
  return(t(Z)%*%Tg.f_q_star%*%eta)
}

# Modified Score Vector - LSMLE
Psi_LSMLE=function(Theta,y,X,Z,alpha,linkobj)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  Beta=Theta[1:ncol(X)]
  Gamma=Theta[1:ncol(Z)+ncol(X)]
  
  psi_beta=Psi_Beta_LSMLE(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha,linkobj=linkobj)
  psi_gamma=Psi_Gamma_LSMLE(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha,linkobj=linkobj)
  
  return(c(psi_beta,psi_gamma))
}

# Robust Point Estimation - LSMLE
Robst.LSMLE.Beta.Reg=function(y,x,z,start_theta,alpha,linkobj,tolerance,maxit)
{
  theta=list()
  link.mu=attributes(linkobj)$name.link.mu
  link.phi=attributes(linkobj)$name.link.phi
  if(missing(tolerance)){tolerance=1e-3}
  if(missing(maxit)){maxit=150}
  if(missing(start_theta)){
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link.mu,link.phi=link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
  }
  theta$x=rep(0,length(start_theta))
  theta$fvec=10
  theta$msg=NULL
  #browser()
  theta=tryCatch(nleqslv(start_theta,Psi_LSMLE_Cpp,jac=Psi_LSMLE_Jacobian_C,y=y,X=x,Z=z,alpha=alpha,link_mu=link.mu,link_phi=link.phi,control=list(ftol=tolerance,maxit=maxit),method="Newton"),error=function(e){
   theta$msg<-e$message
   return(theta)})
  theta$converged=F
  if(all(abs(theta$fvec)<tolerance) & !all(theta$fvec==0)){theta$converged=T}
  
  return(theta)
}

# Auto Selecting tuning parameter algorithm
Opt.Tuning.LSMLE=function(y,x,z,link,link.phi,control)
{
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LSMLE.list=list()
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
    LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
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
  #browser()
  k=k+1
  while(sqv.unstable)
  {
    control$start=Est.param
    LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
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

#Psi_LSMLE Jacobian
Psi_LSMLE_Jacobian=function(Theta,y,X,Z,alpha,linkobj)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  Beta=Theta[1:ncol(X)]
  Gamma=Theta[1:ncol(Z)+ncol(X)]
  q=1-alpha
  mu_hat=linkobj$linkfun.mu$inv.link(X%*%Beta)
  phi_hat=linkobj$linkfun.phi$inv.link(Z%*%Gamma)
  phi_q=phi_hat/q
  a.q=mu_hat*phi_q
  b.q=(1-mu_hat)*phi_q
  y_star=log(y)-log(1-y)
  y_dagger=log(1-y)
  mu_star=suppressWarnings(digamma(a.q)-digamma(b.q))
  mu_dagger=suppressWarnings(digamma(b.q)-digamma(phi_q))
  
  d.linkmu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  d2.linkmu=linkobj$linkfun.mu$d2.linkfun(mu_hat)
  d.linkphi=linkobj$linkfun.phi$d.linkfun(phi_hat)
  d2.linkphi=linkobj$linkfun.phi$d2.linkfun(phi_hat)
  
  Tb=diag(inverse(d.linkmu))
  Tg=diag(inverse(d.linkphi))
  #browser()
  
  u_mu.mu=-(phi_q)^2*(trigamma(a.q)+trigamma(b.q))
  u_mu.phi=((y_star-mu_star)-phi_q*(mu_hat*trigamma(a.q)-(1-mu_hat)*trigamma(b.q)))/q
  u_mu=phi_q*(y_star-mu_star)
  u_phi.phi=(trigamma(phi_q)-trigamma(a.q)*mu_hat^2-trigamma(b.q)*(1-mu_hat)^2)/(q^2)
  u_phi=(mu_hat*(y_star-mu_star)+y_dagger-mu_dagger)/q
  f_q_star=degbeta(y_star,mu_hat,phi_q)^(alpha)
  
  #core1=as.numeric(f_q_star*(((1-2*mu_hat)/(mu_hat-mu_hat^2))*u_mu+u_mu.mu+alpha*u_mu^2))
  core1=as.numeric(f_q_star*(-(d2.linkmu/d.linkmu)*u_mu+u_mu.mu+alpha*u_mu^2))
  core2=as.numeric(f_q_star*(u_mu.phi+alpha*u_mu*u_phi))
  core3=as.numeric(f_q_star*(-(d2.linkphi/d.linkphi)*u_phi+u_phi.phi+alpha*u_phi^2))
  
  J11=t(X)%*%Tb%*%diag(core1)%*%Tb%*%X
  J12=J21=t(X)%*%Tb%*%diag(core2)%*%Tg%*%Z
  J22=t(Z)%*%Tg%*%diag(core3)%*%Tg%*%Z
  
  J=rbind(cbind(J11,J12),cbind(t(J21),J22))
  return(J)
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
