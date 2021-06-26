# Score Function LMDPDE - Beta
Psi_Beta=function(Beta,Gamma,y,X,Z,alpha,linkobj)
{
  mu_hat=linkobj$linkfun.mu$inv.link(X%*%Beta)
  phi_hat=linkobj$linkfun.phi$inv.link(Z%*%Gamma)
  
  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)
  b_alpha=b0*(1+alpha)
  
  y_star=log(y)-log(1-y)
  mu_star=suppressWarnings(digamma(a0)-digamma(b0))
  mu_star_alpha=suppressWarnings(digamma(a_alpha)-digamma(b_alpha))
  #Matrixes
  Phi.Tb=phi_hat*inverse(linkobj$linkfun.mu$d.linkfun(mu_hat))
  Phi.Tb.W_alpha=diag(Phi.Tb*degbeta(y_star,mu_hat,phi_hat)^(alpha))
  Phi.Tb.C_alpha=diag(Phi.Tb*exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0)))

  return(t(X)%*%(Phi.Tb.W_alpha%*%(y_star-mu_star)-Phi.Tb.C_alpha%*%(mu_star_alpha-mu_star)))
}

# Score Function LMDPDE - Gamma
Psi_Gamma=function(Beta,Gamma,y,X,Z,alpha,linkobj)
{
  d.link.phi=linkobj$linkfun.phi$d.linkfun
  inv.link.mu=linkobj$linkfun.mu$inv.link
  inv.link.phi=linkobj$linkfun.phi$inv.link
  
  mu_hat=inv.link.mu(X%*%Beta)
  phi_hat=inv.link.phi(Z%*%Gamma)

  a0=mu_hat*phi_hat
  b0=(1-mu_hat)*phi_hat
  a_alpha=a0*(1+alpha)
  b_alpha=b0*(1+alpha)
  
  y_dagger=log(1-y)
  y_star=log(y)-y_dagger
  mu_dagger=suppressWarnings(digamma(b0)-digamma(phi_hat))
  mu_star=suppressWarnings(digamma(a0)-digamma(b0))
  mu_star_alpha=suppressWarnings(digamma(a_alpha)-digamma(b_alpha))
  kappa_alpha=suppressWarnings(mu_hat*(mu_star_alpha-mu_star)+digamma(b_alpha)-digamma(b0)+digamma(phi_hat)-digamma(phi_hat*(1+alpha)))
  eta=mu_hat*(y_star-mu_star)+y_dagger-mu_dagger
  C_alpha=exp(lbeta(a_alpha,b_alpha)-(1+alpha)*lbeta(a0,b0))
  Tg=inverse(d.link.phi(phi_hat))
  W_alpha=degbeta(y_star,mu_hat,phi_hat)^(alpha)
  #Matrixes
  Tg.W_alpha=diag(Tg*W_alpha)
  Tg.C_alpha=diag(Tg*C_alpha)
  return(t(Z)%*%(Tg.W_alpha%*%eta-Tg.C_alpha%*%kappa_alpha))
}


# Modified Score Vector - LMDPDE
Psi_LMDPDE=function(Theta,y,X,Z,alpha,linkobj)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  Beta=Theta[1:ncol(X)]
  Gamma=Theta[1:ncol(Z)+ncol(X)]
  
  psi_beta=Psi_Beta(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha,linkobj=linkobj)
  psi_gamma=Psi_Gamma(Beta,Gamma,y=y,X=X,Z=Z,alpha=alpha,linkobj=linkobj)
  
  return(c(psi_beta,psi_gamma))
}

# Robust Point Estimation - LMDPDE
Robst.LMDPDE.Beta.Reg=function(y,x,z,start_theta,alpha,linkobj,tolerance,maxit)
{
  theta=list()
  link.mu=attributes(linkobj)$name.link.mu
  link.phi=attributes(linkobj)$name.link.phi
  if(missing(tolerance)){tolerance=1e-3}
  if(missing(maxit)){maxit=150}
  if(missing(start_theta))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link.mu,link.phi=link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
  }
  theta$x=rep(0,length(start_theta))
  theta$fvec=10
  theta$msg=theta$error=NULL
  theta=tryCatch(nleqslv(start_theta,Psi_LMDPDE_Cpp,jac=Psi_LMDPDE_Jacobian_C,y=y,X=x,Z=z,alpha=alpha,link_mu=link.mu,link_phi=link.phi,control=list(ftol=tolerance,maxit=maxit),jacobian=TRUE,method="Newton"),error=function(e){
    theta$msg<-e$message
    return(theta)})
  theta$converged=F
  if(all(abs(theta$fvec)<tolerance) & !all(theta$fvec==0)){theta$converged=T}
  return(theta)
}


# Auto Selecting tuning parameter algorithm
Opt.Tuning.LMDPDE=function(y,x,z,link,link.phi,control)
{
  #browser()
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LMDPDE.list=LMDPDE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.5,0.02)
  K=length(alpha_tuning)
  M=control$M
  L=control$L
  n=length(y)
  unstable=F
  sqv.unstable=T
  est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
  Est.param=do.call("c",est.log.lik$coefficients)
  ponto.inicial.robst=ponto.inicial.temp=Initial.points(y,x,z)
  names(ponto.inicial.robst)=names(ponto.inicial.temp)=names(Est.param)=c(colnames(x),colnames(z))
  p=length(Est.param)
  for(k in 1:(M+1))
  {
    control$start=Est.param
    LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(!LMDPDE.par$converged)
    {
      control$start=ponto.inicial.temp
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    }
    if(!LMDPDE.par$converged)
    {
      control$start=ponto.inicial.robst
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
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
    control$start=Est.param
    LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    if(!LMDPDE.par$converged)
    {
      control$start=ponto.inicial.temp
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
    }
    if(!LMDPDE.par$converged)
    {
      control$start=ponto.inicial.robst
      LMDPDE.par=tryCatch(LMDPDE.Beta.Reg(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LMDPDE.par$converged<-FALSE; return(LMDPDE.par)})
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
    LMDPDE.par.star$Optimal.Tuning=TRUE
  }
  return(LMDPDE.par.star)
}


# Sandwich Matrix - LMDPDE
LMDPDE_Cov_Matrix=function(mu,phi,X,Z,alpha,linkobj)
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
  
  Tb=diag(inverse(linkobj$linkfun.mu$d.linkfun(mu)))
  Tg=diag(inverse(linkobj$linkfun.phi$d.linkfun(phi)))
  
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
  
  #V=n*solve(Lambda)%*%Sigma%*%t(solve(Lambda))
  V=n*MASS::ginv(Lambda)%*%Sigma%*%t(MASS::ginv(Lambda))
  
  result=list()
  result$Lambda=Lambda
  result$Sigma=Sigma
  
  result$Cov=V/n
  result$Std.Error=c(sqrt(diag(V/n)))
  
  return(result)
}

#Psi LMDPDE Jacobian
Psi_LMDPDE_Jacobian=function(Theta,y,X,Z,alpha,linkobj)
{
  #X=as.matrix(X)
  #Z=as.matrix(Z)
  Beta=Theta[1:ncol(X)]
  Gamma=Theta[1:ncol(Z)+ncol(X)]
  q=1-alpha
  mu_hat=linkobj$linkfun.mu$inv.link(X%*%Beta)
  phi_hat=linkobj$linkfun.phi$inv.link(Z%*%Gamma)
  a.0=mu_hat*phi_hat
  b.0=(1-mu_hat)*phi_hat
  a.alpha=a.0*(1+alpha)
  b.alpha=b.0*(1+alpha)
  y_dagger=log(1-y)
  y_star=log(y)-y_dagger
  mu_star=suppressWarnings(digamma(a.0)-digamma(b.0))
  mu_star_alpha=suppressWarnings(digamma(a.alpha)-digamma(b.alpha))
  mu_dagger=suppressWarnings(digamma(b.0)-digamma(phi_hat))
  mu_alpha_dagger=digamma(b.alpha)-digamma(phi_hat*(1+alpha))
  
  d.linkmu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  d2.linkmu=linkobj$linkfun.mu$d2.linkfun(mu_hat)
  d.linkphi=linkobj$linkfun.phi$d.linkfun(phi_hat)
  d2.linkphi=linkobj$linkfun.phi$d2.linkfun(phi_hat)
  
  Tb=diag(inverse(d.linkmu))
  Tg=diag(inverse(d.linkphi))
  #browser()
  
  c_star=exp(lbeta(a.alpha,b.alpha)-(1+alpha)*lbeta(a.0,b.0))
  k_mu.mu=(phi_hat^2)*c_star*((1+alpha)*((mu_star_alpha-mu_star)^2+trigamma(a.alpha)+trigamma(b.alpha))-trigamma(a.0)-trigamma(b.0))
  k_mu=phi_hat*c_star*(mu_star_alpha-mu_star)
  u_mu=phi_hat*(y_star-mu_star)
  u_phi=mu_hat*(y_star-mu_star)+y_dagger-mu_dagger
  u_mu.mu=-(phi_hat^2)*(trigamma(a.0)+trigamma(b.0))
  f_q_star=degbeta(y_star,mu_hat,phi_hat)^(alpha)
  
  core1=as.numeric(f_q_star*(u_mu.mu+alpha*u_mu^2)-(k_mu.mu+(d2.linkmu/d.linkmu)*(f_q_star*u_mu-k_mu)))
  
  D1=mu_hat*(mu_star_alpha-mu_star)+(mu_alpha_dagger-mu_dagger)
  dd=(1+alpha)*(D1)*(mu_star_alpha-mu_star)+((1+alpha)*(mu_hat*trigamma(a.alpha)-(1-mu_hat)*trigamma(b.alpha))-(mu_hat*trigamma(a.0)-(1-mu_hat)*trigamma(b.0)))
  #dd=(1+alpha)*(mu_hat*(mu_star_alpha-mu_star)+(mu_alpha_dagger-mu_dagger))*(mu_star_alpha-mu_star)+((1+alpha)*(mu_hat*trigamma(a.alpha)-(1-mu_hat)*trigamma(b.alpha))-(mu_hat*trigamma(a.0)-(1-mu_hat)*trigamma(b.0)))
  k_mu.phi=c_star*(mu_star_alpha-mu_star+phi_hat*dd)
  u_mu.phi=(y_star-mu_star)-phi_hat*(mu_hat*trigamma(a.0)-(1-mu_hat)*trigamma(b.0))
  
  core2=as.numeric(f_q_star*(u_mu.phi+alpha*u_mu*u_phi)-k_mu.phi)
  
  nu.alpha=trigamma(a.alpha)*mu_hat^2+trigamma(b.alpha)*(1-mu_hat)^2-trigamma(phi_hat*(1+alpha))
  nu.0=trigamma(a.0)*mu_hat^2+trigamma(b.0)*(1-mu_hat)^2-trigamma(phi_hat)
  u_phi.phi=trigamma(phi_hat)-(trigamma(a.0)*mu_hat^2+trigamma(b.0)*(1-mu_hat)^2)
  k_phi=c_star*(D1)
  #k_phi=c_star*(mu_hat*(mu_star_alpha-mu_star)+mu_alpha_dagger-mu_dagger)
  
  k_phi.phi=c_star*((1+alpha)*((mu_hat*(mu_star_alpha-mu_star)+mu_alpha_dagger-mu_dagger)^2+nu.alpha)-nu.0)
  
  core3=as.numeric(f_q_star*(u_phi.phi+alpha*u_phi^2)-(d2.linkphi/d.linkphi)*(u_phi*f_q_star-k_phi)-k_phi.phi)
  
  J11=t(X)%*%Tb%*%diag(core1)%*%Tb%*%X
  J12=J21=t(X)%*%Tb%*%diag(core2)%*%Tg%*%Z
  J22=t(Z)%*%Tg%*%diag(core3)%*%Tg%*%Z
  
  J=rbind(cbind(J11,J12),cbind(t(J21),J22))
  return(J)
}


#Hat matrix
hatvalues.robustbetareg.LMDPDE=function(object)
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

