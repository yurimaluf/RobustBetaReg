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
  theta=tryCatch(nleqslv(start_theta,Psi_LMDPDE_Cpp,jac=Psi_LMDPDE_Jacobian_C,y=y,X=x,Z=z,alpha=alpha,link_mu=link.mu,link_phi=link.phi,control=list(ftol=tolerance,maxit=maxit),method="Newton",jacobian = T),error=function(e){
    theta$msg<-e$message
    return(theta)})
  theta$converged=F
  if(all(abs(theta$fvec)<tolerance) & !all(theta$fvec==0) & all(diag(theta$jac)<0)){theta$converged=T}
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
    LMDPDE.par.star$Optimal.Tuning=TRUE
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
    LMDPDE.par.star$Optimal.Tuning=TRUE
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


#Hat matrix
hatvalues.LMDPDE=function(object)
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

