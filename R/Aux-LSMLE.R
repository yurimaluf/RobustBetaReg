# Auto Selecting tuning parameter algorithm
Opt.Tuning.LSMLE=function(y,x,z,link,link.phi,control)
{
  #browser()
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
  if(is.null(est.log.lik))
  {
    LSMLE.par$coefficients=Initial.points(y,x,z)
    LSMLE.par$Tuning=0
    LSMLE.par$converged=FALSE
    LSMLE.par$sqv=0
    LSMLE.par$Optimal.Tuning=TRUE
    LSMLE.par$message="Error in chol.default(K)"
    return(LSMLE.par)
  }
  Est.param=do.call("c",est.log.lik$coefficients)
  ponto.inicial.robst=ponto.inicial.temp=Initial.points(y,x,z)
  names(ponto.inicial.robst)=names(ponto.inicial.temp)=names(Est.param)=c(colnames(x),colnames(z))
  p=length(Est.param)
  for(k in 1:(M+1))
  {
    #control$start=Est.param
    control$start=ponto.inicial.robst
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      #control$start=ponto.inicial.robst
      control$start=Est.param
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)
    }
    if(is.null(LSMLE.par) || any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LSMLE.list[[k]]<-LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
  }
  
  #browser()
  sqv=as.numeric(SQV_Cpp(zq.t,n,p))
  if(all(sqv<=L))
  {
    LSMLE.par.star<-LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  
  k=k+1
  while(sqv.unstable)
  {
    #control$start=Est.param
    control$start=ponto.inicial.temp
    
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      #control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=Est.param
      #control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)  
    }
    if(any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      unstable=T
      break
    }
    LSMLE.list[[k]]=LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
    sqv=as.numeric(SQV_Cpp(zq.t,n,p))
    sqv.test=sqv[(k-M):(k-1)]
    if(all(sqv.test<=L) || k==K )
    {
      sqv.unstable=F
    }
    k=k+1
  }
  if(k>=K || unstable)
  {
    #LSMLE.par.star=LSMLE.list[[1]]
    LSMLE.par.star=LSMLE.fit(y,x,z,alpha=0,link=link,link.phi=link.phi)
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
  }else{
    LSMLE.par.star=LSMLE.list[[(k-1-M)]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
  }
  return(LSMLE.par.star)
}

# Auto Selecting tuning parameter algorithm 2
Opt.Tuning.LSMLE.2=function(y,x,z,link,link.phi,control)
{
  #browser()
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LSMLE.list=LSMLE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.4,0.02)
  K=length(alpha_tuning)
  M1=10
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
  for(k in 1:M1)
  {
    #browser()
    #control$start=Est.param
    control$start=ponto.inicial.robst
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      #control$start=ponto.inicial.robst
      control$start=Est.param
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)
    }
    if(is.null(LSMLE.par) || any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LSMLE.list[[k]]<-LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
  }
  
  sqv=as.numeric(SQV_Cpp(zq.t,n,p))
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: Todos a baixo de L
  {
    LSMLE.par.star<-LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  if(alpha.ind<8){#Step-2: 
    LSMLE.par.star<-LSMLE.list[[alpha.ind+1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  
  k=11
  while(sqv.unstable)
  {
    #control$start=Est.param
    control$start=ponto.inicial.temp
    
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      #control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=Est.param
      #control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)  
    }
    if(any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      unstable=T
      break
    }
    LSMLE.list[[k]]=LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
    sqv=as.numeric(SQV_Cpp(zq.t,n,p))
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
    LSMLE.par.star=LSMLE.fit(y,x,z,alpha=0,link=link,link.phi=link.phi)
    #LSMLE.par.star=LSMLE.fit(y,x,z,alpha=alpha_tuning[k-1],link=link,link.phi=link.phi)
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
  }else{
    LSMLE.par.star=LSMLE.list[[(k-1-M)]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
  }
  return(LSMLE.par.star)
}

# Auto Selecting tuning parameter algorithm 3
Opt.Tuning.LSMLE.3=function(y,x,z,link,link.phi,control)
{
  #browser()
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LSMLE.list=LSMLE.par=list()
  zq.t=NULL
  alpha_tuning=seq(0,0.4,0.02)
  K=length(alpha_tuning)
  M1=11
  M=control$M
  L=control$L
  n=length(y)
  unstable=F
  sqv.unstable=T
  est.log.lik=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
  # if(is.null(est.log.lik))
  # {
  #   LSMLE.par$coefficients=Initial.points(y,x,z)
  #   LSMLE.par$Tuning=0
  #   LSMLE.par$converged=FALSE
  #   LSMLE.par$sqv=0
  #   LSMLE.par$Optimal.Tuning=TRUE
  #   LSMLE.par$message="Error in chol.default(K)"
  #   return(LSMLE.par)
  # }
  Est.param=do.call("c",est.log.lik$coefficients)
  ponto.inicial.robst=ponto.inicial.temp=Initial.points(y,x,z)
  names(ponto.inicial.robst)=names(ponto.inicial.temp)=names(Est.param)=c(colnames(x),colnames(z))
  p=length(Est.param)
  for(k in 1:M1)
  {
    #browser()
    #control$start=Est.param
    control$start=ponto.inicial.robst
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      #control$start=ponto.inicial.robst
      control$start=Est.param
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)
    }
    if(is.null(LSMLE.par) || any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LSMLE.list[[k]]<-LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
  }
  #browser()
  sqv=as.numeric(SQV_Cpp(zq.t,n,p))
  alpha.ind=max(0,which(sqv>L))
  if(alpha.ind==0)#Step-2: Todos a baixo de L
  {
    LSMLE.par.star<-LSMLE.list[[1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  if(alpha.ind<8){ 
    LSMLE.par.star<-LSMLE.list[[alpha.ind+1]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    rm(LSMLE.list)
    return(LSMLE.par.star)
  }
  
  reached=FALSE
  k=M1+1
  while(sqv.unstable & !reached)
  {
    #browser()
    #control$start=Est.param
    control$start=ponto.inicial.temp
    
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.robst
      #control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      control$start=Est.param
      #control$start=ponto.inicial.robst
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)  
    }
    if(any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      unstable=T
      break
    }
    LSMLE.list[[k]]=LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
    sqv=as.numeric(SQV_Cpp(zq.t,n,p))
    sqv.test=sqv[(k-M):(k-1)]
    if(all(sqv.test<=L) )
    {
      sqv.unstable=F
    }
    k=k+1
    if(k>=K)#Condition of Step 6
    {
      reached=TRUE
    }
  }
  if(reached)
  {
    #browser()
    k=suppressWarnings(max(1,min(which(rollapply(sqv<L,M,sum)==M)))+M+1)
  }
  if(k>=K || unstable)
  {
    LSMLE.par.star=LSMLE.list[[1]]
    #LSMLE.par.star=LSMLE.fit(y,x,z,alpha=0,link=link,link.phi=link.phi)
    #LSMLE.par.star=LSMLE.fit(y,x,z,alpha=alpha_tuning[k-1],link=link,link.phi=link.phi)
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
    LSMLE.par.star$message="Lack of stability"
  }else{
    LSMLE.par.star=LSMLE.list[[(k-1-M)]]
    LSMLE.par.star$sqv=sqv
    LSMLE.par.star$Optimal.Tuning=TRUE
  }
  return(LSMLE.par.star)
}



# Auto Selecting tuning parameter algorithm List
Opt.Tuning.LSMLE.List=function(y,x,z,link,link.phi,B,control)
{
  #browser()
  if(missing(control)){control=robustbetareg.control()}
  control$alpha.optimal=FALSE
  LSMLE.list=LSMLE.par=list()
  zq.t=parametro=stdError=NULL
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
  for(k in 1:B)
  {
    #browser()
    #control$start=Est.param
    control$start=ponto.inicial.robst
    LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e) {LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    if(!LSMLE.par$converged)
    {
      control$start=ponto.inicial.temp
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(!LSMLE.par$converged)
    {
      #control$start=ponto.inicial.robst
      control$start=Est.param
      LSMLE.par=tryCatch(LSMLE.fit(y,x,z,alpha=alpha_tuning[k],link=link,link.phi=link.phi,control = control),error=function(e){LSMLE.par$converged<-FALSE; return(LSMLE.par)})
    }
    if(LSMLE.par$converged)
    {
      ponto.inicial.temp=do.call("c",LSMLE.par$coefficients)
    }
    if(is.null(LSMLE.par) || any(is.na(do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error))) || is.null(do.call("c",LSMLE.par$std.error)))
    {
      sqv.unstable=F
      unstable=T
      break
    }
    LSMLE.list[[k]]<-LSMLE.par
    zq.t=unname(rbind(zq.t,do.call("c",LSMLE.par$coefficients)/do.call("c",LSMLE.par$std.error)))
    
    parametro=unname(rbind(parametro,do.call("c",LSMLE.par$coefficients)))
    stdError=unname(rbind(stdError,do.call("c",LSMLE.par$std.error)))
  }
  
  sqv=as.numeric(SQV_Cpp(zq.t,n,p))
  lista=list(sqv=sqv,par=parametro,std=stdError)
  return(lista)
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
  
  #mean link function 
  Tb=diag((linkobj$linkfun.mu$d.linkfun(mu))^(-1))
  #precision link funtion
  Tg=diag((linkobj$linkfun.phi$d.linkfun(phi))^(-1))
  
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
  #temp
  result$K=t(X)%*%Tb%*%Q.inv%*%(Phi^2)%*%C.k0%*%Lambda_mu_mu%*%Tb%*%X
  result$Std.Error=suppressWarnings(c(sqrt(diag(V/n))))
  
  return(result)
}

#Hat matrix
hatvalues.LSMLE=function(object)
{
  mu_hat=object$fitted.values$mu.predict
  phi_hat=object$fitted.values$phi.predict
  y=object$y
  X=object$model$mean
  linkobj=set.link(link.mu=object$link,link.phi=object$link.phi)
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  
  #y_star=log(y)-log(1-y)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  W.PHI=diag(x=phi_hat*V_star*((d.link.mu)^(-2)))
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  return(diag(H))
}

#Supremum
sup_K_psi_C=function(Theta_0,ind_free,ind_fix,eta_0,Beta,Gamma,X,Z,alpha,linkobj,thrd)
{
  m=ncol(X)
  k=ncol(Z)
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
  return(mean(log(abs(Psi.moment))))
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
  #browser()
  La=La_Cpp(mu, phi, mu_0, phi_0, alpha, Kx, Kz,thrd=thrd)
  return(La)
}


###################################################################################################

sup_K_psi=function(Theta_0,ind_free,ind_fix,eta_0,Beta,Gamma,X,Z,alpha,linkobj,thrd)
{
  m=ncol(X)
  k=ncol(Z)
  #Vector assembling 
  T_0=rep(0,(m+k))
  T_0[ind_free]=Theta_0
  T_0[ind_fix]=eta_0
  B_0=T_0[1:m]
  G_0=T_0[(m+1):(m+k)]
  
  l0=rep(0,(m+k))
  return(suppressWarnings(-(stats::nlminb(l0,K2_psi,Beta=Beta,Gamma=Gamma,Beta_0=B_0,Gamma_0=G_0,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd))$objective))
}


#Kumulant of Psi function
K2_psi=function(l,Beta,Gamma,Beta_0,Gamma_0,X,Z,alpha,linkobj,thrd)
{
  k=ncol(X)
  m=ncol(Z)
  l1=l[1:k]
  l2=l[(k+1):(k+m)]
  Psi.moment=LaplaceApx2(Beta,Gamma,Beta_0,Gamma_0,l1,l2,X,Z,alpha,linkobj)
  return(mean(log(abs(Psi.moment))))
}

LaplaceApx2=function(Beta,Gamma,Beta_0,Gamma_0,l1,l2,X,Z,alpha,linkobj)
{
  #t.max=25
  #t.min=-25
  #browser()
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
  for(i in 1:n)
  {
    opt=root.p(mu=mu[i],phi=phi[i],mu_0=mu_0[i],phi_0=phi_0[i],alpha=alpha,Kx=Kx[i],Kz=Kz[i])
    y.0=opt$root
    p.0=opt$value
    
    D.S=Derivatives.set(y.0,mu=mu[i],phi=phi[i],mu_0=mu_0[i],phi_0=phi_0[i],alpha=alpha,Kx=Kx[i],Kz=Kz[i])
    d2p=D.S$D2
    d3p=D.S$D3
    d4p=D.S$D4
    
    q=1-alpha
    phi_q=phi[i]/q
    D.f=DX.f(y.0,mu=mu[i],phi=phi_q,alpha=alpha)
    df=D.f$df
    d2f=D.f$d2f
    d3f=D.f$d3f
    d4f=D.f$d4f
    
    a1=d4p/(8*d2p^2)-5*(d3p^2)/(24*d2p^3)
    La=c(La,sqrt(2*pi/abs(d2p))*exp(p.0)*(1+a1))
    #La=c(La,a1)
  }
  return(La)
}

# Find Saddlepoint
root.p=function(mu,phi,mu_0,phi_0,alpha,Kx,Kz)
{
  #browser()
  result=list()
  roots=rootSolve::uniroot.all(d.p,lower=-25,upper=25,mu=mu,phi=phi,mu_0=mu_0,phi_0=phi_0,alpha=alpha,Kx=Kx,Kz=Kz)
  if(length(roots)>1)
  {
    p.vector=p(roots,mu,phi,mu_0,phi_0,alpha,Kx,Kz)
    index=which.max(p.vector)
  }else{
    p.vector=p(roots,mu,phi,mu_0,phi_0,alpha,Kx,Kz)
    index=1
  }
  result$root=roots[index]
  result$value=p.vector[index]
  return(result)
}

#First Derivative of expoent of integrand
d.p=function(y_star,mu,phi,mu_0,phi_0,alpha,Kx,Kz)
{
  q=1-alpha
  phi_q=phi/q
  a=mu*phi_q
  b=(1-mu)*phi_q
  muA=digamma(a)-digamma(b)
  muB=digamma(b)-digamma(phi_q)
  EXP.Y=exp(y_star)
  
  u1=phi_q*(y_star-muA)
  u2=(mu*(y_star-muA)-Rmpfr::log1pexp(y_star)-muB)/q
  d.u1=phi_q
  d.u2=(mu-EXP.Y/(1+EXP.Y))/q
  f=degbeta(y_star,mu = mu, phi=phi_q)^alpha
  df=(alpha*phi_q*(mu+(mu-1)*EXP.Y)*f)/(1+EXP.Y)
  dlnf=phi_0*(mu_0+(mu_0-1)*EXP.Y)/(1+EXP.Y)
  
  return((Kx*d.u1+Kz*d.u2)*f+(Kx*u1+Kz*u2)*df+dlnf)
}

#Expoent of Integrand Function
p=function(y_star,mu,phi,mu_0,phi_0,alpha,Kx,Kz)
{
  q=1-alpha
  phi_q=phi/q
  a=mu*phi_q
  b=(1-mu)*phi_q
  muA=digamma(a)-digamma(b)
  muB=digamma(b)-digamma(phi_q)
  u1=Kx*phi_q*(y_star-muA)
  u2=Kz*(mu*(y_star-muA)-Rmpfr::log1pexp(y_star)-muB)/q
  
  return((u2+u1)*degbeta(y_star,mu = mu, phi=phi_q)^alpha+log(degbeta(y_star,mu = mu_0, phi=phi_0)))
  #return(exp(-(lbeta(a,b)+b*y_star+phi_q*log(1+exp(-y_star)))))
}


#Derivatives of expoent integrand
Derivatives.set=function(y_star,mu,phi,mu_0,phi_0,alpha,Kx,Kz)
{
  result=list()
  q=1-alpha
  phi_q=phi/q
  a=mu*phi_q
  b=phi_q-a
  muA=digamma(a)-digamma(b)
  muB=digamma(b)-digamma(phi_q)
  EXP.Y=exp(y_star)
  
  u1=phi_q*(y_star-muA)
  u2=(mu*(y_star-muA)-Rmpfr::log1pexp(y_star)-muB)/q
  d.u1=phi_q
  d.u2=(mu-EXP.Y/(1+EXP.Y))/q
  d2.u2=-(1/q)*EXP.Y/(1+EXP.Y)^2
  d3.u2=(1/q)*EXP.Y*(EXP.Y-1)/(1+EXP.Y)^3
  d4.u2=(-1/q)*EXP.Y*(-4*EXP.Y+exp(2*y_star)+1)/(1+EXP.Y)^4
  
  f=degbeta(y_star,mu = mu, phi=phi_q)^alpha
  D.f=DX.f(y_star,mu,phi_q,alpha)
  df=D.f$df
  d2f=D.f$d2f
  d3f=D.f$d3f
  d4f=D.f$d4f
  
  DX.ln.f=dX.ln.f(y_star,mu_0,phi_0)
  dlnf=DX.ln.f$d.ln.f
  d2lnf=DX.ln.f$d2.ln.f
  d3lnf=DX.ln.f$d3.ln.f
  d4lnf=DX.ln.f$d4.ln.f
  
  D2=(Kz*d2.u2)*f+2*(Kx*d.u1+Kz*d.u2)*df+(Kx*u1+Kz*u2)*d2f+d2lnf
  D3=(Kz*d3.u2)*f+3*(Kz*d2.u2)*df+3*(Kx*d.u1+Kz*d.u2)*d2f+(Kx*u1+Kz*u2)*d3f+d3lnf
  D4=(Kz*d4.u2)*f+4*(Kz*d3.u2)*df+6*(Kz*d2.u2)*d2f+4*(Kx*d.u1+Kz*d.u2)*d3f+(Kx*u1+Kz*u2)*d4f+d4lnf
  
  result$D2=D2
  result$D3=D3
  result$D4=D4
  
  return(result)
}

#Derivative of EGB
DX.f=function(y_star,mu,phi,alpha)
{
  result=list()
  EGB=degbeta(y_star,mu,phi)^alpha
  AP=alpha*phi
  AMP=AP*mu
  AMP.m1=AP*(mu-1)
  EXP.Y=exp(y_star)
  
  df=(AMP+AMP.m1*EXP.Y)*EGB/(1+EXP.Y)
  df2=(EGB/(1+EXP.Y)^2)*(AMP^(2)+EXP.Y^(2)*AMP.m1^2+AP*EXP.Y*(2*AMP.m1*mu-1))
  df3=(EGB/(1+EXP.Y)^3)*(AMP^2*(AMP+EXP.Y*(AMP-AP-2))+(AMP.m1^2)*EXP.Y^(2)*(AMP+AMP.m1*EXP.Y+2)+AP*(2*AMP.m1*mu-1)*EXP.Y*(AMP+EXP.Y*(AMP.m1-1)+1))
  
  fa=EGB/(1+EXP.Y)^4
  a1=AMP^(2)*(AMP-AP-2)*EXP.Y*(AMP+EXP.Y*(AMP-AP-2)+1)
  a12=AMP^(3)*(AMP+EXP.Y*(AMP.m1-3))
  a21=AMP.m1^(2)*(AMP+2)*exp(2*y_star)*(AMP+EXP.Y*(AMP.m1-1)+2)
  a22=(AMP.m1*EXP.Y)^(3)*(AMP+AMP.m1*EXP.Y+3)
  a31=AP*(2*AMP.m1*mu-1)*(AMP+1)*EXP.Y*(AMP+EXP.Y*(AMP.m1-2)+1)
  a32=AP*(2*AMP.m1*mu-1)*(AMP.m1-1)*EXP.Y^(2)*(AMP+EXP.Y*(AMP.m1-1)+2)
  
  df4=fa*(a1+a12+a21+a22+a31+a32)
  result$df=df
  result$d2f=df2
  result$d3f=df3
  result$d4f=df4
  return(result)
}


#Derivative of log-EGB
dX.ln.f=function(y_star,mu,phi,alpha)
{
  if(missing(alpha)){alpha=1}
  EXP.Y=exp(y_star)
  AP=alpha*phi
  result=list()
  
  result$d.ln.f=AP*(mu+(mu-1)*EXP.Y)/(1+EXP.Y)
  result$d2.ln.f=-AP*EXP.Y/(1+EXP.Y)^2
  result$d3.ln.f=AP*EXP.Y*(EXP.Y-1)/(1+EXP.Y)^3
  result$d4.ln.f=-AP*EXP.Y*(-4*EXP.Y+(EXP.Y^2)+1)/(1+EXP.Y)^4
  return(result)
}

#Laplace approximation second order
LaplaceApx2_C_B=function(Beta,Gamma,Beta_0,Gamma_0,l1,l2,X,Z,alpha,linkobj,thrd)
{
  #browser()
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
  La=y0=p0=NULL
  for(i in 1:n)
  {
    opt=root.p_C(mu=mu[i],phi=phi[i],mu_0=mu_0[i],phi_0=phi_0[i],alpha=alpha,Kx=Kx[i],Kz=Kz[i])
    #opt=root.p(mu=mu[i],phi=phi[i],mu_0=mu_0[i],phi_0=phi_0[i],alpha=alpha,Kx=Kx[i],Kz=Kz[i])
    y0=c(y0,opt$root)
    p0=c(p0,opt$value)
  } 
  #browser()
  La=La_CppB(p0,y0, mu, phi, mu_0, phi_0, alpha, Kx, Kz,thrd=8) 
  
  return(La)
}

root.p_C=function(mu,phi,mu_0,phi_0,alpha,Kx,Kz)
{
  result=list()
  #browser()
  #if(is.na(Kx) || is.na(Kz)){browser()}
  roots=rootSolve::uniroot.all(dp_LSMLE_Cpp,lower=-25,upper=25,mu=mu,phi=phi,mu_0=mu_0,phi_0=phi_0,alpha=alpha,Kx=Kx,Kz=Kz)
  #root.teste=uniroot_R(lower=-25, upper=25,mu=mu,phi=phi,mu_0=mu_0,phi_0=phi_0,alpha=alpha,Kx=Kx,Kz=Kz)
  #root.teste2=uniroot_C(lower=-25, upper=25,mu=mu,phi=phi,mu_0=mu_0,phi_0=phi_0,alpha=alpha,Kx=Kx,Kz=Kz)
  if(length(roots)>1)
  {
    p.vector=sapply(1:length(roots),function(k) p_LSMLE_Cpp(roots[k],mu,phi,mu_0,phi_0,alpha,Kx,Kz))
    index=which.max(p.vector)
  }else{
    p.vector=tryCatch(p_LSMLE_Cpp(roots,mu,phi,mu_0,phi_0,alpha,Kx,Kz),error=function(e){problema<-TRUE; return(NULL)})
    index=1
  }
  result$root=roots[index]
  result$value=p.vector[index]
  return(result)
}



########################################################################################################
#' @export
coef.LSMLE=function(object,model=c("full","mean","precision"))
{
  cf <- object$coefficients
  model=match.arg(model)
  switch(model, mean = {
    cf$mean
  }, precision = {
    cf$precision
  }, full = {
    nam1 <- names(cf$mean)
    nam2 <- names(cf$precision)
    cf <- c(cf$mean, cf$precision)
    names(cf) <- c(nam1, if (identical(nam2, "(Phi)")) "(phi)" else paste("(phi)",nam2, sep = "_"))
    cf
  })
}

#' @export
predict.LSMLE = function(object, newdata = NULL, type = c("response", "link", "precision", "variance", "quantile"), at = 0.5) 
{
  #browser()
  type <- match.arg(type)
  if (type == "quantile") {
    qfun <- function(at, mu, phi) {
      rval <- sapply(at, function(p) qbeta(p, mu * phi, (1 - mu) * phi))
      if (length(at) > 1L) {
        if (NCOL(rval) == 1L) 
          rval <- matrix(rval, ncol = length(at), dimnames = list(unique(names(rval)),NULL))
        colnames(rval) <- paste("q_", at, sep = "")
      }
      else {
        rval <- drop(rval)
      }
      rval
    }
  }
  if (missing(newdata)) {
    rval <- switch(type, response = {
      object$fitted.values$mu.predict
    }, link = {
      set.link(object$link,object$link.phi)$linkfun.mu$linkfun(object$fitted.values$mu.predict)
    }, precision = {
      object$fitted.values$phi.predict
    }, variance = {
      object$fitted.values$mu.predict*(1-object$fitted.values$mu.predict)/(1+object$fitted.values$phi.predict)
    }, quantile = {
      qfun(at, object$fitted.values$mu.predict, object$fitted.values$phi.predict)
    })
    return(rval)
  }else{
    #browser()
    newdata=tryCatch(as.data.frame(newdata),error=function(e){newdata})
    x=model.matrix(object$formula,data=newdata,rhs = 1L)
    z=model.matrix(object$formula,data=newdata,rhs = 2L)
    
    rval <- switch (type, response = {
      set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean) 
    }, link = {
      mu_predict=set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
      set.link(object$link,object$link.phi)$linkfun.mu$linkfun(mu_predict)
    }, precision = {
      set.link(object$link,object$link.phi)$linkfun.phi$inv.link(z%*%object$coefficients$precision)
    }, variance = {
      mu_predict=set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
      phi_predict=set.link(object$link,object$link.phi)$linkfun.phi$inv.link(z%*%object$coefficients$precision)
      mu_predict*(1-mu_predict)/phi_predict
    }, quantile={
      mu_predict=set.link(object$link,object$link.phi)$linkfun.mu$inv.link(x%*%object$coefficients$mean)
      phi_predict=set.link(object$link,object$link.phi)$linkfun.phi$inv.link(z%*%object$coefficients$precision)
      qfun(at,mu_predict,phi_predict)
    })
    return(rval)
  }
}



