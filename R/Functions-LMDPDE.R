#' @rdname robustbetareg
#'
# Robust Estimation - LMDPDE
LMDPDE.fit=function(y,x,z,alpha=NULL,link="logit",link.phi="log",control=robustbetareg.control(...), ...)
{
  result=theta=list()
  #Arguments Checking
  alpha.optimal=control$alpha.optimal
  start_theta=control$start
  if(!is.null(alpha)){alpha.optimal=FALSE}
  if(is.null(alpha) & alpha.optimal==FALSE){alpha=0}
  linkobj=make.link(link.mu = link,link.phi = link.phi)
  k=ncol(x)
  m=ncol(z)
  if(alpha.optimal)
  {
    return(Opt.Tuning.LMDPDE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
  }
  #Point Estimation
  theta=tryCatch(nleqslv(start_theta,Psi_LMDPDE_Cpp,jac=Psi_LMDPDE_Jacobian_C,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control=list(ftol=control$tolerance,maxit=control$maxit),method="Newton",jacobian = T),error=function(e){
  theta$msg<-e$message;return(theta)})
  theta$converged=F
  if(all(abs(theta$fvec)<control$tolerance) & !all(theta$fvec==0) & all(diag(theta$jac)<0)){theta$converged=T}
  
  #Predict values
  beta=theta$x[1:k]
  gamma=theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  eta=x%*%beta
  mu_hat=linkobj$linkfun.mu$inv.link(eta)
  phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)
  
  #Expected Standard Error
  MM=tryCatch(LMDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj),error=function(e) NULL)
  if(is.null(MM))
  {
    vcov= std.error.LMDPDE=NULL
  }else{
    vcov=MM$Cov
    std.error.LMDPDE=MM$Std.Error
  }

  #Register of output values 
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients
  result$vcov=vcov
  pseudor2 <- if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LMDPDE)))
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start=start_theta #Alterar caso optar por duas tentativas de ponto inicial
  result$weights=degbeta(y_star,mu_hat,phi_hat)^(alpha)#Weights
  result$Tuning=alpha
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)
  result$link=link
  result$link.phi=link.phi
  result$Optimal.Tuning=alpha.optimal
  result$pseudo.r.squared=pseudor2
  result$control=control
  if(any(is.na(std.error.LMDPDE)))
  {
    str1="Standard-Error is unvailable"
  }
  if(!is.null(theta$msg))
  {
    str2=theta$msg
  }
  if(!is.null(str1)||!is.null(str2))
  {
    result$Message=c(str1,str2)
  }
  if(!any(is.na(std.error.LMDPDE)))
  {
    names(std.error.LMDPDE)<-c(colnames(x),colnames(z))
    se.beta=std.error.LMDPDE[1:k]
    se.gamma=std.error.LMDPDE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
  }
  class(result)="LMDPDE"
  return(result)
}

#' @export
robustbetareg.control.LMDPDE=function(object)
{
  result=object$control
  return(result)
}


#' @export
WaldTypeTest.LMDPDE=function(object,FUN,...)
{
  general=FALSE
  #browser()
  if(missing(FUN)){general=T}
  if(!object$converged){stop(paste("There is no convergence in the model",deparse(substitute(object))))}
  cl = match.call()
  n=length(object$fitted.values$mu.predict)
  V=object$vcov
  b=object$coefficient$mean
  g=object$coefficient$precision
  p=c(b,g)
  k1=length(b)
  k2=length(g)
  k3=length(p)
  if(general)
  {
    result.beta=result.gamma=NULL
    if(ncol(object$model$mean)>1)
    {
      #Beta
      result.beta=list()
      f.b=function(x){x[2:k1]}
      M=numDeriv::jacobian(f.b,p,...)
      m=f.b(p,...)
      r=length(m)
      if(Matrix::rankMatrix(M)[1]!=r)
      {
        stop("The Rank Matrix is not supported")
      }
      W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
      #Beta Register
      result.beta$W.alpha=as.numeric(W_alpha)
      result.beta$df=r
      result.beta$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    }
    if(ncol(object$model$precision)>1)
    {
      #Gamma
      result.gamma=list()
      f.g=function(x){x[(k1+2):k3]}
      M=numDeriv::jacobian(f.g,p,...)
      m=f.g(p,...)
      r=length(m)
      if(Matrix::rankMatrix(M)[1]!=r){stop("The Rank Matrix is not supported")}
      W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
      #Gamma Register
      result.gamma$W.alpha=as.numeric(W_alpha)
      result.gamma$df=r
      result.gamma$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    }
    result=list(beta.wald=result.beta,gamma.wald=result.gamma,general=general,msg="Results based on LMDPDE.")
    class(result)="WaldTest_LMDPDE"
    return(result)
  }else{
    result=list()
    #Hipothesis
    f=FUN
    M=numDeriv::jacobian(f,p,...)
    m=f(p,...)
    r=length(m)
    n=length(object$fitted.values$mu.predict)
    if(Matrix::rankMatrix(M)[1]!=r){stop("The Rank Matrix is not supported")}
    W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
    #Register
    result$W.alpha=as.numeric(W_alpha)
    result$df=r
    result$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    result$general=general
    result$msg="Results based on LMDPDE."
    class(result)="WaldTest_LMDPDE"
    return(result)
  }
}

#' @export
plotenvelope.LMDPDE=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection"),conf=0.95,n.sim=50,PrgBar=T,control=robustbetareg.control(...), ...)
{
  if(missing(control)){control=robustbetareg.control(object)}
  type = match.arg(type)
  ylim.boolean=hasArg(ylim)
  arg=list(...)
  #par(list(...))
  limit=FALSE
  y.sim=ResEnvelop=NULL
  link=object$link
  link.phi=object$link.phi
  y=object$y
  x=as.matrix(object$model$mean)
  z=as.matrix(object$model$precision)
  n=length(object$fitted.values$mu.predict)
  a=object$fitted.values$mu.predict*object$fitted.values$phi.predict
  b=(1-object$fitted.values$mu.predict)*object$fitted.values$phi.predict
  residual=residuals(object,type=type)
  k=1
  if(PrgBar){pb = txtProgressBar(min = 0, max = n.sim, style = 3)}
  while(k<=n.sim & !limit)
  {
    y.sim=pmax(pmin(sapply(seq(1,n,1),function(i) rbeta(1,a[i],b[i])),1-.Machine$double.eps),.Machine$double.eps)
    est.mle=betareg.fit(x,y,z,link=link,link.phi = link.phi)
    start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
    LMDPDE.sim=tryCatch(LMDPDE.Beta.Reg(y=y.sim,x=x,z=z,alpha=object$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){LMDPDE.sim$converged<-FALSE; return(LMDPDE.sim)})
    if(LMDPDE.sim$converged)
    {
      if(type=="sweighted2")
      {
        ResEnvelop=rbind(ResEnvelop,sort(LMDPDE.sim$residuals,decreasing = F))  
      }else{
        ResEnvelop=rbind(ResEnvelop,sort(residuals(LMDPDE.sim,type=type),decreasing = F))
      }
      k=k+1
      # update progress bar
      if(PrgBar){setTxtProgressBar(pb, k)}
    }
    if(k==(2*n.sim)){limit=T}
  }
  Envelope=apply(ResEnvelop,2,quantile,c((1-conf)/2,0.5,1-(1-conf)/2))
  if(!ylim.boolean)
  {
    ylim <- range(Envelope[1,],Envelope[2,],Envelope[3,],residual)
    arg$ylim=ylim
  }
  #browser()
  par(mar=c(5.0,5.0,4.0,2.0),pch=16, cex=1.0, cex.lab=1.0, cex.axis=1.0, cex.main=1.5)
  ARG=append(list(y=residual,main="Envelope Plot", xlab="Normal quantiles",ylab="Residuals"),arg)
  do.call(qqnorm,ARG)
  par(new=T)
  ARG=modifyList(ARG,list(y=Envelope[1,],axes=F,main = "",xlab="",ylab="",type="l",lty=1,lwd=1.0))
  do.call(qqnorm,ARG)
  par(new=T)
  ARG=modifyList(ARG,list(y=Envelope[2,],lty=2))
  do.call(qqnorm,ARG)
  par(new=T)
  ARG=modifyList(ARG,list(y=Envelope[3,],lty=1))
  do.call(qqnorm,ARG)
}

#' @rdname residuals.LSMLE
#'   
#' @export
residuals.LMDPDE=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection"))
{
  type = match.arg(type)
  y=object$y
  x=object$model$mean
  z=object$model$precision
  linkobj=make.link(link.mu=object$link,link.phi=object$link.phi)
  mu.predict=object$fitted.values$mu.predict
  phi.predict=object$fitted.values$phi.predict
  if(type=="sweighted2")
  {
    res=sweighted2_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,X=x,linkobj=linkobj)
  }
  if(type=="sweighted")
  {
    res=sweighted_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  }
  if(type=="pearson")
  {
    res=pearson_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  }
  if(type=="weighted")
  {
    res=weighted_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  }
  if(type=="sweighted.gamma")
  {
    res=sweighted.gamma_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  }
  if(type=="sweighted2.gamma")
  {
    res=sweighted2.gamma_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,Z=z,linkobj=linkobj)
  }
  if(type=="combined")
  {
    res=combined_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  }
  if(type=="combined.projection")
  {
    res=combined.projection_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,X=x,Z=z,linkobj=linkobj)
  }
  return(res)
}

#' @export
cooks.distance.LMDPDE=function(object,...)
{
  h =hatvalues.LMDPDE(object)
  k = length(object$coefficients$mean)
  res=residuals(object,type="pearson")
  return(h*(res^2)/(k*(1-h)^2))
}


  