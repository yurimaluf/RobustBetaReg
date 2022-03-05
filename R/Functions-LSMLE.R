#' @rdname robustbetareg
#' 
#' @param x,z  numeric regressor matrix for mean and precision model respectively, defaulting to an intercept only.
#' 
# Robust Estimation - LSMLE
LSMLE.fit=function(y,x,z,alpha=NULL,link="logit",link.phi="log",control=robustbetareg.control(...),...)
{
  #options(warn = 2) #Convert warnings in errors
  result=theta=list()
  #browser()
  #Arguments Checking
  alpha.optimal=control$alpha.optimal
  start_theta=control$start
  if(!is.null(alpha)){alpha.optimal=FALSE}
  if(is.null(alpha) & alpha.optimal==FALSE){alpha=0}
  linkobj=set.link(link.mu = link,link.phi = link.phi)
  k=ncol(x)
  m=ncol(z)
  if(alpha.optimal)
  {
    return(Opt.Tuning.LSMLE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
    theta$x=start_theta
    theta$converged=mle$converged
  }
  #Point Estimation
  #browser()
  q=1-alpha
  check=TRUE
  theta=tryCatch(nleqslv(start_theta,Psi_LSMLE_Cpp,jac=Psi_LSMLE_Jacobian_C,y=y,X=x,Z=z,alpha=alpha,link_mu=link,link_phi=link.phi,control=list(ftol=control$tolerance,maxit=control$maxit),method="Newton",jacobian = T),error=function(e){
    theta$msg<-e$message;check<<-F;return(theta)})
    if(check){
      if(all(abs(theta$fvec)<control$tolerance) & !all(theta$fvec==0) & all(diag(theta$jac)<0)){
        theta$converged=T
      }else{theta$converged=F}  
    }else{theta$x=start_theta;theta$converged=F}
  
  #Predict values
  beta=theta$x[1:k]
  gamma=theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  eta=x%*%beta
  mu_hat=linkobj$linkfun.mu$inv.link(eta)
  phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)
  
  #Expected Standard Error
  M.LSMLE=tryCatch(LSMLE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj),error=function(e) NULL)
  if(is.null(M.LSMLE))
  {
    vcov=std.error.LSMLE=NaN
  }else{
    vcov=M.LSMLE$Cov
    std.error.LSMLE=M.LSMLE$Std.Error
  }
  
  #Register of output values 
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients
  result$vcov=vcov
  pseudor2 = if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LSMLE)))
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start=start_theta #Alterar caso optar por duas tentativas de ponto inicial
  result$weights=(degbeta(y_star,mu_hat,phi_hat/q))^(alpha) #Pesos 
  result$Tuning=alpha
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)
  result$n=length(mu_hat)
  result$link=link
  result$link.phi=link.phi
  result$Optimal.Tuning=alpha.optimal
  result$pseudo.r.squared=pseudor2
  result$control=control
  if(any(is.na(std.error.LSMLE)))
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
  if(!any(is.na(std.error.LSMLE)))
  {
    names(std.error.LSMLE)<-c(colnames(x),colnames(z))
    se.beta=std.error.LSMLE[1:k]
    se.gamma=std.error.LSMLE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
  }
  class(result)="LSMLE"
  #gc()
  return(result)
}

#' @export
robustbetareg.control.LSMLE=function(object)
{
  result=object$control
  return(result)
}

#' @export
WaldTypeTest.LSMLE=function(object,FUN,...)
{
  general=FALSE
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
    result=list(beta.wald=result.beta,gamma.wald=result.gamma,general=general,msg="Results based on LSMLE.")
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
    result$msg="Results based on LSMLE."
  }
  class(result)="WaldTest_LSMLE"
  return(result)
}

#' @export
SaddlepointTest.LSMLE=function(object,FUN=NULL,...,thrd)
{
  if(!object$converged){stop(paste("There is no convergence in the model",deparse(substitute(object))))}
  result=list()
  msg="Results based on LSMLE."
  beta_hat=object$coefficient$mean
  if(missing(thrd)){thrd=tryCatch(parallel::detectCores(),error=function(e) 1)}
  else{thrd=tryCatch(suppressWarnings(max(1,round(abs(thrd)))),error=function(e) 1)}
  gamma_hat=object$coefficient$precision
  X=object$model$mean
  Z=object$model$precision
  m=ncol(X)
  k=ncol(Z)
  l0=rep(0,m+k)
  B_0=G_0=general=NULL
  N=object$n
  linkobj=set.link(object$link,object$link.phi)
  alpha=object$Tuning
  if(is.null(FUN))
  {
    general=T
    #Beta
    ind_fix=1:m
    ind_free=(m+1):(m+k)
    h_beta=tryCatch(suppressWarnings(stats::optim(par=gamma_hat,sup_K_psi_C,ind_free=ind_free,ind_fix=ind_fix,eta_0=beta_hat,Beta=beta_hat,Gamma=gamma_hat,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd,method="BFGS")$value),error=function(e) NULL)
    if(is.null(h_beta)){msg="The test does not reach the convergence"}
    #Result Register
    result$SaddlePointTest=as.numeric(2*N*h_beta)
    result$df=m
    result$pValue=as.numeric(1-pchisq(2*N*h_beta,df=m))
  }else{
    general=F
    result=list()
    #Hipothesis
    g=FUN
    eta_0=g(c(beta_hat,gamma_hat),...)
    ind_free=which(c(beta_hat,gamma_hat)==eta_0)
    ind_fix=which(c(beta_hat,gamma_hat)!=eta_0)
    if(length(ind_fix)<(m+k))
    {
      h_beta=tryCatch(suppressWarnings(stats::optim(par=c(beta_hat,gamma_hat)[ind_free],sup_K_psi_C,ind_free=ind_free,ind_fix=ind_fix,eta_0=eta_0[ind_fix],Beta=beta_hat,Gamma=gamma_hat,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd,method="BFGS")$value),error=function(e) NULL)
      if(is.null(h_beta) || is.infinite(h_beta)){
        h_beta=tryCatch(suppressWarnings(stats::optim(par=c(beta_hat,gamma_hat)[ind_free],sup_K_psi,ind_free=ind_free,ind_fix=ind_fix,eta_0=eta_0[ind_fix],Beta=beta_hat,Gamma=gamma_hat,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd,method="BFGS")$value),error=function(e) NULL)
      }
      if(is.null(h_beta)){msg="The test does not reach the convergence"}
      df=length(ind_fix)
    }else{
      B_0=eta_0[1:m]
      G_0=eta_0[(m+1):(m+k)]
      df=m+k
      h_beta=tryCatch(suppressWarnings(-(stats::nlminb(l0,K2_psi_C,Beta=beta_hat,Gamma=gamma_hat,Beta_0=B_0,Gamma_0=G_0,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd))$objective),error=function(e) NULL)
      if(is.null(h_beta) || is.infinite(h_beta)){
        h_beta=tryCatch(suppressWarnings(-(stats::nlminb(l0,K2_psi,Beta=beta_hat,Gamma=gamma_hat,Beta_0=B_0,Gamma_0=G_0,X=X,Z=Z,alpha=alpha,linkobj=linkobj,thrd=thrd))$objective),error=function(e) NULL)
      }
      if(is.null(h_beta)){msg="The test does not reach the convergence"}
    }
    #Result Register
    result$SaddlePointTest=as.numeric(2*N*h_beta)
    result$df=df
    result$pValue=as.numeric(1-pchisq(2*N*h_beta,df=df))
  }
  result$msg=msg
  result$thrd=thrd
  result$general=general
  class(result)="SaddlepointTest_LSMLE"
  return(result)
}


#' @export  
plotenvelope.LSMLE=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"),conf=0.95,n.sim=50,PrgBar=T,control=robustbetareg.control(...), ...)
{
  #browser()
  if(missing(control)){control=robustbetareg.control(object)}
  type = match.arg(type)
  ylim.boolean=hasArg(ylim)
  arg=list(...)
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
    #start=c(3.5,-3.4,6)
    LSMLE.sim=tryCatch(LSMLE.fit(y=y.sim,x=x,z=z,alpha=object$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){LSMLE.sim$converged<-FALSE; return(LSMLE.sim)})
    if(LSMLE.sim$converged)
    {
      if(type=="sweighted2")
      {
        ResEnvelop=rbind(ResEnvelop,sort(LSMLE.sim$residuals,decreasing = F))  
      }else{
        #browser()
        LSMLE.sim$y=y.sim
        LSMLE.sim$model$mean=object$model$mean
        LSMLE.sim$model$precision=object$model$precision
        ResEnvelop=rbind(ResEnvelop,sort(residuals(LSMLE.sim,type=type),decreasing = F))
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
  par(mar=c(5.0,5.0,4.0,2.0),pch=16, cex=1.0, cex.lab=1.0, cex.axis=1.0, cex.main=1.5)
  #ARG=append(list(y=residual,main="Envelope Plot", xlab="Normal quantiles",ylab="Residuals"),arg)
  ARG=append(list(y=residual,main="", xlab="Normal quantiles",ylab="Residuals"),arg)
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

#' Residuals Method for robustbetareg Objects
#'  
#'  
#'    
#' 
#' Extract various types of residuals from  robust beta regression models: Pearson residuals (raw residuals scaled by square root of variance function)
#' and different kinds of weighted residuals suggested by Espinheira et al. (2008) and Espinheira et al. (2017).
#' 
#' 
#' 
#' @param object fitted model object of class "LMDPDE" or "LSMLE".
#' @param type character indicating type of residuals.  
#' 
#' @details The definitions of the first four residuals are provided in Espinheira et al. (2008):  Equation 2 for "\code{pearson}", Equation 6 for "\code{weighted}", Equation 7 for "\code{sweighted}", and Equation 8 for "\code{sweighted2}".
#' For the last four residuals the definitions are described in Espinheira et al. (2017): Last equation of Equation 7 and Equation 10 for "\code{sweighted.gamma}" and "\code{sweighted2.gamma}" respectively, Equation 9 for "\code{combined}", and Equation 11 for "\code{combined.projection}".     
#' 
#' @references  \href{https://www.tandfonline.com/doi/abs/10.1080/02664760701834931}{Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008). On Beta Regression Residuals. Journal of Applied Statistics, 35(4), 407–419.}
#' @references  \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201600136}{Espinheira, P.L., Santos, E.G.and Cribari-Neto, F. (2017). On nonlinear beta regression residuals. Biometrical Journal, 59(3), 445-461.} 
#' 
#' @examples 
#' fit=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' residuals(fit,type="sweighted")
#' 
#' @export
residuals.LSMLE=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection","sweighted3"))
{
  type = match.arg(type)
  y=object$y
  x=object$model$mean
  z=object$model$precision
  linkobj=set.link(link.mu=object$link,link.phi=object$link.phi)
  mu.predict=object$fitted.values$mu.predict
  phi.predict=object$fitted.values$phi.predict
  switch(type,sweighted2={
    res=sweighted2_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,X=x,linkobj=linkobj)
  },
  sweighted={
    res=res=sweighted_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  pearson={
    res=pearson_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  weighted={
    res=weighted_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  sweighted.gamma={
    res=sweighted.gamma_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  sweighted2.gamma={
    res=sweighted2.gamma_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,Z=z,linkobj=linkobj)
  },
  combined={
    res=combined_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y)
  },
  combined.projection={
    res=combined.projection_res(mu_hat=mu.predict,phi_hat=phi.predict,y=y,X=x,Z=z,linkobj=linkobj)
  },
  sweighted3={
    res=sweighted3_res(mu_hat=mu.predict,phi_hat=phi.predict,alpha=object$Tuning,y=y,X=x,linkobj=linkobj)
  },
  
  stop(gettextf("%s residual not recognised", sQuote(type)),
       domain = NA))
  
  return(res)
}

#' @export
cooks.distance.LSMLE=function(object,...)
{
  # h =hatvalues.LSMLE(object)
  # k = length(object$coefficients$mean)
  # res=residuals(object,type="pearson")
  # return(h*(res^2)/(k*(1-h)^2))
  #browser()
  p=length(do.call("c",object$coefficients))
  linkobj=set.link(link.mu = object$link, link.phi = object$link.phi)
  y_hat=linkobj$linkfun.mu$inv.link(object$model$mean%*%object$coefficients$mean)
  MSE=as.numeric(t(object$y-y_hat)%*%(object$y-y_hat)/(object$n-p))
  D=NULL
  for(i in 1:object$n){
    fit.temp=robustbetareg(object$formula,data=object$data[-i,],alpha=object$Tuning)
    y_hat_temp=linkobj$linkfun.mu$inv.link(object$model$mean%*%fit.temp$coefficients$mean)
    D=c(D,t(y_hat-y_hat_temp)%*%(y_hat-y_hat_temp)/(MSE*p))
  }
  return(D)
}

