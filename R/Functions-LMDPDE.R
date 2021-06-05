#' Robust Estimation - LMDPDE
#' 
#' Fit robust beta regression models for rates and proportions via LMDPDE.
#' 
#' For more details see:...
#' 
#' @param y The numeric response vector (with values in (0,1)).
#' @param x The numeric regressor matrix for mean model.
#' @param z The numeric regressor matrix for precision model, defaulting to an intercept only.
#' @param alpha The tuning with values (0,1), for robust estimation. When alpha is equal zero is equivalent of MLE. 
#' @param start_theta A numeric vector with an initial guess of the root of estimation equation.
#' @param alpha.optimal A logical value. If TRUE the tuning parameter should be selected automatic.
#' 
#' @return Return a list of components:
#'  \itemize{
#'   \item coefficients - the vector with elements "mean" and "precision" containing the coefficients from the respective models,
#'   \item vcov - the covariance matrix of all parameters in the model,
#'   \item mu_hat - the vector of fitted means,
#'   \item phi_hat - the vector of predicted precision,
#'   \item weight - the weights used,
#'   \item Tuning - the employed tuning parameter,
#'   \item Res.Beta - a vector of standardized weighted residual 2,
#'   \item start - the starting values for the parameters passed to the Newton-Raphson algorithm,
#'   \item std.error - the standard error vector of coefficients,
#'   \item Tab - the summary table of overall result.
#' }
#'
#' @export 
LMDPDE.Beta.Reg=function(y,x,z,alpha,start_theta,alpha.optimal,control)
{
  result=theta=list()
  #browser()
  #Arguments Checking
  if(missing(control)){control=robustbetareg.control()}#Enquanto essa funcao for public, qnd private deletar
  if(missing(alpha.optimal)){alpha.optimal=TRUE}
  if(!missing(alpha)){alpha.optimal=FALSE}
  ocontrol = control
  k=ncol(x)
  m=ncol(z)
  if(alpha.optimal)
  {
    return(Opt.Tuning.LMDPDE(y,x,z,control))
  }
  if(missing(alpha))
  {
    if(m==1)
    {
      mle=suppressWarnings(betareg(y~x[,-1]|1))  
    }
    else{mle=suppressWarnings(betareg(y~x[,-1]|z[,-1]))}
    theta$x=as.numeric(c(mle$coefficients$mean,mle$coefficients$precision))
    alpha=0
    beta=theta$x[1:k]
    gamma=theta$x[seq.int(length.out = m) + k]
    coefficients = list(mean = beta, precision = gamma)
    mu_hat=h1(x,beta)
    phi_hat=h2(z,gamma)
    theta$inter=NULL
    theta$converged=mle$converged
    theta$fvec=Psi_LMDPDE(theta$x,y=y,X=x,Z=z,alpha=0)
    vcov=mle$vcov
    std.error.LMDPDE=sqrt(diag(vcov))
  }else{
    #Point Estimation
    theta=Robst.LMDPDE.Beta.Reg(y,x,z,start_theta=start_theta,alpha=alpha,tolerance=control$tolerance,maxit=control$maxit)
    # if(!theta$converged)
    # {
    #   initial.point.rbst=Initial.points(y,x,z)
    #   theta=Robst.LMDPDE.Beta.Reg(y,x,z,start_theta=initial.point.rbst,alpha=alpha,tolerance=tolerance,maxit=maxit)
    #   initial.point=initial.point.rbst
    # }
    #Predict values
    beta=theta$x[1:k]
    gamma=theta$x[seq.int(length.out = m) + k]
    coefficients = list(mean = beta, precision = gamma)
    mu_hat=h1(x,beta)
    phi_hat=h2(z,gamma)
    #Expected Standard Error
    MM=LMDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha)
    vcov=MM$Cov
    std.error.LMDPDE=MM$Std.Error
  }
  
  #Register of output values 
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients
  result$vcov=vcov
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LMDPDE)))
  result$mu_predict=mu_hat
  if(m==1){result$phi_predict=phi_hat[1]}
  if(m!=1){result$phi_predict=phi_hat}
  result$start=start_theta #Alterar caso optar por duas tentativas de ponto inicial
  result$weights=degbeta(y_star,mu_hat,phi_hat)^(alpha)#Weights
  result$Tuning=alpha
  result$Psi.Value=theta$fvec
  result$Res.Beta=Residual_Beta(mu_hat,phi_hat,y=y,X=x,Z=z)
  result$x=list(mean = x, precision = z)
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
  #Assembly of summary table 
  coef.b=coef.g=b=g=obs.b=obs.g=NULL
  variable=colnames(x)
  variable2=colnames(z)
  variable[1]="Intercept"
  variable2[1]="Intercept Disp"
  if(!any(is.na(std.error.LMDPDE)))
  {
    for(i in 1:k)
    {
      coef=paste0("b",(i-1))
      coef.b=c(coef.b,paste0("b",i))
      p.valor=2-2*pnorm(abs(beta[i]/std.error.LMDPDE[i]))
      obs=star.obs(p.valor)
      if(p.valor<2e-16){p.valor="<2e-16"}
      obs.b=c(obs.b,obs)
      b_=formatC(c(formatC(beta[i]),formatC(std.error.LMDPDE[i]),formatC(beta[i]/std.error.LMDPDE[i]),formatC(p.valor)))
      b=rbind(b,c(variable[i],coef,b_))
    }
    for(i in 1:m)
    {
      coef=paste0("g",(i-1))
      coef.g=c(coef.g,paste0("g",i))
      p.valor=2-2*pnorm(abs(gamma[i]/std.error.LMDPDE[i+k]))
      obs=star.obs(p.valor)
      if(p.valor<2e-16){p.valor="<2e-16"}
      obs.g=c(obs.g,obs)
      g_=formatC(c(formatC(gamma[i]),formatC(std.error.LMDPDE[i+k]),formatC(gamma[i]/std.error.LMDPDE[i+k]),formatC(p.valor)))
      g=rbind(g,c(variable2[i],coef,g_))
    }
    bg=as.data.frame(rbind(b,g))
    bg=cbind(bg,c(obs.b,obs.g))
    bg=format.data.frame(bg,trim=T,width=0.1)
    
    colnames(bg)=c("Variable","Coef","Estimate","Std.Error","z-Value","Pr(>|z|)","")
    
    names(std.error.LMDPDE)<-c(coef.b,coef.g)
    se.beta=std.error.LMDPDE[1:k]
    se.gamma=std.error.LMDPDE[1:m+k]
    coef.std.error = list(se.mean = se.beta, se.precision = se.gamma)
    result$std.error=coef.std.error
    result$Tab=bg 
  }
  class(result)="LMDPDE"
  return(result)
}


#' @export
WaldTypeTest.LMDPDE=function(object,g,...)
{
  result=list()
  #browser()
  p=c(object$coefficient$mean,object$coefficient$precision)
  M=numDeriv::jacobian(g,p,...)
  m=g(p,...)
  r=length(m)
  n=length(object$mu_predict)
  if(Matrix::rankMatrix(M)[1]!=r)
  {
    #return(result$msg="The Rank Matrix is not supported")
  }
  if(!object$converged)
  {
    ###
  }
  V=object$vcov
  W_alpha=t(m)%*%solve(M%*%V%*%t(M))%*%(m)
  result$converged=object$converged
  result$W.alpha=as.numeric(W_alpha)
  result$df=r
  result$pValue=as.numeric(1-pchisq(W_alpha,df=r))
  result$msg="Results based on LMDPDE."
  
  return(result)
}


#' @export
plotenvelope.LMDPDE=function(robustbetareg.obj,ylim,n.sim,index,control)
{
  if(missing(n.sim)){n.sim=100}
  if(missing(control)){control=robustbetareg.control()}
  if(missing(index)){index=NULL}
  if(missing(ylim)){ylim=NULL}
  y.sim=ResEnvelop=NULL
  x=as.matrix(robustbetareg.obj$x$mean)
  z=as.matrix(robustbetareg.obj$x$precision)
  n=length(robustbetareg.obj$mu_predict)
  a=robustbetareg.obj$mu_predict*robustbetareg.obj$phi_predict
  b=(1-robustbetareg.obj$mu_predict)*robustbetareg.obj$phi_predict
  for(i in 1:n)
  {
    y.temp=pmax(pmin(rbeta(n.sim,a[i],b[i]),1-.Machine$double.eps),.Machine$double.eps)
    y.sim=cbind(y.sim,y.temp)
  }
  for(i in 1:n.sim)
  {
    LMDPDE.sim=LMDPDE.Beta.Reg(y=y.sim[i,],x=x,z=z,alpha=robustbetareg.obj$Tuning,start_theta=robustbetareg.obj$start,control=control)
    if(LMDPDE.sim$converged)
    {
      res_beta=sort(LMDPDE.sim$Res.Beta,decreasing = F)
      ResEnvelop=rbind(ResEnvelop,res_beta)
    }
  }
  Envelope=apply(ResEnvelop,2,quantile,c(0.025,0.5,1-0.025))
  if(is.null(ylim))
  {
    ylim <- range(Envelope[1,],Envelope[2,],Envelope[3,],robustbetareg.obj$Res.Beta)  
  }
  par(mar=c(5.0,5.0,4.0,2.0))
  reg=qqnorm(robustbetareg.obj$Res.Beta, main="", xlab="Normal quantiles", ylab="Residuals", ylim=ylim, pch=16, cex=1.0, cex.lab=1.0, cex.axis=1.0, cex.main=1.0)
  if(!is.null(index))
  {
    points(reg$x[index],reg$y[index],pch=20)
    text(reg$x[index],reg$y[index],paste0("(",index,")"),pos = c(1))  
  }
  par(new=T)
  qqnorm(Envelope[1,],axes=F,main = "",xlab="",ylab="",type="l",ylim=ylim,lty=1,lwd=1.0)
  par(new=T)
  qqnorm(Envelope[2,],axes=F,main = "",xlab="",ylab="", type="l",ylim=ylim,lty=2,lwd=1.0)
  par(new=T)
  qqnorm(Envelope[3,],axes=F,xlab="",main = "", ylab="", type="l",ylim=ylim,lty=1,lwd=1.0)
}
