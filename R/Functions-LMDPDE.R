# Robust Estimation - LMDPDE
LMDPDE.Beta.Reg=function(y,x,z,alpha,link,link.phi,control=robustbetareg.control(...), ...)
{
  result=theta=list()
  #Arguments Checking
  alpha.optimal=control$alpha.optimal
  start_theta=control$start
  if(!missing(alpha) & !is.null(alpha)){alpha.optimal=FALSE}
  if(is.null(alpha) & alpha.optimal==FALSE){alpha=0}
  if(missing(link)){link="logit"}
  if(missing(link.phi)){link.phi="log"}
  linkobj=make.link(link.mu = link,link.phi = link.phi)
  k=ncol(x)
  m=ncol(z)
  if(alpha.optimal)
  {
    return(Opt.Tuning.LMDPDE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta) || missing(alpha))
  {
    if(is.null(start_theta) || missing(alpha))
    {
      mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
      start_theta=as.numeric(do.call("c",mle$coefficients))
    }
  }
  if(missing(alpha))
  {
    theta$x=as.numeric(c(mle$coefficients$mean,mle$coefficients$precision))
    alpha=0
    beta=theta$x[1:k]
    gamma=theta$x[seq.int(length.out = m) + k]
    coefficients = list(mean = beta, precision = gamma)
    mu_hat=linkobj$linkfun.mu$inv.link(x%*%beta)
    phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)
    theta$inter=NULL
    theta$converged=mle$converged
    theta$fvec=Psi_LMDPDE(theta$x,y=y,X=x,Z=z,alpha=0,linkobj=linkobj)
    vcov=mle$vcov
    std.error.LMDPDE=sqrt(diag(vcov))
  }else{
    #Point Estimation
    theta=Robst.LMDPDE.Beta.Reg(y,x,z,start_theta=start_theta,alpha=alpha,linkobj=linkobj,tolerance=control$tolerance,maxit=control$maxit)
    #Predict values
    beta=theta$x[1:k]
    gamma=theta$x[seq.int(length.out = m) + k]
    coefficients = list(mean = beta, precision = gamma)
    mu_hat=linkobj$linkfun.mu$inv.link(x%*%beta)
    phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)
    #Expected Standard Error
    MM=LMDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj)
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
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start=start_theta #Alterar caso optar por duas tentativas de ponto inicial
  result$weights=degbeta(y_star,mu_hat,phi_hat)^(alpha)#Weights
  result$Tuning=alpha
  result$Psi.Value=theta$fvec
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)
  result$model=list(mean = x, precision = z)
  result$y=y
  result$link=link
  result$link.phi=link.phi
  result$Optimal.Tuning=alpha.optimal
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
  if(!object$converged){stop(Paste("There is no convergence in the model",deparse(substitute(object))))}
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
plotenvelope.LMDPDE=function(robustbetareg.obj,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection"),conf=0.95,n.sim=50,control=robustbetareg.control(...), ...)
{
  if(missing(control)){control=robustbetareg.control(robustbetareg.obj)}
  type = match.arg(type)
  ylim.boolean=hasArg(ylim)
  limit=FALSE
  y.sim=ResEnvelop=NULL
  link=robustbetareg.obj$link
  link.phi=robustbetareg.obj$link.phi
  y=robustbetareg.obj$y
  x=as.matrix(robustbetareg.obj$model$mean)
  z=as.matrix(robustbetareg.obj$model$precision)
  n=length(robustbetareg.obj$fitted.values$mu.predict)
  a=robustbetareg.obj$fitted.values$mu.predict*robustbetareg.obj$fitted.values$phi.predict
  b=(1-robustbetareg.obj$fitted.values$mu.predict)*robustbetareg.obj$fitted.values$phi.predict
  residual=residuals(robustbetareg.obj,type=type)
  k=1
  pb = txtProgressBar(min = 0, max = n.sim, style = 3)
  #browser()
  while(k<=n.sim & !limit)
  {
    y.sim=pmax(pmin(sapply(seq(1,n,1),function(i) rbeta(1,a[i],b[i])),1-.Machine$double.eps),.Machine$double.eps)
    est.mle=betareg.fit(x,y,z,link=link,link.phi = link.phi)
    start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
    LMDPDE.sim=tryCatch(LMDPDE.Beta.Reg(y=y.sim,x=x,z=z,alpha=robustbetareg.obj$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){LMDPDE.sim$converged<-FALSE; return(LMDPDE.sim)})
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
      setTxtProgressBar(pb, k)
    }
    if(k==(2*n.sim)){limit=T}
  }
  Envelope=apply(ResEnvelop,2,quantile,c((1-conf)/2,0.5,1-(1-conf)/2))
  if(!ylim.boolean)
  {
    ylim <- range(Envelope[1,],Envelope[2,],Envelope[3,],residual)
  }
  par(mar=c(5.0,5.0,4.0,2.0))
  qqnorm(residual, main="Envelope Plot", xlab="Normal quantiles",ylab="Residuals",ylim=ylim, pch=16, cex=1.0, cex.lab=1.0, cex.axis=1.0, cex.main=1.0,...)
  par(new=T)
  qqnorm(Envelope[1,],axes=F,main = "",xlab="",ylab="",ylim=ylim,type="l",lty=1,lwd=1.0,...)
  par(new=T)
  qqnorm(Envelope[2,],axes=F,main = "",xlab="",ylab="",ylim=ylim, type="l",lty=2,lwd=1.0,...)
  par(new=T)
  qqnorm(Envelope[3,],axes=F,xlab="",main = "", ylab="",ylim=ylim, type="l",lty=1,lwd=1.0,...)
}

#' @export
print.LMDPDE=function(obj)
{
  cat("Call: \n")      
  print(obj$call)
  cat("\n")
  cat("Coefficients (mean model with",obj$link,"link):\n")
  print(obj$coefficients$mean)
  cat("\n")
  cat("Coefficients (precision model with",obj$link.phi,"link):\n")
  print(obj$coefficients$precision)
}

#' @export
print.WaldTest_LMDPDE=function(obj)
{
  if(obj$general)
  {
    #browser()
    cat("-- Wald Type Test -- \n")
    if(!is.null(obj$beta.wald))
    {
      p.valor=obj$beta.wald$pValue
      obs=star.obs(p.valor)
      if(p.valor<=2e-16){p.valor="<2e-16"}
      if(p.valor>2e-16){p.valor=paste0("=",obj$beta.wald$pValue)}
      
      cat("Null Hypothesis: all mean coefficients equal to zero \n")  
      cat(paste0("Value=",formatC(obj$beta.wald$W.alpha),", df=",obj$beta.wald$df,", p-Value",p.valor,obs,"\n"))  
    }
    if(!is.null(obj$gamma.wald))
    {
      p.valor=obj$gamma.wald$pValue
      obs=star.obs(p.valor)
      if(p.valor<=2e-16){p.valor="<2e-16"}
      if(p.valor>2e-16){p.valor=paste0("=",obj$gamma.wald$pValue)}
      
      cat("Null Hypothesis: all precision coefficients equal to zero \n")  
      cat(paste0("Value=",formatC(obj$gamma.wald$W.alpha),", df=",obj$gamma.wald$df,", p-Value",p.valor,obs,"\n"))  
    }
  }else{
    cat("-- Wald Type Test -- \n")
    p.valor=obj$pValue
    obs=star.obs(p.valor)
    if(p.valor<=2e-16){p.valor="<2e-16"}
    if(p.valor>2e-16){p.valor=paste0("=",obj$pValue)}
    cat("Null Hypothesis: set by the user \n") 
    cat(paste0("Value=",formatC(obj$W.alpha),", df=",obj$df,", p-Value",p.valor,obs,"\n")) 
  }
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  cat("Results based on LMDPDE \n")
}


#' @export
summary.LMDPDE=function(obj)
{
  b=g=obs.b=obs.g=NULL
  beta=obj$coefficients$mean
  gamma=obj$coefficients$precision
  variable=names(obj$coefficients$mean)
  variable2=names(obj$coefficients$precision)
  std.error.beta=obj$std.error$se.mean
  std.error.gamma=obj$std.error$se.precision
  k=length(beta)
  m=length(gamma)
  for(i in 1:k)
  {
    p.valor=2-2*pnorm(abs(beta[i]/std.error.beta[i]))
    obs=star.obs(p.valor)
    if(p.valor<2e-16){p.valor="<2e-16"}
    obs.b=c(obs.b,obs)
    b_=formatC(c(formatC(beta[i]),formatC(std.error.beta[i]),formatC(beta[i]/std.error.beta[i]),formatC(p.valor)))
    b=rbind(b,c(variable[i],b_))
  }
  b.df=as.data.frame(b)
  b.df=cbind(b.df,obs.b)
  b.df=format.data.frame(b.df,trim=T,width=0.1)
  colnames(b.df)=c("","Estimate","Std. Error", "z value", "Pr(>|z|)","")
  for(i in 1:m)
  {
    p.valor=2-2*pnorm(abs(gamma[i]/std.error.gamma[i]))
    obs=star.obs(p.valor)
    if(p.valor<2e-16){p.valor="<2e-16"}
    obs.g=c(obs.g,obs)
    g_=formatC(c(formatC(gamma[i]),formatC(std.error.gamma[i]),formatC(gamma[i]/std.error.gamma[i]),formatC(p.valor)))
    g=rbind(g,c(variable2[i],g_))
  }
  g.df=as.data.frame(g)
  g.df=cbind(g.df,obs.g)
  g.df=format.data.frame(g.df,trim=T,width=0.1)
  colnames(g.df)=c("","Estimate","Std. Error", "z value", "Pr(>|z|)","")
  
  cat("Call: \n")      
  print(obj$call)
  cat("\n")
  cat("Coefficients (mean model with",obj$link,"link):\n")
  print(b.df,row.names=FALSE)
  cat("\n")
  cat("Phi coefficients (precision model with",obj$link.phi,"link):\n")
  print(g.df,row.names=FALSE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  cat("\n")
  cat("Type of estimator: LMDPDE \n")
  cat(paste0("Tuning value: alpha=",obj$Tuning,"\n"))
  if(obj$Optimal.Tuning)
  {
    cat("Tuning constant generated by the data-driven algorithm")  
  }
  if(!obj$Optimal.Tuning)
  {
    cat("Tuning constant selected by the user")  
  }
}


#' Interactive plots for diagnostic of robust betareg models
#'  
#'  
#'    
#' 
#' Several types of standard diagnostic plots can be produced interactively, involving various kinds of residuals, influence measures, weights etc. 
#' 
#' 
#' 
#' @param object fitted model object of class "LMDPDE" or "LSMLE".
#' @param ask logical. If TRUE the user is asked before each plot.
#' @param ... other parameters to be passed through to plotting functions.
#' 
#' #' @examples 
#' rbr=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' plot(rbr)
#' 
#' @export
plot.LMDPDE=function(object,ask=TRUE,...)
{
  getinfo=Sys.info()
  user=getinfo[which(names(getinfo)=="user")]
  text.main2="the graph number >\n [1] Residuals \n [2] Residuals x Linear predictor \n [3] Cook's Distance \n [4] Weights \n [5] Weigths x Residuals \n [0] Exit \n"
  text.main=paste("Select",text.main2)
  if(!is.na(user))
  {
    user=paste0(toupper(substring(user,1,1)),substring(user,2))
    text.main=paste0(user,", select ",text.main2)
  }
  text.n1="Select the residual number: \n [1] sweighted2 \n [2] sweighted \n [3] pearson \n [4] weighted \n [5] sweighted.gamma \n [6] sweighted2.gamma \n [7] combined \n [8] combined.projection \n [9] back \n [0] exit \n"
  text.n2="Select the residual type to match with linear predictor: \n [1] sweighted2 \n [2] sweighted \n [3] pearson \n [4] weighted \n [5] sweighted.gamma \n [6] sweighted2.gamma \n [7] combined \n [8] combined.projection \n [9] back \n [0] exit \n"
  text.n5="Select the residual type to match with weights: \n [1] sweighted2 \n [2] sweighted \n [3] pearson \n [4] weighted \n [5] sweighted.gamma \n [6] sweighted2.gamma \n [7] combined \n [8] combined.projection \n [9] back \n [0] exit \n"
  show=TRUE
  while(show){
    show.1=show.2=show.5=TRUE
    rstudioapi::sendToConsole("",execute=F,focus=T,echo=T)
    n <- as.numeric(readline(cat(crayon::green(text.main))))
    if(n==1)
    {
      while(show & show.1)
      {
        m <-readline(prompt = cat(crayon::green(text.n1)))
        if(m==1)
        {
          plot(residuals(object,type="sweighted2"),xlab="Obs. number",ylab="Standardized Weighted 2 Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==2)
        {
          plot(residuals(object,type="sweighted"),xlab="Obs. number",ylab="Standardized Weighted Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==3)
        {
          plot(residuals(object,type="pearson"),xlab="Obs. number",ylab="Pearson Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==4)
        {
          plot(residuals(object,type="weighted"),xlab="Obs. number",ylab="Weighted Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==5)
        {
          plot(residuals(object,type="sweighted.gamma"),xlab="Obs. number",ylab="Standardized Weighted Gamma Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==6)
        {
          plot(residuals(object,type="sweighted2.gamma"),xlab="Obs. number",ylab="Standardized Weighted 2 Gamma Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==7)
        {
          plot(residuals(object,type="combined"),xlab="Obs. number",ylab="Combined Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==8)
        {
          plot(residuals(object,type="combined.projection"),xlab="Obs. number",ylab="Combined Projection Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
        }
        if(m==0){show=FALSE}
        if(m==9){show.1=FALSE}
      }
    }
    if(n==2)
    {
      while(show & show.2)
      {
        m <-readline(prompt = cat(crayon::green(text.n2)))
        if(m==1)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted2"),xlab="Linear predictor",ylab="Standardized Weighted 2 Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==2)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted"),xlab="Linear predictor",ylab="Standardized Weighted Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==3)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="pearson"),xlab="Linear predictor",ylab="Pearson Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==4)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="weighted"),xlab="Linear predictor",ylab="Weighted Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==5)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted.gamma"),xlab="Linear predictor",ylab="Standardized Weighted Gamma Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==6)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="sweighted2.gamma"),xlab="Linear predictor",ylab="Standardized Weighted 2 Gamma Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==7)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="combined"),xlab="Linear predictor",ylab="Combined Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==8)
        {
          plot(object$model$mean%*%object$coefficients$mean,residuals(object,type="combined.projection"),xlab="Linear predictor",ylab="Combined Projection Residual",main="Residuals vs linear predictor",...)
          abline(h=0)
        }
        if(m==0){show=FALSE}
        if(m==9){show.2=FALSE}
      }
    }
    if(n==3)
    {
      plot(cooks.distance(object),type="h",xlab="Obs. number",ylab="Cook's distance",main="Cook's distance plot",...)
    }
    if(n==4)
    {
      plot(object$weights,xlab="Obs. number",ylab="Weights",main = "Weights plot",...)
      abline(h=0)
    }
    if(n==5)
    {
      while(show & show.5)
      {
        m <-readline(prompt = cat(crayon::green(text.n5)))
        if(m==1)
        {
          plot(residuals(object,type="sweighted2"),object$weights,ylab="Weights",xlab="Standardized Weighted 2 Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==2)
        {
          plot(residuals(object,type="sweighted"),object$weights,ylab="Weights",xlab="Standardized Weighted Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==3)
        {
          plot(residuals(object,type="pearson"),object$weights,ylab="Weights",xlab="Pearson Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==4)
        {
          plot(residuals(object,type="weighted"),object$weights,ylab="Weights",xlab="Weighted Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==5)
        {
          plot(residuals(object,type="sweighted.gamma"),object$weights,ylab="Weights",xlab="Standardized Weighted Gamma Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==6)
        {
          plot(residuals(object,type="sweighted2.gamma"),object$weights,ylab="Weights",xlab="Standardized Weighted 2 Gamma Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==7)
        {
          plot(residuals(object,type="combined"),object$weights,ylab="Weights",xlab="Combined Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==8)
        {
          plot(residuals(object,type="combined.projection"),object$weights,ylab="Weights",xlab="Combined Projection Residual",main="Weights vs residuals",...)
          abline(h=0)
        }
        if(m==0){show=FALSE}
        if(m==9){show.5=FALSE}
      }
    }
    if(n==0)
    {
      show=FALSE
    }
  }
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
  h =hatvalues.robustbetareg(object)
  k = length(object$coefficients$mean)
  res=residuals(object,type="pearson")
  return(h*(res^2)/(k*(1-h)^2))
}

  