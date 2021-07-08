# Robust Estimation - LSMLE
LSMLE.Beta.Reg=function(y,x,z,alpha,link,link.phi,control=robustbetareg.control(...), ...)
{
  #options(warn = 2) #Convert warnings in errors
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
    return(Opt.Tuning.LSMLE(y,x,z,link,link.phi,control))
  }
  if(is.null(start_theta) || missing(alpha))
  {
    mle=tryCatch(suppressWarnings(betareg.fit(x,y,z,link=link,link.phi = link.phi)),error=function(e) NULL)
    start_theta=as.numeric(do.call("c",mle$coefficients))
  }
  if(missing(alpha))
  {
    theta$x=as.numeric(do.call("c",mle$coefficients))
    alpha=0
    q=1-alpha
    alpha=0
    beta=theta$x[1:k]
    gamma=theta$x[seq.int(length.out = m) + k]
    coefficients = list(mean = beta, precision = gamma)
    eta=x%*%beta
    mu_hat=linkobj$linkfun.mu$inv.link(eta)
    phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)
    theta$inter=NULL
    theta$converged=mle$converged
    theta$fvec=Psi_LSMLE(theta$x,y=y,X=x,Z=z,alpha=0,linkobj=linkobj)
    vcov=mle$vcov
    std.error.LSMLE=sqrt(diag(vcov))
  }else{
    #Point Estimation
    q=1-alpha
    theta=Robst.LSMLE.Beta.Reg(y,x,z,start_theta=start_theta,alpha=alpha,linkobj=linkobj,tolerance=control$tolerance,maxit=control$maxit)
    #Predict values
    beta=theta$x[1:k]
    gamma=theta$x[seq.int(length.out = m) + k]
    coefficients = list(mean = beta, precision = gamma)
    eta=x%*%beta
    mu_hat=linkobj$linkfun.mu$inv.link(eta)
    phi_hat=linkobj$linkfun.phi$inv.link(z%*%gamma)
    #Expected Standard Error
    M.LSMLE=LSMLE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha,linkobj=linkobj)
    vcov=M.LSMLE$Cov
    std.error.LSMLE=M.LSMLE$Std.Error
  }
  
  #Register of output values 
  y_star=log(y)-log(1-y)
  str1=str2=NULL
  result$coefficients=coefficients
  result$vcov=vcov
  pseudor2 <- if (var(eta) * var(qlogis(y)) <= 0){NA}
  else{cor(eta, linkobj$linkfun.mu$linkfun(y))^2}
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LSMLE)))
  if(m==1){phi_predict=phi_hat[1]}
  if(m!=1){phi_predict=phi_hat}
  result$fitted.values=list(mu.predict = mu_hat, phi.predict = phi_predict)
  result$start=start_theta #Alterar caso optar por duas tentativas de ponto inicial
  result$weights=(degbeta(y_star,mu_hat,phi_hat/q))^(alpha) #Pesos 
  result$Tuning=alpha
  result$Psi.Value=theta$fvec
  result$residuals=sweighted2_res(mu_hat,phi_hat,y=y,X=x,linkobj = linkobj)
  result$model=list(mean = x, precision = z)
  result$y=y
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
    result=list(beta.wald=result.beta,gamma.wald=result.gamma,general=general,msg="Results based on LSMLE.")
    class(result)="WaldTest_LSMLE"
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
    # result$H0=f
    # result$arg=formalArgs(f)
    # result$arg.values=formals(f)
    result$pValue=as.numeric(1-pchisq(W_alpha,df=r))
    result$general=general
    result$msg="Results based on LSMLE."
    class(result)="WaldTest_LSMLE"
    return(result)
  }
}


#' @export  
plotenvelope.LSMLE=function(robustbetareg.obj,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection"),conf=0.95,n.sim=50,PrgBar=T,control=robustbetareg.control(...), ...)
{
  if(missing(control)){control=robustbetareg.control(robustbetareg.obj)}
  type = match.arg(type)
  ylim.boolean=hasArg(ylim)
  arg=list(...)
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
  if(PrgBar){pb = txtProgressBar(min = 0, max = n.sim, style = 3)}
  #browser()
  while(k<=n.sim & !limit)
  {
    y.sim=pmax(pmin(sapply(seq(1,n,1),function(i) rbeta(1,a[i],b[i])),1-.Machine$double.eps),.Machine$double.eps)
    est.mle=betareg.fit(x,y,z,link=link,link.phi = link.phi)
    start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
    LSMLE.sim=tryCatch(LSMLE.Beta.Reg(y=y.sim,x=x,z=z,alpha=robustbetareg.obj$Tuning,link=link,link.phi=link.phi,start=start),error=function(e){LSMLE.sim$converged<-FALSE; return(LSMLE.sim)})
    if(LSMLE.sim$converged)
    {
      if(type=="sweighted2")
      {
        ResEnvelop=rbind(ResEnvelop,sort(LSMLE.sim$residuals,decreasing = F))  
      }else{
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

#' @export
print.LSMLE=function(obj)
{
  cat("Call: \n")      
  print(obj$call)
  cat("\n")
  cat("Coefficients (mean model with",obj$link,"link):\n")
  print(obj$coefficients$mean)
  cat("\n")
  cat("Coefficients (precision model with",obj$link.phi,"link):\n")
  print(obj$coefficients$precision)
  if(!obj$converged)
  {
    cat("\n")
    cat("The algorithm did not reach the convergence.")
  }
}

#' @export
print.WaldTest_LSMLE=function(obj)
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
  cat("Results based on LSMLE \n")
}


#' @export
summary.LSMLE=function(obj)
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
  if(obj$Tuning==0){cat("Type of estimator: MLE \n")}
  else{cat("Type of estimator: LSMLE \n")}
  cat(paste0("Pseudo R-squared: ",round(obj$pseudo.r.squared,4)),"\n")
  cat(paste0("Tuning value: alpha=",obj$Tuning),"\n")
  if(obj$Optimal.Tuning)
  {
    cat("Tuning generated by the data-driven algorithm")  
  }
  if(!obj$Optimal.Tuning)
  {
    cat("Tuning selected by the user")  
  }
}


#' @rdname plot.LMDPDE
#'   
#' @export
plot.LSMLE=function(object,ask=TRUE,...)
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
          res=residuals(object,type="sweighted2")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted 2 Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==2)
        {
          res=residuals(object,type="sweighted")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==3)
        {
          res=residuals(object,type="pearson")
          plot(res,xlab="Obs. number",ylab="Pearson Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==4)
        {
          res=residuals(object,type="weighted")
          plot(res,xlab="Obs. number",ylab="Weighted Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==5)
        {
          res=residuals(object,type="sweighted.gamma")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted Gamma Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==6)
        {
          residuals(object,type="sweighted2.gamma")
          plot(res,xlab="Obs. number",ylab="Standardized Weighted 2 Gamma Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==7)
        {
          res=residuals(object,type="combined")
          plot(res,xlab="Obs. number",ylab="Combined Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
        }
        if(m==8)
        {
          res=residuals(object,type="combined.projection")
          plot(res,xlab="Obs. number",ylab="Combined Projection Residual",main="Residuals vs indices of obs.",...)
          abline(h=0)
          q <-readline(prompt = cat(crayon::green("Identify points? \n [1] Yes \n [2] No")))
          if(q==1)
          {
            identify(res,pos=T,plot=T)
          }
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
#' rbr=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' residuals(rbr,type="sweighted")
#' 
#' @export
residuals.LSMLE=function(object,type=c("sweighted2","pearson","weighted","sweighted","sweighted.gamma","sweighted2.gamma","combined","combined.projection"))
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
cooks.distance.LSMLE=function(object,...)
{
  h =hatvalues.LSMLE(object)
  k = length(object$coefficients$mean)
  res=residuals(object,type="pearson")
  return(h*(res^2)/(k*(1-h)^2))
}


