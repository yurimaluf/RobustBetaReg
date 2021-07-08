#' Robust Beta Regression 
#' 
#' Fit robust beta regression models for rates and proportionsvia LMDPDE and LSMLE using a parametrization with mean (depending through a link function on the covariates) and precision parameter (called phi).
#' 
#' For more details see:...
#' 
#' @param formula symbolic description of the model (of type y ~ x or y ~ x | z).
#' @param data arguments controlling formula.
#' @param alpha the tuning with values (0,1), for robust estimation. When alpha is equal to zero is equivalent of MLE. 
#' @param type character specification of the type of estimator. Currently, LSMLE (default) and LMDPDE.
#' @param control a list of control arguments specified via \code{\link[=robustbetareg.control]{robustbetareg.control}}. 
#' 
#' @return Return a list of components:
#'  \itemize{
#'   \item coefficients - the vector with elements "mean" and "precision" containing the coefficients from the respective models
#'   \item vcov - the covariance matrix of all parameters in the model
#'   \item mu_hat - the vector of fitted means
#'   \item phi_hat - the vector of predicted precision
#'   \item weights - the weights used
#'   \item Tuning - the employed tuning parameter
#'   \item Res.Beta - a vector of standardized weighted residual 2
#'   \item start - the starting values for the parameters passed to the Newton-Raphson algorithm
#'   \item std.error - the standard error vector of coefficients
#'   \item Tab - the summary table of overall result.
#' }
#'
#' @export  
robustbetareg = function(formula,data,alpha,type=c("LSMLE","LMDPDE"),link = c("logit", "probit", "cloglog", "cauchit", "loglog"),
link.phi = NULL,control=robustbetareg.control(...), ...)
{
  cl = match.call()
  type = match.arg(type)
  ocontrol=control
  if(missing(data)){data=environment(formula)}
  if(!missing(alpha)){control$alpha.optimal=FALSE}
  if(missing(alpha)){alpha=NULL}
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  formula = Formula::as.Formula(formula)
  oformula <- as.formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~1)
    simple_formula <- TRUE
  }else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf1=model.frame(formula,data=data)
  y=model.response(mf1)
  x=model.matrix(formula,data=mf1,rhs = 1L)
  z=model.matrix(formula,data=mf1,rhs = 2L)
  #Tratamento erro modelo
  if(length(y) < 1){stop("empty model")} 
  if(!(min(y) > 0 & max(y) < 1)){stop("invalid dependent variable, all observations must be in (0, 1)")} 
  if(!is.null(control$start) & (ncol(x)+ncol(z))!=length(control$start) ){stop("Invalid initial starting point")}
  if(!is.null(alpha)){if(alpha < 0 || alpha > 1){stop("invalid tuning constant, the value must be in [0, 1)")}}
   
  link = match.arg(link)
  if(is.null(link.phi))
    {
    link.phi <- if(simple_formula){"identity"}
    else "log" 
  } 
  if(is.null(control$start))
  {
    est.mle=betareg(oformula,data,link=link,link.phi = link.phi)
    control$start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
  }
  if(type=="LMDPDE")
  {
    result=LMDPDE.Beta.Reg(y,x,z,alpha=alpha,link=link,link.phi=link.phi,control=control)
  }
  if(type=="LSMLE")
  {
    result=LSMLE.Beta.Reg(y,x,z,alpha=alpha,link=link,link.phi=link.phi,control=control)
  }
  result$call <- cl
  result$formula=as.formula(formula)
  return(result)
}


#' Control Parameter for Robust Beta Regression 
#' 
#' Various parameters that control fitting of robust beta regression models using \code{\link[=robustbetareg]{robustbetareg.}}
#' 
#' For more details see:...
#' 
#' @param object fitted model object of class "LMDPDE" or "LSMLE".
#' @param start a numeric vector with an initial guess of the root of estimation equation.
#' @param alpha.optimal a logical value. If TRUE the tuning parameter should be selected automatic.
#' @param tolerance numeric tolerance for convergence.
#' @param maxit integer specifying the maxit argument of iterations used by the Newton-Raphson algorithm.
#' @param L a parameter of auto selecting algorithm of tuning parameter (default L=0.02).
#' @param M a integer parameter value of auto selecting algorithm of tuning parameter (default M=3).
#' @param ... currently not used.
#' 
#' @return A list with the arguments specified.
#'
#' @export  
robustbetareg.control=function(object,start=NULL,alpha.optimal=TRUE,tolerance=1e-3,maxit=250,L=0.02,M=3,...)
{
  UseMethod("robustbetareg.control")
}

#' @export  
robustbetareg.control.default=function(start=NULL,alpha.optimal=TRUE,tolerance=1e-3,maxit=250,L=0.02,M=3,...)
{
  result <- list(start=start,alpha.optimal=alpha.optimal,tolerance = tolerance, maxit = maxit, L=L, M=M)
  result <- c(result, list(...))
  return(result)
}


#' The EGB of the second type
#' 
#' Density and random generation for the exponential generalized beta (EGB) of the second type, 
#' 
#' For more details see 
#' 
#' @param y_star logit transformation of original data.
#' @param mu mu parameter.
#' @param phi phi parameter.
#' @param log a logical value. If TRUE return the log of density function.
#' 
#' @return Return value of density function or a random sample.
#' 
#'@export 
degbeta=function(y_star,mu,phi,log=FALSE)
{
  #options(warn = 2) #Converte warnings em erros 
  a0=mu*phi
  b0=(1-mu)*phi
  if(log==F){
    #k=exp(-(lbeta(a0,b0)+b0*y_star+phi*log(1+exp(-y_star))))
    k=tryCatch(exp(-(lbeta(a0,b0)+b0*y_star+phi*log(1+exp(-y_star)))),error=function(e){stop("Error")})
  }
  if(log==T){
    #k=-(lbeta(a0,b0)+b0*y_star+phi*log(1+exp(-y_star)))
    k=tryCatch(-(lbeta(a0,b0)+b0*y_star+phi*log(1+exp(-y_star))),error=function(e){stop("Error")})
  }
  return(k)
}

#' @rdname degbeta
#' 
#' @param n number of observations.
#' 
#'@export 
regbeta=function(n,mu,phi)
{
  h=rbeta(mu*phi,(1-mu)*phi)
  return(log(h)-log(1-h))
}

#' 
Newton.Raphson=function(x,M,tol,details,FUN=f, ...)
{
  if(missing(M)){M=100}
  if(missing(details)){details=F}
  if(missing(tol)){tol=1e-3}
  result=list()
  convergence=F
  f=FUN
  x.h=p=x
  y.h=f(x, ...)
  #browser()
  msg1=msg2=NULL
  for(i in 1:M)
  {
    A=numDeriv::jacobian(f,p, ...)
    #A=jacobian(f,p)#teste
    #A#teste
    if(any(is.na(A)))
    {
      msg1="NA values in Jacobian matrix."
    }
    b=-f(p, ...)
    #b=-f(p)#teste
    #b#teste
    if(any(is.na(b)))
    {
      msg2="NA values in fn."
    }
    if(!is.null(msg1) || !is.null(msg2))
    {
      result$msg=paste0(msg1," ",msg2)
      break
    }
    #Solve(A,b)
    #w=solve(A)%*%b
    w=MASS::ginv(A)%*%b
    p=p+w
    x.h=rbind(as.vector(p),x.h)
    y.h=rbind(f(p,...),y.h)  
    if(sqrt(sum(w^2))<tol & max(abs(y.h[1,]))<tol)
    {
      convergence=T
      break
    }
  }
  result$sol=as.vector(p)
  if(details)
  {
    result$x.h=x.h
    result$f.vec=y.h
    result$iter=i
  }
  result$converged=convergence
  return(result)
}

#' 
SQV=function(zq,n,p)
{
  return(sqrt(rowSums(diff(zq)^2))/(sqrt(n)*p))
}

#'
pearson_res=function(mu_hat,phi_hat,y)
{
  var.y=mu_hat*(1-mu_hat)/(1+phi_hat)
  ri=(y-mu_hat)/sqrt(var.y)
  return(ri)
}

#'
sweighted_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  nu=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  ri=(y_star-mu_star)/sqrt(nu)#standardized weighted residuals
  
  return(ri)
}

#'
sweighted2_res=function(mu_hat,phi_hat,y,X,linkobj)
{
  n=length(mu_hat)
  m=length(phi_hat)
  if(m==1){phi_hat=rep(phi_hat,n)}
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  y_star=log(y)-log(1-y)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  
  #W=diag(x=phi_hat*V_star*(inverse(d.link.mu))^2)#Depende da funcao de ligacao. Aqui funcao ligacao logit
  #PHI=diag(phi_hat)#
  W.PHI=diag(x=phi_hat*V_star*(inverse(d.link.mu))^2)
  
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  
  nu=V_star*(1-diag(H))
  diff=(y_star-mu_star)
  ri=diff/sqrt(nu)#standardized weighted residuals
  
  return(ri)
}


#'
weighted_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  
  nu=V_star*phi_hat
  ri=(y_star-mu_star)/sqrt(nu)#standardized weighted residuals
  
  return(ri)
}

#'
sweighted.gamma_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  xi=trigamma((1-mu_hat)*phi_hat)-trigamma(phi_hat)  
  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  ri=a/sqrt(xi) #standardized weighted residuals gamma
  return(ri)
}

#'
sweighted2.gamma_res=function(mu_hat,phi_hat,y,Z,linkobj)
{
  y_star=log(y/(1-y))
  n=length(mu_hat)
  m=length(phi_hat)
  if(m==1){phi_hat=rep(phi_hat,n)}
  d.link.phi=linkobj$linkfun.phi$d.linkfun(phi_hat)
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  xi=(trigamma(mu_hat*phi_hat)*mu_hat^2+trigamma((1-mu_hat)*phi_hat)*(1-mu_hat)^2-trigamma(phi_hat))
  Z=as.matrix(Z)
  D=diag(xi*(inverse(d.link.phi))^2)##Para funcao ligacao de phi
  G=sqrt(D)%*%Z%*%solve(t(Z)%*%D%*%Z)%*%t(Z)%*%sqrt(D)
  
  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  nu=xi*(1-diag(G))
  ri=a/sqrt(nu) #standardized weighted residuals
  
  return(ri)
}


#'
combined_res=function(mu_hat,phi_hat,y)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  zeta=trigamma(mu_hat*phi_hat)*(1+mu_hat)^2+trigamma((1-mu_hat)*phi_hat)*mu_hat^2-trigamma(phi_hat)
  ri=(y_star-mu_star+a)/sqrt(zeta) #combined residuals
  return(ri)
}


#'
combined.projection_res=function(mu_hat,phi_hat,y,X,Z,linkobj)
{
  n=length(mu_hat)
  m=length(phi_hat)
  if(m==1){phi_hat=rep(phi_hat,n)}
  #derivada da inversa da funcao ligacao
  d.link.mu=linkobj$linkfun.mu$d.linkfun(mu_hat)
  d.link.phi=linkobj$linkfun.phi$d.linkfun(phi_hat)
  
  d.inv.g.beta=mu_hat-mu_hat^2
  d.inv.g.gamma=phi_hat

  y_star=log(y/(1-y))
  res.beta=sweighted_res(mu_hat,phi_hat,y)
  res.gamma=sweighted.gamma_res(mu_hat,phi_hat,y)
  v=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  b=mu_hat*v-trigamma((1-mu_hat)*phi_hat)
  xi=trigamma(mu_hat*phi_hat)*mu_hat^2+trigamma((1-mu_hat)*phi_hat)*(1-mu_hat)^2-trigamma(phi_hat)
  w=phi_hat*v*(inverse(d.link.mu))^2
  TT=diag(inverse(d.link.mu))
  HH=diag(inverse(d.link.phi))
  B=diag(b)
  
  W.PHI=diag(x=phi_hat*w)
  W_PHI=diag(x=phi_hat/w)
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  
  D=diag(xi*(inverse(d.link.phi))^2)
  G=sqrt(D)%*%Z%*%solve(t(Z)%*%D%*%Z)%*%t(Z)%*%sqrt(D)
  M=sqrt(W_PHI)%*%TT%*%B%*%HH%*%sqrt(D)
  
  h=1-diag(H)
  g=1-diag(G)
  m=diag(M)
  nu=pmax(h+g+2*m,.Machine$double.eps)
  ri=(res.beta+res.gamma)/sqrt(nu)
  return(ri)
}

#'
Initial.points=function(y,X,Z)
{
  X=as.matrix(X)
  Z=as.matrix(Z)
  x2=X[,-1]
  ystar=log(y/(1-y))
  lmbeta=suppressWarnings(robustbase::lmrob(ystar~x2))
  betaini=as.numeric(lmbeta$coefficients)
  muini=exp(lmbeta$fitted.values-Rmpfr::log1pexp(lmbeta$fitted.values))
  #muini=exp(lmbeta$fitted.values)/(1+exp(lmbeta$fitted.values))
  sigma2ini=sigma(lmbeta)^2*muini^2*(1-muini)^2
  phiini=as.numeric(muini*(1-muini)/sigma2ini)
  #phiini=inverse(sigma(lmbeta)^2*muini*(1-muini))
  #muini<<-muini
  #sigma2ini<<-sigma2ini
  if(dim(Z)[2]>1)
  {
    phistar=log(phiini)
    gammaini=solve(t(Z)%*%Z)%*%t(Z)%*%phistar
  }else{gammaini=log(mean(phiini,na.rm=T))
  }
  gammaini=as.numeric(c(gammaini))
  Est.param.Rbst=c(betaini,gammaini)
  return(Est.param.Rbst)
}

#'
star.obs=function(p.valor)
{
  obs=NULL
  if(p.valor<0.001){
    obs="***"
  } else if(p.valor<0.01){
    obs="**"
  } else if(p.valor<0.05)
  {
    obs="*"
  } else if(p.valor<0.1)
  {
    obs="."
  } else {
    obs=" "
  }
  return(obs)
}

#'
Z.q=function(ObjRbst)
{
  return(c(ObjRbst$coefficients$mean,ObjRbst$coefficients$precision)/c(ObjRbst$std.error$se.mean,ObjRbst$std.error$se.precision))
}

#'
Par.q=function(ObjRbst)
{
  return(c(ObjRbst$coefficients$mean,ObjRbst$coefficients$precision))
}

#'
SE.q=function(ObjRbst)
{
  return(c(ObjRbst$std.error$se.mean,ObjRbst$std.error$se.precision))
}

#' Simulated Envelope of Residuals
#'  
#'  
#'    
#' 
#' Plot a simulated envelope of beta residuals, from LMDPDE and LSMLE objects.
#' 
#' 
#' 
#' 
#' @param robustbetareg.obj Fitted model object of class "LMDPDE" or "LSMLE" (see \code{\link[RobustBetaReg:robustbetareg]{robustbetareg}}). 
#' @param n.sim the number of simulation sample. Deafault n.sim=50. 
#' @param conf the confidence level of the envelopes required. The default is to find 95% confidence envelopes.
#' @param control a list of control arguments specified via \code{\link[RobustBetaReg:robustbetareg.control]{robustbetareg.control}}.
#' @param ... other parameters to be passed through to plotting functions.
#' 
#' @return Return a simulated envelope graphic.
#' 
#' @examples 
#' rbr=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' plotenvelope(rbr,n.sim=100)
#' 
#' @export
plotenvelope=function(robustbetareg.obj,n.sim,conf,control=robustbetareg.control(...),...)
{
 UseMethod("plotenvelope")
}

#' Robust Wald-type Tests
#'  
#'  
#'    
#' 
#' Wald-type tests for both simple and composite hypothesis for independent but non-homogeneous observations, based on LMDPDE and LSMLE.
#' 
#' 
#' 
#' @param object fitted model object of class "LMDPDE" or "LSMLE"
#' @param FUN the function representing the null hypothesis to be tested  
#' @param ... Further arguments to be passed 
#' 
#' @references \href{https://www.tandfonline.com/doi/abs/10.1080/02664760701834931}{Basu, A., Ghosh, A., Martin, N. et al. Robust Wald-type tests for non-homogeneous observations based on the minimum density power divergence estimator. Metrika 81, 493–522 (2018)}
#' 
#' @examples 
#' rbr=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' g=function(theta,B){theta[2:3]-B}#Hiphothesis to be tested
#' WaldTypeTest(rbr,FUN=g,B=c(0,0))#Testing income=persons=0
#' 
#' @export
WaldTypeTest=function(object,FUN,...)
{
  UseMethod("WaldTypeTest")
}


#'
hatvalues=function(object)
{
  UseMethod("hatvalues")
}


#'
inverse=function(x)
{
  return(x^(-1))
}

#' @noRd
ddnorm=function(x)
{
  #-(EQL::hermite(x,1))*exp(-x^(2)/2)/sqrt(2*pi)
  return(-x*exp(-x^(2)/2)/sqrt(2*pi))
}

#'
make.link=function(link.mu,link.phi)
{
  #c("logit", "probit", "cloglog", "cauchit", "loglog")
  if(missing(link.mu)){link.mu="logit"}
  if(missing(link.phi)){link.phi="log"}
  #Mean Links
  if(link.mu=="logit")
  {
    linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(log(mu)-log(1-mu))
    }
    d.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(inverse(mu-mu^2))
    }
    d2.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return((2*mu-1)/(mu*(1-mu))^2)
    }
    inv.link=function(eta)
    {
      return(as.numeric(pmax(pmin(exp(eta-Rmpfr::log1pexp(eta)),1-.Machine$double.eps),.Machine$double.eps)))
    }
  }
  if(link.mu=="probit")
  {
    linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(qnorm(mu))
    }
    d.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(inverse(dnorm(qnorm(mu))))
    }
    d2.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(-ddnorm(qnorm(mu))/(dnorm(qnorm(mu)))^3)
    }
    inv.link=function(eta)
    {
      return(as.numeric(pmax(pmin(pnorm(eta),1-.Machine$double.eps),.Machine$double.eps)))
    }
  }
  if(link.mu=="cloglog")
  {
    linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(log(-log(1-mu)))
    }
    d.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(inverse(-(1-mu)*log(1-mu)))
    }
    d2.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(-(log(1-mu)+1)/((1-mu)*log(1-mu))^2)
    }
    inv.link=function(eta)
    {
      return(as.numeric(pmax(pmin(1-exp(-exp(-eta)),1-.Machine$double.eps),.Machine$double.eps)))
    }
  }
  if(link.mu=="cauchit")
  {
    linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(tan(pi*(mu-0.5)))
    }
    d.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(pi*pracma::sec(pi*(mu-0.5))^2)
    }
    d2.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(2*pi*tan(pi*(mu-0.5))*sec(pi*(mu-0.5))^2)
    }
    inv.link=function(eta)
    {
      return(as.numeric(pmax(pmin(0.5+atan(eta)/pi,1-.Machine$double.eps),.Machine$double.eps)))
    }
  }
  if(link.mu=="loglog")
  {
    linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(-log(-log(mu)))
    }
    d.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(-inverse(mu*log(mu)))
    }
    d2.linkfun=function(mu)
    {
      mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return((log(mu)+1)/(mu*log(mu))^2)
    }
    inv.link=function(eta)
    {
      return(as.numeric(pmax(pmin(exp(-exp(-eta)),1-.Machine$double.eps),.Machine$double.eps)))
    }
  }
  Linkfun.Mu=list(linkfun=linkfun,d.linkfun=d.linkfun,d2.linkfun=d2.linkfun,inv.link=inv.link)
  ### Precision Links
  if(link.phi=="log")
  {
    linkfun.phi=function(phi)
    {
      phi=pmax(phi,.Machine$double.eps)
      return(log(phi))
    }
    d.linkfun.phi=function(phi)
    {
      phi=pmax(phi,.Machine$double.eps)
      return(inverse(phi))
    }
    d2.linkfun.phi=function(phi)
    {
      phi=pmax(phi,.Machine$double.eps)
      return(-inverse(phi^2))
    }
    inv.link.phi=function(eta)
    {
      return(as.numeric(exp(eta)))
    }
  }
  if(link.phi=="identity")
  {
    linkfun.phi=function(phi)
    {
      phi=pmax(phi,.Machine$double.eps)
      return(phi)
    }
    d.linkfun.phi=function(phi)
    {
      return(rep(1,length(phi)))
    }
    d2.linkfun.phi=function(phi)
    {
      return(rep(0,length(phi)))
    }
    inv.link.phi=function(eta)
    {
      return(as.numeric(eta))
    }
  }
  if(link.phi=="sqrt")
  {
    linkfun.phi=function(phi)
    {
      phi=pmax(phi,.Machine$double.eps)
      return(sqrt(phi))
    }
    d.linkfun.phi=function(phi)
    {
      return(inverse(2*sqrt(phi)))
    }
    d2.linkfun.phi=function(phi)
    {
      return(-inverse(4*phi^(3/2)))
    }
    inv.link.phi=function(eta)
    {
      return(as.numeric(eta^2))
    }
  }
  Linkfun.Phi=list(linkfun=linkfun.phi,d.linkfun=d.linkfun.phi,d2.linkfun=d2.linkfun.phi,inv.link=inv.link.phi)
  ###
  linkobj=structure(list(linkfun.mu=Linkfun.Mu,linkfun.phi=Linkfun.Phi),name.link.mu=link.mu,name.link.phi=link.phi,class="link-rbr")
  return(linkobj)
}

#' @export
soma=function(x1,x2)
{
  times4(x1+x2)
}

#' @export
soma2=function(x1,x2)
{
  times3(x1+x2)
}
