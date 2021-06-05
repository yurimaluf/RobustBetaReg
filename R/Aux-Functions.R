#' Robust Beta Regression 
#' 
#' Fit robust beta regression models for rates and proportionsvia LMDPDE and LSMLE using a parametrization with mean (depending through a link function on the covariates) and precision parameter (called phi).
#' 
#' For more details see:...
#' 
#' @param formula symbolic description of the model (of type y ~ x or y ~ x | z).
#' @param data arguments controlling formula.
#' @param alpha the tuning with values (0,1), for robust estimation. When alpha is equal zero is equivalent of MLE. 
#' @param type character specification of the type of estimator. Currently, LMDPDE (default) and LSMLE.
#' @param start a numeric vector with an initial guess of the root of estimation equation.
#' @param alpha.optimal a logical value. If TRUE the tuning parameter should be selected automatic.
#' @param control a list of control arguments specified via \code{\link[RobustBetaReg:robustbetareg.control]{robustbetareg.control}}. 
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
robustbetareg = function(formula,data,alpha,type=c("LMDPDE","LSMLE"),alpha.optimal, start, control=robustbetareg.control(...), ...)
{
  cl = match.call()
  type = match.arg(type)
  ocontrol=control
  #browser()
  if(missing(data)){data=environment(formula)}
  if(missing(alpha.optimal)){alpha.optimal=TRUE}
  if(!missing(alpha)){alpha.optimal=FALSE}
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  formula = Formula::as.Formula(formula)
  mf1=model.frame(formula,data=data)
  y=model.response(mf1)
  x=model.matrix(formula,data=mf1,rhs = 1L)
  z=model.matrix(formula,data=mf1,rhs = 2L)
  if(missing(start))
  {
    est.mle=betareg(formula,data)
    start=c(est.mle$coefficients$mean,est.mle$coefficients$precision)
  }
  start_theta=start
  if(type=="LMDPDE")
  {
    result=LMDPDE.Beta.Reg(y,x,z,alpha=alpha,start_theta,alpha.optimal,control=control)
  }
  if(type=="LSMLE")
  {
    result=LSMLE.Beta.Reg(y,x,z,alpha=alpha,start_theta,alpha.optimal,control=control)
  }
  result$call <- cl
  result$formula=as.formula(formula)
  return(result)
}


#' Control Parameter for Robust Beta Regression 
#' 
#' Various parameters that control fitting of robust beta regression models using robustbetareg
#' 
#' For more details see:...
#' 
#' @param tolerance numeric tolerance for convergence.
#' @param maxit integer specifying the maxit argument of iterations used by the Newton-Raphson algorithm.
#' @param L A parameter of auto selecting algorithm of tuning parameter (default L=0.02).
#' @param M A integer parameter value of auto selecting algorithm of tuning parameter (default M=3).
#' 
#' @return A list with the arguments specified.
#'
#' @export  
robustbetareg.control=function(tolerance=1e-3,maxit=250,L=0.02,M=3,...)
{
  #browser()
  result <- list(tolerance = tolerance, maxit = maxit, L=L, M=M)
  result <- c(result, list(...))
  return(result)
}


#' Inv-Logit Link function
#' 
#' Inverse logit link function for parameter mu  
#' 
#' This function... - RobustBetaReg
#' 
#' @param x Design matrix
#' @param B Beta parameter vector 
#' 
#' @return Return mu parameter for each observation
#'
#' @export 
h1=function(x,B)
{
  #x=x.intercept(x)
  x=as.matrix(x)
  k=as.numeric(x%*%B)
  #h=exp(k-log1pexp(k))
  h=pmax(pmin(exp(k-Rmpfr::log1pexp(k)),1-.Machine$double.eps),.Machine$double.eps)
  return(h)
}

#' Inv-log Link function
#' 
#' Inverse log link function for parameter phi  
#' 
#' This function... - RobustBetaReg
#' 
#' @param z Design matrix
#' @param G Gamma parameter vector 
#' 
#' @return Return mu parameter for each observation
#' 
#' @export
h2=function(z,G)
{
  z=as.matrix(z)
  k=as.numeric(z%*%G)
  return(exp(k))
}

#' The EGB of the second type
#' 
#' Exponential Generalized Beta (EGB) of the second type density function.
#' 
#' For more details see 
#' 
#' @param y_star Logit transformation of original data
#' @param mu mu parameter
#' @param phi phi parameter
#' @param log LOGICAL parameter. If TRUE return the log of density function
#' 
#' @return Return value of density function.
#' 
#'@export 
degbeta=function(y_star,mu,phi,log)
{
  if(missing(log)){log=F}
  a0=mu*phi
  b0=(1-mu)*phi
  if(log==F){
    k=exp(-(lbeta(a0,b0)+b0*y_star+phi*log(1+exp(-y_star))))  
  }
  if(log==T){
    k=-(lbeta(a0,b0)+b0*y_star+phi*log(1+exp(-y_star)))
  }
  return(k)
}

#' The EGB of the second type
#' 
#' Random generation for Exponential Generalized Beta (EGB) of the second type.
#' 
#' For more details see 
#' 
#' @param n number of observations
#' @param mu mu parameter
#' @param phi phi parameter
#' 
#' @return Return value of density function.
#' 
#'@export 
regbeta=function(n,mu,phi)
{
  h=rbeta(mu*phi,(1-mu)*phi)
  return(log(h)-log(1-h))
}


#' Newton-Raphson
#' 
#' Newton-Raphson algorithm. The function solves a system of nonlinear equations. 
#' 
#' @param p0 Initial point 
#' @param M Max number of iteration
#' @param tol Tolerance 
#' @param details A logical indicating if the each iteration should be return.
#' @param FUN A function of x returning a vector of function values with the same length as the vector x.
#' 
#' @return Return solution of a system of nonlinear equations.
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
    w=ginv(A)%*%b
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

#' SQV
#' 
#' Standardized Quadratic Variations  
#' 
#' This function... - RobustBetaReg
#' 
#' @param zq Standardized parameter vector
#' @param n Sample size 
#' @param p Parameter vector length
#' 
#' @return Return SQV vector.
#' 
SQV=function(zq,n,p)
{
  return(sqrt(rowSums(diff(zq)^2))/(sqrt(n)*p))
}

#' Residual Beta
#' 
#' Residual Type Standardized Weighted Residual 2.
#' 
#' For more details see: On beta regression residuals: Patrícia L. Espinheira a , Silvia L.P. Ferrari and Francisco Cribari-Neto; Journal of Applied Statistics Vol. 35, No. 4, April 2008 
#' 
#' @param mu_hat ...
#' @param phi_hat ...
#' @param y ...
#' @param X ...
#' @param Z ...
#' 
#' @return Return value ...
#'
Residual_Beta=function(mu_hat,phi_hat,y,X,Z)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  V_star=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  
  W=diag(x=phi_hat*V_star*(mu_hat-mu_hat^2)^2)
  PHI=diag(phi_hat)
  H=sqrt(W%*%PHI)%*%X%*%solve(t(X)%*%PHI%*%W%*%X)%*%t(X)%*%sqrt(PHI%*%W)
  
  nu=V_star*(1-diag(H))
  diff=(y_star-mu_star)
  ri=diff/sqrt(nu)#standardized weighted residuals
  
  return(ri)
}

#' Residual Gamma
#' 
#' Residual Type Gamma
#' 
#' For more details see: Espinheira, P. L., Santos, E. G., and Cribari Neto, F. (2017). On nonlinear beta regression residuals. Biometrical Journal, 59:445–461. 
#' 
#' @param mu_hat ...
#' @param phi_hat ...
#' @param y ...
#' @param X ...
#' @param Z ...
#' 
#' @return Return value ...
#'
Residual_Gamma=function(mu_hat,phi_hat,y,X,Z)
{
  y_star=log(y/(1-y))
  mu_star=digamma(mu_hat*phi_hat)-digamma((1-mu_hat)*phi_hat)
  xi=(trigamma(mu_hat*phi_hat)*mu_hat^2+trigamma((1-mu_hat)*phi_hat)*(1-mu_hat)^2-trigamma(phi_hat))
  Z=as.matrix(Z)
  D=diag(xi*(phi_hat)^2)
  G=sqrt(D)%*%Z%*%solve(t(Z)%*%D%*%Z)%*%t(Z)%*%sqrt(D)
  
  a=mu_hat*(y_star-mu_star)+log(1-y)-(digamma((1-mu_hat)*phi_hat)-digamma(phi_hat))
  nu=xi*(1-diag(G))
  ri=a/sqrt(nu) #standardized weighted residuals
  
  return(ri)
}


#' Robust Initial Points
#' 
#' Robust Initial Points
#' 
#' @param y ...
#' @param X ...
#' @param Z ...
#' 
#' @return Return value ...
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

#' Teste ***
#' 
#' star observation
#' 
#' For more details see 
#' 
#' @param p.value 
#' 
#' @return Return number stars ***.
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


#' Teste2 ***
#' 
#' star observation
#' 
#' For more details see 
#' 
#' @param ObjRbst 
#' 
#' @return Return number stars ***.
#'
Z.q=function(ObjRbst)
{
  return(c(ObjRbst$coefficients$mean,ObjRbst$coefficients$precision)/c(ObjRbst$std.error$se.mean,ObjRbst$std.error$se.precision))
}

#' Teste3 ***
#' 
#' star observation
#' 
#' For more details see 
#' 
#' @param ObjRbst 
#' 
#' @return Return number stars ***.
#'
Par.q=function(ObjRbst)
{
  return(c(ObjRbst$coefficients$mean,ObjRbst$coefficients$precision))
}

#' Teste4 ***
#' 
#' star observation
#' 
#' For more details see 
#' 
#' @param ObjRbst 
#' 
#' @return Return number stars ***.
#'
SE.q=function(ObjRbst)
{
  return(c(ObjRbst$std.error$se.mean,ObjRbst$std.error$se.precision))
}

#' Envelope Bla
#' 
#' Funcao TESTE.  
#' 
#' This function.
#' 
#' @param x1 Numero
#' @param x2 Numero 
#' 
#' @return Return a soma 
#' 
#' @export
soma=function(x1,x2)
{
  h=x1+x2
  return(h)
}

#' Simulated Envelope of Residuals
#'  
#'  
#'    
#' 
#' Plot a simulated envelope of beta residuals, from LMDPDE and LSMLE.
#' 
#' 
#' 
#' 
#' @param robustbetareg.obj Object of robust beta regression (see \code{\link[RobustBetaReg:robustbetareg]{robustbetareg}}). 
#' @param n.sim the number of simulation sample 
#' @param ylim the y limits of the plot
#' @param index index of points to draw 
#' @param control a list of control arguments specified via \code{\link[RobustBetaReg:robustbetareg.control]{robustbetareg.control}}. 
#' 
#' @return Return a simulated envelope graphic.
#' 
#' @examples 
#' rbr.obj=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' plotenvelope(rbr.obj,n.sim=100)
#' 
#' @export
plotenvelope=function(robustbetareg.obj,n.sim,ylim,index,control)
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
#' @param object Object of robust beta regression estimator
#' @param g A funtion representing the null hypothesis to be tested  
#' @param ... Further arguments to be passed 
#' 
#' @references \url{https://link.springer.com/article/10.1007%2Fs00184-018-0653-4}
#' 
#' @examples 
#' rbr.obj=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' g=function(theta,B){theta[2:3]-B}#Hiphothesis to be tested
#' WaldTypeTest(rbr.obj,g,B=c(0,0))#Testing income=persons=0
#' 
#' @export
WaldTypeTest=function(object,g,...)
{
  UseMethod("WaldTypeTest")
}
