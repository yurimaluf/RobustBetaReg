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
h2=function(z,G)
{
  z=as.matrix(z)
  k=as.numeric(z%*%G)
  return(exp(k))
}

#' EGB of the second type
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
Gen_Beta=function(y_star,mu,phi,log)
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
  #y.h=f(p0)#teste
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
    w=solve(A)%*%b
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
  #X=x.intercept(X)
  #Z=x.intercept(Z)
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

