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
    est.mle=suppressWarnings(betareg(oformula,data,link=link,link.phi = link.phi))
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

#' Simulated Envelope of Residuals
#'  
#'  
#'    
#' 
#' Plot a simulated envelope of beta residuals, from LSMLE and LMDPDE objects.
#' 
#' 
#' 
#' 
#' @param robustbetareg.obj Fitted model object of class "LSMLE" or "LMDPDE" (see \code{\link[RobustBetaReg:robustbetareg]{robustbetareg}}). 
#' @param n.sim the number of simulation sample. Deafault n.sim=50. 
#' @param conf the confidence level of the envelopes required. The default is to find 95% confidence envelopes.
#' @param control a list of control arguments specified via \code{\link[RobustBetaReg:robustbetareg.control]{robustbetareg.control}}.
#' @param ... other parameters to be passed through to plotting functions.
#' 
#' @return Return a simulated envelope graphic.
#' 
#' @examples 
#' fit=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' plotenvelope(fit,n.sim=100)
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
#' Wald-type tests for both simple and composite hypothesis for independent but non-homogeneous observations, based on LSMLE and LMDPDE.
#' 
#' 
#' 
#' @param object fitted model object of class "LSMLE" or "LMDPDE" 
#' @param FUN the function representing the null hypothesis to be tested  
#' @param ... Further arguments to be passed 
#' 
#' @references \href{https://www.tandfonline.com/doi/abs/10.1080/02664760701834931}{Basu, A., Ghosh, A., Martin, N. et al. Robust Wald-type tests for non-homogeneous observations based on the minimum density power divergence estimator. Metrika 81, 493–522 (2018)}
#' 
#' @examples 
#' fit=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' h0=function(theta,B){theta[2:3]-B}#Hiphothesis to be tested
#' WaldTypeTest(fit,h0,B=c(0,0))#Testing income=persons=0
#' 
#' @export
WaldTypeTest=function(object,FUN,...)
{
  UseMethod("WaldTypeTest")
}


#' Robust Saddlepoint Test
#'  
#'  
#'    
#' 
#' Saddlepoint tests for both simple and composite hypothesis for independent but non-homogeneous observations, based on LSMLE and LMDPDE. 
#' 
#' 
#' 
#' @param object fitted model object of class "LSMLE" or "LMDPDE"  
#' @param FUN the function representing the null hypothesis to be tested  
#' @param ... further arguments to be passed 
#' @param thrd number (integer) of threads to speed up the process. If missing the value is autodetected by the available number of multi-core processor
#' 
#' @references \href{https://www.sciencedirect.com/science/article/pii/S0047259X09001183}{Lo, S. N., Ronchetti, E. Robust and accurate inference for generalized linear models. Journal of Multivariate Analysis, 100, 2126–2136 (2009)}
#' 
#' @examples 
#' fit=robustbetareg(I(food/income)~income+persons|1,data=FoodExpenditure,alpha=0.08)
#' h0=function(theta,B){c(theta[1],B)}#Hiphothesis to be tested
#' SaddlepointTest(fit,h0,B=c(0,0))#H0: income=persons=0
#' 
#' @export
SaddlepointTest=function(object,FUN,...,thrd)
{
  UseMethod("SaddlepointTest")
}

