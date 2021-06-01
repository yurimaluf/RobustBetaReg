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
#' @param L A parameter of auto selecting algorithm of tuning parameter (default L=0.02).
#' @param M A integer parameter value of auto selecting algorithm of tuning parameter (default M=3).
#' @param tolerance The function value tolerance.
#' @param maxit The maximum number of iterations used by the algorithm.
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
LMDPDE.Beta.Reg=function(y,x,z,alpha,start_theta,alpha.optimal,L,M,tolerance,maxit)
{
  result=theta=list()
  if(missing(alpha.optimal)){alpha.optimal=TRUE}
  if(!missing(alpha)){alpha.optimal=FALSE}
  if(missing(z)){z=rep(1,length(y))}
  if(missing(tolerance)){tolerance=1e-3}
  if(missing(maxit)){maxit=250}
  if(missing(L)){L=0.02}
  if(missing(M)){M=3}
  x=as.matrix(x)
  z=as.matrix(z)
  k=ncol(x)
  m=ncol(z)
  if(alpha.optimal)
  {
    return(Opt.Tuning.LMDPDE(y,x,z,tolerance,L=L,M=M))
  }
  if(missing(alpha))
  {
    if(dim(z)[2]==1)
    {
      mle=suppressWarnings(betareg(y~x[,-1]|1))  
    }
    else{mle=suppressWarnings(betareg(y~x[,-1]|z[,-1]))}
    theta$x=as.numeric(c(mle$coefficients$mean,mle$coefficients$precision))
    alpha=0
    theta$inter=NULL
    theta$converged=mle$converged
    theta$fvec=rep(0,length(theta$x))
    vcov=mle$vcov
    std.error.LMDPDE=sqrt(diag(vcov))
    #return(mle)
  }
  if(missing(start_theta))
  {
    #start_theta=Initial.points(y,x,z)
    est.log.lik=tryCatch(suppressWarnings(betareg(y~x[,-1]|1)),error=function(e) NULL)
    start_theta=as.numeric(c(est.log.lik$coefficients$mean,est.log.lik$coefficients$precision))
  }
  initial.point=start_theta
  #Point Estimation
  theta=Robst.LMDPDE.Beta.Reg(y,x,z,start_theta=start_theta,alpha=alpha,tolerance=tolerance,maxit=maxit)
  if(!theta$converged)
  {
    initial.point.rbst=Initial.points(y,x,z)
    theta=Robst.LMDPDE.Beta.Reg(y,x,z,start_theta=initial.point.rbst,alpha=alpha,tolerance=tolerance,maxit=maxit)
    if(!theta$converged){initial.point=initial.point.rbst}
  }
  #Predict values
  beta=theta$x[1:ncol(x)]
  gamma=theta$x[seq.int(length.out = m) + k]
  coefficients = list(mean = beta, precision = gamma)
  mu_hat=h1(x,beta)
  phi_hat=h2(z,gamma)
  y_star=log(y)-log(1-y)
  
  #Expected Standard Error
  MM=LMDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha)
  vcov=MM$Cov
  std.error.LMDPDE=MM$Std.Error
  
  #Register of output values 
  str1=str2=NULL
  result$coefficients=coefficients
  result$vcov=vcov
  if(!is.null(theta$iter)){result$iter=theta$iter}
  result$converged=(theta$converged & all(!is.na(std.error.LMDPDE)))
  result$mu_predict=mu_hat
  if(m==1){result$phi_hat=phi_hat[1]}
  if(m!=1){result$phi_hat=phi_hat}
  result$start=initial.point
  result$weights=Gen_Beta(y_star,mu_hat,phi_hat)^(alpha)#Weights
  result$Tuning=alpha
  result$Psi.Value=theta$fvec
  result$Res.Beta=Residual_Beta(mu_hat,phi_hat,y=y,X=x,Z=z)
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
  return(result)
}

#' Robust Wald type Test - LMDPDE
#' 
#' Wald type test statistic based on LMDPDE
#' 
#' For more details see: 
#' 
#' @param M teste
#' @param m teste
#' @param y The numeric response vector (with values in (0,1)).
#' @param x The numeric regressor matrix for mean model.
#' @param z The numeric regressor matrix for precision model, defaulting to an intercept only.
#' @param alpha The tuning with values (0,1), for robust estimation. When alpha is equal zero is equivalent of MLE. 
#' @param start_theta A numeric vector with an initial guess of the root of estimation equation.
#' @param tolerance The function value tolerance.
#' 
#' @return Return a list of components:
#'  \itemize{
#'   \item coefficients - a vector with elements "mean" and "precision" containing the coefficients from the respective models.
#'   \item vcov - covariance matrix of all parameters in the model.
#' }
#'
LMDPDE_Wald_Type=function(M,m,y,x,z,start_theta,alpha,tolerance)
{
  result=list()
  if(rankMatrix(M)[1]!=dim(M)[1]){return(c("The Rank Matrix is not supported"))}
  if(missing(tolerance)){tolerance=1e-4}
  if(missing(start_theta))
  {
    start_theta=Initial.points(y,x,z)
  }
  if(alpha==0)
  {
    if(dim(z)[2]==1)
    {
      mle=suppressWarnings(betareg(y~x[,-1]|1))  
    }
    else{mle=suppressWarnings(betareg(y~x[,-1]|z[,-1]))}
  }
  theta=Robst.LMDPDE.Beta.Reg(y,x,z,start_theta,alpha,tolerance)
  n=length(y)
  #Valores Preditos
  beta=theta$x[1:dim(as.matrix(x))[2]]
  gama=theta$x[(dim(as.matrix(x))[2]+1):(dim(as.matrix(x))[2]+dim(as.matrix(z))[2])]
  mu_hat=h1(x,beta)
  phi_hat=h2(z,gama)
  
  V=LMDPDE_Cov_Matrix(mu_hat,phi_hat,x,z,alpha=alpha)$Cov/n
  df=rankMatrix(M)[1]
  W_alpha=n*t(m)%*%solve(M%*%V%*%t(M))%*%(m)
  pValue=1-pchisq(W_alpha,df=df)
  
  result$W_alpha=as.numeric(W_alpha)
  result$df=df
  result$pValue=as.numeric(pValue)
  result$theta=theta$x
  result$fvec=theta$fvec
  
  return(result)
}

