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
  W.PHI=diag(x=phi_hat*V_star*((d.link.mu)^(-2)))
  
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
  D=diag(xi*(d.link.phi)^(-2))##Para funcao ligacao de phi
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
  
  y_star=log(y/(1-y))
  res.beta=sweighted_res(mu_hat,phi_hat,y)
  res.gamma=sweighted.gamma_res(mu_hat,phi_hat,y)
  v=trigamma(mu_hat*phi_hat)+trigamma((1-mu_hat)*phi_hat)
  b=mu_hat*v-trigamma((1-mu_hat)*phi_hat)
  xi=trigamma(mu_hat*phi_hat)*mu_hat^2+trigamma((1-mu_hat)*phi_hat)*(1-mu_hat)^2-trigamma(phi_hat)
  w=phi_hat*v*((d.link.mu)^(-2))
  TT=diag((d.link.mu)^(-1))
  HH=diag((d.link.phi)^(-1))
  B=diag(b)
  
  W.PHI=diag(x=phi_hat*w)
  W_PHI=diag(x=phi_hat/w)
  H=sqrt(W.PHI)%*%X%*%solve(t(X)%*%W.PHI%*%X)%*%t(X)%*%sqrt(W.PHI)
  
  D=diag(xi*(d.link.phi)^(-2))
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
hatvalues=function(object)
{
  UseMethod("hatvalues")
}

#' Create a Link for beta regression model
#' 
#' This function provides several link functions for beta regression models
#' 
#' For more details see:...
#' 
#' @param link.mu character; one of "logit"(default), "probit", "cauchit", "cloglog" and loglog
#' @param link.phi character; one of "log"(default) "indentity" and "sqrt"
#' @param thrd number (integer) of threads to speed up the process. If missing, the value is autodetected by the available number of multi-core processor
#' 
#' @return A structure with link function, inverse link function, derivative deta/dmu and d2eta/dmu2, for both models: mean and precision.
#'
#' @examples 
#' links=set.link(link.mu="cauchit",link.phi="sqrt")
#' attributes(links)
#' 
set.link=function(link.mu="logit",link.phi="log",thrd)
{
  if(missing(thrd)){thrd=tryCatch(parallel::detectCores(),error=function(e) 1)}
  thrd=tryCatch(suppressWarnings(max(1,round(abs(thrd)))),error=function(e) 1)
  #browser()
  switch(link.mu, logit = {
    linkfun <- function(mu) .Call(C_logit_link, mu, thrd )
    linkinv <- function(eta) .Call(C_logit_linkinv, eta, thrd)
    d.linkfun <- function(mu) .Call(C_logit_deta_dmu, mu, thrd)
    d2.linkfun <- function(mu) .Call(C_logit_d2eta_dmu2, mu, thrd)
    #valideta <- function(eta) TRUE
  }, probit = {
    linkfun <- function(mu){mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return(qnorm(mu))}
    linkinv <- function(eta) {thresh <- -qnorm(.Machine$double.eps)
      eta <- pmin(pmax(eta, -thresh), thresh)
      return(pnorm(eta))}
    d.linkfun <- function(mu){mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      return((dnorm(qnorm(mu)))^(-1))}
    d2.linkfun=function(mu){mu=pmax(pmin(mu,1-.Machine$double.eps),.Machine$double.eps)
      qmu=qnorm(mu)
      ddnorm=-qmu*exp(-qmu^(2)/2)/sqrt(2*pi)
      return(-ddnorm/(dnorm(qmu))^3)}
  }, cauchit = {
    linkfun <- function(mu) .Call(C_cauchit_link, mu, thrd)
    linkinv <- function(eta) .Call(C_cauchit_linkinv, eta, thrd)
    d.linkfun <- function(mu).Call(C_cauchit_deta_dmu, mu, thrd)
    d2.linkfun <- function(mu).Call(C_cauchit_d2eta_dmu2, mu, thrd)
    #valideta <- function(eta) TRUE
  }, cloglog = {
    linkfun <- function(mu) .Call(C_cloglog_link, mu, thrd)
    linkinv <- function(eta) .Call(C_cloglog_linkinv, eta, thrd)
    d.linkfun <- function(mu).Call(C_cloglog_deta_dmu, mu, thrd)
    d2.linkfun <- function(mu).Call(C_cloglog_d2eta_dmu2, mu, thrd)
    #valideta <- function(eta) TRUE
  }, loglog = {
    linkfun <- function(mu) .Call(C_loglog_link, mu, thrd)
    linkinv <- function(eta) .Call(C_loglog_linkinv, eta, thrd)
    d.linkfun <- function(mu).Call(C_loglog_deta_dmu, mu, thrd)
    d2.linkfun <- function(mu).Call(C_loglog_d2eta_dmu2, mu, thrd)
    #valideta <- function(eta) TRUE
  }, stop(gettextf("%s link.mu not recognised", sQuote(link.mu)),
          domain = NA))

  Linkfun.Mu=list(linkfun=linkfun,d.linkfun=d.linkfun,d2.linkfun=d2.linkfun,inv.link=linkinv)

  switch(link.phi, log = {
    linkfun.phi <- function(phi) .Call(C_log_link, phi, thrd)
    linkinv.phi <- function(eta) .Call(C_log_linkinv, eta, thrd)
    d.linkfun.phi <- function(phi) .Call(C_log_deta_dphi, phi, thrd)
    d2.linkfun.phi <- function(phi) .Call(C_log_d2eta_dphi2, phi, thrd)
    #linkfun.phi <- function(phi) as.numeric(log(pmax(phi,.Machine$double.eps)))
    #linkinv.phi <- function(eta) as.numeric(exp(eta))
    #d.linkfun.phi <- function(phi)  as.numeric((pmax(phi,.Machine$double.eps))^(-1))
    #d2.linkfun.phi <- function(phi) as.numeric(-(pmax(phi,.Machine$double.eps))^(-2))
    #valideta <- function(eta) TRUE
  }, identity = {
    linkfun.phi <- function(phi) pmax(phi,.Machine$double.eps)
    linkinv.phi <- function(eta) as.numeric(eta)
    d.linkfun.phi <- function(phi) rep(1,length(phi))
    d2.linkfun.phi <- function(phi) rep(0,length(phi))
    #valideta <- function(eta) TRUE
  }, sqrt = {
    linkfun.phi <- function(phi) .Call(C_sqrt_link, phi, thrd)
    linkinv.phi <- function(eta) .Call(C_sqrt_linkinv, eta, thrd)
    d.linkfun.phi <- function(phi) .Call(C_sqrt_deta_dphi, phi, thrd)
    d2.linkfun.phi <- function(phi) .Call(C_sqrt_d2eta_dphi2, phi, thrd)
    #valideta <- function(eta) TRUE
  }, stop(gettextf("%s link.phi not recognised", sQuote(link.phi)), domain = NA))

  Linkfun.Phi=list(linkfun=linkfun.phi,d.linkfun=d.linkfun.phi,d2.linkfun=d2.linkfun.phi,inv.link=linkinv.phi)

  return(structure(list(linkfun.mu=Linkfun.Mu,linkfun.phi=Linkfun.Phi),name.link.mu=link.mu,name.link.phi=link.phi,thrd=thrd,class="link-rbr"))
}

#'
set.link2=function(link.mu="logit",link.phi="log")
{
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
      return((mu-mu^2)^(-1))
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
      return((dnorm(qnorm(mu)))^(-1))
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
      return((-(1-mu)*log(1-mu))^(-1))
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
      return(2*pi*tan(pi*(mu-0.5))*pracma::sec(pi*(mu-0.5))^2)
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
      return(-(mu*log(mu))^(-1))
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
      return((phi)^(-1))
    }
    d2.linkfun.phi=function(phi)
    {
      phi=pmax(phi,.Machine$double.eps)
      return(-(phi^(-2)))
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
      return((2*sqrt(phi))^(-1))
    }
    d2.linkfun.phi=function(phi)
    {
      #return(-(0.25*phi^(-3/2)))
      return(-0.25*(phi^(-3/2)))
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


