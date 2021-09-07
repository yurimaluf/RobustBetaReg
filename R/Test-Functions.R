#' @export
somaC <- function(x, y) .Call(testeC, x, y)

#' @export
rbeta.cont=function(eta.mu,eta.phi,link = c("logit", "probit", "cloglog", "cauchit", "loglog"),link.phi = c("log","sqrt"),ep,seed,tau=1)
{
  #browser()
  if(missing(seed)){seed=as.numeric(Sys.time())}
  if(missing(ep)){ep=0}
  link = match.arg(link)
  link.phi = match.arg(link.phi)
  ep=round(ep)
  linkobj=set.link(link,link.phi)
  mu=linkobj$linkfun.mu$inv.link(eta.mu)
  phi=linkobj$linkfun.phi$inv.link(eta.phi)
  result=list()
  n=length(mu)
  set.seed(seed)
  y=rbeta(n,mu*phi,(1-mu)*phi)
  ind.c=NULL
  if(ep!=0)
  {
    ind.c=sample(seq(1,n),ep)
    for(i in ind.c)
    {
      mu_l=(1-mu[ind.c])^(tau)/((1-mu[ind.c])^(tau)+mu[ind.c]^tau)
      phi_l=(mu[ind.c]^(tau)+(1-mu[ind.c]^(tau)))/(phi[ind.c]*mu[ind.c]*(1-mu[ind.c]))^(tau)
      y[ind.c]=rbeta(1,mu_l*phi_l,(1-mu_l)*phi_l)  
    }
  }
  result$sample=pmax(pmin(y,1-.Machine$double.eps),.Machine$double.eps)
  result$ind=ind.c
  return(result)
}



#' 
rbeta_cont=function(mu,phi,ep,seed,tau)
{
  if(missing(seed)){seed=as.numeric(Sys.time())}
  if(missing(tau)){tau=1}
  limit=1e-5
  result=list()
  n=length(mu)
  set.seed(seed)
  y=rbeta(n,mu*phi,(1-mu)*phi)
  ind.c=NULL
  if(ep!=0)
  {
    ind.c=sample(seq(1,n),ep)
    for(i in ind.c)
    {
      mu_l=(1-mu[ind.c])^(tau)/((1-mu[ind.c])^(tau)+mu[ind.c]^tau)
      phi_l=(mu[ind.c]^(tau)+(1-mu[ind.c]^(tau)))/(phi[ind.c]*mu[ind.c]*(1-mu[ind.c]))^(tau)
      y[ind.c]=rbeta(1,mu_l*phi_l,(1-mu_l)*phi_l)  
    }
  }
  y[y<limit]=limit
  y[y>(1-limit)]=(1-limit)
  
  result$sample=y
  result$ind=ind.c
  return(result)
}

