#library(SQUAREM)
#library(gaussquad)


# If x is a n-column vector, turn it into n by 1 matrix
# If x is a matrix, keep it
#' Title
tomatrix = function(x){
  if(is.vector(x)){
    x = as.matrix(x)
  }
  return(x)
}

# To simplify computation, compute a part of post_pi_vash in advance
#' Title
getA = function(n,k,v,alpha.vec,modalpha.vec,sehat){
  A = v/2*log(v/2)-lgamma(v/2)+(v/2-1)*outer(rep(1,k),2*log(sehat))+outer(alpha.vec*log(modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2),rep(1,n))
  return(A)
}

#estimate mixture proportions of se's prior by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#' Title 
EMest_se = function(sehat,g,prior,maxiter=5000, v,unimodal, singlecomp, estpriormode){ 
  
  pi.init = g$pi
  k = ncomp(g)
  n = length(sehat)
  tol = min(0.1/n,1e-4) # set convergence criteria to be more stringent for larger samples
  
  if(unimodal=='variance'){
    c.init = g$beta[1]/(g$alpha[1]+1)
  }else if(unimodal=='precision'){
    c.init = g$beta[1]/(g$alpha[1]-1)
  }
  c.init = max(c.init,1e-5)
  
  EMfit = IGmixEM(sehat, v, c.init, g$alpha, pi.init, prior, unimodal,singlecomp, estpriormode, tol, maxiter)
  
  loglik = EMfit$B # actually return log lower bound not log-likelihood! 
  converged = EMfit$converged
  niter = EMfit$niter
  loglik.final = EMfit$B[length(EMfit$B)]
  
  g$pi=EMfit$pihat 
  g$c=EMfit$chat 
  if(singlecomp==TRUE){
    g$alpha=EMfit$alphahat
  }
  if(unimodal=='variance'){
    g$beta=g$c*(g$alpha+1)
  }else if(unimodal=='precision'){
    g$beta=g$c*(g$alpha-1)
  }
  
  return(list(loglik=loglik.final,converged=converged,g=g,niter=niter))
}


#' Title
IGmixEM = function(sehat, v, c.init, alpha.vec, pi.init, prior, unimodal,singlecomp, estpriormode, tol, maxiter){
  q = length(pi.init)
  n = length(sehat)
  
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  
  
  if(singlecomp==FALSE){
    if(estpriormode==TRUE){
      params.init=c(log(c.init),pi.init)
      A=getA(n=n,k=q,v,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,sehat=sehat)
      
      res = squarem(par=params.init,fixptfn=fixpoint_se, objfn=penloglik_se, 
                    A=A,n=n,k=q,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sehat=sehat,prior=prior,
                    control=list(maxiter=maxiter,tol=tol))
      return(list(chat = exp(res$par[1]), pihat=res$par[2:(length(res$par))], B=-res$value.objfn, 
                  niter = res$iter, converged=res$convergence))
    }else{
      A=getA(n=n,k=q,v,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,sehat=sehat)
      
      res = squarem(par=pi.init,fixptfn=fixpoint_pi, objfn=penloglik_pi, 
                    A=A,n=n,k=q,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sehat=sehat,prior=prior,c=c.init,
                    control=list(maxiter=maxiter,tol=tol))
      return(list(chat = c.init, pihat=res$par, B=-res$value.objfn, 
                  niter = res$iter, converged=res$convergence))
    }
    
  }else{
    params.init=c(log(c.init),log(alpha.vec))
    res=optim(params.init,loglike.se.ac,gr=gradloglike.se.ac,method='L-BFGS-B',
              lower=c(NA,0), upper=c(NA,log(100)),
              n=n,k=1,v=v,sehat=sehat,pi=pi,unimodal=unimodal,
              control=list(maxit=maxiter,pgtol=tol))
    return(list(chat = exp(res$par[1]), pihat=1, B=-res$value, 
                niter = res$counts[1], converged=res$convergence,
                alphahat=exp(res$par[2])))
  }
  
}


# Estimate the single inv-gamma prior distn params (moments matching)
# Prior: s^2~IG(a,b)
# sehat^2~s^2*Gamma(df/2,df/2)
#' Title
momentm = function(sehat,df){
  n = length(sehat)
  e = 2*log(sehat)-digamma(df/2)+log(df/2)
  ehat = mean(e)
  a = solve_trigamma(mean((e-ehat)^2*n/(n-1)-trigamma(df/2)))
  b = a*exp(ehat+digamma(df/2)-log(df/2))
  return(list(a=a,b=b))
}

# Solve trigamma(y)=x
#' Title
solve_trigamma = function(x){
  if(x > 1e7){
    y.new = 1/sqrt(x)
  }else if (x < 1e-6){
    y.new = 1/x
  }else{    
    y.old = 0.5+1/x
    delta = trigamma(y.old)*(1-trigamma(y.old)/x)/psigamma(y.old,deriv=2)
    y.new = y.old+delta
    while(-delta/y.new <= 1e-8){
      y.old = y.new
      delta = trigamma(y.old)*(1-trigamma(y.old)/x)/psigamma(y.old,deriv=2)
      y.new = y.old+delta
    }
  }
  return(y.new)
}

# prior of se: se|pi,alpha.vec,c ~ pi*IG(alpha_i,c*(alpha_i-1))
# Likelihood: sehat^2|se^2 ~ sj*Gamma(v/2,v/2)
# pi, alpha.vec, c: known
# Posterior weight of P(se|sehat) (IG mixture distn)
#' Title 
post_pi_vash = function(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi){
  post.pi.mat=t(pi*exp(A+alpha.vec*log(c)-(alpha.vec+v/2)*log(outer(c*modalpha.vec,v/2*sehat^2,FUN="+"))))
  return(pimat=post.pi.mat)
}

#' Title
fixpoint_se = function(params,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior){
  logc=params[1]
  pi=params[2:(length(params))]

  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,exp(logc),pi)
  m.rowsum = rowSums(mm)
  classprob = mm/m.rowsum
  newpi = colSums(classprob)+prior-1
  newpi = ifelse(newpi<1e-5,1e-5,newpi)
  newpi = newpi/sum(newpi);
  est=optim(logc,loglike.se,gr=gradloglike.se,method='BFGS',n=n,k=k,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sehat=sehat,pi=newpi)
  newc=exp(est$par[1])
  newlogc=est$par[1]
  params = c(newlogc,newpi)
  return(params)
}

#' Title
fixpoint_pi = function(pi,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior,c){  
  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  classprob = mm/m.rowsum
  newpi = colSums(classprob)+prior-1
  newpi = ifelse(newpi<1e-5,1e-5,newpi)
  newpi = newpi/sum(newpi);
  return(newpi)
}

#' Title
penloglik_se = function(params,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior){
  c=exp(params[1])
  pi=params[2:(length(params))]
  priordens = sum((prior-1)*log(pi))
  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik = sum(log(m.rowsum))
  return(-(loglik+priordens))
}

#' Title
penloglik_pi = function(pi,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior,c){
  priordens = sum((prior-1)*log(pi))
  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik = sum(log(m.rowsum))
  return(-(loglik+priordens))
}

# Log-likelihood: L(sehat^2|c,pi,alpha.vec)
#' Title 
loglike.se = function(logc,n,k,alpha.vec,modalpha.vec,v,sehat,pi){  
  c=exp(logc)
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike.se (w.r.t logc)
#' Title
gradloglike.se = function(logc,n,k,alpha.vec,modalpha.vec,v,sehat,pi){
  c=exp(logc)
  
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  classprob = pimat/rowSums(pimat)
  gradmat = c*classprob*(outer(rep(1,n),alpha.vec/c)
                         -outer(rep(1,n),(alpha.vec+v/2)*modalpha.vec)/
                           (outer(v/2*sehat^2,c*modalpha.vec,FUN='+')))
  grad = sum(-gradmat)
  return(grad)
}

# Log-likelihood: L(sehat|c,pi,alpha.vec)
#' Title
loglike.se.a = function(logalpha.vec,c,n,k,v,sehat,pi,unimodal){  
  alpha.vec=exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike.se for single component prior (w.r.t logalpha)
#' Title
gradloglike.se.a = function(logalpha.vec,c,n,k,v,sehat,pi,unimodal){
  alpha.vec=exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  grad=-alpha.vec*sum(log(c)+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+v/2)
                      -c*(alpha.vec+v/2)/(c*modalpha.vec+v/2*sehat^2)-log(c*modalpha.vec+v/2*sehat^2))
  return(grad)
}

# Log-likelihood: L(sehat|c,pi,alpha.vec)
#' Title 
loglike.se.ac = function(params,n,k,v,sehat,pi,unimodal){
  c=exp(params[1])
  alpha.vec=exp(params[2:length(params)])
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 +outer(rep(1,n),alpha.vec)*(log(outer(rep(1,n),c*modalpha.vec))-log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
                                 -(v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  #classprob=pimat/rowSums(pimat)
  logl = sum(log(rowSums(pimat)))
  logl = min(logl,1e200)
  logl = max(logl,-1e200)
  return(-logl)
}

# Gradient of funtion loglike.se for single component prior (w.r.t logc and logalpha)
#' Title
gradloglike.se.ac=function(params,n,k,v,sehat,pi,unimodal){
  c=exp(params[1])
  alpha.vec=exp(params[2:(length(params))])
  if(unimodal=='variance'){
    modalpha.vec=alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec=alpha.vec-1
  }
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  
  classprob = pimat/rowSums(pimat)
  gradmat.c = c*classprob*(outer(rep(1,n),alpha.vec/c)
                           -outer(rep(1,n),(alpha.vec+v/2)*modalpha.vec)/
                             (outer(v/2*sehat^2,c*modalpha.vec,FUN='+')))
  grad.c = sum(-gradmat.c)
  
  grad.a=-alpha.vec*sum(log(c)+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+v/2)
                        -c*(alpha.vec+v/2)/(c*modalpha.vec+v/2*sehat^2)-log(c*modalpha.vec+v/2*sehat^2))
  res = c(grad.c,grad.a)
  res = pmin(res,1e200)
  res = pmax(res,-1e200)
  return(res)
}

#compute posterior shape (alpha1) and rate (beta1)
#' Title
post.igmix = function(m,betahat,sebetahat,v){
  n = length(sebetahat)
  alpha1 = outer(rep(1,n),m$alpha+v/2)
  beta1 = outer(m$beta,v/2*sebetahat^2,FUN="+")
  ismissing = is.na(sebetahat)
  beta1[,ismissing]=m$beta
  return(list(alpha=alpha1,beta=t(beta1)))
}

# Moderated t test
#' Title
mod_t_test=function(betahat,se,pi,v){
  n=length(betahat)
  k=length(pi)/n
  pvalue=rep(NA,n)
  completeobs=(!is.na(betahat) & !is.na(apply(se,1,sum)))
  temppvalue=pt(outer(betahat[completeobs],rep(1,k))/se[completeobs,],df=v,lower.tail=TRUE)
  temppvalue=pmin(temppvalue,1-temppvalue)*2
  pvalue[completeobs]=apply(pi[completeobs,]*temppvalue,1,sum)
  return(pvalue)
}


vash.core = function(sehat,df,betahat,randomstart,singlecomp,unimodal,
                     prior,g,maxiter,estpriormode,priormode,completeobs){
  if(!is.null(g)){
    maxiter = 1 # if g is specified, don't iterate the EM
    prior = rep(1,length(g$pi)) #prior is not actually used if g specified, but required to make sure EM doesn't produce warning
    l = length(g$pi)
  } else {   
    mm = momentm(sehat[completeobs],df)
    mm$a = max(mm$a,1e-5)
    if(singlecomp==TRUE){
      alpha = mm$a
    }else{
      if(mm$a>1){
        alpha = 1+((64/mm$a)^(1/6))^seq(-3,6)*(mm$a-1)
      }else{
        alpha = mm$a*2^seq(0,13)
      }
      #alpha=alpha[alpha<=70] # avoid numerical error
    }
    mm$b = max(mm$b, 1e-5)
    if(unimodal=='precision'){
      alpha = unique(pmax(alpha,1+1e-5)) # alpha<=1 not allowed
      beta = mm$b/(mm$a-1)*(alpha-1)
    }else if(unimodal=='variance'){
      beta = mm$b/(mm$a+1)*(alpha+1)
    }
    
    if(estpriormode==FALSE & !is.null(priormode)){
      if(unimodal=='precision'){
        beta = priormode*(alpha-1)
      }else if(unimodal=='variance'){
        beta = priormode*(alpha+1)
      }
    }else if (estpriormode==TRUE & !is.null(priormode)){
      warning('Flag estpriormode=TRUE, vash will still estimate the prior mode instead of using the input value!')
    }       
    l = length(alpha)
  }
  
  
  if(is.null(prior)){
    prior = rep(1,l)
  }
  
  if(randomstart){
    pi.se = rgamma(l,1,1)
  } else {   
    pi.se=rep(1,l)/l
  }
  pi.se=pi.se/sum(pi.se)
  
  if (!is.null(g)){
    g=igmix(g$pi,g$alpha,g$beta)
  }else{
    g=igmix(pi.se,alpha,beta)
  }
  
  if(length(prior)!=l){
    stop("invalid prior specification")
  }
  
  pi.fit.se = EMest_se(sehat[completeobs],g,prior,maxiter,df,unimodal,singlecomp,estpriormode)
  return(pi.fit.se)
}



#' 
#' @title  Main Variance Adaptive SHrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sehat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param sehat, a p vector of standard errors
#' @param df: appropriate degrees of freedom for (chi-square) distribution of sehat
#' @param betahat: a p vector of estimates
#' @param randomstart: bool, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param singlecomp: bool, indicating whether to use a single inverse-gamma distribution as the prior distribution for the variances
#' @param unimodal: unimodal constraint for the prior distribution of the variances ("variance") or the precisions ("precision")
#' @param prior: string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param g: the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param maxiter: maximum number of iterations of the EM algorithm
#' @param estpriormode: bool, indicating whether to estimate the mode of the unimodal prior
#' @param priormode: specified prior mode (only works when estpriormode=FALSE).
#'
#' @return a list with elements fitted.g is fitted mixture
#' 
#' @export
#' 
#' @examples 
#' se = rigamma(100,1,1)
#' sehat = se*rchisq(100,5)/5
#' se.vash = vash(sehat,5)
#' plot(sehat,se.vash$PosteriorMean,xlim=c(0,10),ylim=c(0,10))
vash = function(sehat,df,
                betahat=NULL,
                randomstart=FALSE,   
                singlecomp = FALSE,
                unimodal = c("auto","variance","precision"),
                prior=NULL,
                g=NULL,
                maxiter = 5000,
                estpriormode = TRUE,
                priormode = NULL){
  
  #method provides a convenient interface to set a particular combinations of parameters for prior an
  #If method is supplied, use it to set up specific values for these parameters; provide warning if values
  #are also specified by user
  #If method is not supplied use the user-supplied values (or defaults if user does not specify them)
  
  
  
  if(missing(unimodal)){
    unimodal = match.arg(unimodal) 
  }
  if(!is.element(unimodal,c("auto","variance","precision"))) stop("Error: invalid type of unimodal.")  
  
  completeobs = (!is.na(sehat))
  n=sum(completeobs)
  
  # If some standard errors are almost 0, add a small pseudocount to prevent numerical errors
  sehat[sehat==0] = min(min(sehat[sehat>0]),1e-6)
  
  if(n==0){
    stop("Error: all input values are missing.")
  }  
  
  if(unimodal=='auto' & !singlecomp){
    pifit.prec = vash.core(sehat,df,betahat,randomstart,singlecomp,unimodal='precision',
                          prior,g,maxiter,estpriormode,priormode,completeobs)
    pifit.var = vash.core(sehat,df,betahat,randomstart,singlecomp,unimodal='variance',
                          prior,g,maxiter,estpriormode,priormode,completeobs)
    if (pifit.prec$loglik >= pifit.var$loglik){
      pi.fit.se = pifit.prec
      opt.unimodal = 'precision'
    }else{
      pi.fit.se = pifit.var
      opt.unimodal = 'variance'
    }
  }else if(unimodal=='auto' & singlecomp){
    pi.fit.se = vash.core(sehat,df,betahat,randomstart,singlecomp,unimodal='variance',
                          prior,g,maxiter,estpriormode,priormode,completeobs)
    opt.unimodal = NA
  }else{
    pi.fit.se = vash.core(sehat,df,betahat,randomstart,singlecomp,unimodal,
                          prior,g,maxiter,estpriormode,priormode,completeobs)
    opt.unimodal = NA
  }
  
  post.se = post.igmix(pi.fit.se$g,rep(numeric(0),n),sehat[completeobs],df)
  postpi.se = t(matrix(rep(pi.fit.se$g$pi,length(sehat)),ncol=length(sehat)))
  postpi.se[completeobs,] = t(comppostprob(pi.fit.se$g,rep(numeric(0),n),sehat[completeobs],df))
  
  PosteriorMean.se = rep(pi.fit.se$g$c,length=length(sehat))
  #PosteriorSD.se = rep(0,length=n)
  
  
  #PosteriorSD.se[completeobs] = postsd(pi.fit.se$g,NULL,sehat[completeobs],df) 
  PosteriorShape.se = t(matrix(rep(pi.fit.se$g$alpha,length(sehat)),ncol=length(sehat)))
  PosteriorShape.se[completeobs,] = post.se$alpha
  PosteriorRate.se = t(matrix(rep(pi.fit.se$g$beta,length(sehat)),ncol=length(sehat)))
  PosteriorRate.se[completeobs,] = post.se$beta
  
  #PosteriorMean.se[completeobs] = sqrt(postmean(pi.fit.se$g,rep(numeric(0),n),sehat[completeobs],df))
  #PosteriorMean.se = sqrt(apply(postpi.se*PosteriorRate.se/PosteriorShape.se,1,sum))
  PosteriorMean.se[completeobs] = sqrt(1/apply(postpi.se*PosteriorShape.se/PosteriorRate.se,1,sum))

  
  if(length(betahat)==n){
    pvalue=mod_t_test(betahat,sqrt(PosteriorRate.se/PosteriorShape.se),
                      postpi.se,PosteriorShape.se*2)
    qvalue=qvalue(pvalue)$qval
  }else if(length(betahat)==0){
    pvalue=NULL
    qvalue=NULL
  }else{
    warning("betahat has different length as sehat, cannot compute moderated t-tests")
    pvalue=NULL
    qvalue=NULL
  }
  
  result = list(fitted.g=pi.fit.se$g,
                PosteriorMean=PosteriorMean.se,
                #PosteriorSD=PosteriorSD.se,
                PosteriorShape=PosteriorShape.se,
                PosteriorRate=PosteriorRate.se,
                PosteriorPi=postpi.se,
                pvalue=pvalue,
                qvalue=qvalue,
                unimodal=unimodal,
                opt.unimodal=opt.unimodal,
                fit=pi.fit.se,call=match.call(),data=list(sehat=sehat))
  class(result)= "vash"
  return(result)
}