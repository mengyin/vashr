#' @import ashr qvalue SQUAREM
#' @importFrom stats optim pt rgamma
#' @title  Main Variance Adaptive SHrinkage function
#'
#' @description Takes vectors of standard errors (sehat), and applies shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for variances.
#'
#' @details See vignette for more details.
#' 
#' @param sehat a p vector of observed standard errors
#' @param df scalar, appropriate degree of freedom for (chi-square) distribution of sehat
#' @param betahat a p vector of estimates (optional)
#' @param randomstart logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param singlecomp logical, indicating whether to use a single inverse-gamma distribution as the prior distribution for the variances
#' @param unimodal put unimodal constraint on the prior distribution of variances ("variance") or precisions ("precision"). This also can be automatically chosen ("auto") by comparing the likelihoods.
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param g the prior distribution for variances (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param estpriormode logical, indicating whether to estimate the mode of the unimodal prior
#' @param priormode specified prior mode (only works when estpriormode=FALSE).
#' @param scale a scalar or a p vector, such that sehat/scale are exchangeable, i.e, can be modeled by a common unimodal prior. 
#' @param maxiter maximum number of iterations of the EM algorithm
#'
#' @return vash returns an object of \code{\link[base]{class}} "vash", a list with the following elements
#' \item{fitted.g}{Fitted mixture prior}
#' \item{sd.post}{A vector consisting the posterior estimate of standard deviations (square root of variances) from the mixture}
#' \item{PosteriorPi}{A p*k matrix consisting the mixture proportions of the posterior distribution of variance (k is number of mixture components)}
#' \item{PosteriorShape}{A p*k matrix consisting the shape parameters of the inverse-gamma posterior distribution of variance (k is number of mixture components)}
#' \item{PosteriorRate}{A p*k matrix consisting the rate parameters of the inverse-gamma posterior distribution of variance (k is number of mixture components)}
#' \item{pvalue}{A vector of p-values}
#' \item{qvalue}{A vector of q-values}
#' \item{fit}{The fitted mixture object}
#' \item{unimodal}{Denoting whether unimodal variance prior or unimodal precision prior has been used}
#' \item{opt.unimodal}{Denoting whether unimodal variance prior or unimodal precision prior model gets higher likelihood}
#' \item{call}{A call in which all of the specified arguments are specified by their full names}
#' \item{data}{A list consisting the input sehat, df and betahat}
#' 
#' @export
#' 
#' @examples 
#' ##An simple example
#' #generate true variances (sd^2) from an inverse-gamma prior
#' sd = sqrt(1/rgamma(100,5,5)) 
#' #observed standard errors are estimates of true sd's
#' sehat = sqrt(sd^2*rchisq(100,7)/7) 
#' #run the vash function
#' fit = vash(sehat,df=7) 
#' #plot the shrunk sd estimates against the observed standard errors
#' plot(sehat, fit$sd.post, xlim=c(0,10), ylim=c(0,10)) 
#' 
#' ##Running vash with a pre-specified g, rather than estimating it
#' sd = sqrt(c(1/rgamma(100,5,4),1/rgamma(100,10,9)))
#' sehat = sqrt(sd^2*rchisq(100,7)/7)
#' true_g = igmix(c(0.5,0.5),c(5,10),c(4,9)) # define true g 
#' #Passing this g into vash causes it to i) take the shape and the rate for each component from this g, 
#' #and ii) initialize pi to the value from this g.
#' se.vash = vash(sehat, df=7, g=true_g)
#' 
#' @references Lu, M., & Stephens, M. (2016). Variance Adaptive Shrinkage (vash): Flexible Empirical Bayes estimation of variances. bioRxiv, 048660.
#' 
vash = function(sehat, df,
                betahat = NULL,
                randomstart = FALSE,   
                singlecomp = FALSE,
                unimodal = c("auto","variance","precision"),
                prior = NULL,
                g = NULL,
                estpriormode = TRUE,
                priormode = NULL,
                scale = 1,
                maxiter = 5000){
  
  ## 1.Handling Input Parameters
  
  # Set unimodal constraint on variance or precision prior
  # if no unimodal mode specified, select a default "auto" mode,
  # which chooses the model with higher likelihood
  if(missing(unimodal)){
    unimodal = match.arg(unimodal) 
  }
  if(!is.element(unimodal,c("auto","variance","precision"))) stop("Error: invalid type of unimodal.")  
  if(!missing(betahat) & length(betahat)!=length(sehat)) stop("Error: sehat and betahat must have same lengths.")
  if(!missing(g) & class(g)!="igmix") stop("Error: invalid type of g.")
  if(!is.numeric(df) | length(df)>1) stop("Error: invalid type of df.")
  
  # Set observations with infinite standard errors to missing
  # later these missing observations will be ignored in EM, and posterior will be same as prior.
  sehat[sehat==Inf] = NA
  sehat = sehat/scale
  if(min(sehat)<0) stop("Error: sehat/scale must be non-negative.")
  completeobs = (!is.na(sehat))
  n = sum(completeobs)
  
  # If some standard errors are almost 0, add a small pseudocount to prevent numerical errors
  sehat[sehat==0] = min(min(sehat[sehat>0]),1e-6)
  
  if(n==0){
    stop("Error: all input values are missing.")
  }  
  
  ## 2. Fitting the mixture prior
  
  # If singlecomp==TRUE, then both variance and precision priors are unimodal
  # otherwise need decide whether using unimodal variance prior model or 
  # unimodal precision prior model by likelihood
  if(unimodal=='auto' & !singlecomp){
    pifit.prec = est_prior(sehat,df,betahat,randomstart,singlecomp,unimodal='precision',
                           prior,g,maxiter,estpriormode,priormode,completeobs)
    pifit.var = est_prior(sehat,df,betahat,randomstart,singlecomp,unimodal='variance',
                          prior,g,maxiter,estpriormode,priormode,completeobs)
    if (pifit.prec$loglik >= pifit.var$loglik){
      mix.fit = pifit.prec
      opt.unimodal = 'precision'
    }else{
      mix.fit = pifit.var
      opt.unimodal = 'variance'
    }
  }else if(unimodal=='auto' & singlecomp){
    mix.fit = est_prior(sehat,df,betahat,randomstart,singlecomp,unimodal='variance',
                        prior,g,maxiter,estpriormode,priormode,completeobs)
    opt.unimodal = NA
  }else{
    mix.fit = est_prior(sehat,df,betahat,randomstart,singlecomp,unimodal,
                        prior,g,maxiter,estpriormode,priormode,completeobs)
    opt.unimodal = NA
  }
  
  ## 3. Posterior inference
  # compute posterior distribution
  post.se = post.igmix(mix.fit$g,rep(numeric(0),n),sehat[completeobs],df)
  postpi.se = t(matrix(rep(mix.fit$g$pi,length(sehat)),ncol=length(sehat)))
  #postpi.se[completeobs,] = t(comppostprob(mix.fit$g,rep(numeric(0),n),sehat[completeobs],df))
  data = list(x=rep(numeric(0),n), s=sehat[completeobs], v=df)
  postpi.se[completeobs,] = t(comp_postprob(mix.fit$g, data = data))
  
  PosteriorMean.se = rep(mix.fit$g$c,length=length(sehat))
  
  PosteriorShape.se = t(matrix(rep(mix.fit$g$alpha,length(sehat)),ncol=length(sehat)))
  PosteriorShape.se[completeobs,] = post.se$alpha
  PosteriorRate.se = t(matrix(rep(mix.fit$g$beta,length(sehat)),ncol=length(sehat)))
  PosteriorRate.se[completeobs,] = post.se$beta
  
  PosteriorMean.se[completeobs] = sqrt(1/apply(postpi.se*PosteriorShape.se/PosteriorRate.se,1,sum))
  
  # obtain p-values by moderated t-test
  # and then compute q-values from the p-values
  if(length(betahat)==n){
    pvalue = mod_t_test(betahat,scale*sqrt(PosteriorRate.se/PosteriorShape.se),
                        postpi.se,PosteriorShape.se*2)
    qvalue = qvalue(pvalue)$qval
  }else if(length(betahat)==0){
    pvalue = NULL
    qvalue = NULL
  }else{
    warning("betahat has different length as sehat, cannot compute moderated t-tests")
    pvalue = NULL
    qvalue = NULL
  }
  
  result = list(fitted.g=mix.fit$g,
                sd.post=PosteriorMean.se,
                PosteriorPi=postpi.se,
                PosteriorShape=PosteriorShape.se,
                PosteriorRate=PosteriorRate.se,
                pvalue=pvalue,
                qvalue=qvalue,
                unimodal=unimodal,
                opt.unimodal=opt.unimodal,
                fit=mix.fit,call=match.call(),data=list(sehat=sehat,betahat=betahat,df=df))
  class(result) = "vash"
  return(result)
}


# If x is a n-column vector, turn it into n by 1 matrix
# If x is a matrix, keep it
tomatrix = function(x){
  if(is.vector(x)){
    x = as.matrix(x)
  }
  return(x)
}

# To simplify computation, compute a part of post_pi_vash in advance
getA = function(n,k,v,alpha.vec,modalpha.vec,sehat){
  A = v/2*log(v/2)-lgamma(v/2)+(v/2-1)*outer(rep(1,k),2*log(sehat))+outer(alpha.vec*log(modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2),rep(1,n))
  return(A)
}

#' 
#' @title Compute log-likelihood of the variance model with mixture inverse-gamma prior
#' 
#' @description Suppose we observe standard errors sehat, where \eqn{sehat^2~s^2 \times \chi^2_{df}/df}, and $s^2$ comes from an inverse-gamma mixture prior. This function computes the log-likelihood of the model.
#' 
#' @param sehat a p vector of observed standard errors
#' @param df the degree of freedom
#' @param pi a k vector, the mixture proportions of the k component inverse-gamma mixture prior
#' @param alpha a k vector, the shape parameters of the k component inverse-gamma mixture prior
#' @param beta a k vector, the rate parameters of the k component inverse-gamma mixture prior
#'
#' @return log-likehihood
#' 
#' @examples 
#' sehat = abs(rnorm(10))
#' loglike(sehat, df=10, pi=c(0.1,0.3,0.5), alpha=c(5,8,10), beta=c(2,4,9))
#' 
#' @export
#' 
loglike = function(sehat,df,pi,alpha,beta){
  k = length(pi)
  n = length(sehat)
  
  pimat = outer(rep(1,n),pi)*exp(df/2*log(df/2)-lgamma(df/2)
                                 +(df/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha*log(beta)-lgamma(alpha)+lgamma(alpha+df/2))
                                 -outer(rep(1,n),alpha+df/2)*log(outer(df/2*sehat^2,beta,FUN="+")))
  logl = sum(log(rowSums(pimat)))
  return(logl)
}

#' @title Estimate mixture proportions and mode of the unimodal inverse-gaama mixture variance prior.
#' @description Estimate mixture proportions and mode of the unimodal inverse-gaama mixture variance prior.
#'
#' @param sehat n vector of standard errors of observations
#' @param g the initial prior distribution for variances
#' @param prior numeric vector indicating Dirichlet prior on mixture proportions
#' @param df appropriate degrees of freedom for chi-square distribution of sehat^2
#' @param unimodal put unimodal constraint on the prior distribution of variances ("variance") or precisions ("precision")
#' @param singlecomp logical, indicating whether to use a single inverse-gamma distribution as the prior distribution for the variances
#' @param estpriormode logical, indicating whether to estimate the mode of the unimodal prior
#' @param maxiter maximum number of iterations of the EM algorithm
#' 
#' @return A list, including the final loglikelihood, the fitted prior g, number of iterations and a flag to indicate convergence.
#' 
#' @examples 
#' fitted.prior = est_mixprop_mode(sehat=abs(rnorm(100)),g=igmix(c(.5,.5),c(1,3),c(1,3)),
#' prior=c(1,1),df=10,unimodal="variance",singlecomp=FALSE,estpriormode=TRUE)
#' 
#' @export
#' 
est_mixprop_mode = function(sehat, g, prior, df, unimodal, singlecomp, estpriormode, maxiter=5000){ 
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
  
  EMfit = IGmixEM(sehat, df, c.init, g$alpha, pi.init, prior, unimodal,singlecomp, estpriormode, tol, maxiter)
  
  converged = EMfit$converged
  niter = EMfit$niter
  
  g$pi = EMfit$pihat 
  g$c = EMfit$chat 
  if(singlecomp==TRUE){
    g$alpha = EMfit$alphahat
  }
  if(unimodal=='variance'){
    g$beta = g$c*(g$alpha+1)
  }else if(unimodal=='precision'){
    g$beta = g$c*(g$alpha-1)
  }
  loglik = loglike(sehat,df,g$pi,g$alpha,g$beta)
  
  return(list(loglik=loglik, converged=converged, g=g, niter=niter))
}

# Estimate the mixture inverse-gamma prior by EM algorithm
IGmixEM = function(sehat, v, c.init, alpha.vec, pi.init, prior, unimodal,singlecomp, estpriormode, tol, maxiter){
  q = length(pi.init)
  n = length(sehat)
  
  if(unimodal=='variance'){
    modalpha.vec = alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec = alpha.vec-1
  }
  
  if(singlecomp==FALSE){
    if(estpriormode==TRUE){
      params.init = c(log(c.init),pi.init)
      A = getA(n=n,k=q,v,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,sehat=sehat)
      
      res = squarem(par=params.init,fixptfn=fixpoint_se, objfn=penloglik_se, 
                    A=A,n=n,k=q,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sehat=sehat,prior=prior,
                    control=list(maxiter=maxiter,tol=tol))
      return(list(chat = exp(res$par[1]), pihat=res$par[2:(length(res$par))], B=-res$value.objfn, 
                  niter = res$iter, converged=res$convergence))
    }else{
      A = getA(n=n,k=q,v,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,sehat=sehat)
      
      res = squarem(par=pi.init, fixptfn=fixpoint_pi, objfn=penloglik_pi, 
                    A=A,n=n,k=q,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sehat=sehat,prior=prior,c=c.init,
                    control=list(maxiter=maxiter,tol=tol))
      return(list(chat=c.init, pihat=res$par, B=-res$value.objfn, 
                  niter=res$iter, converged=res$convergence))
    }
    
  }else{
    params.init = c(log(c.init),log(alpha.vec))
    res = optim(params.init,loglike.se.ac,gr=gradloglike.se.ac,method='L-BFGS-B',
              lower=c(NA,-3), upper=c(NA,log(100)),
              n=n,k=1,v=v,sehat=sehat,pi=pi,unimodal=unimodal,
              control=list(maxit=maxiter,pgtol=tol))
    return(list(chat=exp(res$par[1]), pihat=1, B=-res$value, 
                niter=res$counts[1], converged=res$convergence,
                alphahat=exp(res$par[2])))
  }
  
}


#' @title Estimate the single component inverse-gamma prior of variances by matching the moments
#' @description Suppose the true variances come from a single component inverse-gamma prior, and 
#' standard errors (noisy estimates of the square root of true variances) are observed.
#' This function estimates the single component inverse-gamma prior of variances by matching the moments. 
#'
#' @param sehat n vector of standard errors of observations
#' @param df degrees of freedom for chi-square distribution of sehat^2
#' 
#' @return The shape (a) and rate (b) of the fitted inverse-gamma prior 
#' @examples 
#' sd = sqrt(1/rgamma(100,5,5)) #true prior: IG(5,5)
#' sehat = sqrt(sd^2*rchisq(100,7)/7) 
#' est_singlecomp_mm(sehat=sehat,df=7) 
#' 
#' @export
#' 
est_singlecomp_mm = function(sehat,df){
  n = length(sehat)
  e = 2*log(sehat)-digamma(df/2)+log(df/2)
  ehat = mean(e)
  a = solve_trigamma(mean((e-ehat)^2*n/(n-1)-trigamma(df/2)))
  b = a*exp(ehat+digamma(df/2)-log(df/2))
  return(list(a=a,b=b))
}

# Solve trigamma(y)=x
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
post_pi_vash = function(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi){
  post.pi.mat = t(pi*exp(A+alpha.vec*log(c)-(alpha.vec+v/2)*log(outer(c*modalpha.vec,v/2*sehat^2,FUN="+"))))
  return(pimat=post.pi.mat)
}

# fix point function of updating both pi and c
fixpoint_se = function(params,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior){
  logc = params[1]
  pi = params[2:(length(params))]

  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,exp(logc),pi)
  m.rowsum = rowSums(mm)
  classprob = mm/m.rowsum
  newpi = colSums(classprob)+prior-1
  newpi = ifelse(newpi<1e-5,1e-5,newpi)
  newpi = newpi/sum(newpi);
  est = optim(logc,loglike.se,gr=gradloglike.se,method='BFGS',n=n,k=k,alpha.vec=alpha.vec,modalpha.vec=modalpha.vec,v=v,sehat=sehat,pi=newpi)
  newc = exp(est$par[1])
  newlogc = est$par[1]
  params = c(newlogc,newpi)
  return(params)
}

# fix point function of updating pi (c known)
fixpoint_pi = function(pi,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior,c){  
  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  classprob = mm/m.rowsum
  newpi = colSums(classprob)+prior-1
  newpi = ifelse(newpi<1e-5,1e-5,newpi)
  newpi = newpi/sum(newpi);
  return(newpi)
}

# penalized log-likelihood (c unknown)
penloglik_se = function(params,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior){
  c = exp(params[1])
  pi = params[2:(length(params))]
  priordens = sum((prior-1)*log(pi))
  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik = sum(log(m.rowsum))
  return(-(loglik+priordens))
}

# penalized log-likelihood (c known)
penloglik_pi = function(pi,A,n,k,alpha.vec,modalpha.vec,v,sehat,prior,c){
  priordens = sum((prior-1)*log(pi))
  mm = post_pi_vash(A,n,k,v,sehat,alpha.vec,modalpha.vec,c,pi)
  m.rowsum = rowSums(mm)
  loglik = sum(log(m.rowsum))
  return(-(loglik+priordens))
}

# Log-likelihood: L(sehat^2|c,pi,alpha.vec)
loglike.se = function(logc,n,k,alpha.vec,modalpha.vec,v,sehat,pi){  
  c = exp(logc)
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike.se (w.r.t logc)
gradloglike.se = function(logc,n,k,alpha.vec,modalpha.vec,v,sehat,pi){
  c = exp(logc)
  
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
loglike.se.a = function(logalpha.vec,c,n,k,v,sehat,pi,unimodal){  
  alpha.vec = exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec = alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec = alpha.vec-1
  }
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),alpha.vec*log(c*modalpha.vec)-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 -outer(rep(1,n),alpha.vec+v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
  logl = sum(log(rowSums(pimat)))
  return(-logl)
}

# Gradient of funtion loglike.se for single component prior (w.r.t logalpha)
gradloglike.se.a = function(logalpha.vec,c,n,k,v,sehat,pi,unimodal){
  alpha.vec = exp(logalpha.vec)
  if(unimodal=='variance'){
    modalpha.vec = alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec = alpha.vec-1
  }
  grad = -alpha.vec*sum(log(c)+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+v/2)
                      -c*(alpha.vec+v/2)/(c*modalpha.vec+v/2*sehat^2)-log(c*modalpha.vec+v/2*sehat^2))
  return(grad)
}

# Log-likelihood: L(sehat|c,pi,alpha.vec)
loglike.se.ac = function(params,n,k,v,sehat,pi,unimodal){
  c = exp(params[1])
  alpha.vec = exp(params[2:length(params)])
  if(unimodal=='variance'){
    modalpha.vec = alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec = alpha.vec-1
  }
  pimat = outer(rep(1,n),pi)*exp(v/2*log(v/2)-lgamma(v/2)
                                 +(v/2-1)*outer(2*log(sehat),rep(1,k))
                                 +outer(rep(1,n),-lgamma(alpha.vec)+lgamma(alpha.vec+v/2))
                                 +outer(rep(1,n),alpha.vec)*(log(outer(rep(1,n),c*modalpha.vec))-log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))
                                 -(v/2)*log(outer(v/2*sehat^2,c*modalpha.vec,FUN="+")))

  logl = sum(log(rowSums(pimat)))
  logl = min(logl,1e200) # avoid numerical overflow
  logl = max(logl,-1e200)
  return(-logl)
}

# Gradient of funtion loglike.se for single component prior (w.r.t logc and logalpha)
gradloglike.se.ac=function(params,n,k,v,sehat,pi,unimodal){
  c = exp(params[1])
  alpha.vec = exp(params[2:(length(params))])
  if(unimodal=='variance'){
    modalpha.vec = alpha.vec+1
  }else if(unimodal=='precision'){
    modalpha.vec = alpha.vec-1
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
  
  grad.a = -alpha.vec*sum(log(c)+log(modalpha.vec)+alpha.vec/modalpha.vec-digamma(alpha.vec)+digamma(alpha.vec+v/2)
                        -c*(alpha.vec+v/2)/(c*modalpha.vec+v/2*sehat^2)-log(c*modalpha.vec+v/2*sehat^2))
  res = c(grad.c,grad.a)
  res = pmin(res,1e200)
  res = pmax(res,-1e200)
  return(res)
}

# compute posterior shape (alpha1) and rate (beta1)
post.igmix = function(m,betahat,sebetahat,v){
  n = length(sebetahat)
  alpha1 = outer(rep(1,n),m$alpha+v/2)
  beta1 = outer(m$beta,v/2*sebetahat^2,FUN="+")
  ismissing = is.na(sebetahat)
  beta1[,ismissing] = m$beta
  return(list(alpha=alpha1,beta=t(beta1)))
}

#' 
#' @title Moderated t-test with standard errors moderated by mixture prior
#' @description Moderated t-test with standard errors moderated by mixture prior. 
#' 
#' @param betahat a p vector of observed effect sizes
#' @param se a p*k matrix of moderated standard errors
#' @param pi a p*k matrix of mixture proportions, row sums are 1.
#' @param df a p*k matrix of degrees of freedom
#'
#' @return A p vector of p-values testing if the effect sizes are 0. 
#' 
#' @examples 
#' # p=10, k=5
#' betahat = rnorm(10)
#' se = matrix(abs(rnorm(10*5)), ncol=5)
#' pi = matrix(rep(0.2,10*5), ncol=5)
#' df = matrix(rep(11:15,each=10), ncol=5)
#' mod_t_test(betahat,se,pi,df)
#' 
#' @export
#' 
mod_t_test=function(betahat,se,pi,df){
  n = length(betahat)
  k = length(pi)/n
  pvalue = rep(NA,n)
  completeobs = (!is.na(betahat) & !is.na(apply(se,1,sum)))
  temppvalue = pt(outer(betahat[completeobs],rep(1,k))/se[completeobs,],df=df,lower.tail=TRUE)
  temppvalue = pmin(temppvalue,1-temppvalue)*2
  pvalue[completeobs] = apply(pi[completeobs,]*temppvalue,1,sum)
  return(pvalue)
}

#' 
#' @title Fit the mixture inverse-gamma prior of variance
#' @description Fit the mixture inverse-gamma prior of variance, given the variance estimates (sehat^2). 
#' 
#' @param sehat a p vector of observed standard errors
#' @param df appropriate degrees of freedom for (chi-square) distribution of sehat
#' @param betahat a p vector of estimates (optional)
#' @param randomstart logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param singlecomp logical, indicating whether to use a single inverse-gamma distribution as the prior distribution for the variances
#' @param unimodal put unimodal constraint on the prior distribution of variances ("variance") or precisions ("precision")
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param g the prior distribution for variances (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param maxiter maximum number of iterations of the EM algorithm
#' @param estpriormode logical, indicating whether to estimate the mode of the unimodal prior
#' @param priormode specified prior mode (only works when estpriormode=FALSE).
#' @param completeobs a p vector of non-missing flags
#'
#' @return The fitted mixture prior (g) and convergence info
#' 
#' @export
#' 
est_prior = function(sehat, df, betahat, randomstart, singlecomp, unimodal,
                     prior, g, maxiter, estpriormode, priormode, completeobs){
  
  # Set up initial parameters
  if(!is.null(g)){
    maxiter = 1 # if g is specified, don't iterate the EM
    estpriormode = FALSE
    prior = rep(1,length(g$pi)) #prior is not actually used if g specified, but required to make sure EM doesn't produce warning
    l = length(g$pi)
  } else {   
    # use moment matching to initialize IG shape and rate
    mm = est_singlecomp_mm(sehat[completeobs],df)
    mm$a = max(mm$a,1e-5)
    if(singlecomp==TRUE){
      alpha = mm$a
    }else{
      # let alpha be a grid of values ranging from small to large
      if(mm$a>1){
        alpha = 1+((64/mm$a)^(1/6))^seq(-3,6)*(mm$a-1)
      }else{
        alpha = mm$a*2^seq(0,13)
      }
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
  }else{
    if(!is.numeric(prior)){
      if(prior=="nullbiased"){ # set up prior to favour "null"
        prior = rep(1,l)
        prior[1] = l-1 #prior 10-1 in favour of null by default
      }else if(prior=="uniform"){
        prior = rep(1,l)
      }
    }
    if(length(prior)!=l | !is.numeric(prior)){
      stop("Error: invalid prior specification")
    }
  }
  
  if(randomstart){
    pi.se = rgamma(l,1,1)
  } else {   
    pi.se = rep(1,l)/l
  }
  pi.se = pi.se/sum(pi.se)
  
  if (!is.null(g)){
    g = igmix(g$pi,g$alpha,g$beta)
  }else{
    g = igmix(pi.se,alpha,beta)
  }
  
  # fit the mixture prior
  mix.fit = est_mixprop_mode(sehat[completeobs],g,prior,df,unimodal,singlecomp,estpriormode,maxiter)
  return(mix.fit)
}



