#' ------------------------------------------------------------------------------
#' MH_gaussian:  Metropolis sampling with a Gaussian proposal
#' ------------------------------------------------------------------------------



MH_gaussian <- function(log.f, 
                        p = 1, 
                        curr, 
                        iters = 1,
                        step.sizes,...){
  
  #' ----------
  #' Arguments
  #' ----------
  #'  
  #' log.f            : the log of the density function
  #' 
  #' p                : the number of variables being sampled
  #'
  #' curr             : the current value obtained
  #'  
  #' iters            : for each transition, 
  #'                    the Markov chain sampling is run 'iters_per.iter' times
  #'                    
  #' step.sizes       : the standard deviations of Gaussian proposal
  #'                    stepsizes[i] for the `i'th  variable 
  #'                  
  #' ...              : other arguments passed by log_f
  #' ------------------------------------------------------------------------
  
  ## Check the dimension
  
  if(p != length(curr))
    stop("The number of variables in initial values doesn't match the dimension of the samples")
  
  # if(!is.finite(log.f(curr,...))) 
  #   stop("The log of the density is zero given initial values") 
  
  ## Evaluate the log density at the current value
  
  y <- rep(0, p)
  
  old_log.f <- log.f(curr,...)
  
  
  ## One-step transition
  
  y <- curr 
  
  iters.inside <- 0
  
  repeat{
    
    iters.inside <- iters.inside + 1
    
    # 1) propose a value
    y.prop <- rnorm(p) * step.sizes + y
    
    # 2) accept or reject
    new_log.f <- log.f(y.prop,...)
    
    accept.prob <- exp(new_log.f - old_log.f)
    #' acceptance probability
    #alpha <- min(accept.prob, 1)
    # jumping with probability alpha
    
    alpha_bin <- (runif(p) <= accept.prob)
    
    y[alpha_bin] <- y.prop[alpha_bin]
    
    old_log.f[alpha_bin] <- new_log.f[alpha_bin]
    
    if (iters.inside == iters) break
  }
  
  ## Return the value 
  return(y)
}


#' ------------------------------------------------------------------------------
#' MH_weight:  Metropolis sampling with a Gaussian proposal
#' ------------------------------------------------------------------------------



MH_weight <- function(log.f, 
                      p = 1, 
                      curr, 
                      iters = 1,
                      step.sizes,...){
  
  #' ----------
  #' Arguments
  #' ----------
  #'  
  #' log.f            : the log of the density function
  #' 
  #' p                : the number of variables being sampled
  #'
  #' curr             : the current value obtained
  #'  
  #' iters            : for each transition, 
  #'                    the Markov chain sampling is run 'iters_per.iter' times
  #'                    
  #' step.sizes       : the standard deviations of Gaussian proposal
  #'                    stepsizes[i] for the `i'th  variable 
  #'                  
  #' ...              : other arguments passed by log_f
  #' ------------------------------------------------------------------------
  
  ## Check the dimension
  
  if(p != length(curr))
    stop("The number of variables in initial values doesn't match the dimension of the samples")
  
  # if(!is.finite(log.f(curr,...))) 
  #   stop("The log of the density is zero given initial values") 
  
  ## Evaluate the log density at the current value
  
  w <- rep(0, p)
  
  old_log.f <- log.f(curr,...)
  
  
  ## One-step transition
  
  w <- curr 
  
  iters.inside <- 0
  
  repeat{
    
    iters.inside <- iters.inside + 1
    
    # 1) propose a value
    #w.prop <- exp(rnorm(p) * step.sizes + log(w))
    repeat{
      w.prop <- rnorm(p) * step.sizes + w
      if ((w.prop>0) & (w.prop<1)) break
    }
    
    # 2) accept or reject
    new_log.f <- log.f(w.prop,...)
    
    accept.prob <- exp(new_log.f - old_log.f)
    #' acceptance probability
    #alpha <- min(accept.prob, 1)
    # jumping with probability alpha
    
    alpha_bin <- (runif(p) <= accept.prob)
    
    w[alpha_bin] <- w.prop[alpha_bin]
    
    old_log.f[alpha_bin] <- new_log.f[alpha_bin]
    
    if (iters.inside == iters) break
  }
  
  ## Return the value 
  #pi <- w/sum(w)
  return(w)
}




n = 10000
log_q = function(theta, y=3, n=10) {
  if (theta<0 | theta>1) return(-Inf)
  (y-0.5)*log(theta)+(n-y-0.5)*log(1-theta)
}
current = 0.5 # Initial value
samps = rep(NA,n)
for (i in 1:n) {
  proposed = rnorm(1, current, 0.4) # tuning parameter is 0.4
  logr = log_q(proposed)-log_q(current)
  if (log(runif(1)) < logr) current = proposed
  samps[i] = current
}
length(unique(samps))/n # acceptance rate



