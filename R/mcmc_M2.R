#' ------------------------------------------------------------------------
#' irt.mcmc: 
#' ---------
#' 
#' ---   Return  : 
#' ----------------------------------------------------------------------------


M2.tmcmc <- function(Y, t_resp,
                             par1,
                             par2,
                             M,
                             L,
                             step.sizes) {
  
  #' Arguments
  #' ---------------------------------------------
  #' Y           : N x J data matrix
  #' 
  #' par1        : a list of parameters on the 1st level: theta, xi, eta, beta, ...
  #' par2        : a list of parameters on the 2nd level
  #' hyp_par     :
  #' nburnin     : 
  #' M           : the total No. of iterations
  #' step.sizes  :
  #' L           : No. of chains
  #' --------------------------------------------
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
  alpha_theta = hyp_par$alpha_theta
  beta_theta = hyp_par$beta_theta
  
  alpha_beta = hyp_par$alpha_beta
  beta_beta = hyp_par$beta_beta
  
  mu_b = hyp_par$mu_b
  s2_b = hyp_par$s2_b
  
  alpha_xi = hyp_par$alpha_xi
  beta_xi = hyp_par$beta_xi
  
  alpha_eta = hyp_par$alpha_eta
  beta_eta = hyp_par$beta_eta
  
  alpha_delta = hyp_par$alpha_delta
  beta_delta = hyp_par$beta_delta
  
  
  
  #' Initialisation
  #'---------------
  
  #'================================================
  
  
  pi_xi <- matrix(nrow = M, ncol = L)
  pi_xi[1,] <- runif(L, 0.01, 0.1) #rep(par2$pi_xi, L)
  
  pi_eta <- matrix(nrow = M, ncol = L)
  pi_eta[1,] <- runif(L, 0.35, 0.45)
  
  
  Sigma <- list()
  Sigma[[1]] <- replicate(L, matrix(c(runif(1, 0.1, 0.5), 0, 0, runif(L, 0.1, 0.5)), 
                                    nrow = 2, ncol = 2), simplify=FALSE)
  
  Omega <- list()
  Omega[[1]] <- replicate(L, matrix(c(runif(1, 0.1, 0.5), 0, 0, runif(L, 0.1, 0.5)), 
                                    nrow = 2, ncol = 2), simplify=FALSE)
  
  Mu <- array(dim = c(M, 2, L))
  Mu[1,,] <- rep(c(0,0), L)
  
  delta <- matrix(nrow = M, ncol = L)
  delta[1,] <- runif(L, 0.5, 1.5)
  
  gamma <- matrix(nrow = M, ncol = L)
  gamma[1,] <- runif(L, 0.5, 1.5)
  
  kappa <- matrix(nrow = M, ncol = L)
  kappa[1,] <- runif(L, 0.5, 2)
  
  #'================================================
  
  
  theta <- array(dim = c(M, N, L))
  theta[1,,] <- rep(rnorm(N, 0, 0.1), L)
  
  tau <- array(dim = c(M, N, L))
  tau[1,,] <- rep(rnorm(N, 0, 0.1), L)
  
  xi <- array(dim = c(M, N, L))
  xi[1,,] <- rep(rbinom(N, 1, pi_xi[1,]))
  
  
  beta <- array(dim = c(M, J, L))
  beta[1,,] <- rep(rnorm(J, 0, 0.1), L)
  
  alpha <- array(dim = c(M, J, L))
  alpha[1,,] <- rep(rnorm(J, 0, 0.1), L)
  
  eta <- array(dim = c(M, J, L))
  eta[1,,] <- rep(rbinom(J, 1, pi_eta[1,]))
  
  
  #' -------------------------------------------------------------------------
  #' log posterior density functions of theta, tau, beta, alpha, delta, gamma
  #' -------------------------------------------------------------------------
  
  f.theta <- function(theta, s2_theta, xi, eta, beta, delta){
    
    temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta #'N x J
    temp <- rowSums(Y * temp1 - log(1 + exp(temp1))) + log(dnorm(theta, 0, sqrt(s2_theta)))
    return(temp)
  }
  
  f.tau <- function(tau, s2_tau, xi, eta, alpha, gamma, kappa){
    
    temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
    temp <- rowSums(log(dlnorm(t_resp, temp2, sqrt(kappa)))) + log(dnorm(tau, 0, sqrt(s2_tau)))
    return(temp)
  }
  
  
  f.beta <- function(beta, mu_beta, s2_beta, xi, theta, eta, delta){
    
    temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
    temp <- colSums(Y * temp1 - log(1 + exp(temp1))) + log(dnorm(beta, mu_beta, sqrt(s2_beta)))
    return(temp)
  }
  
  f.alpha <- function(alpha, mu_alpha, s2_alpha, xi, tau, eta, gamma, kappa){
    
    temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
    temp <- colSums(log(dlnorm(t_resp, temp2, sqrt(kappa)))) + log(dnorm(tau, 0, sqrt(s2_tau)))
    return(temp)
  }
  
  f.delta <- function(delta, alpha_delta, beta_delta, xi, theta, eta, beta){
    
    temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
    temp <- sum(Y * temp1 - log(1 + exp(temp1))) + log(dgamma(delta, alpha_delta, scale = beta_delta))
    return(temp)
  }
  
  f.gamma <- function(gamma, alpha_gamma, beta_gamma, xi, tau, eta, alpha){
    
    temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
    temp <- sum(log(dlnorm(t_resp, temp2, sqrt(kappa)))) + log(dgamma(gamma, alpha_gamma, scale = beta_gamma))
  }
  
  f.kappa <- function(kappa, alpha_kappa, beta_kappa, xi, tau, eta, alpha){
    
    temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
    temp <- sum(log(dlnorm(t_resp, temp2, sqrt(kappa)))) 
    + log(dgamma(kappa, alpha_kappa, scale = beta_kappa))
    
  }
  
  
  #' -------------------------------------
  #'  Posterior samplers for pi_xi, pi_eta
  #' -------------------------------------
  
  f.pi_xi <- function(pi_xi, xi, alpha_xi, beta_xi) {
    log(dbeta(pi_xi, alpha_xi, beta_xi)) + sum(xi*log(pi_xi) + (1-xi)*log(1-pi_xi))
  }
  
  f.pi_eta <- function(pi_eta, eta, alpha_eta, beta_eta) {
    log(dbeta(pi_eta, alpha_eta, beta_eta)) + sum(eta*log(pi_eta) + (1-eta)*log(1-pi_eta))
  }
  
  
  #'---------------
  #' Main loop
  #'---------------
  
  for (m in 2:M) {
    
    swap1_index <- sample(1:(L-1), L-1)
    
    ## Updata LVs
    for (l in 1:L){
      
      
      theta[m,,l] <- MH_gaussian(log.f = f.theta,
                                 p = N,
                                 curr = theta[m-1,,l],
                                 step.sizes = step.sizes[1],
                                 s2_theta = Sigma[[m-1]][[l]][1,], xi = xi[m-1,,l], 
                                 eta = eta[m-1,,l], beta = beta[m-1,,l], 
                                 delta = delta[m-1,l])
      
      tau[m,,l] <- MH_gaussian(log.f = f.tau,
                               p = N,
                               curr = tau[m-1,,l],
                               step.sizes = step.sizes[1],
                               s2_tau = Sigma[[m-1]][[l]][,2], xi = xi[m-1,,l], 
                               eta = eta[m-1,,l], alpha = alpha[m-1,,l], 
                               gamma = gamma[m-1,l], kappa = kappa[m-1,l])
      
      xi[m,,l] <- sample.xi(xi[m-1,,l], pi_xi[m-1,l], theta[m,,l], 
                            eta[m-1,,l], beta[m-1,,l], delta[m-1,l])
      
      
      beta[m,,l] <- MH_gaussian(log.f = f.beta,
                                p = J,
                                curr = beta[m-1,,l],
                                step.sizes = step.sizes[2],
                                mu_beta = Mu[m-1, 1, l], 
                                s2_beta = Omega[[m-1]][[l]][1,], 
                                xi = xi[m,,l], theta = theta[m,,l], 
                                eta = eta[m-1,,l], delta = delta[m-1,l])
      
      alpha[m,,l] <- MH_gaussian(log.f = f.alpha,
                                 p = J,
                                 curr = alpha[m-1,,l],
                                 step.sizes = step.sizes[2],
                                 mu_alpha = Mu[m-1, 2, l], 
                                 s2_alpha = Omega[[m-1]][[l]][,2], 
                                 xi = xi[m,,l], tau = tau[m,,l], eta = eta[m-1,,l], 
                                 gamma = gamma[m-1,l], kappa = kappa[m-1,l])
      
      eta[m,,l] <- sample.eta(eta[m-1,,l], pi_eta[m-1,l], xi[m,,l], theta[m,,l], 
                              beta[m,,l], delta[m-1,l])
      
      delta[m,l] <- MH_gaussian(log.f = f.delta,
                                p = 1,
                                curr = delta[m-1,l],
                                step.sizes = step.sizes[3],
                                alpha_delta = alpha_delta, beta_delta = beta_delta, 
                                xi = xi[m,,l], theta = theta[m,,l], 
                                eta = eta[m,,l], beta = beta[m,,l])
      
      gamma[m,l] <- MH_gaussian(log.f = f.gamma,
                                p = 1,
                                curr = gamma[m-1,l],
                                step.sizes = step.sizes[3],
                                alpha_gamma = alpha_gamma, beta_gamma = beta_gamma, 
                                xi = xi[m,,l], tau = tau[m,,l], 
                                eta = eta[m,,l], alpha = alpha[m,,l])
      
      kappa[m,l] <- MH_gaussian(log.f = f.kappa,
                                p = 1,
                                curr = kappa[m-1,l],
                                step.sizes = step.sizes[3],
                                alpha_kappa = alpha_kappa, beta_kappa = beta_kappa, 
                                xi = xi[m,,l], tau = tau[m,,l], 
                                eta = eta[m,,l], alpha = alpha[m,,l])
      
      
      ## Update Pars
      
      pi_xi[m,l] <- rbeta(1, (alpha_xi + sum(xi[m,,l])), (beta_xi + sum(xi[m,,l]==0)))
      
      pi_eta[m,l] <- rbeta(1, (alpha_eta + sum(eta[m,,l])), (beta_eta + sum(eta[m,,l]==0)))
      
      Sigma[[m]][[l]] <- MCMCpack::riwish(v = N + 3,
                                          S = matrix(c(2,0,0,2), ncol = 2) 
                                          + matrix(c(sum(theta[m,,l]^2), 
                                                     sum(theta[m,,l]*tau[m,,l]), 
                                                     sum(theta[m,,l]*tau[m,,l]), 
                                                     sum(tau[m,,l]^2)), ncol = 2))
      
      Mu[m,,l] <- solve(matrix(c(1/25, 0, 0, 1/25), ncol = 2) + solve(J*Omega[[m-1]][[l]]))*
        solve(Omega[[m-1]][[l]])*
        
        Omega[[m]][[l]] <- MCMCpack::riwish(v = J + 3,
                                            S = matrix(c(2,0,0,2), ncol = 2) 
                                            + matrix(c(sum((beta[m,,l]-Mu[m,1,l])^2), 
                                                       sum((beta[m,,l]-Mu[m,1,l])*(alpha[m,,l]-Mu[m,2,l])), 
                                                       sum((beta[m,,l]-Mu[m,1,l])*(alpha[m,,l]-Mu[m,2,l])), 
                                                       sum((alpha[m,,l]-Mu[m,2,l])^2)), ncol = 2))
      
      
    } #' end of loop l
    
    
    ## start the coupling update
    
    swap1 <- c(swap1_index[l-1], swap1_index[l-1] + 1)   
    
    log_A1 <- f.pi_xi(pi_xi[m,swap1[1]], xi[m,,swap1[1]], 
                      alpha_xi = alpha_xi, beta_xi = beta_xi) + 
      f.pi_xi(pi_xi[m,swap1[2]], xi[m,,swap1[2]], alpha_xi = alpha_xi, beta_xi = beta_xi) -
      f.pi_xi(pi_xi[m,swap1[2]], xi[m,,swap1[1]], alpha_xi = alpha_xi, beta_xi = beta_xi) - 
      f.pi_xi(pi_xi[m,swap1[1]], xi[m,,swap1[2]], alpha_xi = alpha_xi, beta_xi = beta_xi) 
    
    if (log(runif(1)) < log_A1) {
      pi_xi[m,swap1] = rev(pi_xi[m,swap1])
      xi[m,,swap1] = rev(xi[m,,swap1])
    }
    
    
    swap2 <- swap1 
    
    log_A2 <- f.pi_eta(pi_eta[m,swap2[1]], eta[m,,swap2[1]], 
                       alpha_eta = alpha_eta, beta_eta = beta_eta) + 
      f.pi_eta(pi_eta[m,swap2[2]], eta[m,,swap2[2]], alpha_eta = alpha_eta, beta_eta = beta_eta) -
      f.pi_eta(pi_eta[m,swap2[2]], eta[m,,swap2[1]], alpha_eta = alpha_eta, beta_eta = beta_eta) - 
      f.pi_eta(pi_eta[m,swap2[1]], eta[m,,swap2[2]], alpha_eta = alpha_eta, beta_eta = beta_eta) 
    
    
    if (log(runif(1)) < log_A2) {
      pi_eta[m,swap2] = rev(pi_eta[m,swap2])
      eta[m,,swap2] = rev(eta[m,,swap2])
    }
    
    ## end of the coupling update
    
  } #' end of loop m
  
  
  result <- list()
  
  result$par1 <- list(xi = xi, theta = theta, tau = tau,
                      eta = eta, beta = beta, alpha = alpha)
  result$par2 <- list(delta = delta, gamma = gamma, kappa = kappa,
                      pi_xi = pi_xi, Sigma = Sigma,
                      pi_eta = pi_eta, Mu = Mu, Omega = Omega)
  result$M <- M
  result$nburnin <- nburnin  
  result$data <- Y
  
  return(result)
}