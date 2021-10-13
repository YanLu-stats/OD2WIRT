
#' -----------------------------------------------------------------------
#' A Parallel-Tempered MCMC Algorithm for Estimating the Full Model
#' -----------------------------------------------------------------------
#' 
#' ---   Return  : Posterior samples etc.
#' -----------------------------------------------------------------------


M2.tmcmc <- function(Y, t_resp,
                             par1,
                             par2,
                             M,
                             L,
                             K,
                             step.sizes) {
  
  #' Arguments
  #' ---------------------------------------------
  #' Y           : N x J data matrix
  #' t_resp      : N x J response time matrix
  #' par1        : a list of parameters on the 1st level: latent variables
  #' par2        : a list of parameters on the 2nd level
  #' nburnin     : 
  #' M           : the total No. of iterations
  #' L           : for calculating mDIC
  #' K           : No. of parallel chains 
  #' step.sizes  :
  #' --------------------------------------------
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
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
                                 delta = delta[m-1,l], 
                                 N = N, J = J, Y = Y)
      
      tau[m,,l] <- MH_gaussian(log.f = f.tau,
                               p = N,
                               curr = tau[m-1,,l],
                               step.sizes = step.sizes[1],
                               s2_tau = Sigma[[m-1]][[l]][,2], xi = xi[m-1,,l], 
                               eta = eta[m-1,,l], alpha = alpha[m-1,,l], 
                               gamma = gamma[m-1,l], kappa = kappa[m-1,l], 
                               N = N, J = J, t_resp = t_resp)
      
      xi[m,,l] <- sample.xi(xi = xi[m-1,,l], pi_xi = pi_xi[m-1,l], theta = theta[m,,l], 
                            eta = eta[m-1,,l], beta = beta[m-1,,l], delta = delta[m-1,l], 
                            tau = tau[m-1,,l], alpha = alpha[m-1,,l], 
                            gamma = gammadelta[m-1,l], kappa = kappa[m-1,l],
                            N = N, J = J, Y = Y, t_resp = t_resp)
      
      beta[m,,l] <- MH_gaussian(log.f = f.beta,
                                p = J,
                                curr = beta[m-1,,l],
                                step.sizes = step.sizes[2],
                                mu_beta = Mu[m-1, 1, l], 
                                s2_beta = Omega[[m-1]][[l]][1,], 
                                xi = xi[m,,l], theta = theta[m,,l], 
                                eta = eta[m-1,,l], delta = delta[m-1,l],
                                N = N, J = J, Y = Y)
      
      alpha[m,,l] <- MH_gaussian(log.f = f.alpha,
                                 p = J,
                                 curr = alpha[m-1,,l],
                                 step.sizes = step.sizes[2],
                                 mu_alpha = Mu[m-1, 2, l], 
                                 s2_alpha = Omega[[m-1]][[l]][,2], 
                                 xi = xi[m,,l], tau = tau[m,,l], eta = eta[m-1,,l], 
                                 gamma = gamma[m-1,l], kappa = kappa[m-1,l],
                                 N = N, J = J, t_resp = t_resp)
      
      eta[m,,l] <- sample.eta(eta = eta[m-1,,l], pi_eta = pi_eta[m-1,l], theta = theta[m,,l], 
                              xi = xi[m,,l], beta = beta[m,,l], delta = delta[m-1,l], 
                              tau = tau[m,,l], alpha = alpha[m,,l], 
                              gamma = gammadelta[m-1,l], kappa = kappa[m-1,l],
                              N = N, J = J, Y = Y, t_resp = t_resp)
      
      delta[m,l] <- MH_gaussian(log.f = f.delta,
                                p = 1,
                                curr = delta[m-1,l],
                                step.sizes = step.sizes[3],
                                sigma_delta = 2.5, 
                                xi = xi[m,,l], theta = theta[m,,l], 
                                eta = eta[m,,l], beta = beta[m,,l],
                                N = N, J = J, Y = Y)
      
      gamma[m,l] <- MH_gaussian(log.f = f.gamma,
                                p = 1,
                                curr = gamma[m-1,l],
                                step.sizes = step.sizes[3],
                                sigma_gamma = 2.5,
                                xi = xi[m,,l], tau = tau[m,,l], 
                                eta = eta[m,,l], alpha = alpha[m,,l],
                                N = N, J = J, t_resp = t_resp)
      
      kappa[m,l] <- MH_gaussian(log.f = f.kappa,
                                p = 1,
                                curr = kappa[m-1,l],
                                step.sizes = step.sizes[3],
                                alpha_kappa = 0.001, beta_kappa = 0.001, 
                                xi = xi[m,,l], tau = tau[m,,l], 
                                eta = eta[m,,l], alpha = alpha[m,,l],
                                N = N, J = J, t_resp = t_resp)
      
      
      pi_xi[m,l] <- MH_weight(log.f = f.pi_xi,
                              p = 1,
                              curr = pi_xi[m-1],
                              step.size = step.sizes[3],
                              xi = xi[m,,l], alpha_xi = 2, beta_xi = 2)
      
      
      Sigma[[m]][[l]] <- MH_gaussian(log.f = f.Sigma,
                                 p = 1,
                                 curr = Sigma[[m-1]][[l]],
                                 step.sizes = step.sizes[3],
                                 theta = theta[m,,l], tau = tau[m,,l],
                                 scale_Sigma = matrix(c(2,0,0,2), ncol = 2), 
                                 df_Sigma = 2, N = N)
      
      
      pi_eta[m,l] <- MH_weight(log.f = f.pi_eta,
                               p = 1,
                               curr = pi_eta[m-1],
                               step.size = step.sizes[3],
                               eta = eta[m,,l], alpha_eta = 2, beta_eta = 2)
      
      Mu[m,,l] <- MH_gaussian(log.f = f.mu_beta,
                                p = 1,
                                curr = Mu[m-1,,l],
                                step.sizes = step.sizes[3],
                                beta = beta[m,,l], alpha = alpha[m,,l],
                                Omega = Omega[[m-1]][[l]],
                                J = J)
      
      Omega[[m]][[l]] <- MH_gaussian(log.f = f.Omega,
                                     p = 1,
                                     curr = Omega[[m-1]][[l]],
                                     step.sizes = step.sizes[3],
                                     Mu = Mu[m,,l],
                                     beta = beta[m,,l], alpha = alpha[m,,l],
                                     scale_Omega = matrix(c(2,0,0,2), ncol = 2), 
                                     df_Omega = 2, J = J)
      

    } #' end of loop l
    
    
    ## coupling update
    
    swap1 <- c(swap1_index[l-1], swap1_index[l-1] + 1)   
    
    log_A1 <- f.pi_xi(pi_xi[m,swap1[1]], xi[m,,swap1[1]], 
                      alpha_xi = 2, beta_xi = 2) + 
      f.pi_xi(pi_xi[m,swap1[2]], xi[m,,swap1[2]], alpha_xi = 2, beta_xi = 2) -
      f.pi_xi(pi_xi[m,swap1[2]], xi[m,,swap1[1]], alpha_xi = 2, beta_xi = 2) - 
      f.pi_xi(pi_xi[m,swap1[1]], xi[m,,swap1[2]], alpha_xi = 2, beta_xi = 2) 
    
    if (log(runif(1)) < log_A1) {
      pi_xi[m,swap1] = rev(pi_xi[m,swap1])
      xi[m,,swap1] = rev(xi[m,,swap1])
    }
    
    
    swap2 <- swap1 
    
    log_A2 <- f.pi_eta(pi_eta[m,swap2[1]], eta[m,,swap2[1]], 
                       alpha_eta = 2, beta_eta = 2) + 
      f.pi_eta(pi_eta[m,swap2[2]], eta[m,,swap2[2]], alpha_eta = 2, beta_eta = 2) -
      f.pi_eta(pi_eta[m,swap2[2]], eta[m,,swap2[1]], alpha_eta = 2, beta_eta = 2) - 
      f.pi_eta(pi_eta[m,swap2[1]], eta[m,,swap2[2]], alpha_eta = 2, beta_eta = 2) 
    
    
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
