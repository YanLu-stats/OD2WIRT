#' ------------------------------------------------------------------------
#' irt.mcmc: 
#' ---------
#' 
#' ---   Return  : 
#' ----------------------------------------------------------------------------


M1.mcmc <- function(Y, 
                    par1,
                    par2,
                    M,
                    L,
                    temps = 2^(0:3),
                    step.sizes,
                    seed = 123) {
  
  #' Arguments
  #' ---------------------------------------------
  #' Y           : N x J data matrix
  #' par1        : a list of parameters on the 1st level: theta, xi, eta, beta, ...
  #' par2        : a list of parameters on the 2nd level
  #' M           : the total No. of iterations
  #' L           : replicates of posterior samples for calculating mDIC
  #' temps       : temperature ladders for tempered MCMC
  #' step.sizes  :
  #' --------------------------------------------
  
  N <- nrow(Y)
  J <- ncol(Y)
  
  
  #' Initialisation
  #'---------------
  
  pi_xi <- c()
  pi_xi[1] <- par2$pi_xi
  
  s2_theta <- c()
  s2_theta[1] <- par2$s2_theta
  
  
  pi_eta <- c()
  pi_eta[1] <- par2$pi_eta
  
  mu_beta <- c()
  mu_beta[1] <- par2$mu_beta
  
  s2_beta <- c()
  s2_beta[1] <- par2$s2_beta
  
  
  delta <- c()
  delta[1] <- par2$delta
  
  
  #'================================================
  
  theta <- matrix(nrow = M, ncol = N)
  theta[1,] <- par1$theta
  
  xi <- matrix(nrow = M, ncol = N)
  xi[1,] <- par1$xi
  
  
  beta <- matrix(nrow = M, ncol = J)
  beta[1,] <- par1$beta
  
  eta <- matrix(nrow = M, ncol = J)
  eta[1,] <- par1$eta
  
  
  #'================================================  
  
  theta_rep <- list()
  theta_rep[[1]] <- rnorm(L, 0, sqrt(s2_theta[1]))
  
  beta_rep <- list()
  beta_rep[[1]] <- rnorm(L, mu_beta[1], sqrt(s2_beta[1]))
  
  xi_rep <- list()
  xi_rep[[1]] <- rbinom(L, 1, pi_xi[1])
  
  eta_rep <- list()
  eta_rep[[1]] <- rbinom(L, 1, pi_eta[1])
  
  
  #'---------------
  #' Main loop
  #'---------------
  set.seed(seed)
  
  for(m in 2:M) {
    
    ## Updata LVs
    
    theta[m,] <- MH_gaussian(log.f = f.theta,
                             p = N,
                             curr = theta[m-1,],
                             step.sizes = step.sizes[1],
                             s2_theta = s2_theta[m-1], xi = xi[m-1,], eta = eta[m-1,], 
                             beta = beta[m-1,], delta = delta[m-1], N = N, J = J, Y = Y)
    
    xi[m,] <- sample.xi(xi[m-1,], pi_xi[m-1], theta[m,], eta[m-1,], beta[m-1,], delta[m-1], N = N, J = J, Y = Y)
    
    beta[m,] <- MH_gaussian(log.f = f.beta,
                            p = J,
                            curr = beta[m-1,],
                            step.sizes = step.sizes[1],
                            mu_beta = mu_beta[m-1], s2_beta = s2_beta[m-1], xi = xi[m,],
                            theta = theta[m,], eta = eta[m-1,], delta = delta[m-1], N = N, J = J, Y = Y)
    
    eta[m,] <- sample.eta(eta[m-1,], pi_eta[m-1], xi[m,], theta[m,], beta[m,], delta[m-1], N = N, J = J, Y = Y)
    
    
    delta[m] <- MH_gaussian(log.f = f.delta,
                            p = 1,
                            curr = delta[m-1],
                            step.sizes = step.sizes[3],
                            sigma_delta = 2.5, 
                            xi = xi[m,], theta = theta[m,], eta = eta[m,], beta = beta[m,], N = N, J = J, Y = Y)
    
    
    pi_xi[m] <- MH_weight(log.f = f.pi_xi,
                          p = 1,
                          curr = pi_xi[m-1],
                          step.size = step.sizes[3],
                          xi = xi[m,], alpha_xi = 2, beta_xi = 2)
    
    s2_theta[m] <- MH_gaussian(log.f = f.s2_theta,
                               p = 1,
                               curr = s2_theta[m-1],
                               step.sizes = step.sizes[3],
                               alpha_theta = 1/2, beta_theta = 1,
                               theta = theta[m,], N = N)
    
    pi_eta[m] <- MH_weight(log.f = f.pi_eta,
                           p = 1,
                           curr = pi_eta[m-1],
                           step.size = step.sizes[3],
                           eta = eta[m,], alpha_eta = 2, beta_eta = 2)
    
    mu_beta[m] <- MH_gaussian(log.f = f.mu_beta,
                              p = 1,
                              curr = mu_beta[m-1],
                              step.sizes = step.sizes[3],
                              mu_b = 0, s2_b = 25,
                              beta = beta[m,], s2_beta = s2_beta[m-1], J = J)
    
    s2_beta[m] <- MH_gaussian(log.f = f.s2_beta,
                              p = 1,
                              curr = s2_beta[m-1],
                              step.sizes = step.sizes[3],
                              alpha_beta = 1/2, beta_beta = 1,
                              beta = beta[m,], mu_beta = mu_beta[m], J = J)
    
    ## Replicated LVs
    
    theta_rep[[m]] <- rnorm(L, 0, sqrt(s2_theta[m]))
    
    beta_rep[[m]] <- rnorm(L, mu_beta[m], sqrt(s2_beta[m]))
    
    xi_rep[[m]] <- rbinom(L, 1, pi_xi[m])
    
    eta_rep[[m]] <- rbinom(L, 1, pi_eta[m])
    
    
  } #' end of loop m
  
  
  result <- list()
  
  result$par1 <- list(xi = xi, theta = theta, eta = eta, beta = beta)
  result$par2 <- list(delta = delta, 
                      pi_xi = pi_xi, s2_theta = s2_theta, 
                      pi_eta = pi_eta, mu_beta = mu_beta, s2_beta = s2_beta)
  
  result$rep <- list(theta = theta_rep, beta = beta_rep, xi = xi_rep, eta = eta_rep)
  
  result$M <- M
  result$L <- L
  result$data <- Y
  
  return(result)
}
