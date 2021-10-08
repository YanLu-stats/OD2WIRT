
#' ------------------------------------------------------------------------------
#' simdata_M2: Data Generation from The Full Model M2
#' ------------------------------------------------------------------------------


simData_M2 <- function(#myseed, 
  N, J,
  pi_xi, pi_eta,
  s2_theta, cov_theta_tau, s2_tau,
  s2_beta, cov_beta_alpha, s2_alpha,
  mu_beta, mu_alpha,
  delta, gamma, kappa){
  #'--- Person ----------------------------
  #set.seed(myseed)
  
  Sigma <- matrix(c(s2_theta, cov_theta_tau, cov_theta_tau, s2_tau), 2, 2)
  
  Theta <- MASS::mvrnorm(N, c(0,0), Sigma)
  theta <- Theta[,1]
  tau <- Theta[,2]
  
  xi.rand.samples <- rbinom(N, 1, pi_xi)
  xi <- sort(xi.rand.samples)
  
  #'--- Item ----------------------------
  
  Mu <- c(mu_beta, mu_alpha)
  Omega <- matrix(c(s2_beta, cov_beta_alpha, cov_beta_alpha, s2_alpha), 2, 2)
  
  Beta <- MASS::mvrnorm(J, Mu, Omega)
  beta <- Beta[,1]
  alpha <- Beta[,2]
  
  eta.rand.samples <- rbinom(J, 1, pi_eta)
  eta <- sort(eta.rand.samples)
  
  #'--- Global --------------------------
  
  delta <- runif(1, 0.5, 1) #rgamma(1, 2, 2)
  gamma <- runif(1, 0.5, 1)
  kappa <- runif(1, 0.4, 0.8)
  
  #'-------------------------------------------------------------------------------------
  logits <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  probs <- 1 / (1 + exp(-logits))
  resp.prob <- matrix(probs, ncol = J)
  obs.resp <- matrix(sapply(c(resp.prob), rbinom, n = 1, size = 1), ncol = J)
  
  mu_t <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma
  time.resp <- matrix(sapply(c(mu_t), rlnorm, n = 1, sdlog = sqrt(kappa)), ncol = J)
  #'-------------------------------------------------------------------------------------
  
  output        <- list()
  
  output$resp.prob <-  resp.prob # List 1
  #' response probs
  output$mean.log.time <- mu_t   # List 2
  
  output$resp   <- obs.resp      # List 3
  output$time <- time.resp       # List 4
  #' data
  
  output$par1 <- list(theta = theta, tau = tau, xi = xi, 
                      eta = eta, beta = beta, alpha = alpha) # List 5
  
  output$par2 <- list(pi_xi = pi_xi, pi_eta = pi_eta,
                      s2_theta = s2_theta, cov_theta_tau = cov_theta_tau, s2_tau = s2_tau,
                      s2_beta = s2_beta, cov_beta_alpha = cov_beta_alpha, s2_alpha = s2_alpha,
                      mu_beta = mu_beta, mu_alpha = mu_alpha,
                      delta = delta, gamma = gamma, kappa = kappa) # List 6
  #output$seed <- myseed
  return(output)
}

