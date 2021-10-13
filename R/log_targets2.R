#' --------------------------------------------------------
#' log posterior density functions for the Full Model M2
#' --------------------------------------------------------

#' ----------------------------------------------------------
#' log posterior density functions of 1st level parameters
#' ----------------------------------------------------------


f.theta <- function(theta, s2_theta, xi, eta, beta, delta, N, J, Y){
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta #'N x J
  temp <- rowSums(Y * temp1 - log(1 + exp(temp1))) + log(dnorm(theta, 0, sqrt(s2_theta)))
  return(temp)
}


f.tau <- function(tau, s2_tau, xi, eta, alpha, gamma, kappa, N, J, t_resp){
  temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
  temp <- rowSums(log(dlnorm(t_resp, temp2, sqrt(kappa)))) + log(dnorm(tau, 0, sqrt(s2_tau)))
  return(temp)
}



f.beta <- function(beta, mu_beta, s2_beta, xi, theta, eta, delta, N, J, Y){
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  temp <- colSums(Y * temp1 - log(1 + exp(temp1))) + log(dnorm(beta, mu_beta, sqrt(s2_beta)))
  return(temp)
}



f.alpha <- function(alpha, mu_alpha, s2_alpha, xi, tau, eta, gamma, kappa, N, J, t_resp){
  temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
  temp <- colSums(log(dlnorm(t_resp, temp2, sqrt(kappa)))) + log(dnorm(tau, 0, sqrt(s2_tau)))
  return(temp)
}



f.xi <- function(xi, pi_xi, theta, eta, beta, delta, tau, alpha, gamma, kappa, N, J, Y, t_resp){
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta #'N x J
  temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma
  temp <- rowSums(Y * temp1 - log(1 + exp(temp1))) + xi*log(pi_xi)+(1-xi)*log(1-pi_xi) 
  + rowSums(log(dlnorm(t_resp, temp2, sqrt(kappa))))
  return(temp)
}

f.eta <- function(eta, pi_eta, theta, xi, beta, delta, tau, alpha, gamma, kappa, N, J, Y, t_resp){
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma
  temp <- colSums(Y * temp1 - log(1 + exp(temp1))) + eta*log(pi_eta)+(1-eta)*log(1-pi_eta) +
    colSums(log(dlnorm(t_resp, temp2, sqrt(kappa))))
  return(temp)
}




#' ---------------------------------------------------------
#'  log posterior density functions for 2nd level parameters
#' --------------------------------------------------------- 

f.delta <- function(delta, sigma_delta, xi, theta, eta, beta, N, J, Y){
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  temp <- sum(Y * temp1 - log(1 + exp(temp1))) + log(extraDistr::dhcauchy(delta, sigma_delta))
  return(temp)
}

f.gamma <- function(gamma, sigma_gamma, xi, tau, eta, alpha, N, J, t_resp){
  temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
  temp <- sum(log(dlnorm(t_resp, temp2, sqrt(kappa)))) + log(extraDistr::dhcauchy(gamma, sigma_gamma))
}

f.kappa <- function(kappa, alpha_kappa, beta_kappa, xi, tau, eta, alpha, N, J, t_resp){
  temp2 <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(eta) * gamma #'N x J
  temp <- sum(log(dlnorm(t_resp, temp2, sqrt(kappa)))) 
  + log(nimble::dinvgamma(kappa, alpha_kappa, scale = beta_kappa))
}


f.pi_xi <- function(pi_xi, xi, alpha_xi, beta_xi) {
  log(dbeta(pi_xi, alpha_xi, beta_xi)) + sum(xi*log(pi_xi) + (1-xi)*log(1-pi_xi))
}

f.pi_eta <- function(pi_eta, eta, alpha_eta, beta_eta) {
  log(dbeta(pi_eta, alpha_eta, beta_eta)) + sum(eta*log(pi_eta) + (1-eta)*log(1-pi_eta))
}

f.Sigma <- function(Sigma, theta, tau, scale_Sigma, df_Sigma){
  log (MCMCpack::riwish(v = df_Sigma, S = scale_Sigma))
  + sum(log(mvtnorm::dmvnorm(c(theta,tau), c(0,0), Sigma)))
}

f.Omega <- function(Omega, beta, alpha, Mu, scale_Omega, df_Omega){
  log (MCMCpack::riwish(v = df_Omega, S = scale_Omega))
  + sum(log(mvtnorm::dmvnorm(c(beta,alpha), Mu, Omega)))
}

f.Mu <- function(Mu, beta, alpha, Omega){
  log(mvtnorm::dmvnorm(Mu, c(0,0), scale_Mu)) 
  + sum(log(mvtnorm::dmvnorm(c(beta,alpha), Mu, Omega)))

}


f.xi <- function(xi, pi_xi, theta, eta, beta, delta, N, J, Y,
                 tau, alpha, gamma, kappa, t_resp){
  
  temp <- matrix(0, N, 2)
  
  temp1a <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta)
  temp1b <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J))
  temp[,1] <- rowSums(Y * temp1a - log(1 + exp(temp1a))) + log(1-pi_xi)
              + rowSums(log(dlnorm(t_resp, temp1b, sqrt(kappa))))
  
  temp2a <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + rep(1, N) %*% t(eta)*delta
  temp2b <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - rep(1, N) %*% t(eta)*gamma
  temp[,2] <- rowSums(Y * temp2a - log(1 + exp(temp2a))) + log(pi_xi)
              + rowSums(log(dlnorm(t_resp, temp2b, sqrt(kappa))))
  
  temp <- temp - apply(temp, 1, max)
  temp <- exp(temp)
  temp.sum <- rowSums(temp)
  prob <- temp/temp.sum
  return(prob)
}


sample.xi <- function(xi, pi_xi, theta, eta, beta, delta, N, J, Y, tau, alpha, gamma, kappa, t_resp){
  prob_xi <- f.xi(xi, pi_xi, theta, eta, beta, delta, N, J, Y, tau, alpha, gamma, kappa, t_resp)
  post <- rbinom(N, 1, prob_xi[,2])
  return(post)
}



f.eta <- function(eta, pi_eta, xi, theta, beta, delta, N, J, Y,
                  tau, alpha, gamma, kappa, N, J, Y, t_resp){
  
  temp <- matrix(0, J, 2)
  
  temp1a <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta)
  temp1b <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J))
  temp[,1] <- colSums(Y * temp1a - log(1 + exp(temp1a))) + log(1-pi_eta) 
  + colSums(log(dlnorm(t_resp, temp1b, sqrt(kappa))))
  
  temp2a <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(rep(1, J))*delta
  temp2b <- rep(1, N) %*% t(alpha) - tau %*% t(rep(1, J)) - xi %*% t(rep(1, J))*gamma
  temp[,2] <- colSums(Y * temp2a - log(1 + exp(temp2a))) + log(pi_eta) 
  + colSums(log(dlnorm(t_resp, temp2b, sqrt(kappa))))
  
  temp <- temp - apply(temp, 1, max)
  temp <- exp(temp)
  temp.sum <- rowSums(temp)
  prob <- temp/temp.sum
  return(prob)
}


sample.eta <- function(eta, pi_eta, xi, theta, beta, delta, N, J, Y, tau, alpha, gamma, kappa, t_resp){
  prob_eta <- f.eta(eta, pi_eta, xi, theta, beta, delta, N, J, Y)
  post <- rbinom(J, 1, prob_eta[,2])
  return(post)
}


