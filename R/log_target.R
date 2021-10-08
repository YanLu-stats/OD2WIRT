
#' -----------------------------------------------------
#' log posterior density functions of theta, beta, delta
#' -----------------------------------------------------


f.theta <- function(theta, s2_theta, xi, eta, beta, delta, N, J, Y){
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta #'N x J
  temp <- rowSums(Y * temp1 - log(1 + exp(temp1))) + log(dnorm(theta, 0, sqrt(s2_theta)))
  return(temp)
}


f.xi <- function(xi, pi_xi, theta, eta, beta, delta, N, J, Y){
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta #'N x J
  temp <- rowSums(Y * temp1 - log(1 + exp(temp1))) + xi*log(pi_xi)+(1-xi)*log(1-pi_xi)
  return(temp)
}



f.eta <- function(eta, pi_eta, theta, xi, beta, delta, N, J, Y){
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  temp <- colSums(Y * temp1 - log(1 + exp(temp1))) + eta*log(pi_eta)+(1-eta)*log(1-pi_eta)
  return(temp)
}




f.beta <- function(beta, mu_beta, s2_beta, xi, theta, eta, delta, N, J, Y){
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  temp <- colSums(Y * temp1 - log(1 + exp(temp1))) + log(dnorm(beta, mu_beta, sqrt(s2_beta)))
  return(temp)
}



f.delta <- function(delta, sigma_delta, xi, theta, eta, beta, N, J, Y){
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(eta) * delta
  #temp <- sum(Y * temp1 - log(1 + exp(temp1))) + log(dgamma(delta, alpha_delta, scale = beta_delta))
  temp <- sum(Y * temp1 - log(1 + exp(temp1))) + log(extraDistr::dhcauchy(delta, sigma_delta))
  return(temp)
}


f.s2_theta <- function(s2_theta, theta, alpha_theta, beta_theta, N){
  log(nimble::dinvgamma(s2_theta, alpha_theta, beta_theta)) 
  + sum(log(dnorm(theta, 0, sqrt(s2_theta)))) 
}
#log(nimble::dinvgamma(s2_theta, (1 + N/2), (1 + sum(theta^2)/2))) 

f.s2_beta <- function(s2_beta, mu_beta, beta, alpha_beta, beta_beta, J){
  log(nimble::dinvgamma(s2_beta, alpha_beta, beta_beta))
  + sum(log(dnorm(beta, mu_beta, sqrt(s2_beta))))
}
#log(nimble::dinvgamma(s2_theta, (1 + J/2), (1 + sum((beta - mu_beta)^2)/2)))

f.mu_beta <- function(mu_beta, s2_beta, beta, mu_b, s2_b, J){
  log(dnorm(mu_beta, mu_b, s2_b^(-1/2)))
    sum(log(dnorm(beta, mu_beta, sqrt(s2_beta))))
}
#log(dnorm(mu_beta, sum(beta)/(J + s2_beta/s2_b), (J/s2_beta + 1/s2_b)^(-1/2)))



f.pi_xi <- function(pi_xi, xi, alpha_xi, beta_xi) {
  log(dbeta(pi_xi, alpha_xi, beta_xi)) + sum(xi*log(pi_xi) + (1-xi)*log(1-pi_xi))
}

f.pi_eta <- function(pi_eta, eta, alpha_eta, beta_eta) {
  log(dbeta(pi_eta, alpha_eta, beta_eta)) + sum(eta*log(pi_eta) + (1-eta)*log(1-pi_eta))
}


f.xi <- function(xi, pi_xi, theta, eta, beta, delta, N, J, Y){
  
  temp <- matrix(0, N, 2)
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta)
  temp[,1] <- rowSums(Y * temp1 - log(1 + exp(temp1))) + log(1-pi_xi)
  
  temp2 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + rep(1, N) %*% t(eta)*delta
  temp[,2] <- rowSums(Y * temp2 - log(1 + exp(temp2))) + log(pi_xi)
  
  temp <- temp - apply(temp, 1, max)
  temp <- exp(temp)
  temp.sum <- rowSums(temp)
  prob <- temp/temp.sum
  return(prob)
}

sample.xi <- function(xi, pi_xi, theta, eta, beta, delta, N, J, Y){
  
  prob_xi <- f.xi(xi, pi_xi, theta, eta, beta, delta, N, J, Y)
  post <- rbinom(N, 1, prob_xi[,2])
  return(post)
}


f.eta <- function(eta, pi_eta, xi, theta, beta, delta, N, J, Y){
  
  temp <- matrix(0, J, 2)
  
  temp1 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta)
  temp[,1] <- colSums(Y * temp1 - log(1 + exp(temp1))) + log(1-pi_eta)
  
  temp2 <- theta %*% t(rep(1, J)) - rep(1, N) %*% t(beta) + xi %*% t(rep(1, J))*delta
  temp[,2] <- colSums(Y * temp2 - log(1 + exp(temp2))) + log(pi_eta)
  
  temp <- temp - apply(temp, 1, max)
  temp <- exp(temp)
  temp.sum <- rowSums(temp)
  prob <- temp/temp.sum
  return(prob)
}


sample.eta <- function(eta, pi_eta, xi, theta, beta, delta, N, J, Y){
  
  prob_eta <- f.eta(eta, pi_eta, xi, theta, beta, delta, N, J, Y)
  post <- rbinom(J, 1, prob_eta[,2])
  return(post)
}



