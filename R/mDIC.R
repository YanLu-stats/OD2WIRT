
#'-----------------------------------------
#' M1: Reduced Model
#'-----------------------------------------

## Posterior mean of the deviance OR posterior expected deviance
dev_rep.fun <- function(Y, theta, beta, xi, eta, delta){
  J <- ncol(Y)
  N <- nrow(Y)
  logits <- theta - beta + xi*eta*delta
  prob <- exp(logits)/(1+exp(logits))
  llk <- log(mean(prob))
  dev <- -2*J*N*llk
  return(dev)
}

## Deviance based on the posterior means of model parameters
D_est.fun <- function(n.rep, n.burnin, res){
  Y <- res$data
  J <- ncol(Y)
  N <- nrow(Y)
  M <- res$M
  post.delta <- res$par2$delta[(n.burnin+1):M]
  post.delta.mean <- mean(post.delta)
  post.pi_xi <- res$par2$pi_xi[(n.burnin+1):M]
  post.pi_xi.mean <- mean(post.pi_xi)
  post.pi_eta <- res$par2$pi_eta[(n.burnin+1):M]
  post.pi_eta.mean <- mean(post.pi_eta)
  post.s2_theta <- res$par2$s2_theta[(n.burnin+1):M]
  post.s2_theta.mean <- mean(post.s2_theta)
  post.mu_beta <- res$par2$mu_beta[(n.burnin+1):M]
  post.mu_beta.mean <- mean(post.mu_beta)
  post.s2_beta <- res$par2$s2_beta[(n.burnin+1):M]
  post.s2_beta.mean <- mean(post.s2_beta)
  
  theta <- rnorm(n.rep, 0, sqrt(post.s2_theta.mean))
  beta <- rnorm(n.rep, post.mu_beta.mean, sqrt(post.s2_beta.mean))
  xi <- rbinom(n.rep, 1, post.pi_xi)
  eta <- rbinom(n.rep, 1, post.pi_eta)
  
  logits <- theta - beta + xi*eta*post.delta.mean
  prob <- exp(logits)/(1+exp(logits))
  llk <- log(mean(prob))
  dev <- -2*J*N*llk
  return(dev)
}


#'-----------------------------------------
#' M0: Model without the Cheating Effect
#'-----------------------------------------

dev_rep.fun0 <- function(Y, theta, beta){
  J <- ncol(Y)
  N <- nrow(Y)
  logits <- theta - beta 
  prob <- exp(logits)/(1+exp(logits))
  llk <- log(mean(prob))
  dev <- -2*J*N*llk
  return(dev)
}


D_est.fun0 <- function(n.rep, n.burnin, res){
  Y <- res$data
  J <- ncol(Y)
  N <- nrow(Y)
  M <- res$M
  
  post.s2_theta <- res$par2$s2_theta[(n.burnin+1):M]
  post.s2_theta.mean <- mean(post.s2_theta)
  post.mu_beta <- res$par2$mu_beta[(n.burnin+1):M]
  post.mu_beta.mean <- mean(post.mu_beta)
  post.s2_beta <- res$par2$s2_beta[(n.burnin+1):M]
  post.s2_beta.mean <- mean(post.s2_beta)
  
  theta <- rnorm(n.rep, 0, sqrt(post.s2_theta.mean))
  beta <- rnorm(n.rep, post.mu_beta.mean, sqrt(post.s2_beta.mean))
  
  logits <- theta - beta
  prob <- exp(logits)/(1+exp(logits))
  llk <- log(mean(prob))
  dev <- -2*J*N*llk
  return(dev)
}



#'-----------------------------------------
#' M2: Full Model
#'-----------------------------------------


alpha <- rnorm(170, -0.4, sqrt(0.4))
tau <- rnorm(1624, 0, sqrt(0.25))

temp2 <- rep(1, 1624) %*% t(alpha) - tau %*% t(rep(1, 170)) - xi.real%*%t(eta.real)*0.6 #'N x J
temp <- log(dlnorm(t.real, temp2, sqrt(0.8)))

prob <- exp(logits)/(1+exp(logits))
llk <- log(mean(prob))
dev <- -2*J*N*llk


dev_rep.fun1 <- function(Y, theta, beta, xi, eta, delta){
  J <- ncol(Y)
  N <- nrow(Y)
  logits <- theta - beta + xi*eta*delta
  prob <- exp(logits)/(1+exp(logits))
  llk <- log(mean(prob))
  dev <- -2*J*N*llk
  return(dev)
}


D_est.fun1 <- function(n.rep, n.burnin, res){
  Y <- res$data
  J <- ncol(Y)
  N <- nrow(Y)
  M <- res$M
  post.delta <- res$par2$delta[(n.burnin+1):M]
  post.delta.mean <- mean(post.delta)
  post.pi_xi <- res$par2$pi_xi[(n.burnin+1):M]
  post.pi_xi.mean <- mean(post.pi_xi)
  post.pi_eta <- res$par2$pi_eta[(n.burnin+1):M]
  post.pi_eta.mean <- mean(post.pi_eta)
  post.s2_theta <- res$par2$s2_theta[(n.burnin+1):M]
  post.s2_theta.mean <- mean(post.s2_theta)
  post.mu_beta <- res$par2$mu_beta[(n.burnin+1):M]
  post.mu_beta.mean <- mean(post.mu_beta)
  post.s2_beta <- res$par2$s2_beta[(n.burnin+1):M]
  post.s2_beta.mean <- mean(post.s2_beta)
  
  theta <- rnorm(n.rep, 0, sqrt(post.s2_theta.mean))
  beta <- rnorm(n.rep, post.mu_beta.mean, sqrt(post.s2_beta.mean))
  xi <- rbinom(n.rep, 1, post.pi_xi)
  eta <- rbinom(n.rep, 1, post.pi_eta)
  
  logits <- theta - beta + xi*eta*post.delta.mean
  prob <- exp(logits)/(1+exp(logits))
  llk <- log(mean(prob))
  dev <- -2*J*N*llk
  return(dev)
}

