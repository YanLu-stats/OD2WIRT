# OD2WIRT
Two-way Outlier Detection using IRT models
---
title: "ReadMe"
output: html_document
---


<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# OD2WIRT

<!-- badges: start -->
<!-- badges: end -->

The goal of OD2WIRT is to provide a MCMC routine for the Two-way Outlier Detection Model.


## Example: Data Generation

This is an example showing how to simulate binary item response data and response time data.

```{r example}
source("simData.R")

Y.sim <- simdata_M2(N = 1000, J = 50, pi_xi = 0.2, pi_eta = 0.5,
                    s2_theta = 0.5, cov_theta_tau = 0.2, s2_tau = 0.5,
                    s2_beta = 0.8, cov_beta_alpha = 0.2, s2_alpha = 0.5,
                    mu_beta = 0, mu_alpha = 0, delta = 1, gamma = 1, kappa = 0.5)

```

To get an item response data matrix Y, use
```{r example}
Y <- Y.sim$obs.resp
```

To get a response time matrix resp.T, use
```{r example}
resp.T <- Y.sim$obs.time
```

To get a list of parameter values used in the simulation, use
```{r example}
par1.sim <- Y.sim$par1
par2.sim <- Y.sim$par2
```


## Example: Model Estimation

### Reduced Model

The example below shows you how to estimate the reduced model using an Metropolis-Hastings (MH) based {M}C{M}C algorithm.

Source the log density functions and {M}C{M}C samplers first:
```{r example}
source("MH_rw.R") ## Random-walk MH sampler
source("log_target.R") ## 
source("mcmc_M1.R")  ## Main function of MCMC for the reduced model
source("mcmc_M2.R")  ## Main function of MCMC for the full model
#source("tmcmc_M1.R") 
#source("tmcmc_M2.R") 

#' functions from other R files may be used as well
```

Alternatively, you can source R files needed to implement the algorithm all at once:
```{r example}
file.sources <- list.files(pattern="*.R")
sapply(file.sources,source, .GlobalEnv)
```

Then specify the data, initial values for parameters, the lengths of burn-in and {M}C{M}C chains, step sizes, and temperature spacing. Run the algorithm:
```{r example}
## Initials

res <- M1.mcmc(Y = Y, 
                par1 = Y.sim$par1,
                par2 = Y.sim$par2,
                M = 5000,
                L = 1000,
                step.sizes = c(0.05, 0.02, 0.01))

```
Obtain posterior samples and summary statistics like the mean:
```{r example}

post.xi <- res$par1$xi
post.xi.mean <- colMeans(post.xi)
post.eta <- res$par1$eta
post.eta.mean <- colMeans(post.eta)

post.theta <- res$par1$theta
post.theta.mean <- colMeans(post.theta)
post.beta <- res$par1$beta
post.beta.mean <- colMeans(post.beta)

post.delta <- res$par2$delta
post.delta.mean <- mean(post.delta)
post.pi_xi <- res$par2$pi_xi
post.pi_xi.mean <- mean(post.pi_xi)
post.pi_eta <- res$par2$pi_eta
post.pi_eta.mean <- mean(post.pi_eta)
post.s2_theta <- res$par2$s2_theta
post.s2_theta.mean <- mean(post.s2_theta)
post.mu_beta <- res$par2$mu_beta
post.mu_beta.mean <- mean(post.mu_beta)
post.s2_beta <- res$par2$s2_beta
post.s2_beta.mean <- mean(post.s2_beta)
```

Calculate the marginal DIC for model comparison:
```{r example}
source("mDIC.R")

```

### Full Model

Estimate the full model using the tempered {M}C{M}C algorithm. Let $L=1$ if the untempered version of {M}C{M}C is used:
```{r example}
source("mcmc_M2.R") ## Main function of MCMC for the full model
```


```{r example}
## Initials

res2 <- M2.tmcmc(Y = Y,
                 t_resp = resp.T,
                 par1 = par1.sim,
                 par2 = par1.sim,
                 M = 5000,
                 L = 1,
                 step.sizes = c(0.05, 0.02, 0.01))
  

```
