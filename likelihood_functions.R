
## This file contains all the negative log likelihood functions used in the case studies.

# For each model, we provide two versions of the same function, one only using base R functions
# computing the forward algorithm in R, and one using our R package LaMa to speed up computation and shortening the code.


## brief overview of the objects used throughout this script:

# theta.star: vector containing (unconstrained) initial values for the model parameters
# X: observed data sequence
# N: number of states


# install.packages("LaMa")
library(LaMa) # just to be save, this will be loaded when the functions are sourced.
# install.packages("expm")
library(expm)


# HMM case study ----------------------------------------------------------

## homogeneous HMM

### long version
mllk_slow = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  mu.turn = theta.star[2 * N + 1:N]
  kappa = exp(theta.star[3 * N + 1:N])
  ## transition probability matrix
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[4 * N + 1:(N * (N - 1))])
  Gamma = Gamma / rowSums(Gamma)
  ## initial distribution -> stationary
  delta = solve(t(diag(N) - Gamma + 1), rep(1, N))
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind, j] = dgamma(X$step[ind], shape = mu[j]^2 / sigma[j]^2, 
                              scale = sigma[j]^2 / mu[j]) *
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  foo = delta %*% diag(allprobs[1, ])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for(t in 2:nrow(X)){
    foo = phi %*% Gamma %*% diag(allprobs[t, ])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  return(-l)
}

### short version using the package LaMa

# with the building-block functions from the R package LaMa, we can write the 
# negative log-likelihood function much shorter and make it much faster (by a factor of 10-20).
mllk_fast = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  mu.turn = theta.star[2 * N + 1:N]
  kappa = exp(theta.star[3 * N + 1:N])
  ## transition probability matrix
  Gamma = LaMa::tpm(theta.star[4 * N + 1:(N * (N - 1))]) # convenient LaMa function
  ## initial distribution -> stationary
  delta = LaMa::stationary(Gamma) # convenient LaMa function
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind, j] = dgamma(X$step[ind], shape = mu[j]^2 / sigma[j]^2, 
                             scale = sigma[j]^2 / mu[j]) *
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  -LaMa::forward(delta, Gamma, allprobs) # convenient LaMa function
}


## inhomogeneous HMM with covariate time of day

### long version
mllk_HMM_slow = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  mu.turn = theta.star[2 * N + 1:N]
  kappa = exp(theta.star[3 * N + 1:N])
  ## regression coefs for tpm
  beta = matrix(theta.star[4 * N + 1:(3 * N * (N - 1))], ncol = 3)
  ## 24 unique tpms
  Gamma = array(dim = c(N, N, 24))
  tod = unique(data$timeOfDay)
  for(t in 1:24){
    gamma = diag(N)
    gamma[!gamma] = exp(beta[, 1] + beta[, 2] * sin(2 * pi * tod[t] / 24) +
                          beta[, 3] * cos(2 * pi * tod[t] / 24))
    Gamma[, , t] = gamma / rowSums(gamma)
  }
  ## initial distribution -> periodically stationary (Koslik et al., 2023)
  GammaT = Gamma[, , X$timeOfDay2[1]]
  for(t in 2:24) GammaT = GammaT %*% Gamma[, , X$timeOfDay2[t]]
  delta = solve(t(diag(N) - GammaT + 1), rep(1, N))
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind, j] = dgamma(X$step[ind], shape = mu[j]^2 / sigma[j]^2, 
                              scale = sigma[j]^2 / mu[j]) *
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  foo = delta %*% diag(allprobs[1, ])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for(t in 2:nrow(X)){
    foo = phi %*% Gamma[, , X$timeOfDay2[t]] %*% diag(allprobs[t, ])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  return(-l)
}

### short version using the package LaMa
mllk_HMM_fast = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  mu.turn = theta.star[2 * N + 1:N]
  kappa = exp(theta.star[3 * N + 1:N])
  ## regression coefs for tpm
  beta = matrix(theta.star[4 * N + 1:(3 * N * (N - 1))], ncol = 3)
  ## 24 unique tpms
  Gamma = tpm_p(tod = unique(X$timeOfDay), L = 24, beta = beta) 
  # tpm_p() does the sine and cosine expansion by default. For even higher efficiency, 
  # a design matrix should be calculated before model fitting, which can be supplied to tpm_p(). 
  # But this way is more convenient, when speed does not matter too much.
  ## initial distribution -> periodically stationary
  delta = stationary_p(Gamma, X$timeOfDay2[1])
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind, j] = dgamma(X$step[ind], shape = mu[j]^2 / sigma[j]^2, 
                              scale = sigma[j]^2 / mu[j]) *
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  -forward_p(delta, Gamma, allprobs, X$timeOfDay2)
}


# Continuos-time HMM case study -------------------------------------------

### long version
mllk_ct_slow = function(theta.star, X){
  beta = matrix(theta.star[1:4], ncol = 2)
  sigma = exp(theta.star[4+1:2])
  # structured generator matrix
  Q = matrix(0, 3, 3)
  Q[1,2:3] = exp(theta.star[6+1:2])
  Q[2,3] = exp(theta.star[9])
  diag(Q) = -rowSums(Q)
  delta = c(exp(theta.star[10:11]),0)
  delta = delta/sum(delta)
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = array(dim = c(3,3,n_p-1))
    for(t in 1:(n_p-1)){
      Qube[,,t] = expm::expm(Q*timediff[t])
    }
    allprobs = matrix(1, nrow = n_p, ncol = 3)
    ind = which(!is.na(X_p$fev))
    for(j in 1:2){
      allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
    }
    allprobs[,3] = 0
    allprobs[which(X_p$fev==999),] = c(0,0,1)
    # forward algorithm to calculate the log-likelihood recursively
    foo = delta%*%diag(allprobs[1,])
    phi = foo/sum(foo)
    l_p = log(sum(foo))
    for(t in 2:n_p){
      foo = phi%*%Qube[,,t-1]%*%diag(allprobs[t,])
      phi = foo/sum(foo)
      l_p = l_p + log(sum(foo))
    }
    l = l + l_p
  }
  return(-l)
}

### short version using the package LaMa
mllk_ct_fast = function(theta.star, X){
  beta = matrix(theta.star[1:4], ncol = 2)
  sigma = exp(theta.star[4+1:2])
  # structured generator matrix
  Q = matrix(0, 3, 3)
  Q[1,2:3] = exp(theta.star[6+1:2])
  Q[2,3] = exp(theta.star[9])
  diag(Q) = -rowSums(Q)
  delta = c(exp(theta.star[10:11]),1)
  delta = delta/sum(delta)
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = LaMa::tpm_cont(Q, timediff) # exp(Q*dt)
    allprobs = matrix(1, nrow = n_p, ncol = 3)
    ind = which(!is.na(X_p$fev))
    for(j in 1:2){
      allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
    }
    allprobs[,3] = 0
    allprobs[which(X_p$fev==999),] = c(0,0,1)
    # forward algorithm to calculate the log-likelihood recursively
    l = l + LaMa::forward_g(delta, Qube, allprobs)
  }
  return(-l)
}


# State-space model case study --------------------------------------------

### long version in base R
mllk_ssm_slow = function(theta.star, y, bm, m){
  # transforming parameters to natural
  phi = plogis(theta.star[1])
  sigma = exp(theta.star[2])
  beta = exp(theta.star[3])
  # defining intervals and gridpoints for numerical integration
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximating tpm resulting from midpoint quadrature
  Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
  delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
  # approximating state-dependent density based on midpoints
  allprobs = matrix(1, length(y), m)
  ind = which(!is.na(y))
  allprobs[ind] = t(sapply(y[ind], dnorm, mean = 0, sd = beta * exp(bstar/2)))
  # forward algorithm to calculate the approximate log-likelihood recursively
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for(t in 2:nrow(X)){
    foo = phi%*%Gamma%*%diag(allprobs[t,])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  return(-l)
}

### short version using the package LaMa
mllk_ssm_fast = function(theta.star, y, bm, m){
  # transforming parameters to natural
  phi = plogis(theta.star[1])
  sigma = exp(theta.star[2])
  beta = exp(theta.star[3])
  mu = theta.star[4]
  # defining intervals and gridpoints for numerical integration
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximating tpm resulting from midpoint quadrature
  Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
  delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
  # approximating state-dependent density based on midpoints
  allprobs = matrix(1, length(y), m)
  ind = which(!is.na(y))
  allprobs[ind,] = t(sapply(y[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))
  # forward algorithm to calculate the approximate log-likelihood recursively
  -forward(delta, Gamma, allprobs)
}


# Continuous-time SSM case study ------------------------------------------

### longer version in base R, in this case only the forward algorithm differs
mllk_ctSSM_slow = function(theta.star, X, deltat, bm, m){
  # OU process parameters
  theta = exp(theta.star[1])
  sigma = exp(theta.star[2])
  # intercept of the state-dependent process
  beta = theta.star[3] # intercept
  # construction of intervals for numerical integration
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximate transition densities with time-dependent variance
  Gamma = array(0, dim = c(m, m, length(deltat))) # transition probability matrices for unique time differences
  for (t in 1:length(deltat)) {
    Dt = deltat[t]
    Gamma[,,t] = sapply(bstar, dnorm, mean = exp(-theta * Dt) * bstar,
                        sd = sqrt((1 - exp(-2 * theta * Dt)) * sigma^2 / (2 * theta))) * h
  }
  # initial distribution = stationary distribution of OU process
  delta = dnorm(bstar, 0, sqrt(sigma^2 / (2 * theta))) * h
  # approximating state-dependent density based on interval midpoints
  allprobs = t(sapply(X$seven_meter_success, dbinom, size = 1, p = plogis(beta + bstar)))
  
  # forward algorithm to calculate the approximate log-likelihood recursively
  uIDs = unique(X$uID)
  l = 0
  for(i in 1:length(uIDs)){
    ind = which(X$uID == uIDs[i])
    Xi = X[ind,]
    
    # forward algorithm
    foo = delta%*%diag(allprobs[1,])
    li = log(sum(foo))
    phi = foo / sum(foo)
    for(t in 2:nrow(Xi)){
      foo = phi %*% Gamma[,,Xi$match2Array[t]] %*% diag(allprobs[t,])
      li = li + log(sum(foo))
      phi = foo / sum(foo)
    }
    l = l + li
  }
  -sum(l)
}

### short version using the package LaMa
mllk_ctSSM_fast = function(theta.star, X, deltat, bm, m){
  # OU process parameters
  theta = exp(theta.star[1])
  sigma = exp(theta.star[2])
  # intercept of the state-dependent process
  beta = theta.star[3] # intercept
  # construction of intervals for numerical integration
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximate transition densities with time-dependent variance
  Gamma = array(0, dim = c(m, m, length(deltat))) # transition probability matrices for unique time differences
  for (t in 1:length(deltat)) {
    Dt = deltat[t]
    Gamma[,,t] = sapply(bstar, dnorm, mean = exp(-theta * Dt) * bstar, 
                        sd = sqrt((1 - exp(-2 * theta * Dt)) * sigma^2 / (2 * theta))) * h
  }
  # initial distribution = stationary distribution of OU process
  delta = dnorm(bstar, 0, sqrt(sigma^2 / (2 * theta))) * h 
  # approximating state-dependent density based on interval midpoints
  allprobs = t(sapply(X$seven_meter_success, dbinom, size = 1, p = plogis(beta + bstar)))
  
  # forward algorithm to calculate the approximate log-likelihood recursively
  -LaMa::forward_g(delta, Gamma[,,X$match2Array], allprobs, trackInds)
}


# MMPP case study ---------------------------------------------------------

### longer version in base R
mllk_mmpp_slow = function(theta.star, timediff, N=2){
  lambda = exp(theta.star[1:N]) # state specific rates
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
  # we split the Omega matrix into (Q-Lambda)*dt and Lambda
  Qube = array(dim = c(N,N,length(timediff)))
  for(t in 1:length(timediff)){
    Qube[,,t] = expm::expm((Q-diag(lambda))*timediff[t])
  }
  allprobs = matrix(lambda, nrow = length(timediff)+1, ncol = N, byrow = T) # Lambda
  allprobs[1,] = 1 # no Lambda matrix at the start of the likelihood
  # forward algorithm in R
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for(t in 1:length(timediff)){
    foo = phi %*% Qube[,,t] %*% diag(allprobs[t+1,])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  -l
}

### short version using the package LaMa
mllk_mmpp_fast = function(theta.star, timediff, N=2){
  lambda = exp(theta.star[1:N]) # state specific rates
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary initial distribution
  # we split the Omega matrix into (Q-Lambda)*dt and Lambda for efficiency
  Qube = LaMa::tpm_cont(Q-diag(lambda), timediff) # (Q-Lambda)*dt
  allprobs = matrix(lambda, nrow = length(timediff)+1, ncol = N, byrow = T) # Lambda
  allprobs[1,] = 1 # no Lambda matrix at the start of the likelihood
  # forward algorithm in C++ using LaMa
  -LaMa::forward_g(delta, Qube, allprobs)
}


