# Data --------------------------------------------------------------------

library(msm)
data = fev

# getting rid of the outliers at 1000

# this is wrong: 999 signals death
# how to deal with this properly?

library(dplyr)
# data = data %>% filter(fev < 800)
data$fev[which(data$fev == 999)] = NA
 
## little EDA
plot(data$days, data$fev, xlab = "time", ylab = "fav", bty = "n")
hist(data$fev, prob = T, border = "white", main = "Histogram of fev", xlab = "fev")

# Likelihood functions ----------------------------------------------------

### long version

mllk_ct_long = function(theta.star, N, X){
  mu = theta.star[1:N]
  sigma = exp(theta.star[N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[2*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary distribution of the
  # continuous-time Markov chain
  
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = array(dim = c(N,N,n_p-1))
    for(t in 1:(n_p-1)){
      Qube[,,t] = expm::expm(Q*timediff[t])
    }
    allprobs = matrix(1, nrow = n_p, ncol = N)
    ind = which(!is.na(X_p$fev))
    for(j in 1:N){
      allprobs[ind,j] = dnorm(X_p$fev[ind], mu[j], sigma[j])
    }
    # forward algorithm to calculate the log-likelihood recursively
    foo = delta%*%diag(allprobs[1,])
    l_p = log(sum(foo))
    phi = foo / sum(foo)
    for(t in 2:n_p){
      foo = phi%*%Qube[,,t-1]%*%diag(allprobs[t,])
      l_p = l_p + log(sum(foo))
      phi = foo / sum(foo)
    }
    l = l + l_p # adding patient contribution to overall log likelihood
  }
  return(-l)
}

### shorter and faster version
mllk_ct_short = function(theta.star, N, X){
  mu = theta.star[1:N]
  sigma = exp(theta.star[N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[2*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary distribution of the
  # continuous-time Markov chain
  
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = Lcpp::tpm_cont(Q, timediff)
    allprobs = matrix(1, nrow = n_p, ncol = N)
    ind = which(!is.na(X_p$fev))
    for(j in 1:N){
      allprobs[ind,j] = dnorm(X_p$fev[ind], mu[j], sigma[j])
    }
    # forward algorithm to calculate the log-likelihood recursively
    l = l + Lcpp::forward_g(delta, Qube, allprobs)
  }
  return(-l)
}


# Fitting a continuous-time HMM -------------------------------------------

## 2-state model

# initial parameters for numerical optimisation
mean0 = c(50, 100)
sd0 = c(20,20)
qs0 = c(1/10, 1/10)

# working scale initial parameter vector
theta.star0 = c(mean0, log(sd0), log(qs0))

# fitting the model by numerically optimising the log-likelihood with nlm
mod_fev1 = nlm(mllk_ct_short, theta.star0, N = 2, X = data,
               iterlim = 1000, print.level = 2, stepmax = 10)

# obtaining the estimated parameters
N = 2
## state-dependent distributions
(mu = mod_fev1$estimate[1:N])
(sigma = exp(mod_fev1$estimate[N+1:N]))
# generator matrix
Q = diag(N)
Q[!Q] = exp(mod_fev1$estimate[2*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
round(Q,4)
(delta = solve(t(Q+1), rep(1,N), tol = 1e-20))

# visualising the fitted state-dependent and marginal distributions
color = c("orange", "deepskyblue") # color vector

hist(data$fev, prob = T, border = "white", main = "Histogram of fev", xlab = "fev")
for(j in 1:N) curve(delta[j]*dnorm(x, mu[j], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dnorm(x, mu[1], sigma[1])+delta[2]*dnorm(x, mu[2], sigma[2]),
        add = T, lwd = 2, lty = 2)

##################################################################

## 3-state model

# initial parameters for numerical optimisation
mean0 = c(50, 90, 100)
sd0 = c(25, 20, 20)
qs0 = rep(1/100, 6)

# working scale initial parameter vector
theta.star0 = c(mean0, log(sd0), log(qs0))

# fitting the model by numerically optimising the log-likelihood with nlm
mod_fev2 = nlm(mllk_ct_short, theta.star0, N = 3, X = data,
               iterlim = 1000, print.level = 2, stepmax = 10)

# obtaining the estimated parameters
N = 3
## state-dependent distributions
(mu = mod_fev2$estimate[1:N])
(sigma = exp(mod_fev2$estimate[N+1:N]))
# generator matrix
Q = diag(N)
Q[!Q] = exp(mod_fev2$estimate[2*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
Q
(delta = solve(t(Q+1), rep(1,N), tol = 1e-20))

# visualising the fitted state-dependent and marginal distributions
color = c("orange", "deepskyblue", "seagreen2") # color vector

hist(data$fev, prob = T, border = "white", main = "Histogram of fev", xlab = "fev", breaks = 80)
for(j in 1:N) curve(delta[j]*dnorm(x, mu[j], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dnorm(x, mu[1], sigma[1])+delta[2]*dnorm(x, mu[2], sigma[2])+
        delta[3]*dnorm(x, mu[3], sigma[3]),
      add = T, lwd = 2, lty = 2)



# Adding covariate accute -------------------------------------------------

# Likelihood functions ----------------------------------------------------

### long version

mllk_ct_long_cov = function(theta.star, N, X){
  beta = matrix(theta.star[1:(2*N)], ncol = 2)
  sigma = exp(theta.star[2*N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[3*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary distribution of the
  # continuous-time Markov chain
  
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = array(dim = c(N,N,n_p-1))
    for(t in 1:(n_p-1)){
      Qube[,,t] = expm::expm(Q*timediff[t])
    }
    allprobs = matrix(1, nrow = n_p, ncol = N)
    ind = which(!is.na(X_p$fev))
    for(j in 1:N){
      allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
    }
    # forward algorithm to calculate the log-likelihood recursively
    foo = delta%*%diag(allprobs[1,])
    l_p = log(sum(foo))
    phi = foo / sum(foo)
    for(t in 2:n_p){
      foo = phi%*%Qube[,,t-1]%*%diag(allprobs[t,])
      l_p = l_p + log(sum(foo))
      phi = foo / sum(foo)
    }
    l = l + l_p # adding patient contribution to overall log likelihood
  }
  return(-l)
}

### shorter and faster version
mllk_ct_short_cov = function(theta.star, N, X){
  beta = matrix(theta.star[1:(2*N)], ncol = 2)
  sigma = exp(theta.star[2*N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[3*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary distribution of the
  # continuous-time Markov chain
  
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = Lcpp::tpm_cont(Q, timediff) # exp(Q*dt)
    allprobs = matrix(1, nrow = n_p, ncol = N)
    ind = which(!is.na(X_p$fev))
    for(j in 1:N){
      allprobs[ind,j] = dnorm(X_p$fev[ind], beta[j,1]+beta[j,2]*X_p$acute[ind], sigma[j])
    }
    # forward algorithm to calculate the log-likelihood recursively
    l = l + Lcpp::forward_g(delta, Qube, allprobs)
  }
  return(-l)
}


## 2-state model

# initial parameters for numerical optimisation
beta0 = c(50, 100) # intercept
sd0 = c(20,20)
qs0 = c(1/10, 1/10)

# working scale initial parameter vector
theta.star0 = c(beta0, rep(0,2), log(sd0), log(qs0))

# fitting the model by numerically optimising the log-likelihood with nlm
mod_fev3 = nlm(mllk_ct_short_cov, theta.star0, N = 2, X = data,
               iterlim = 1000, print.level = 2, stepmax = 10, hessian = TRUE)

I_fev3 = solve(mod_fev3$hessian)
# obtaining the estimated parameters
N = 2
## state-dependent distributions
(beta = matrix(mod_fev3$estimate[1:(2*N)], ncol = 2))
beta_sds = matrix(sqrt(diag(I_fev3)[1:(2*N)]), ncol = 2)
beta+2*beta_sds
beta-2*beta_sds
# covariate effects are significant

(sigma = exp(mod_fev3$estimate[2*N+1:N]))
Q = diag(N) # generator matrix
Q[!Q] = exp(mod_fev3$estimate[3*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
Q
(delta = solve(t(Q+1), rep(1,N), tol = 1e-20))


par(mfrow = c(2,1))
hist(data$fev[which(data$acute==0)], prob = T, border = "white", 
     main = "accute = 0", xlab = "fev", breaks = 50, ylim = c(0,0.02), xlim = c(0,150))
for(j in 1:N) curve(delta[j]*dnorm(x, beta[j,1], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dnorm(x, beta[1,1], sigma[1])+delta[2]*dnorm(x, beta[2,1], sigma[2]),
      add = T, lwd = 2, lty = 2)

hist(data$fev[which(data$acute==1)], prob = T, border = "white", 
     main = "acute = 1", xlab = "fev", breaks = 50, ylim = c(0,0.02), xlim = c(0,150))
for(j in 1:N) curve(delta[j]*dnorm(x, beta[j,1]+beta[j,2], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dnorm(x, beta[1,1]+beta[1,2], sigma[1])+delta[2]*dnorm(x, beta[2,1]+beta[2,2], sigma[2]),
      add = T, lwd = 2, lty = 2)



# same analysis but with gamma state-dependent distributions --------------

### shorter and faster version
mllk_ct_short_cov_gamma = function(theta.star, N, X){
  beta = matrix(theta.star[1:(2*N)], ncol = 2)
  sigma = exp(theta.star[2*N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[3*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary distribution of the
  # continuous-time Markov chain
  
  # for the rest: loop over patients
  l = 0 # initialize log likelihood
  ptnums = unique(X$ptnum) # patient numbers
  for(i in ptnums){
    X_p = X[which(X$ptnum==i),]
    n_p = nrow(X_p)
    timediff = diff(X_p$days)
    Qube = Lcpp::tpm_cont(Q, timediff) # exp(Q*dt)
    allprobs = matrix(1, nrow = n_p, ncol = N)
    ind = which(!is.na(X_p$fev))
    for(j in 1:N){
      mu = exp(beta[j,1]+beta[j,2]*X_p$acute[ind])
      allprobs[ind,j] = dgamma(X_p$fev[ind], shape=mu^2/sigma[j]^2, scale=sigma[j]^2/mu)
    }
    # forward algorithm to calculate the log-likelihood recursively
    l = l + Lcpp::forward_g(delta, Qube, allprobs)
  }
  return(-l)
}


# initial parameters for numerical optimisation
beta0 = c(50, 100) # intercept
sd0 = c(20,20)
qs0 = c(1/10, 1/10)

# working scale initial parameter vector
theta.star0 = c(log(beta0), rep(0,2), log(sd0), log(qs0))

# fitting the model by numerically optimising the log-likelihood with nlm
mod_fev4 = nlm(mllk_ct_short_cov_gamma, theta.star0, N = 2, X = data,
               iterlim = 1000, print.level = 2, stepmax = 10, hessian = TRUE)
I_fev4 = solve(mod_fev4$hessian)

N = 2
## state-dependent distributions
(beta = matrix(mod_fev4$estimate[1:(2*N)], ncol = 2))
beta_sds = matrix(sqrt(diag(I_fev4)[1:(2*N)]), ncol = 2)
beta+2*beta_sds
beta-2*beta_sds
# covariate effects are significant

(sigma = exp(mod_fev4$estimate[2*N+1:N]))
Q = diag(N) # generator matrix
Q[!Q] = exp(mod_fev4$estimate[3*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
Q
(delta = solve(t(Q+1), rep(1,N), tol = 1e-20))


par(mfrow = c(2,1))
hist(data$fev[which(data$acute==0)], prob = T, border = "white", 
     main = "accute = 0", xlab = "fev", breaks = 50, ylim = c(0,0.02), xlim = c(0,150))
mu = exp(beta[,1]); shape = mu^2/sigma^2; scale = sigma^2/mu
for(j in 1:N) curve(delta[j]*dgamma(x, shape=shape[j], scale = scale[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dgamma(x, shape=shape[1], scale=scale[1])+delta[2]*dgamma(x, shape=shape[2], scale=scale[2]),
      add = T, lwd = 2, lty = 2)

hist(data$fev[which(data$acute==1)], prob = T, border = "white", 
     main = "acute = 1", xlab = "fev", breaks = 50, ylim = c(0,0.02), xlim = c(0,150))
mu = exp(beta[,1]+beta[,2]); shape = mu^2/sigma^2; scale = sigma^2/mu
for(j in 1:N) curve(delta[j]*dgamma(x, shape=shape[j], scale = scale[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dgamma(x, shape=shape[1], scale=scale[1])+delta[2]*dgamma(x, shape=shape[2], scale=scale[2]),
      add = T, lwd = 2, lty = 2)
