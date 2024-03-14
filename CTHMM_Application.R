# Data --------------------------------------------------------------------

library(msm)
data = fev

# getting rid of the outliers at 1000
library(dplyr)
data = data %>% filter(fev < 800)
 
## little EDA
plot(data$days, data$fev, xlab = "time", ylab = "fav", bty = "n")
hist(data$fev, prob = T, border = "white", main = "Histogram of fev", xlab = "fev")

order(data$days)

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
