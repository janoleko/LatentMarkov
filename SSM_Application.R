library(fHMM)

# data = download_data("^GSPC")
data = download_data("BTC-USD")
data$return = c(NA, diff(log(data$Close)))

hist(data$return, xlim = c(-0.1,0.1), prob = TRUE, border = "white", breaks = 100)
# heavier tails than a normal distribution

plot(data$return, type = "l")

# Likelihood functions ----------------------------------------------------

### long version
mllk_ssm_long = function(theta.star, y, bm, m){
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
  Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
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

### short version using the package Lcpp
library(Lcpp)
# with the building-block functions from the R package Lcpp, we can write the 
# negative log-likelihood function much shorter and make it fuch faster (by a factor of 10-20).
mllk_ssm_short = function(theta.star, y, bm, m){
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
  Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
  delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
  # approximating state-dependent density based on midpoints
  allprobs = matrix(1, length(y), m)
  ind = which(!is.na(y))
  allprobs[ind,] = t(sapply(y[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))
  # forward algorithm to calculate the approximate log-likelihood recursively
  -forward(delta, Gamma, allprobs)
}


# Model fitting -----------------------------------------------------------

# initial parameters for numerical optimisation
hist(data$return, xlim = c(-0.1,0.1), prob = TRUE, border = "white", breaks = 100)
curve(dnorm(x, 0, sd(data$return[-1])), add = TRUE)
curve(dnorm(x, 0, 0.005), add = TRUE)
curve(dnorm(x, 0, 0.05), add = TRUE)

phi0 = 0.9
sigma0 = 0.6
beta0 = 0.02
mu0 = 0.001
qnorm(0.025, 0, sigma0/sqrt(1-phi0^2)) # 0.025 quantiale of stationary distribution

# working scale initial parameter vector
theta.star0 = c(qlogis(phi0), log(sigma0), log(beta0), mu0)

mod_BC1 = nlm(mllk_ssm_short, theta.star0, y = data$return, bm = 3.5, m = 300,
                    print.level = 2, hessian = TRUE, stepmax = 10)

# obtaining the estimated parameters
(phi = plogis(mod_SP1$estimate[1]))
(sigma = exp(mod_SP1$estimate[2]))
(beta = exp(mod_SP1$estimate[3]))
(mu = mod_SP1$estimate[4]) # expected return slightly positive
qnorm(0.025, 0, sigma/sqrt(1-phi^2))

## state decoding
bm = 3; m = 300
# defining intervals and gridpoints for numerical integration
b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
h = b[2]-b[1] # interval width
bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
# approximating tpm resulting from midpoint quadrature
Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution of AR(1) process
# approximating state-dependent density based on midpoints
allprobs = matrix(1, length(data$return), m)
ind = which(!is.na(data$return))
allprobs[ind,] = t(sapply(data$return[ind], dnorm, mean = mu, sd = beta * exp(bstar/2)))

mod_BC1$rawstates = viterbi(delta, Gamma, allprobs)
mod_BC1$states = bstar[mod_BC1$rawstates]

par(mfrow = c(1,1))
a = beta*exp(mod_BC1$states/2)
library(scales)
plot(data$return, type = "l", bty = "n", ylab = "return")
lines(2*beta*exp(mod_BC1$states/2)-0.45, type = "l", col = "orange")

