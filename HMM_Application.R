
# Data --------------------------------------------------------------------

elephant_data = read.csv("http://www.rolandlangrock.com/elephant_data.csv")

# Originally, 15 elephants from the Etosha National Park where fitted with GPS 
# collars to collect hourly movement data (@tsalyuk2019temporal). 
# For simplicity, the raw data have already been preprocessed and we consider 
# the track of only one individual.

## Defining time of day variable (could also be done with lubridate)
elephant_data$timestamp = strptime(elephant_data$timestamp, 
                                   "%Y-%m-%d %H:%M:%S", tz = "GMT")
hours = as.numeric(format(elephant_data$timestamp, "%H"))
minutes = as.numeric(format(elephant_data$timestamp, "%M"))
elephant_data$timeOfDay = hours + (minutes/60)

## calculating step lengths and turning angles
data = moveHMM::prepData(trackData = elephant_data, 
                         coordNames = c("location.long", "location.lat"), 
                         type = "LL")
data$timeOfDay2 = round(data$timeOfDay)

# we only have 2 step lengths that are exactly zero. 
sum(na.omit(data$step) == 0)
# This does not really justify using a zero inflated step-length distribution
# and we just replace these zeros with very small values
data$step[which(data$step == 0)] = 1e-5


# Likelihood functions ----------------------------------------------------

## homogeneous HMM

### long version
mllk_long = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N])
  ## transition probability matrix
  Gamma = diag(N)
  Gamma[!Gamma] = exp(theta.star[4*N+1:(N*(N-1))])
  Gamma = Gamma / rowSums(Gamma)
  ## initial distribution -> stationary
  delta = solve(t(diag(N)-Gamma+1), rep(1,N))
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
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
mllk_short = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N])
  ## transition probability matrix
  Gamma = tpm(theta.star[4*N+1:(N*(N-1))])
  ## initial distribution -> stationary
  delta = stationary(Gamma)
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  -forward(delta, Gamma, allprobs)
}


# Fitting a homogeneous HMM -----------------------------------------------

# initial parameters for numerical optimisation
stepMean0 = c(0.1, 0.4, 1)
stepSd0 = c(0.1, 0.4, 1)
angleMean0 = rep(0,3)
angleConcentration0 = c(0.1, 0.9, 2)

# working scale initial parameter vector
theta.star0 = c(log(stepMean0), log(stepSd0), # steppars
                angleMean0, log(angleConcentration0), # anglepars
                rep(-2, 6)) # tpmpars

# fitting the model by numerically optimising the log-likelihood with nlm
mod_elephant1 = nlm(mllk_short, theta.star0, X = data, N = 3, 
                    print.level = 2, hessian = TRUE)
# or slow
# mod_elephant1 = nlm(mllk_long, theta.star0, X = data, N = 3, 
#                     print.level = 2, hessian = TRUE)

# obtaining the estimated parameters
## parameters for state-dependent distributions
N = 3
(mu = exp(mod_elephant1$estimate[1:N]))
(sigma = exp(mod_elephant1$estimate[N+1:N]))
(mu.turn = mod_elephant1$estimate[2*N+1:N])
(kappa = exp(mod_elephant1$estimate[3*N+1:N]))
## transition probability matrix
Gamma = tpm(mod_elephant1$estimate[4*N+1:(N*(N-1))])
round(Gamma, 3)
## stationary distribution (proportion of time spent in each state)
(delta = stationary(Gamma))

# visualising the fitted state-dependent and marginal distributions
color = c("orange", "deepskyblue", "seagreen2") # color vector

# histogram and the marginal distribution
pdf("./figs/hmm_elephant.pdf", width = 8, height = 3)
par(mfrow = c(1,3))
## step length
hist(data$step, prob = TRUE, border = "white", xlim = c(0,4), breaks = 50,
     main = "", xlab = "step length", ylab = "density", cex = 0.9)
shape = mu^2/sigma^2; scale = sigma^2/mu
# component distributions
for(j in 1:N){
  curve(delta[j]*dgamma(x, shape=shape[j], scale=scale[j]), 
        add = T, col = color[j], lwd = 2, n = 300)
}
# marginal as sum
curve(delta[1]*dgamma(x, shape=shape[1], scale=scale[1])+
        delta[2]*dgamma(x, shape=shape[2], scale=scale[2])+
        delta[3]*dgamma(x, shape=shape[3], scale=scale[3]), 
      add = T, lty = 2, lwd = 2, n = 300)
legend("topright", col = color, lwd = 2, legend = paste("state", 1:3), bty = "n")

## turning angle
hist(data$angle, prob = TRUE, border = "white", 
     main = "", xlab = "turning angle", ylab = "density", cex = 0.9)
# component distributions
for(j in 1:N){
  curve(delta[j]*CircStats::dvm(x, mu.turn[j], kappa[j]), 
        add = T, col = color[j], lwd = 2, n = 300)
}
# marginal as sum
curve(delta[1]*CircStats::dvm(x, mu.turn[1], kappa[1])+
        delta[2]*CircStats::dvm(x, mu.turn[2], kappa[2])+
        delta[3]*CircStats::dvm(x, mu.turn[3], kappa[3]), 
      add = T, lty = 2, lwd = 2, n = 300)
legend("topright", col = color, lwd = 2, legend = paste("state", 1:3), bty = "n")
dev.off()


# state decoding
# for general state decoding we can also use Lcpp. We must however first calculate
# the allprobs matrix
allprobs = matrix(1, nrow(data), N)
ind = which(!is.na(data$step) & !is.na(data$angle))
for(j in 1:N){
  allprobs[ind,j] = dgamma(data$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
    CircStats::dvm(data$angle[ind], mu.turn[j], kappa[j])
}
mod_elephant1$states = viterbi(delta, Gamma, allprobs)

library(scales) # for color transparency
plot(data$x, data$y, type = "l", xlab = "x", ylab = "y")
points(data$x, data$y, pch = 20, col = alpha(color[mod_elephant1$states], 0.4))

plot(data$timestamp, data$step, col = color[mod_elephant1$states], bty = "n", 
     ylab = "step length", xlab = "time")
plot(data$timestamp, data$angle, col = color[mod_elephant1$states], bty = "n", 
     ylab = "turning angle", xlab = "time")


# Covariate effects: Time of day ------------------------------------------

# we assume the animal's behaviour varies with the time of day. 
# Thus we include the time of day as a covariate in our model.

# likelihood functions

### long version
mllk_long_t = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N])
  ## regression coefs for tpm
  beta = matrix(theta.star[4*N+1:(3*N*(N-1))], ncol = 3)
  ## 24 unique tpms
  Gamma = array(dim = c(N,N,24))
  tod = unique(data$timeOfDay)
  for(t in 1:24){
    gamma = diag(N)
    gamma[!gamma] = exp(beta[,1] + beta[,2]*sin(2*pi*tod[t]/24) +
                          beta[,3]*cos(2*pi*tod[t]/24))
    Gamma[,,t] = gamma / rowSums(gamma)
  }
  ## initial distribution -> periodically stationary
  GammaT = Gamma[,,X$timeOfDay2[1]]
  for(t in 2:24){
    GammaT = GammaT%*%Gamma[,,X$timeOfDay2[t]]
  }
  delta = solve(t(diag(N)-GammaT+1), rep(1,N))
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  foo = delta%*%diag(allprobs[1,])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for(t in 2:nrow(X)){
    foo = phi%*%Gamma[,,X$timeOfDay2[t]]%*%diag(allprobs[t,])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  return(-l)
}

### short version using the package Lcpp
mllk_short_t = function(theta.star, X, N){
  # transforming working parameters to natural
  ## parameters for state-dependent distributions
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N])
  ## regression coefs for tpm
  beta = matrix(theta.star[4*N+1:(3*N*(N-1))], ncol = 3)
  ## 24 unique tpms
  Gamma = tpm_p(tod = unique(X$timeOfDay), L=24, beta=beta)
  ## initial distribution -> periodically stationary
  delta = stationary_p(Gamma, X$timeOfDay2[1])
  # allprobs matrix of state-dependent densities
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  # forward algorithm to calculate the log-likelihood recursively
  -forward_p(delta, Gamma, allprobs, X$timeOfDay2)
}


# initial parameters for numerical optimisation
stepMean0 = c(0.1, 0.4, 1)
stepSd0 = c(0.1, 0.4, 1)
angleMean0 = rep(0,3)
angleConcentration0 = c(0.1, 0.9, 2)

# working scale initial parameter vector
theta.star0 = c(log(stepMean0), log(stepSd0), # steppars
                angleMean0, log(angleConcentration0), # anglepars
                rep(-2, 6), runif(6*2,-0.1,0.1)) # tpmpars

# fitting the model by numerically optimising the log-likelihood with nlm
mod_elephant2 = nlm(mllk_short_t, theta.star0, X = data, N = 3, 
                    print.level = 2, iterlim = 1000, hessian = TRUE)

# getting state process parameters
N = 3
beta = matrix(mod_elephant2$estimate[4*N+1:(3*N*(N-1))], ncol = 3)

# Visualising periodically stationary distribution
len = 300
todseq = seq(0,24,length=len)
Delta = matrix(NA, len, N)
for(t in 1:len){
  Gamma = tpm_p(tod = (todseq[t]+1:24) %% 24, beta = beta)
  Delta[t,] = stationary_p(Gamma, t = 1)
}

# Obtaining confidence bands via Monte Carlo
# drawing from the approximate normal distribution of the MLE
thetas = mvtnorm::rmvnorm(1000, mean=mod_elephant2$estimate, 
                          sigma=solve(mod_elephant2$hessian))
Deltas = array(dim = c(len, N, 1000))
for(i in 1:1000){
  cat("\n", i)
  for(t in 1:len){
    beta = matrix(thetas[i,4*N+1:(3*N*(N-1))], ncol = 3)
    Gamma = tpm_p(tod = (todseq[t]+1:24) %% 24, beta = beta)
    Deltas[t,,i] = stationary_p(Gamma, t = 1)
  }
}
DeltaCI = array(dim = c(len, N, 2))
for(j in 1:N){
  for(t in 1:len){
    DeltaCI[t,j,] = quantile(Deltas[t,j,], probs = c(0.025, 0.975))
  }
}

# pdf("./figs/hmm_timeOfDay.pdf", width = 6, height = 4)
# par(mfrow = c(1,1))
plot(NA, xlim = c(0,24), ylim = c(0,1), bty = "n", xaxt = "n",
     xlab = "time of day", ylab = "state occupancy probabilities")
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
for(j in 1:N){ 
  polygon(c(todseq, rev(todseq)), c(DeltaCI[,j,1], rev(DeltaCI[,j,2])), col = alpha(color[j], 0.2), border = F)
  lines(todseq, Delta[,j], lwd = 2, col = color[j]) 
}
legend("top", col = color, lwd = 2, legend = paste("state", 1:3), bty = "n")
dev.off()


